"""
Huber regression estimation for
MIMO model with TOD mechanism
"""
import pandas as pd
import numpy as np
from scipy import sparse
from datetime import datetime

from sklearn.linear_model import LinearRegression
import pymc as pm
import pytensor

from ga.filters.tib import TIBSystem


class MISOTODBAYESIAN(object):
    def __init__(self,
                 p1: int,
                 p2: int,
                 q: int,
                 model: str = "inverse",
                 phi: float = 1.0,
                 tod_index: np.ndarray = None) -> None:
        self.p1 = p1
        self.p2 = p2
        self.q = q
        self.model = model
        self.phi = phi
        self.tod_index = tod_index

    def generate_states_offline(self, sys: TIBSystem,
                                y: np.ndarray):
        if not isinstance(y, np.ndarray):
            raise ValueError('Input should be an ndarray.')

        z = np.empty((len(y)+1, sys.p))
        z[0, :] = np.zeros(sys.p)
        Minv = sparse.linalg.inv(sys.M)
        MBy = y.dot(sys.U.T)        # if it's inverse system

        for ii, _ in enumerate(y):
            zii = sys.N.dot(sys.Q.dot(z[ii, :]))
            zii += MBy[ii, :]       # propagate z_t+1 with z_t and y_t
            z[ii+1, :] = Minv.dot(zii)
        return z        # output z with the same timestamp of y

    def bayesian_fit(self, C: np.ndarray, z: pd.DataFrame,
                     y: pd.DataFrame):
        N = len(y)
        sections = int(np.floor(len(y) / 60))    # one section per quarter
        t_section = np.zeros(len(y)).astype(int)
        tmp_section = np.repeat(np.arange(1, sections+1), N / sections)
        t_section[-len(tmp_section):] = tmp_section
        t_section_t = pytensor.shared(t_section)

        z_t = pytensor.shared(z.values)
        dim = C.shape[1]
        mus = C.ravel()
        sigmas = np.ones(dim)

        with pm.Model() as model:
            # define priors
            chol_beta, *_ = pm.LKJCholeskyCov(
                "chol_cov_alpha", n=dim, eta=2, sd_dist=pm.HalfCauchy.dist(2.5), compute_corr=True
            )

            beta = pm.MvGaussianRandomWalk(
                "beta", mu=np.zeros(dim), chol=chol_beta, shape=(sections+1, dim)
            )

            beta_r = beta[t_section_t]

            sigma = pm.HalfCauchy(name="sigma", beta=10)
            nu = pm.Gamma(name="nu", alpha=8, beta=2)

            mu = (z_t * beta_r).sum(axis=1)
            y = y.values.ravel()

            likelihood = pm.StudentT(name="likelihood", mu=mu, nu=nu, sigma=sigma, observed=y)

            with model:
                idata = pm.sample(1000, cores=4)

        last_coef = idata.posterior['beta'].values[:, :, -1, :]
        coef = last_coef.mean(axis=0).mean(axis=0).reshape(1, -1)
        return coef

    def _fit(self, sys1: TIBSystem, sys2: TIBSystem,
             x: pd.DataFrame,
             y: pd.DataFrame):
        z1 = self.generate_states_offline(sys1, x.values)
        z2 = self.generate_states_offline(sys2, x.values)

        z1_df = pd.DataFrame(z1[:-1], index=y.index)
        z2_df = pd.DataFrame(z2[:-1], index=y.index)

        tod_time = datetime.strptime(self.tod_index[0], '%H:%M:%S').time()
        filt = (y.index.time == tod_time)     # separate over night return estimation with the rest

        r1_z1_df = z1_df.loc[filt]
        r1_y_df = y.loc[filt]

        r2_z2_df = z2_df.loc[~filt]
        r2_y_df = y.loc[~filt]

        C1 = np.zeros((self.q, self.p1))
        C2 = np.zeros((self.q, self.p2))

        yhat = np.zeros(y.shape)
        yres = np.zeros(y.shape)
        ybar = np.zeros(y.shape)
        C1_norm = np.zeros(y.shape[0])
        C2_norm = np.zeros(y.shape[0])

        lm1 = LinearRegression(fit_intercept=False)
        lm2 = LinearRegression(fit_intercept=False)

        lm1.fit(r1_z1_df.values, r1_y_df.values)
        lm2.fit(r2_z2_df.values, r2_y_df.values)

        C1 = self.bayesian_fit(lm1.coef_.reshape(1, -1),
                               r1_z1_df, r1_y_df)
        C2 = lm2.coef_.reshape(1, -1)

        self.sys1 = sys1
        self.sys2 = sys2
        self.C1 = C1
        self.C2 = C2

        self.zhat1 = z1[-1]
        self.zhat2 = z2[-1]
        return yhat, yres, ybar, C1, C2, C1_norm, C2_norm

    def _predict(self, x: pd.DataFrame,
                 z1: np.ndarray,
                 z2: np.ndarray) -> np.ndarray:

        yhat = np.zeros((x.shape[0], 1))
        Minv1 = sparse.linalg.inv(self.sys1.M)
        Minv2 = sparse.linalg.inv(self.sys2.M)

        i = 0
        for ii, xii in x.iterrows():
            tod = xii.name.time()
            if str(tod) == self.tod_index[0]:
                yhat[i] = (self.C1.dot(z1)).ravel()
            else:
                yhat[i] = (self.C2.dot(z2)).ravel()

            eps = x.loc[[ii]].values.reshape(-1, 1)

            zii = self.sys1.N.dot(self.sys1.Q.dot(z1))
            zii += self.sys1.U.dot(eps)
            z1 = Minv1.dot(zii)

            zii = self.sys2.N.dot(self.sys2.Q.dot(z2))
            zii += self.sys2.U.dot(eps)
            z2 = Minv2.dot(zii)
            i += 1

        return yhat
