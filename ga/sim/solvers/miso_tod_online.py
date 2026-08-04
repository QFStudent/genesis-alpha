"""
Online two-system MISO TOD solver (prod miso_tod_online.py).
"""
import pandas as pd
import numpy as np
from scipy import sparse
from datetime import datetime

from ga.filters.tib import TIBSystem


def _rls_update(z: np.ndarray, y: float, C: np.ndarray, P: np.ndarray, lam: float):
    zcol = np.asarray(z, dtype=float).reshape(-1, 1)
    P = P / lam
    denom = 1.0 + float(zcol.T @ P @ zcol)
    k = (P @ zcol) / denom
    err = y - float(C @ zcol)
    C = C + err * k.T
    P = P - k @ zcol.T @ P
    return C, P


class MISO_TWOSYS_Online(object):
    def __init__(self,
                 p1: int,
                 p2: int,
                 q: int,
                 model: str = "inverse",
                 phi1: float = 1.0,
                 phi2: float = 1.0,
                 phi1_step: int = 1,
                 phi2_step: int = 1,
                 R01: np.ndarray = None,
                 R02: np.ndarray = None,
                 tod_index: np.ndarray = None) -> None:
        self.p1 = p1
        self.p2 = p2
        self.q = q
        self.model = model
        self.phi1 = phi1
        self.phi2 = phi2
        self.phi1_step = phi1_step
        self.phi2_step = phi2_step
        self.R01 = R01
        self.R02 = R02
        self.tod_index = tod_index

    def generate_states_offline(self, sys: TIBSystem, y: np.ndarray):
        if not isinstance(y, np.ndarray):
            raise ValueError('Input should be an ndarray.')

        z = np.empty((len(y)+1, sys.p))
        z[0, :] = np.zeros(sys.p)
        Minv = sparse.linalg.inv(sys.M)
        MBy = y.dot(sys.U.T)

        for ii, _ in enumerate(y):
            zii = sys.N.dot(sys.Q.dot(z[ii, :]))
            zii += MBy[ii, :]
            z[ii+1, :] = Minv.dot(zii)
        return z

    def _fit(self, sys1: TIBSystem, sys2: TIBSystem,
             x: pd.DataFrame,
             y: pd.DataFrame):
        z1_path = self.generate_states_offline(sys1, x.values)
        z2_path = self.generate_states_offline(sys2, x.values)
        z1 = z1_path[0].copy()
        z2 = z2_path[0].copy()

        C1 = np.zeros((self.q, self.p1))
        C2 = np.zeros((self.q, self.p2))
        P1 = np.array(self.R01, dtype=float, copy=True)
        P2 = np.array(self.R02, dtype=float, copy=True)

        yhat = np.zeros(y.shape)
        yres = np.zeros(y.shape)
        ybar = np.zeros(y.shape)
        C1_norm = np.zeros(y.shape[0])
        C2_norm = np.zeros(y.shape[0])

        Minv1 = sparse.linalg.inv(sys1.M)
        Minv2 = sparse.linalg.inv(sys2.M)

        i = 0
        for ii, _ in x.iterrows():
            tod = ii.time() if hasattr(ii, 'time') else x.index[i].time()
            yt = float(y.loc[ii].values.ravel()[0])

            if str(tod) == self.tod_index[0]:
                yhat[i] = (C1 @ z1.reshape(-1, 1)).ravel()
                lam = self.phi1 if (i % self.phi1_step == 0) else 1.0
                C1, P1 = _rls_update(z1, yt, C1, P1, lam)
            else:
                yhat[i] = (C2 @ z2.reshape(-1, 1)).ravel()
                lam = self.phi2 if (i % self.phi2_step == 0) else 1.0
                C2, P2 = _rls_update(z2, yt, C2, P2, lam)

            eps = x.loc[[ii]].values.reshape(-1, 1)

            zii = sys1.N.dot(sys1.Q.dot(z1.reshape(-1, 1)))
            zii += sys1.U.dot(eps)
            z1 = Minv1.dot(zii).ravel()

            zii = sys2.N.dot(sys2.Q.dot(z2.reshape(-1, 1)))
            zii += sys2.U.dot(eps)
            z2 = Minv2.dot(zii).ravel()
            i += 1

        self.sys1 = sys1
        self.sys2 = sys2
        self.C1 = C1
        self.C2 = C2
        self.zhat1 = z1_path
        self.zhat2 = z2_path
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
