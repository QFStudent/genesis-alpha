from typing import List

import numpy as np
import pandas as pd

from ga.filters import tib
from ga.sim.tib_model_helpers import mimo_system_matrices_biased
from ga.sim.solvers.miso_tod_online import MISO_TWOSYS_Online

try:
    from fracdiff.sklearn import Fracdiff
except ImportError:  # optional prod dependency
    Fracdiff = None


class TwoSysMISO:
    def __init__(
        self,
        numPoles1: int,
        numPoles2: int,
        penalty: float,
        a: float,
        b: float,
        numInputs: int,
        tod_idx: List[str],
        phi1: float,
        phi2: float,
        phi1_step: int,
        phi2_step: int,
    ):
        self.numPoles1 = numPoles1
        self.numPoles2 = numPoles2
        self.penalty = penalty
        self.a = a
        self.b = b
        self.numInputs = numInputs
        self.tod_idx = tod_idx
        self.phi1 = phi1
        self.phi2 = phi2
        self.phi1_step = phi1_step
        self.phi2_step = phi2_step

    def extract_tod(self, df: pd.DataFrame) -> np.ndarray:
        return pd.unique(df.index.time)

    def transform_feature(self, feature: pd.DataFrame) -> pd.DataFrame:
        if Fracdiff is None:
            raise ImportError("two_sys_miso_fracdiff requires the fracdiff package")
        cum_feature = feature.cumsum(axis=0)
        f = Fracdiff(0.9)
        for col in feature.columns:
            if col == "FEATURE":
                cum_feature[col] = feature[col].values
            else:
                cum_feature[col] = f.fit_transform(
                    cum_feature[col].values.reshape(-1, 1)
                ).ravel()
        return cum_feature

    def fit_train(self, train: pd.DataFrame, test: pd.DataFrame) -> tuple:
        lambd1 = tib.poles_chebyshev_roots(self.numPoles1, a=self.a, b=self.b)
        lambd2 = tib.poles_chebyshev_roots(self.numPoles2, a=self.a, b=self.b)
        nulls1 = tib.mimo_standard_null_vectors(self.numPoles1, self.numInputs)
        nulls2 = tib.mimo_standard_null_vectors(self.numPoles2, self.numInputs)
        sys1 = mimo_system_matrices_biased(lambd1, nulls1.T)
        sys2 = mimo_system_matrices_biased(lambd2, nulls2.T)

        solver_inst = MISO_TWOSYS_Online(
            p1=self.numPoles1,
            p2=self.numPoles2,
            q=1,
            model="inverse",
            phi1=self.phi1,
            phi2=self.phi2,
            phi1_step=self.phi1_step,
            phi2_step=self.phi2_step,
            R01=np.eye(self.numPoles1),
            R02=np.eye(self.numPoles2),
            tod_index=self.tod_idx,
        )

        ext_df = pd.concat([train, test], axis=0)
        frac_df = self.transform_feature(ext_df)
        x = frac_df.loc[train.index, :]
        x = x.iloc[:, 1:]
        y = ext_df.loc[train.index, ["Target"]]
        test_feature = frac_df.loc[test.index, :]
        test_feature = test_feature.iloc[:, 1:]

        yhat, yres, ybar, C1, C2, C1_norm, C2_norm = solver_inst._fit(sys1, sys2, x, y)
        self.C1 = C1
        self.C2 = C2
        self.sys1 = sys1
        self.sys2 = sys2
        self.zhat01 = solver_inst.zhat1[-1]
        self.zhat02 = solver_inst.zhat2[-1]
        self.C = C1.ravel()  # place holder
        self.zhat0 = self.zhat01.copy()  # place holder
        self.sys = self.sys1
        return solver_inst, test_feature

    def predict(self, solver_inst: MISO_TWOSYS_Online, test_feature: pd.DataFrame):
        yhat = solver_inst._predict(
            test_feature,
            self.zhat01.reshape(-1, 1),
            self.zhat02.reshape(-1, 1),
        )
        return yhat

    def run_epoch(self, train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
        solver_inst, test_feature = self.fit_train(train, test)
        yhat = self.predict(solver_inst, test_feature)
        return yhat
