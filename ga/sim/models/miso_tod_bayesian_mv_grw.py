import logging
from typing import List

import numpy as np
import pandas as pd

from ga.filters import tib
from ga.filters.tib import TIBSystem
from ga.sim.tib_model_helpers import mimo_system_matrices_biased
from ga.sim.solvers import miso_tod_bayesian_mv_grw

_logger = logging.getLogger(__name__)


class BAYESIAN_TOD_MISO:
    """
    TOD MIMO model estimated using Huber regression
    """

    def __init__(self, numPoles1: int, numPoles2: int,
                 penalty: float, a: float, b: float,
                 numInputs: int,
                 tod_idx: List[str]):
        self.numPoles1 = numPoles1
        self.numPoles2 = numPoles2
        self.penalty = penalty
        self.a = a
        self.b = b
        self.numInputs = numInputs
        self.tod_idx = tod_idx

    def extract_tod(self, df: pd.DataFrame) -> np.ndarray:
        return pd.unique(df.index.time)

    def fit_train(self, train: pd.DataFrame) -> TIBSystem:
        assert train.columns[0] == 'Target', 'First column should be Target.'

        x = train.iloc[:, 1:]
        y = train[["Target"]]

        lambd1 = tib.poles_chebyshev_roots(self.numPoles1, a=self.a, b=self.b)
        lambd2 = tib.poles_chebyshev_roots(self.numPoles2, a=self.a, b=self.b)
        nulls1 = tib.mimo_standard_null_vectors(self.numPoles1, self.numInputs)
        nulls2 = tib.mimo_standard_null_vectors(self.numPoles2, self.numInputs)
        sys1 = mimo_system_matrices_biased(lambd1, nulls1.T)
        sys2 = mimo_system_matrices_biased(lambd2, nulls2.T)

        solver_inst = miso_tod_bayesian_mv_grw.MISOTODBAYESIAN(p1=self.numPoles1,
                                                              p2=self.numPoles2,
                                                              q=1,
                                                              model="inverse",
                                                              phi=1,
                                                              tod_index=self.tod_idx)

        yhat, yres, ybar, C1, C2, C1_norm, C2_norm = solver_inst._fit(sys1, sys2, x, y)
        self.C1 = C1
        self.C2 = C2
        self.sys1 = sys1
        self.sys2 = sys2
        self.zhat01 = solver_inst.zhat1
        self.zhat02 = solver_inst.zhat2
        self.C = C1.ravel()        # place holder
        self.zhat0 = self.zhat01.copy()  # place holder
        self.sys = self.sys1
        return solver_inst

    def predict(self, solver, raw_test: pd.DataFrame):
        x = raw_test.iloc[:, 1:]
        yhat = solver._predict(x,
                               self.zhat01.reshape(-1, 1),
                               self.zhat02.reshape(-1, 1))
        return yhat

    def run_epoch(self, train: pd.DataFrame,
                  test: pd.DataFrame) -> np.ndarray:
        assert test.index[0] > train.index[-1], 'Look ahead issue.'
        solver = self.fit_train(train)
        yhat = self.predict(solver, test)
        return yhat
