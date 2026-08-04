"""
Schlicht (1989) / Ludsteck VC estimator for random-walk time-varying coefficients.

Reference implementation aligned with VC.m (Mathematica) and Schlicht (2020) IZA DP 12920.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Union

import numpy as np
from scipy import optimize, sparse


@dataclass
class VCEstimateResult:
    """Schlicht VC estimation output (cf. VC.m return {coeff, sdb, sdu, sdi})."""

    coeff: np.ndarray            # (T, n) coefficient paths
    sdb: np.ndarray              # (T, n) coefficient standard deviations
    sdu: float                   # sqrt estimated Var(u)
    sdi: np.ndarray              # (n,) sqrt estimated Var(v_i)
    variance_ratios: np.ndarray  # (n,) r_i = sigma_i^2 / sigma_u^2
    residual_var: float          # Q / (T - n)
    criterion: float             # ML objective at optimum


def make_block_X(x: np.ndarray) -> sparse.csc_matrix:
    """Map (T, n) regressors to block-diagonal X (T, T*n)."""
    x = np.asarray(x, dtype=float)
    T, n = x.shape
    rows, cols, data = [], [], []
    for i in range(T):
        base = i * n
        for j in range(n):
            rows.append(i)
            cols.append(base + j)
            data.append(x[i, j])
    return sparse.csc_matrix((data, (rows, cols)), shape=(T, T * n))


def make_P(n: int, T: int) -> sparse.csc_matrix:
    """Random-walk operator P: (T-1)*n x T*n (VC.m makeP)."""
    rows, cols, data = [], [], []
    for i in range((T - 1) * n):
        coeff = i % n
        t = i // n
        rows.extend([i, i])
        cols.extend([t * n + coeff, (t + 1) * n + coeff])
        data.extend([-1.0, 1.0])
    return sparse.csc_matrix((data, (rows, cols)), shape=((T - 1) * n, T * n))


def make_S(inv_variance_ratios: np.ndarray, n: int, T: int) -> sparse.csc_matrix:
    """Diagonal S with entries 1/r_i repeated over time (VC.m makeS)."""
    diag = np.tile(inv_variance_ratios, T - 1)
    return sparse.diags(diag, format="csc")


def _parse_variance_ratios(
    variance_ratios: Union[float, Sequence[float], np.ndarray],
    n: int,
) -> np.ndarray:
    if np.isscalar(variance_ratios):
        return np.full(n, float(variance_ratios))
    r = np.asarray(variance_ratios, dtype=float).ravel()
    if r.size != n:
        raise ValueError(f"variance_ratios must have length {n}, got {r.size}")
    return r


def solve_coefficients(
    x: np.ndarray,
    y: np.ndarray,
    variance_ratios: np.ndarray,
) -> tuple[np.ndarray, sparse.csc_matrix, sparse.csc_matrix, sparse.csc_matrix, sparse.csc_matrix]:
    """
    Penalized LS solution for fixed variance ratios (paper eq. 4.9).

    Returns
    -------
    a : (T*n,) stacked coefficient vector
    M, X, P, S : system matrices
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    T, n = x.shape
    X = make_block_X(x)
    P = make_P(n, T)
    inv_r = 1.0 / variance_ratios
    S = make_S(inv_r, n, T)
    M = X.T @ X + P.T @ S @ P
    a = sparse.linalg.spsolve(M, X.T @ y)
    return a, M, X, P, S


def ml_criterion(
    variance_ratios: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
) -> float:
    """VC.m critfun: log|M| + (T-n) log Q + (T-1) sum log r."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    T, n = x.shape
    r = np.asarray(variance_ratios, dtype=float)
    a, M, X, P, S = solve_coefficients(x, y, r)
    Pa = P @ a
    u = y - X @ a
    Q = float(u @ u + Pa @ S @ Pa)
    if Q <= 0:
        return np.inf
    sign, logdet = np.linalg.slogdet(M.toarray())
    if sign <= 0:
        return np.inf
    return logdet + (T - n) * np.log(Q) + (T - 1) * np.sum(np.log(r))


def vc_estimate(
    x: np.ndarray,
    y: np.ndarray,
    variance_ratios: Union[float, Sequence[float], np.ndarray] = 1.0,
    optimize_ratios: bool = True,
    maxiter: int = 500,
) -> VCEstimateResult:
    """
    Estimate Schlicht VC model y_t = x_t' b_t + u_t, b_t = b_{t-1} + v_t.

    Parameters
    ----------
    x : (T, n) regressor matrix (no implicit intercept).
    y : (T,) or (T, 1) dependent variable.
    variance_ratios : initial r_i = Var(v_i)/Var(u); scalar expands to all coefs.
    optimize_ratios : if True, minimize ML criterion over r (VC.m NMinimize).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    if x.ndim != 2:
        raise ValueError("x must be a 2D array")
    T, n = x.shape
    if y.shape[0] != T:
        raise ValueError("x and y must have the same number of rows")
    if T <= n:
        raise ValueError("Need T > n for VC estimation")

    r0 = _parse_variance_ratios(variance_ratios, n)
    if optimize_ratios:
        bounds = [(1e-10, 1e6)] * n
        starts = []
        for mult in (0.5, 1.0, 1.5):
            starts.append(np.clip(r0 * mult, 1e-10, None))
        best = None
        for x0 in starts:
            res = optimize.minimize(
                ml_criterion,
                x0,
                args=(x, y),
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": maxiter},
            )
            if best is None or res.fun < best.fun:
                best = res
        r_opt = best.x
        crit = float(best.fun)
    else:
        r_opt = r0.copy()
        crit = ml_criterion(r_opt, x, y)

    a, M, X, P, S = solve_coefficients(x, y, r_opt)
    Pa = P @ a
    u = y - X @ a
    Q = float(u @ u + Pa @ S @ Pa)
    vu = Q / (T - n)

    Minv = np.linalg.inv(M.toarray())
    sdb_flat = np.sqrt(np.clip(vu * np.diag(Minv), 0.0, None))
    coeff = a.reshape(T, n)
    sdb = sdb_flat.reshape(T, n)

    return VCEstimateResult(
        coeff=coeff,
        sdb=sdb,
        sdu=float(np.sqrt(vu)),
        sdi=np.sqrt(vu * r_opt),
        variance_ratios=r_opt,
        residual_var=vu,
        criterion=crit,
    )
