"""Simulate data from a Schlicht VC random-walk coefficient model."""

from __future__ import annotations

import numpy as np

from ga.vc.schlicht import make_block_X, make_P


def simulate_vc_model(
    T: int,
    n: int,
    sigma_u: float,
    sigma_v: np.ndarray,
    *,
    x: np.ndarray | None = None,
    b0: np.ndarray | None = None,
    seed: int = 123,
) -> dict:
    """
    Simulate y_t = x_t' b_t + u_t with b_t = b_{t-1} + v_t.

    If x is omitted, uses x_t = [1, z_t] where z_t ~ N(0, I_{n-1}).
    Returns dict with x, y, b (true paths), u, v, and true variance ratios.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    sigma_v = np.broadcast_to(np.asarray(sigma_v, dtype=float), (n,))

    rng = np.random.default_rng(seed)
    if x is None:
        if n == 1:
            x = np.ones((T, 1))
        else:
            z = rng.standard_normal((T, n - 1))
            x = np.column_stack([np.ones(T), z])

    x = np.asarray(x, dtype=float)
    if x.shape != (T, n):
        raise ValueError(f"x must have shape ({T}, {n})")

    if b0 is None:
        b0 = np.ones(n)
    b = np.zeros((T, n))
    b[0] = b0
    v = np.zeros((T - 1, n))
    for t in range(1, T):
        v[t - 1] = rng.normal(0.0, sigma_v)
        b[t] = b[t - 1] + v[t - 1]

    u = rng.normal(0.0, sigma_u, size=T)
    y = (x * b).sum(axis=1) + u

    return {
        "x": x,
        "y": y,
        "b": b,
        "u": u,
        "v": v,
        "sigma_u": sigma_u,
        "sigma_v": sigma_v,
        "variance_ratios": (sigma_v ** 2) / (sigma_u ** 2),
    }


def fitted_values(coeff: np.ndarray, x: np.ndarray) -> np.ndarray:
    """In-sample fitted y from coefficient path (T, n) and x (T, n)."""
    return (x * coeff).sum(axis=1)


def stacked_state(coeff: np.ndarray) -> np.ndarray:
    """Flatten (T, n) coefficient path to VC stacked vector a."""
    return coeff.reshape(-1, order="C")


def random_walk_residuals(coeff: np.ndarray) -> np.ndarray:
    """Pa where P is the VC random-walk operator."""
    T, n = coeff.shape
    a = stacked_state(coeff)
    P = make_P(n, T)
    return (P @ a).reshape(T - 1, n)
