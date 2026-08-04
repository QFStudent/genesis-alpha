import numpy as np


def calc_tv_coefficients(X: np.ndarray, y: np.ndarray, sigma_i: float) -> np.ndarray:
    """Gaussian-kernel weighted least squares time-varying coefficients."""
    n, k = X.shape
    y = np.ravel(y)
    coefs = np.zeros((n, k))
    t = np.arange(n, dtype=float)
    for i in range(n):
        weights = np.exp(-0.5 * ((t - i) / sigma_i) ** 2)
        sw = np.sqrt(weights)
        Xw = X * sw[:, None]
        yw = y * sw
        coefs[i], _, _, _ = np.linalg.lstsq(Xw, yw, rcond=None)
    return coefs
