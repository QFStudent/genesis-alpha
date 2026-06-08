import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import LinearOperator, eigsh
from numpy.fft import fft, ifft


class FastHankelProduct(LinearOperator):
    """
    Parameters
    ----------
    h: impulse response without the zero lag term, np.ndarray, shape = (n, )
    """

    def __init__(self, h: np.ndarray, dtype='float32'):
        self.h = h
        self.shape = (len(h), len(h))
        self.dtype = np.dtype(dtype)
        self.n = len(self.h)

    def _matvec(self, x):
        """
        Parameters
        x: np.ndarray, shape = (n, )
        """
        chat = np.hstack([self.h[-1], np.zeros(self.n-1), self.h[:-1]])
        Pi = np.eye(self.n)[::-1]    # permutation matrix
        xhat = np.hstack([Pi.dot(x), np.zeros(self.n-1)])
        y = ifft(fft(chat)*fft(xhat))
        return y[:self.n].real


def reduce_fft_truncate(h: np.ndarray, p: int):
    """
    Model reduction using fast hankel matrix vector multiplication
    that can be calculated by FFT

    Parameters
    ----------
    h: array-like, shape = (n, )
    p:  int, the reduced dimension of our model

    Returns
    -------
    lambd: array-like, shape = (p, )
    """
    fft_operator = FastHankelProduct(h)
    evals_large, evecs_large = eigsh(fft_operator, p, which='LM')
    j = np.argsort(np.abs(evals_large))[::-1]
    v, r = la.qr(evecs_large[:, j])
    A = np.dot(v.T[:p, 1:], v[:-1, :p])
    lambd, _ = la.eig(A)
    return lambd


def reduce_svd_truncate(h: np.ndarray, order: int):
    """
    Reduce the system based on the given number of poles
    using svd method
    """
    if h.ndim == 2:
        raise ValueError('h needs to be raveled.')

    n = len(h)
    Gamma = la.hankel(h)
    u, s, v = la.svd(Gamma)

    A = np.dot(v[:order, 1:], v.T[:-1, :order])
    A = np.real(A)

    lambd, _ = la.eig(A)
    return lambd


def reduce_svd_truncate_fixed(h: np.ndarray, order: int):
    """
    Corrected SISO Ho-Kalman / ERA pole recovery (companion to reduce_svd_truncate).

    reduce_svd_truncate has two bugs; this version fixes both:
      1. it builds ``la.hankel(h)`` -- a *triangular* Hankel (lower-right triangle
         zeroed) that silently discards the late impulse-response samples. Harmless
         for a fast-decaying IR, but it corrupts pole recovery when the tail carries
         energy. Here we populate the full square Hankel ``Gamma[i, j] = h[i + j]``.
      2. it estimates the shift with a transpose, ``W[:, 1:] @ W[:, :-1].T``; the
         truncated ``W`` has lost column-orthonormality, so the right inverse of
         ``W[:, :-1]`` is the pseudo-inverse, not the transpose. Here
         ``A = W[:, 1:] @ pinv(W[:, :-1])``.

    Requires at least ``2*order - 1`` impulse-response samples.

    Parameters
    ----------
    h: array-like, shape (N,)   positive-lag impulse response (raveled)
    order: int                  number of poles to recover

    Returns
    -------
    lambd: array-like, shape (order,)   recovered poles (may be complex)
    """
    if h.ndim == 2:
        raise ValueError('h needs to be raveled.')

    h = np.asarray(h, dtype=float).ravel()
    N = len(h)
    if N < 2 * order - 1:
        raise ValueError(f"need >= {2 * order - 1} samples for order {order}, got {N}.")

    m = (N + 1) // 2
    Gamma = la.hankel(h[:m], h[m - 1:2 * m - 1])      # proper square Hankel
    _, _, Vh = la.svd(Gamma)
    W = Vh[:order, :]    # input-balanced controllability factor (orthonormal rows)

    A = np.real(W[:, 1:].dot(np.linalg.pinv(W[:, :-1])))
    lambd, _ = la.eig(A)
    return lambd


def toeplog(h: np.ndarray):
    """
    Convert impulse response to the power series coefficients
    of log transfer function.

    Note: in this calculation, h has h_0 which is 1, and a_0
    should be ignored in the reduction step.
    """
    num_ir = len(h)
    d = np.arange(num_ir-1) + 1
    b = d * h[1:]

    c = h[:-1]
    r = np.zeros(num_ir-1)
    r[0] = h[0]
    x = la.solve_toeplitz((c, r), b)
    a_ = x / d

    a = np.insert(a_, 0, np.log(h[0]), axis=0)
    return a


def info_svd_reduce(h: np.ndarray, order: int, rho: float) -> np.ndarray:
    """
    System reduction by minimizing the log transfer function

    Here h is the positive lag impulse response
    """
    if h.ndim == 2:
        raise ValueError('h needs to be raveled.')

    h_ = np.insert(h, 0, 1, axis=0)    # convert to unity zero lag response
    hr = h_ * (rho ** (np.arange(0, len(h_))))
    a = toeplog(h_)
    lambd = reduce_svd_truncate(a[1:], order) / rho
    return lambd


def info_svd_reduce_fixed(h: np.ndarray, order: int, rho: float) -> np.ndarray:
    """
    Corrected cepstrum / information-distance pole estimation (companion to
    info_svd_reduce).

    Fixes two problems in info_svd_reduce:
      1. it computes a rho-damped IR ``hr`` but then runs ``toeplog`` on the
         *un-damped* ``h_`` -- the damping is silently dropped. Here we use
         ``toeplog(hr)``.
      2. it feeds the raw cepstrum ``a_k`` to the pole-finder. For a rational system
         ``a_k = (1/k) * sum_i lam_i**k``; the ``1/k`` makes ``a_k`` NOT a low-rank
         Hankel sequence, so Ho-Kalman cannot recover the poles from it. Multiplying
         by ``k`` gives the power sums ``g_k = k * a_k = sum_i lam_i**k`` -- a sum of
         geometrics from which the shift recovers the poles (Prony).

    The damping ``rho**k`` (``rho`` in (0, 1]) pulls near-unit-root poles inward for
    conditioning; recovered poles are rescaled by ``1/rho``.

    LIMITATION: reliable for *identification* (``order`` ~ true system order). It does
    NOT give numerically stable *order reduction* (``order`` << true order): the
    power-sum Hankel has no energy ordering (every pole enters with residue 1), so
    aggressive truncation can return spurious / out-of-disk poles. For genuine
    reduction use the impulse-response-domain ``msvdreduce``.

    Parameters
    ----------
    h: array-like, shape (N,)   positive-lag impulse response (raveled)
    order: int                  number of poles
    rho: float                  damping factor in (0, 1]

    Returns
    -------
    lambd: array-like, shape (order,)   recovered poles
    """
    if h.ndim == 2:
        raise ValueError('h needs to be raveled.')

    h_ = np.insert(h, 0, 1, axis=0)                # unity zero-lag term
    hr = h_ * (rho ** np.arange(len(h_)))          # geometric damping (now applied)
    a = toeplog(hr)                                # cepstrum of the damped IR
    g = np.arange(len(a)) * a                      # g_k = k * a_k = power sums sum lam**k
    lambd = reduce_svd_truncate_fixed(g[1:], order) / rho
    return lambd


def blkhankel(ir: np.ndarray): 
    """
    Construct block hankel matrix from given 
    impulse response. 
    """
    di, do, l = ir.shape 
    k = int(np.floor((l + 1) / 2))
    H = np.zeros((di*k, do*k))
    for i in range(k): 
        for j in range(k): 
            H[i*di:(i+1)*di, j*do:(j+1)*do] = ir[:, :, i+j]
    return H 


def lschur(A: np.ndarray): 
    # Shur triangularization and convert A to lower 
    T, Z = la.schur(A, output='real')
    L = T[::-1, ::-1]
    Q = Z[:, ::-1]
    return L, Q 


def msvdreduce(ir: np.ndarray, order: int): 
    do, di, num_ir = ir.shape 
    H = blkhankel(ir)

    U, S, Vh = la.svd(H)    # U.dot(S).dot(VT)
    Vh = Vh[:order, :]      # order reduction by truncation 
    # A = Vh[:, di:].dot(Vh.T[:-di, :])
    A = Vh[:, di:].dot(np.linalg.pinv(Vh[:, :-di]))
    
    A, Q = lschur(A)
    B = Q.T.dot(Vh[:, :di])
    w = np.diag(A)
    u = np.zeros((di, order, 1))

    for i in range(order):
        ui = (B[0] / np.linalg.norm(B[0], ord=2)).reshape(-1, 1)
        u[:, i] = ui
        tmp = np.eye(di) - (1 + np.conj(w[i]))*ui.dot(ui.T)
        B = B[1:, :]
    return w, u[:, :, 0]


def msvdreduce_fixed(ir: np.ndarray, order: int):
    """
    Corrected Hankel-SVD / ERA recovery (companion to msvdreduce).

    msvdreduce returns the correct poles but its null-vector loop is unfinished: it
    computes the inverse deflation and then discards it (``tmp`` is never applied,
    ``B = B[1:, :]`` takes the raw rows), so the returned null vectors are the
    un-deflated rows of B. This version applies the deflation
    ``J_i^{-1} = I - (1 + 1/conj(w_i)) u_i u_i*`` to the remaining rows (Yu eqs
    452-453), so it returns the orthonormal null vectors y_k.

    To rebuild a TIB pair from the result, pass ``(w, u)`` to
    ``ga.filters.tib.build_tib_from_directions`` -- NOT to ``null_basis_realization``,
    which treats its input as *raw* null vectors and would re-apply the Blaschke
    deflation to the already-orthonormalized u_k (double-deflation -> wrong system).

    Recovery is exact at full order (``order`` = system order). For ``order`` < system
    order the poles remain meaningful but the recovered null vectors are approximate.
    Real poles only: like ``msvdreduce`` the poles are read off the (real Schur)
    diagonal; for complex-conjugate pairs use ``balanced_truncation`` +
    ``tib_from_state_space`` instead.

    Parameters:
        ir: impulse response, shape (p, q, L) = (outputs, inputs, lags).
        order: number of poles to recover.

    Returns:
        w: poles, shape (order,).
        u: orthonormal null vectors, shape (q, order) -- columns are unit vectors.
    """
    do, di, num_ir = ir.shape
    H = blkhankel(ir)

    U, S, Vh = la.svd(H)
    Vh = Vh[:order, :]
    A = Vh[:, di:].dot(np.linalg.pinv(Vh[:, :-di]))   # shift: pseudo-inverse, not transpose
    A, Q = lschur(A)
    B = Q.T.dot(Vh[:, :di]).astype(complex)
    w = np.diag(A).astype(complex)
    u = np.zeros((di, order), dtype=complex)

    for i in range(order):
        ui = B[0] / np.linalg.norm(B[0], ord=2)
        u[:, i] = ui
        if B.shape[0] > 1:                            # apply J_i^{-1} to the remaining rows
            Jinv = np.eye(di) - (1.0 + 1.0/np.conj(w[i])) * np.outer(ui, np.conj(ui))
            B = B[1:, :].dot(Jinv)
        else:
            B = B[1:, :]
    return w, u


