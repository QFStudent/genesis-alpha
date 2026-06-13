import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import LinearOperator, eigsh, svds
from numpy.fft import fft, ifft


class FastHankelProduct(LinearOperator):
    """
    Parameters
    ----------
    h: impulse response without the zero lag term, np.ndarray, shape = (n, )
    dtype: declared operator dtype. Defaults to float64: the FFT matvec computes in
        float64 (it returns ``y.real``), so declaring float32 made iterative solvers
        such as eigsh iterate in single precision and lose ~6 digits -- the historical
        reduce_fft_truncate accuracy bug.
    """

    def __init__(self, h: np.ndarray, dtype='float64'):
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


def _scalar_hankel_matvec_fft(g, x, k):
    """
    Fast scalar Hankel matvec via FFT: returns H @ x where H[i, j] = g[i + j]
    (i, j = 0..k-1), using a length-(3k-2) linear convolution. O(k log k).

    g must hold the 2k-1 entries g[0..2k-2] (fully populated, not triangular).
    See docs/derivations/05-fast-hankel-matvec.md.
    """
    g = np.asarray(g, dtype=float)[:2 * k - 1]
    xr = np.asarray(x, dtype=float)[::-1]
    n = (2 * k - 1) + k - 1                      # length of the linear convolution
    nfft = 1 << (max(n, 1) - 1).bit_length()     # next power of two
    conv = np.fft.irfft(np.fft.rfft(g, nfft) * np.fft.rfft(xr, nfft), nfft)
    return conv[k - 1:2 * k - 1]


class BlockHankelOperator(LinearOperator):
    """
    Fully-populated block Hankel matrix as a fast LinearOperator (no dense storage).

    Represents H of shape (p*k, q*k) with p x q blocks  H[i, j] = ir[:, :, i+j],
    matching blkhankel. matvec/rmatvec run in O(p q k log k) via one FFT convolution
    per (output, input) channel pair (Yu eq 441); see
    docs/derivations/05-fast-hankel-matvec.md.

    Unlike FastHankelProduct (which is SISO and encodes the *triangular* single-arg
    Hankel), this uses the full Markov sequence ir[a, b, 0..2k-2].

    Parameters
    ----------
    ir: impulse response, shape (p, q, L) = (outputs, inputs, lags).
    k:  number of block rows/columns (needs 2k-1 <= L).
    """
    def __init__(self, ir, k):
        self.ir = np.asarray(ir, dtype=float)
        self.p, self.q, self.L = self.ir.shape
        self.k = k
        super().__init__(dtype=np.float64, shape=(self.p * k, self.q * k))

    def _matvec(self, x):
        k, p, q = self.k, self.p, self.q
        xb = np.asarray(x, dtype=float).reshape(k, q)
        yb = np.zeros((k, p))
        for a in range(p):
            for b in range(q):
                yb[:, a] += _scalar_hankel_matvec_fft(self.ir[a, b], xb[:, b], k)
        return yb.reshape(-1)

    def _rmatvec(self, y):
        k, p, q = self.k, self.p, self.q
        yb = np.asarray(y, dtype=float).reshape(k, p)
        xb = np.zeros((k, q))
        for a in range(p):
            for b in range(q):
                xb[:, b] += _scalar_hankel_matvec_fft(self.ir[a, b], yb[:, a], k)
        return xb.reshape(-1)


def msvdreduce_fast(ir: np.ndarray, order: int):
    """
    Fast (FFT/Lanczos) Hankel-SVD recovery -- companion to msvdreduce_fixed.

    Same algorithm as msvdreduce_fixed (pseudo-inverse shift + finished null-vector
    deflation), but the rank-`order` partial SVD of the block Hankel is computed by
    Lanczos (scipy.sparse.linalg.svds) on a BlockHankelOperator -- the matrix is never
    formed densely. This is Yu's "Fast Partial Block Hankel SVD" (sec 6.2); it scales
    to long impulse responses where the dense blkhankel + la.svd of msvdreduce would
    not. See docs/derivations/05-fast-hankel-matvec.md.

    Returns (w, u) exactly like msvdreduce_fixed: poles and orthonormal null vectors.
    Rebuild a TIB pair with ga.filters.tib.build_tib_from_directions (NOT
    null_basis_realization). Real poles (diagonal read-off), as in msvdreduce.

    Parameters
    ----------
    ir: impulse response, shape (p, q, L) = (outputs, inputs, lags).
    order: number of poles to recover (must be < min(p*k, q*k), k = floor((L+1)/2)).

    Returns
    -------
    w: poles, shape (order,).
    u: orthonormal null vectors, shape (q, order).
    """
    do, di, num = ir.shape
    k = int(np.floor((num + 1) / 2))
    op = BlockHankelOperator(ir, k)
    U, S, Vh = svds(op, k=order)                  # Lanczos partial SVD, no dense Hankel
    idx = np.argsort(S)[::-1]
    Vh = Vh[idx, :]
    A = Vh[:, di:].dot(np.linalg.pinv(Vh[:, :-di]))   # shift: pseudo-inverse
    A, Q = lschur(A)
    B = Q.T.dot(Vh[:, :di]).astype(complex)
    w = np.diag(A).astype(complex)
    u = np.zeros((di, order), dtype=complex)
    for i in range(order):
        ui = B[0] / np.linalg.norm(B[0], ord=2)
        u[:, i] = ui
        if B.shape[0] > 1:
            Jinv = np.eye(di) - (1.0 + 1.0 / np.conj(w[i])) * np.outer(ui, np.conj(ui))
            B = B[1:, :].dot(Jinv)
        else:
            B = B[1:, :]
    return w, u


def msvdreduce_complex(ir: np.ndarray, order: int, fit_C: bool = True, num_lags: int = None):
    """
    Hankel-SVD reduction to a real TIB realization, for REAL **or COMPLEX** poles.

    Companion to msvdreduce_fixed that removes its real-pole-only limitation. Rather than
    reading poles off the real-Schur diagonal and peeling scalar null vectors -- which is
    wrong for complex-conjugate poles, because the real Schur form carries them in 2x2
    blocks, not as scalar diagonal entries -- this forms the dense shift realization
    (A, B) from the truncated SVD (eig(A) = the poles, real or complex) and re-coordinates
    it into real TIB form via tib_from_state_space, which represents complex-conjugate
    pairs as 2x2 real diagonal blocks. Optionally least-squares fits C so the realization
    reproduces the impulse response.

    Unlike msvdreduce / msvdreduce_fixed (which return (poles, null_vectors)), this returns
    a *realization*: complex poles live inside 2x2 real blocks, not as scalar
    (pole, null vector) pairs. Read the poles via ``la.eigvals(sys.A(dense=True))``
    (complex-aware) -- NOT np.diag, which would give the real parts of the 2x2 blocks.
    (NB: TIBStateSpace.poles() currently fails on a dense-backed A -- it runs la.eig on
    the sparse form; use la.eigvals(sys.A(dense=True)) instead.)

    Parameters
    ----------
    ir : np.ndarray, shape (p, q, L) = (outputs, inputs, lags).
    order : reduced order; needs order <= min(p, q) * floor((L+1)/2), and L >~ 2*order.
    fit_C : if True, also least-squares fit C (p x order) so C A^{k-1} B matches ir.
    num_lags : lags used for the C fit (default: all L).

    Returns
    -------
    sys : TIBStateSpace
        Real, (block-)lower-triangular A_tib with A_tib A_tib^T + B_tib B_tib^T = I.
        Real poles are 1x1 diagonal entries; complex-conjugate pairs are 2x2 real diagonal
        blocks. Poles: la.eigvals(sys.A(dense=True)).
    C : np.ndarray (p, order) or None
        Fitted output matrix (only if fit_C).

    Notes
    -----
    tib_from_state_space solves a discrete Lyapunov equation, so the shift A must be stable;
    on clean IR it is, but a noisy IR can push the truncated shift's eigenvalues onto or
    outside the unit circle.
    """
    from ga.reducers.balanced_truncation import tib_from_state_space
    from ga.filters.tib import krylov_basis

    do, di, num = ir.shape
    k = int(np.floor((num + 1) / 2))
    max_order = min(do, di) * k
    if not (0 < order <= max_order):
        raise ValueError(
            f"order={order} exceeds the recoverable max {max_order} for L={num} "
            f"(p={do}, q={di}); need a longer IR -- roughly L >= 2*order."
        )

    H = blkhankel(ir)
    U, S, Vh = la.svd(H)
    Vh = Vh[:order, :]
    A = Vh[:, di:].dot(np.linalg.pinv(Vh[:, :-di]))   # dense shift; eig(A) = poles (real/complex)
    B = Vh[:, :di]
    sys, _ = tib_from_state_space(A, B)               # complex pairs -> 2x2 real diagonal blocks

    C = None
    if fit_C:
        L = num if num_lags is None else num_lags
        At, Bt = sys.A(dense=True), sys.B(dense=True)
        K = np.asarray(krylov_basis(At, Bt, L)).reshape(At.shape[0], -1)   # (order, q*L)
        C = np.real(ir[:, :, :L].reshape(do, -1) @ np.linalg.pinv(K))      # (p, order)

    return sys, C


