import numpy as np 
from numpy.typing import NDArray
import scipy.linalg as la
import scipy.sparse as sp
from dataclasses import dataclass


@dataclass
class TIBSystem:
    """

    Q is a rotation matrix and is used to keep matrices real
    (even when poles are complex). rows of U are null vectors
    and U[:q, :] is unitary.
    """
    M: np.ndarray
    N: np.ndarray
    Q: np.ndarray
    U: np.ndarray
    C: np.ndarray = None

    def __str__(self):
        return "TIBSystem({}, {})".format(self.p, self.q)

    @property
    def p(self):
        """number of poles"""
        return self.M.shape[0]

    @property
    def q(self):
        """Output/innovation size"""
        return self.U.shape[1]

    @staticmethod
    def _ensure_sparse(A):
        return A if sparse.issparse(A) else sparse.csc_matrix(A)

    @staticmethod
    def _ensure_dense(A):
        return A.toarray() if sparse.issparse(A) else A

    def A(self, dense=False):
        A = sparse.linalg.spsolve(self.M, self.N).dot(self.Q)
        check = self._ensure_dense if dense else self._ensure_sparse
        return check(A)

    def B(self, dense=False):
        B = sparse.linalg.spsolve(self.M, self.U)
        B = B.reshape((self.p, self.q))
        check = self._ensure_dense if dense else self._ensure_sparse
        return check(B)

    def poles(self):
        A = self.A()

        if self.q == 1:
            ev, _ = np.linalg.eig(A.todense())
        else:
            ev = A.diagonal(0)

        return ev


@dataclass(frozen=True)
class TIBStateSpace: 
    A_: np.ndarray
    B_: np.ndarray

    def __str__(self): 
        return "TIBStateSpace({}, {})".format(self.p, self.q)
    
    @property
    def p(self): 
        """number of poles"""
        return self.A_.shape[0]
    
    @property
    def q(self): 
        """Output/innovation size"""
        return self.B_.shape[1]
    
    @staticmethod
    def _ensure_sparse(A): 
        return A if sp.issparse(A) else sp.csr_matrix(A)
    
    @staticmethod
    def _ensure_dense(A): 
        return A.toarray() if sp.issparse(A) else A
    
    def A(self, dense=False): 
        A = self.A_ 
        check = self._ensure_dense if dense else self._ensure_sparse
        return check(A)
    
    def B(self, dense=False): 
        B = self.B_ 
        check = self._ensure_dense if dense else self._ensure_sparse
        return check(B)
    
    def poles(self):
        # la.eig needs a dense array (self.A() returns sparse) and always returns complex;
        # real_if_close drops negligible imaginary parts (< ~1e-8) so a real system gives
        # real poles, while genuine complex-conjugate pairs stay complex.
        ev = la.eig(self.A(dense=True))[0]
        return np.real_if_close(ev, tol=1e8)


def poles_chebyshev_roots(num_poles: int, a: float = 0.0, b: float = 1.0): 
    """
    Compute the Chebyshev roots of the polynomial
    """
    poles = np.cos(np.pi * np.arange(1, num_poles * 2, 2) / (num_poles * 2))
    poles *= 0.5 * (b - a)
    poles += 0.5 * (b + a)
    return poles[::-1]


def siso_system_matrices_real(lambd: np.ndarray):
    """
    Transform prescribed poles to a _real_ banded matrices.
    The poles can have non-zero imaginary part.

    It is expected that non-real poles are in complex conjugate
    pair. To avoid having complex number on the diagonal, this
    returns a rotation matrix Q such all matrix entries are real.
    """
    rho = lambd[np.isreal(lambd)]
    m = len(rho)
    cc = lambd[np.imag(lambd) > 0]
    acc = np.abs(cc)
    psi = np.arccos(2. * acc / (1. + acc**2)) - np.arccos(2. * np.real(cc)/
        (1. + acc**2))
    cc = np.array([cc, np.conj(cc)]).conj().T.flatten()
    mu = np.concatenate([rho, np.abs(cc)])
    mu = np.real(mu)

    if len(mu) != len(lambd):
        raise ValueError('Complex eigenvalues should be conjugate pairs.')

    M, N = siso_system_matrices(mu)
    Q = [np.eye(m)]
    for k in range(len(psi)):
        G = np.array([[np.cos(psi[k]), -np.sin(psi[k])],
                      [np.sin(psi[k]), np.cos(psi[k])]])
        Q.append(G)

    Q = la.block_diag(*Q)
    Q = sparse.csc_matrix(Q)
    U = np.eye(M.shape[0], 1)
    return TIBSystem(M=M, N=N, Q=Q, U=U)


def siso_system_matrices(lambd: np.ndarray):
    """
    Compute bidiagonal matrix representation of scalar TIB
    system with eigenvalues lambda. Repeated eigenvalues are
    accepted.
    """
    c = 1/np.sqrt(1 - np.abs(lambd)**2)
    s = lambd*c
    n = len(lambd)
    M = sparse.spdiags(1.0 / c, 0, n, n) * sparse.spdiags(
        [s.conj(), c], np.array([-1, 0]), n, n)
    N = sparse.spdiags(1.0 / c, 0, n, n) * sparse.spdiags(
        [c, s], np.array([-1, 0]), n, n)
    M = M / np.sqrt(1 - np.abs(lambd[0])**2)
    N = N / np.sqrt(1 - np.abs(lambd[0])**2)
    # Apply signature to conform to Meixner function definition
    S = sparse.spdiags((-1)**np.arange(n), 0, n, n)
    M = M * S          # sparse matrix overload matrix multiplication using *
    N = N * S
    # spsolve is more efficient when sparse matrix is
    # in the CSC matrix format
    M = sparse.csc_matrix(M)
    N = sparse.csc_matrix(N)
    return M, N


def mimo_standard_null_vectors(p: int, q: int) -> np.ndarray: 
    """
    Compute the standard null vectors of the MIMO system
    p: number of poles
    q: number of inputs 
    """
    c = np.array(np.arange(p) % q == 0, dtype=float)
    r = np.zeros(q)
    return la.toeplitz(c, r)


def blaschke_factor(w: complex, z: complex) -> complex: 
    """
    Blaschke factor B_{w}(z)
    w: scalar pole (|w|<1)
    z: complex evaluation point 
    """
    denom = (1.0 - np.conj(w) * z)
    if abs(denom) < 1e-15: 
        # avoid blow-up at pole in evaluation 
        denom = denom + 1e-15 
        print(f"Warning: Blaschke factor is numerically unstable at pole {w}")
    return (z - w) / denom


def blaschke_potapov_factor(w: complex, z: complex, u: np.ndarray) -> np.ndarray: 
    f"""
    Blaschke-Potapov factor β_{w,u}(z)
    w: scalar pole (|w|<1)
    z: complex evaluation point 
    u: unit vector in C^m (shape m×1)
    
    returns m×m vector
    """
    u = u.reshape(-1, 1)
    if abs(np.linalg.norm(u) - 1.0) > 1e-9:
        raise ValueError("u is not a unit vector")
    b = blaschke_factor(w, z) - 1
    return np.eye(u.shape[0], dtype=complex) + b * (u @ u.conj().T) 


def normalize(x: np.ndarray, eps: float = 1e-15) -> np.ndarray: 
    """
    Normalize a vector. 
    """
    n = np.linalg.norm(x)
    if n < eps: 
        raise ValueError("Cannot normalize: vector norm is ~0.")
    return x / n 


def _real_if_close_ndarray(A: np.ndarray, tol: float = 1e-12) -> np.ndarray: 
    """
    If A is complex but has negligible imaginary parts, drop them. 

    Some realizations are computed in complex algebra for generality even when 
    the result is mathematically real. Downstream code (e.g. real-only RLS / 
    Cholesky updates) may require float dtypes. 
    """
    A = np.asarray(A)
    if np.iscomplexobj(A): 
        imag_max = np.max(np.abs(np.imag(A))) if A.size else 0.0 
        if imag_max <= tol: 
            return np.real(A)
    return A 


def null_basis_realization(lambd: NDArray[complex],
                           v: NDArray[complex],
                           return_all: bool = False) -> "TIBStateSpace":
    """
    Null basis realization of a MIMO system. 
    
    Parameters: 
        lambd: poles of the system. 
        v: null basis of the system. shape (m, n) where n is the number of poles and m is the number of inputs
    
    Returns: 
        A: system matrix. 
        B: input matrix. 
    """
    lambd = np.asarray(lambd, dtype=complex).reshape(-1)
    n = len(lambd)

    v = np.asarray(v, dtype=complex)
    if v.ndim != 2 or v.shape[1] != n: 
        raise ValueError(f"v must have shape (m, n). Got {v.shape}, n={n}.")
    m = v.shape[0]

    # Sanity checks 
    if np.any(np.abs(lambd) >= 1.0 + 1e-12):
        raise ValueError("Poles must be inside the unit circle.")
    for k in range(n):
        if np.linalg.norm(v[:, k]) < 1e-12: 
            raise ValueError(f"Null vector {k} is zero.")
    
    # Storage 
    ys = np.zeros((m, n, 1), dtype=complex)
    Js = [] 

    # Initialize k = 1 
    y1 = normalize(v[:, 0].reshape(-1, 1))
    ys[:, 0] = y1 
    t = np.sqrt(max(0.0, 1.0 - np.abs(lambd[0])**2))
    A = np.array([[lambd[0]]], dtype=complex)         # (1, 1)
    B = t * y1.conj().reshape(1, -1)         # B_1 = t_1 y_1* (y_1 normalized)
    L = np.zeros((m, 0), dtype=complex)          
    P = np.eye(m, dtype=complex)

    # Define J_1 (needed starting at k = 2)
    J1 = np.eye(m, dtype=complex) - (1.0 + np.conj(lambd[0])) * y1 @ y1.conj().T
    Js.append(J1)

    A_list = [A.copy()]
    B_list = [B.copy()]

    # Iterate k = 2, ..., n
    for k in range(1, n): 
        L = np.hstack([Js[k-1] @ L, t * ys[:, k-1]])
        t = np.sqrt(max(0.0, 1.0 - np.abs(lambd[k])**2))

        # for every pole lambd[k], we need to build the matrix M_k = Π_{i=1}^{k-1} β_{λ_i, y_i}{λ_k^*}
        z = np.conj(lambd[k])
        M = np.eye(m, dtype=complex)
        for i in range(0, k):
            beta = blaschke_potapov_factor(lambd[i], z, ys[:, i])   # β_{λ_i, y_i}(λ_k*): pole=λ_i, eval=z=λ_k*
            M = M @ beta.conj().T
        ys[:, k] = normalize(M @ v[:, k].reshape(-1, 1))
        Jk = np.eye(m, dtype=complex) - (1.0 + np.conj(lambd[k])) * ys[:, k] @ ys[:, k].conj().T
        Js.append(Jk)

        # Update product P_k = J_{k-1} ... J_1 (for this k, we need P = J_{k-1}...J_1)
        P = Js[k-1] @ P 

        # sanity check 
        tmp = ys[:, k].conj().T @ (L @ L.conj().T + P @ P.conj().T) @ ys[:, k]
        assert np.abs(np.linalg.norm(tmp) - 1) < 1e-10, "Sanity check failed."

        # Build A_k and B_k 
        # Bottom-left row: t_k y_k^* L_k (1 x (k-1))
        a_row = (t * np.conj(ys[:, k]).reshape(1, m)) @ L    # (1, k-1)
        # New B row: t_k y_k^* P (1 x m)
        b_row = (t * np.conj(ys[:, k]).reshape(1, m)) @ P    # (1, m)

        # Expand A and B 
        A_new = np.zeros((k+1, k+1), dtype=complex)
        A_new[:-1, :-1] = A 
        A_new[-1, :-1] = a_row 
        A_new[-1, -1] = lambd[k]
        A = A_new 

        B = np.vstack([B, b_row])

        A_list.append(A.copy())
        B_list.append(B)

    A = _real_if_close_ndarray(A)
    B = _real_if_close_ndarray(B)
    return TIBStateSpace(A_=A, B_=B)


def build_tib_from_directions(lambd: NDArray[complex],
                              y: NDArray[complex]) -> "TIBStateSpace":
    """
    Build a TIB realization from poles and ALREADY-ORTHONORMAL input directions.

    Companion to null_basis_realization for the *recovery* direction.
    null_basis_realization takes raw null vectors v and Blaschke-orthonormalizes them
    into y; this function takes the orthonormal y directly and runs the same forward
    recursion *without* the orthogonalization step. Use it to reconstruct a TIB pair
    from directions recovered by ga.reducers.hankel.msvdreduce_fixed (whose output is
    already the orthonormalized y_k).

    WARNING: do NOT pass recovered y_k to null_basis_realization -- it treats its
    input as raw v and re-applies the Blaschke deflation to already-deflated vectors
    (double-deflation), realizing the wrong transfer function. Pass them here instead.

    Parameters:
        lambd: poles, shape (n,).
        y: orthonormal input directions, shape (m, n) -- columns are unit vectors.

    Returns:
        TIBStateSpace with A (n x n, block-lower-triangular) and B (n x m).
    """
    lambd = np.asarray(lambd, dtype=complex).reshape(-1)
    n = len(lambd)

    y = np.asarray(y, dtype=complex)
    if y.ndim != 2 or y.shape[1] != n:
        raise ValueError(f"y must have shape (m, n). Got {y.shape}, n={n}.")
    m = y.shape[0]

    if np.any(np.abs(lambd) >= 1.0 + 1e-12):
        raise ValueError("Poles must be inside the unit circle.")

    y0 = y[:, 0].reshape(-1, 1)
    t = np.sqrt(max(0.0, 1.0 - np.abs(lambd[0])**2))
    A = np.array([[lambd[0]]], dtype=complex)
    B = t * y0.conj().reshape(1, -1)
    L = np.zeros((m, 0), dtype=complex)
    P = np.eye(m, dtype=complex)
    Js = [np.eye(m, dtype=complex) - (1.0 + np.conj(lambd[0])) * y0 @ y0.conj().T]

    for k in range(1, n):
        L = np.hstack([Js[k-1] @ L, t * y[:, k-1].reshape(-1, 1)])
        t = np.sqrt(max(0.0, 1.0 - np.abs(lambd[k])**2))
        yk = y[:, k].reshape(-1, 1)
        Js.append(np.eye(m, dtype=complex) - (1.0 + np.conj(lambd[k])) * yk @ yk.conj().T)
        P = Js[k-1] @ P

        a_row = (t * yk.conj().reshape(1, m)) @ L
        b_row = (t * yk.conj().reshape(1, m)) @ P

        A_new = np.zeros((k+1, k+1), dtype=complex)
        A_new[:-1, :-1] = A
        A_new[-1, :-1] = a_row
        A_new[-1, -1] = lambd[k]
        A = A_new
        B = np.vstack([B, b_row])

    A = _real_if_close_ndarray(A)
    B = _real_if_close_ndarray(B)
    return TIBStateSpace(A_=A, B_=B)


def null_basis_realization_real(lambd: NDArray[complex],
                                v: NDArray[complex]) -> "TIBStateSpace":
    """
    Real TIB realization for poles that may include complex-conjugate pairs.

    Companion to null_basis_realization that keeps (A, B) -- and hence the impulse
    response and any fitted C -- real when some poles are complex. The construction:

      1. build the realization on the pole *magnitudes* (signed value for real poles)
         via null_basis_realization. Magnitudes are real, so (A_mag, B_mag) is real,
         input-balanced (A A* + B B* = I) and lower-triangular, with each conjugate
         pair sitting as a 2x2 diagonal block [[r, 0], [x, r]];
      2. inject each pair's phase with a 2x2 rotation Q. Because Q is orthogonal the
         input-balance is preserved, and A = A_mag @ Q has the pair as eigenvalues.

    The rotation angle is read off the *realized* block (not a pole-only formula),
    since null_basis_realization's deflation couples each pair to the poles before it:
    for r = A_mag[j, j], x = A_mag[j+1, j], theta = arg(lambda),

        psi = -atan2(x, 2r) + arccos( 2 r cos(theta) / sqrt(4 r^2 + x^2) )

    makes the 2x2 block A_mag[j:j+2, j:j+2] @ G(psi) have eigenvalues r e^{+-i theta}.

    Parameters:
        lambd: poles. Non-real poles must occur in complex-conjugate pairs.
        v: null vectors, shape (q, n) -- one column per pole; the real part is used.

    Returns:
        TIBStateSpace with real, block-lower-triangular A and real B; eig(A) = lambd
        (the state ordering is real poles first, then each conjugate pair).
    """
    lambd = np.asarray(lambd).reshape(-1)
    v = np.asarray(v)
    n = lambd.shape[0]
    if v.ndim != 2 or v.shape[1] != n:
        raise ValueError(f"v must have shape (q, n) with n=len(lambd). Got {v.shape}, n={n}.")

    real_mask = np.abs(np.imag(lambd)) < 1e-12
    pos_idx = np.where(np.imag(lambd) > 1e-12)[0]
    neg_idx = list(np.where(np.imag(lambd) < -1e-12)[0])
    if len(pos_idx) != len(neg_idx):
        raise ValueError("Non-real poles must occur in complex-conjugate pairs.")

    # processing order: real poles first, then each conjugate pair adjacent
    order = list(np.where(real_mask)[0])
    pairs = []                         # (slot of the pair's first state, angle theta)
    used = set()
    for i in pos_idx:
        j = next((k for k in neg_idx
                  if k not in used and abs(lambd[k] - np.conj(lambd[i])) < 1e-9), None)
        if j is None:
            raise ValueError(f"No conjugate partner for pole {lambd[i]}.")
        used.add(j)
        pairs.append((len(order), float(np.angle(lambd[i]))))
        order += [i, j]
    order = np.asarray(order, dtype=int)

    # radial poles: signed value for real poles, magnitude (repeated) for pairs
    radial = np.where(real_mask[order], np.real(lambd[order]), np.abs(lambd[order])).astype(float)
    v_ord = np.real(v[:, order])

    sys_mag = null_basis_realization(radial, v_ord)
    A_mag = np.real(sys_mag.A(dense=True))
    B_mag = np.real(sys_mag.B(dense=True))

    n_real = n - 2 * len(pairs)
    blocks = [np.eye(n_real)] if n_real else []
    for j, theta in pairs:
        r, x = A_mag[j, j], A_mag[j + 1, j]
        rad = np.hypot(2.0 * r, x)
        psi = -np.arctan2(x, 2.0 * r) + np.arccos(np.clip(2.0 * r * np.cos(theta) / rad, -1.0, 1.0))
        c, s = np.cos(psi), np.sin(psi)
        blocks.append(np.array([[c, -s], [s, c]]))
    Q = la.block_diag(*blocks) if blocks else np.eye(n)

    return TIBStateSpace(A_=np.real(A_mag @ Q), B_=B_mag)


def build_tib_from_directions_real(lambd: NDArray[complex],
                                   y: NDArray[complex]) -> "TIBStateSpace":
    """
    Real TIB realization from poles + ALREADY-ORTHONORMAL directions (may be complex pairs).

    Real counterpart of build_tib_from_directions, exactly as null_basis_realization_real
    is to null_basis_realization: build on the pole *magnitudes* via
    build_tib_from_directions (no re-deflation -- the directions are used as-is), then
    inject each conjugate pair's phase with a 2x2 rotation Q. So A = A_mag @ Q is real
    with the pair as eigenvalues, B = B_mag is real, and the input-balance is preserved
    (Q orthogonal). The rotation angle is read off the realized block; see
    null_basis_realization_real for the derivation.

    Use this -- not build_tib_from_directions -- when the poles include complex pairs;
    build_tib_from_directions would place them on a complex diagonal and go complex.

    Parameters:
        lambd: poles. Non-real poles must occur in complex-conjugate pairs.
        y: orthonormal directions, shape (q, n) -- one column per pole; real part is used.

    Returns:
        TIBStateSpace with real, block-lower-triangular A and real B; eig(A) = lambd
        (the state ordering is real poles first, then each conjugate pair).
    """
    lambd = np.asarray(lambd).reshape(-1)
    y = np.asarray(y)
    n = lambd.shape[0]
    if y.ndim != 2 or y.shape[1] != n:
        raise ValueError(f"y must have shape (q, n) with n=len(lambd). Got {y.shape}, n={n}.")

    real_mask = np.abs(np.imag(lambd)) < 1e-12
    pos_idx = np.where(np.imag(lambd) > 1e-12)[0]
    neg_idx = list(np.where(np.imag(lambd) < -1e-12)[0])
    if len(pos_idx) != len(neg_idx):
        raise ValueError("Non-real poles must occur in complex-conjugate pairs.")

    # processing order: real poles first, then each conjugate pair adjacent
    order = list(np.where(real_mask)[0])
    pairs = []                         # (slot of the pair's first state, angle theta)
    used = set()
    for i in pos_idx:
        j = next((k for k in neg_idx
                  if k not in used and abs(lambd[k] - np.conj(lambd[i])) < 1e-9), None)
        if j is None:
            raise ValueError(f"No conjugate partner for pole {lambd[i]}.")
        used.add(j)
        pairs.append((len(order), float(np.angle(lambd[i]))))
        order += [i, j]
    order = np.asarray(order, dtype=int)

    # radial poles: signed value for real poles, magnitude (repeated) for pairs
    radial = np.where(real_mask[order], np.real(lambd[order]), np.abs(lambd[order])).astype(float)
    y_ord = np.real(y[:, order])

    sys_mag = build_tib_from_directions(radial, y_ord)   # directions used as-is (no re-deflation)
    A_mag = np.real(sys_mag.A(dense=True))
    B_mag = np.real(sys_mag.B(dense=True))

    n_real = n - 2 * len(pairs)
    blocks = [np.eye(n_real)] if n_real else []
    for j, theta in pairs:
        r, x = A_mag[j, j], A_mag[j + 1, j]
        rad = np.hypot(2.0 * r, x)
        psi = -np.arctan2(x, 2.0 * r) + np.arccos(np.clip(2.0 * r * np.cos(theta) / rad, -1.0, 1.0))
        c, s = np.cos(psi), np.sin(psi)
        blocks.append(np.array([[c, -s], [s, c]]))
    Q = la.block_diag(*blocks) if blocks else np.eye(n)

    return TIBStateSpace(A_=np.real(A_mag @ Q), B_=B_mag)


def krylov_basis(A: np.ndarray, B: np.ndarray, n: int) -> NDArray[complex]:
    """
    Use the naive method to generate the Krylov basis. 
    """
    if hasattr(A, 'toarray'): 
        A = A.toarray() 
    else: 
        A = np.asarray(A)
    if hasattr(B, 'toarray'): 
        B = B.toarray() 
    else: 
        B = np.asarray(B)

    if B.ndim < 2: 
        B = B.reshape((-1, 1))
    
    nCol = int(2**np.ceil(np.log2(n)))
    p, q = B.shape 
    K = np.empty((p, q * nCol), dtype=A.dtype)
    K[:, :q] = B 

    width = q 
    while width < K.shape[1]: 
        newCols = np.arange(width, 2 * width)
        K[:, newCols] = A @ K[:, :width]
        width = 2 * width 
        A = A @ A 
    
    K = K[:, :q * n]
    if q > 1: 
        K = K.reshape((p, q, n), order='F')
    
    return K 


def mimo_poles_to_ir(lambd: NDArray[complex], v: NDArray[complex], C: NDArray[float], m: int) -> NDArray[float]:
    """
    Convert poles to impulse response.

    Parameters:
        lambd: poles of the system. shape (n,) where n is the number of poles
        v: null basis of the system. shape (n_inputs, n) where n is the number of poles
        C: output matrix. shape (p, n) where p is the number of outputs and n is the number of poles
        m: number of lags (length of the impulse response to generate).
           NOTE: this is *not* the number of inputs despite the name.

    Returns:
        ir: impulse response. shape (p, n_inputs, m) where p is the number of
            outputs, n_inputs the number of inputs, and m the number of lags.
    """
    n = len(lambd)
    tib_sys = null_basis_realization(lambd, v)
    K = krylov_basis(tib_sys.A(), tib_sys.B(), m)
    ir = np.tensordot(C, K, (-1, 0))
    return ir 

