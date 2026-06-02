import numpy as np 
from numpy.typing import NDArray
import scipy.linalg as la
import scipy.sparse as sp
from dataclasses import dataclass


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
        A = self.A() 
        ev = la.eig(A)[0]
        return ev 


def poles_chebyshev_roots(num_poles: int, a: float = 0.0, b: float = 1.0): 
    """
    Compute the Chebyshev roots of the polynomial
    """
    poles = np.cos(np.pi * np.arange(1, num_poles * 2, 2) / (num_poles * 2))
    poles *= 0.5 * (b - a)
    poles += 0.5 * (b + a)
    return poles[::-1]


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

