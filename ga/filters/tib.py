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
        return f"TIBStateSpace({}, {})".format(self.p, self.q)
    
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
    if np.linalg.norm(u) < 1 - 1e-12: 
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
                           return_all: bool = False) -> tuple[NDArray[complex], NDArray[complex]]: 
    """
    Null basis realization of a MIMO system. 
    
    Parameters: 
        lambd: poles of the system. 
        v: null basis of the system. shape (m, n) where n is the number of poles and m is the number of inputs
    
    Returns: 
        A: system matrix. 
        B: input matrix. 
    """
    pass 



