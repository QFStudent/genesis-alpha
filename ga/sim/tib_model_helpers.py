"""TIB MIMO matrix helpers used by ga/sim model and solver modules."""

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp

from ga.filters.tib import TIBSystem


def mimo_system_germs(U: np.ndarray):
    p, q = U.shape
    M = np.eye(p)
    if p > q:
        _, R = np.linalg.qr(U.T)
        M[q:, :q] = -la.solve_triangular(R[:, :q], R[:, q:]).T
        _, M[q:, :] = la.lu(M[q:, :], permute_l=True)

    Omega = np.tril(U.dot(U.T))
    N = M.dot(Omega)
    return map(sp.csc_matrix, (M, N))


def mimo_system_matrices_biased(lambd: np.ndarray, U: np.ndarray,
                                normalize=True):
    assert U.shape[0] == lambd.size, "Incompatible sizes"
    if normalize:
        mag = np.linalg.norm(U, axis=1, keepdims=True)
        U = U / mag

    Mu, Nu = mimo_system_germs(U)
    dw = sp.diags(lambd)
    dcw = sp.diags(1.0 / np.sqrt(1 - lambd * np.conj(lambd)))

    Du = Nu - Mu
    M = (Nu + Du * dw.conj()) * dcw
    N = (Nu * dw + Du) * dcw

    p, q = U.shape
    U = U.copy()
    U[q:, :] = 0
    Q = sp.eye(p)
    return TIBSystem(M=M, N=N, Q=Q, U=U)
