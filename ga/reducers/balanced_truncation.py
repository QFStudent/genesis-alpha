"""
Balanced Truncation for Model Order Reduction

Implements the square-root (SR) method for balanced truncation, which is
numerically more robust than directly computing Grammians.

References:
- Yu Mu, "Notes for Model Reduction", Jan 2026
- Laub, Heath, Paige, Ward, "Computation of System Balancing Transformations"

The algorithm:
1. Solve Lyapunov equations for controllability (P) and observability (Q) Grammians
2. SVD of Hankel factors: P = Z_P.T @ Z_P, Q = Z_Q.T @ Z_Q
3. Construct transformation matrices W, V for balancing
   bridge matrix: Z_P.T @ Z_Q = U @ Σ @ V.T, where Z_P.T is upper triangular
   A_r = W.T @ A @ V, B_r = W.T @ B, C_r = C @ V
4. Compute reduced system: A_r = W.T @ A @ V, B_r = W.T @ B, C_r = C @ V
5. Compute reduced system: A_r = W.T @ A @ V, B_r = W.T @ B, C_r = C @ V
"""

import numpy as np
import scipy.linalg as la
from dataclasses import dataclass
from typing import Tuple, Optional

from ga.filters.tib import TIBStateSpace


@dataclass
class BalancedReductionResult:
    """Result container for balanced truncation."""
    A_r: np.ndarray   # Reduced state matrix (r x r)
    B_r: np.ndarray   # Reduced input matrix (r x m)
    C_r: np.ndarray   # Reduced output matrix (p x r)
    D: np.ndarray     # Feedthrough matrix (unchanged, p x m)
    hsv: np.ndarray   # Hankel singular values
    W: np.ndarray     # Left transformation matrix
    V: np.ndarray     # Right transformation matrix
    info: dict        # Diagnostic information

    @property
    def order(self) -> int:
        """Reduced system order."""
        return self.A_r.shape[0]

    def error_bound(self) -> float:
        """Upper bound on H-infinity error: ||G - G_r||_He <= 2 * sum(neglected HSV)"""
        return 2.0 * np.sum(self.hsv[self.order:])


def solve_lyapunov_discrete(A: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """
    Solve discrete-time Lyapunov equation: A @ X @ A.T - X + Q = 0
    Equivalently: X = A @ X @ A.T + Q

    Parameters:
    A: State matrix (n x n)
    Q: Right-hand side matrix (n x n), typically B @ B.T or C.T @ C

    Returns:
    X: Solution matrix (n x n)
    """
    return la.solve_discrete_lyapunov(A, Q)


def compute_grammians(A: np.ndarray, B: np.ndarray, C: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute controllability and observability Grammians.
    Controllability Grammian P satisfies: P = A @ P @ A.T + B @ B.T
    Observability Grammian Q satisfies: Q = A.T @ Q @ A + C.T @ C

    Parameters:
    A: State matrix (n x n)
    B: Input matrix (n x m)
    C: Output matrix (p x n)

    Returns:
    P: Controllability Grammian (n x n)
    Q: Observability Grammian (n x n)
    """
    P = solve_lyapunov_discrete(A, B @ B.T)
    Q = solve_lyapunov_discrete(A.T, C.T @ C)
    return P, Q


def cholesky_with_fallback(M: np.ndarray, name: str = "matrix") -> np.ndarray:
    """
    Compute Cholesky factor with fallback for near-singular matrices.

    Parameters:
    M: Positive semi-definite matrix
    name: Name for error messages

    Returns:
    Z: Lower Cholesky factor such that M = Z @ Z.T
    """
    # Ensure symmetry
    M = (M + M.T) / 2

    # Check for negative eigenvalues (numerical issues)
    eigvals = la.eigvalsh(M)
    min_eig = np.min(eigvals)

    if min_eig < -1e-10:
        raise ValueError(f"{name} has negative eigenvalue: {min_eig}")

    # Add small regularization if needed
    if min_eig < 1e-12:
        reg = max(1e-12, -min_eig + 1e-12)
        M = M + reg * np.eye(M.shape[0])

    try:
        Z = la.cholesky(M, lower=True)
    except la.LinAlgError:
        # Fallback: use eigendecomposition
        eigvals, eigvecs = la.eigh(M)
        eigvals = np.maximum(eigvals, 1e-12)
        Z = eigvecs @ np.diag(np.sqrt(eigvals))

    return Z


def balanced_truncation(
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    D: Optional[np.ndarray] = None,
    order: Optional[int] = None,
    tol: float = 1e-12,
    return_transformation: bool = True,
) -> BalancedReductionResult:
    """
    Balanced truncation using the square-root method.

    Given a stable discrete-time system (A, B, C, D), compute a reduced-order
    system (A_r, B_r, C_r, D) that minimizes the Hankel norm approximation error.

    Algorithm (Square-Root Method):
    1. Solve Lyapunov equations for P and Q
    2. Cholesky: P = Z_P.T @ Z_P, Q = Z_Q.T @ Z_Q
    3. SVD: Z_P.T @ Z_Q = U @ Σ @ V.T
    4. W = Z_P @ U_1 @ Σ_1^{-1/2}, V = Z_Q @ V_1 @ Σ_1^{-1/2}
    5. A_r = W.T @ A @ V, B_r = W.T @ B, C_r = C @ V
    6. A_r = W.T @ A @ V, B_r = W.T @ B  # V from Z_Q @ V from Σ^{-1/2}

    Parameters:
    A: State matrix (n x n), must be stable (all eigenvalues inside unit circle)
    B: Input matrix (n x m)
    C: Output matrix (p x n)
    D: Feedthrough matrix (p x m), optional (passed through unchanged)
    order: Reduced order r. If None, determined by Hankel singular values
    tol: Tolerance for determining order r
    return_transformation: If True, include W and V in result

    Raises:
    ValueError: If system is unstable or Grammians are ill-conditioned

    Returns:
    BalancedReductionResult containing reduced system and diagnostics
    """
    n = A.shape[0]
    m = B.shape[1]
    p = C.shape[0]

    if D is None:
        D = np.zeros((p, m))

    info = {}

    # Check stability
    eigvals_A = la.eigvals(A)
    max_eig = np.max(np.abs(eigvals_A))
    info['max_eigenvalue'] = max_eig

    if max_eig >= 1.0 - 1e-10:
        raise ValueError(f"System is unstable: max |eigenvalue| = {max_eig}")

    # Step 1: Compute Grammians
    P, Q = compute_grammians(A, B, C)
    info['cond_P'] = np.linalg.cond(P)
    info['cond_Q'] = np.linalg.cond(Q)

    # Step 2: Cholesky factorization
    Z_P = cholesky_with_fallback(P, "Controllability Grammian P")
    Z_Q = cholesky_with_fallback(Q, "Observability Grammian Q")

    # Step 3: SVD of bridge matrix
    # Z_P.T is upper triangular, Z_P.T @ Z_Q is the bridge matrix
    bridge = Z_P.T @ Z_Q
    U, Sigma, Vh = la.svd(bridge, full_matrices=True)
    hsv = Sigma
    info['hankel_singular_values'] = hsv

    # Determine order if not specified
    if order is None:
        order = np.sum(hsv > tol * hsv[0])
        order = max(1, order)  # At least order 1

    # Step 4: Partition and construct transformation matrices
    info['original_order'] = n
    info['reduced_order'] = order

    U_1 = U[:, :order]
    V_1 = Vh[:order, :].T

    Sigma_1 = Sigma[:order]
    Sigma_1_inv_sqrt = np.diag(1.0 / np.sqrt(Sigma_1))

    W = Z_P @ U_1 @ Sigma_1_inv_sqrt
    V_trans = Z_Q @ V_1 @ Sigma_1_inv_sqrt

    # Step 5: Compute reduced system
    A_r = W.T @ A @ V_trans
    B_r = W.T @ B
    C_r = C @ V_trans

    # Sanity checks
    info['biorthogonality_error'] = np.linalg.norm(W.T @ V_trans - np.eye(order))
    info['balanced_P_error'] = np.linalg.norm(W.T @ P @ W - np.diag(Sigma_1))
    info['balanced_Q_error'] = np.linalg.norm(V_trans.T @ Q @ V_trans - np.diag(Sigma_1))

    if not return_transformation:
        W = None
        V_trans = None

    return BalancedReductionResult(
        A_r=A_r,
        B_r=B_r,
        C_r=C_r,
        D=D,
        hsv=hsv,
        W=W,
        V=V_trans,
        info=info,
    )


def extract_poles_from_bt(result: BalancedReductionResult) -> np.ndarray:
    """
    Extract poles from balanced truncation result.

    Returns:
    poles: Eigenvalues of A_r (may be complex)
    """
    return la.eigvals(result.A_r)


def extract_poles_and_nullvecs_from_bt(
    result: BalancedReductionResult,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract poles and approximate null vectors from balanced truncation result.

    This converts the BT output to a form compatible with TIB realization.
    Note: The null vectors are approximate since BT doesn't preserve TIB structure.

    Parameters:
    result: BalancedReductionResult from balanced_truncation()

    Returns:
    -> Tuple[np.ndarray, np.ndarray]:
        poles: Eigenvalues of A_r (may be complex)
        null_vectors: Normalized rows of transformed B matrix, shape (r, m)
    """
    A_r = result.A_r
    B_r = result.B_r

    # Schur decomposition: A_r = Q @ T @ Q.T where T is quasi-upper triangular
    # Convert to lower triangular by flipping
    T, Q = la.schur(A_r, output='real')
    L = T[::-1, ::-1]
    Q_flip = Q[:, ::-1]

    # Poles are on the diagonal of L
    poles = np.diag(L)

    # Transform B to Schur coordinates
    B_tilde = Q_flip.T @ B_r
    r = A_r.shape[0]
    m = B_r.shape[1]
    null_vectors = np.zeros((r, m))

    # Extract null vectors by normalizing each row
    for i in range(r):
        row = B_tilde[i, :]
        norm = np.linalg.norm(row)
        if norm > 1e-12:
            null_vectors[i, :] = row / norm
        else:
            # Fallback: use canonical basis
            null_vectors[i, i % m] = 1.0

    return poles, null_vectors


def tib_from_state_space(
    A: np.ndarray,
    B: np.ndarray,
) -> Tuple[TIBStateSpace, np.ndarray]:
    """
    Re-coordinate a stable, real state-space (A, B) into Triangular Input
    Balanced (TIB) form, exactly preserving the transfer function.

    This is the correct counterpart to extract_poles_and_nullvecs_from_bt().
    That function returns the normalized rows of the *Schur* realization's B,
    which is balanced in the Grammian sense but NOT input-balanced, so feeding
    it to null_basis_realization yields the wrong system (see docs TIBForm).
    Here we instead:

      1. input-balance: similarity by the Cholesky factor of the controllability
         Grammian P (P = A P A* + B B*) so the new controllability Grammian is I
         -- i.e. A1 A1* + B1 B1* = I, the TIB balance identity;
      2. apply an orthogonal real-Schur rotation to make A (block-)lower-
         triangular. Orthogonal transforms preserve the Grammian, so balance is
         retained.

    The result is real, exactly TIB, and reproduces the original impulse
    response. Real poles appear as 1x1 entries on the diagonal; complex-
    conjugate pairs appear as 2x2 real diagonal blocks (so A is block-lower-
    triangular / lower-Hessenberg). Recover the poles via la.eigvals or
    TIBStateSpace.poles().

    Parameters:
    A: state matrix (n x n), stable (all eigenvalues inside the unit circle)
    B: input matrix (n x m)

    Returns:
    Tuple[TIBStateSpace, np.ndarray]:
        sys: TIB realization (A_tib, B_tib)
        transform: (n x n); the output matrix maps as C_tib = C @ transform
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)

    # 1. Controllability Grammian P = A P A* + B B*, factor P = Z Z*
    P = solve_lyapunov_discrete(A, B @ B.T)
    Z = cholesky_with_fallback(P, "Controllability Grammian P")

    # 2. Input-balance: similarity by Z drives the controllability Grammian to I
    A1 = la.solve(Z, A @ Z)
    B1 = la.solve(Z, B)

    # 3. Orthogonal real-Schur: upper quasi-triangular U* A1 U (2x2 blocks for
    #    complex-conjugate pairs). Flipping reverses it to (block-)lower form.
    _, U = la.schur(A1, output="real")
    Uf = U[:, ::-1]
    A_tib = Uf.T @ A1 @ Uf
    B_tib = Uf.T @ B1
    transform = Z @ Uf

    return TIBStateSpace(A_=A_tib, B_=B_tib), transform


def bt_impulse_response(
    result: BalancedReductionResult,
    num_lags: int,
) -> np.ndarray:
    """
    Compute impulse response from balanced truncation result.

    IR[k] = C_r @ A_r^k @ B_r

    Parameters:
    result: BalancedReductionResult
    num_lags: Number of lags to compute

    Returns:
    ir: Impulse response, shape (p, m, num_lags)
    """
    A_r = result.A_r
    B_r = result.B_r
    C_r = result.C_r
    p = C_r.shape[0]
    m = B_r.shape[1]

    A_r_power = np.eye(A_r.shape[0])
    ir = np.zeros((p, m, num_lags))

    for k in range(num_lags):
        ir[:, :, k] = C_r @ A_r_power @ B_r
        A_r_power = A_r_power @ A_r

    return ir


def compare_impulse_responses(
    ir1: np.ndarray,
    ir2: np.ndarray,
    normalize: bool = True,
) -> dict:
    """
    Compare two impulse responses.

    Parameters:
    ir1: First impulse response, shape (p, m, num_lags)
    ir2: Second impulse response, shape (p, m, num_lags)
    normalize: If True, normalize before comparison

    Returns:
    Dictionary with comparison metrics
    """
    # Ensure same length
    min_lags = min(ir1.shape[-1], ir2.shape[-1])
    ir1 = ir1[..., :min_lags]
    ir2 = ir2[..., :min_lags]

    if normalize:
        norm1 = np.max(np.abs(ir1))
        norm2 = np.max(np.abs(ir2))
        if norm1 > 1e-12:
            ir1 = ir1 / norm1
        if norm2 > 1e-12:
            ir2 = ir2 / norm2

    diff = ir1 - ir2

    return {
        'max_error': np.max(np.abs(diff)),
        'mean_error': np.mean(np.abs(diff)),
        'relative_error': np.mean(np.abs(diff)) / max(np.linalg.norm(ir1), 1e-12),
        'early_max_error': np.max(np.abs(diff[..., :min(20, min_lags)])),
    }


def sanity_check_bt(
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    result: BalancedReductionResult,
    verbose: bool = True,
) -> dict:
    """
    Run sanity checks on balanced truncation result.

    Checks from the notes:
    1. U_1.T @ U_1 = I
    2. V_1.T @ V_1 = I
    3. ||Z_P.T @ U_1 @ U_1.T @ Z_P - I||
    4. W.T @ V = Σ_1
    5. W.T @ P @ W = Σ_1
    6. Condition numbers

    Parameters:
    A, B, C: Original system matrices
    result: BalancedReductionResult
    verbose: If True, print results

    Returns:
    Dictionary with check results
    """
    checks = {}

    # Get stored info
    checks['biorthogonality_error'] = result.info.get('biorthogonality_error', np.nan)
    checks['balanced_P_error'] = result.info.get('balanced_P_error', np.nan)
    checks['balanced_Q_error'] = result.info.get('balanced_Q_error', np.nan)
    checks['cond_P'] = result.info.get('cond_P', np.nan)
    checks['cond_Q'] = result.info.get('cond_Q', np.nan)

    # Check reduced system stability
    eigvals_Ar = la.eigvals(result.A_r)
    checks['max_eig_Ar'] = np.max(np.abs(eigvals_Ar))
    checks['is_stable'] = checks['max_eig_Ar'] < 1.0

    # Transformation condition numbers
    if result.W is not None:
        checks['cond_W'] = np.linalg.cond(result.W)
        checks['cond_V'] = np.linalg.cond(result.V)

    if verbose:
        print("Balanced Truncation Sanity Checks:")
        print(f"  Biorthogonality error (W.T @ V - I): {checks['biorthogonality_error']:.2e}")
        print(f"  Balanced P error: {checks['balanced_P_error']:.2e}")
        print(f"  Balanced Q error: {checks['balanced_Q_error']:.2e}")
        print(f"  Condition(P): {checks['cond_P']:.2e}")
        print(f"  Condition(Q): {checks['cond_Q']:.2e}")
        print(f"  Max |eigenvalue(A_r)|: {checks['max_eig_Ar']:.6f}")
        print(f"  Reduced system stable: {checks['is_stable']}")

        if result.W is not None:
            print(f"  Condition(W): {checks['cond_W']:.2e}")
            print(f"  Condition(V): {checks['cond_V']:.2e}")

    return checks
