"""
Tests for ga.filters.tib.mimo_ar2ir -- the multivariate AR-coefficients -> impulse-response
map (forward IR of the AR / whitening inverse A(z)^{-1}, A(z) = I - sum_k A_k z^{-k}).

The defining correctness property is A(z) H(z) = I (the IR inverts the AR polynomial); since
the inverse of a square invertible transfer function is unique, a machine-precision residual is
conclusive. The remaining tests cross-check against independent methods (matrix power, a dense
block-Toeplitz solve, a round-trip ir->ar, and the left/right and scalar recursions).
See docs/derivations/07-mimo-ar-to-ir.md.
"""

import numpy as np
import pytest

from ga.filters.tib import mimo_ar2ir


def _conv_resid(A, h, p, d):
    """max_n || sum_k cal_A_k h_{n-k} - delta_{n0} I ||  with cal_A_0=I, cal_A_k=-A_k."""
    cal = [np.eye(d)] + [-A[k] for k in range(p)]
    return max(np.max(np.abs(sum(cal[k] @ h[n - k] for k in range(0, min(n, p) + 1))
                             - (np.eye(d) if n == 0 else 0.0))) for n in range(len(h)))


@pytest.mark.parametrize("seed", range(20))
def test_inverse_property(seed):
    # A(z) H(z) = I to machine precision, over random shapes incl. p >= m (the correctness test)
    rng = np.random.default_rng(seed)
    d = int(rng.integers(1, 6)); p = int(rng.integers(1, 6)); m = int(rng.integers(2, 50))
    A = rng.standard_normal((p, d, d)) * (0.5 / p)
    h = mimo_ar2ir(A, m)
    assert h.shape == (m, d, d)
    assert _conv_resid(A, h, p, d) < 1e-10


def test_h0_identity_and_h1_equals_A1():
    rng = np.random.default_rng(0)
    A = rng.standard_normal((2, 3, 3)) * 0.3
    h = mimo_ar2ir(A, 10)
    np.testing.assert_allclose(h[0], np.eye(3), atol=1e-12)
    np.testing.assert_allclose(h[1], A[0], atol=1e-12)


def test_var1_matches_matrix_power():
    # VAR(1): h_n = A1^n exactly
    rng = np.random.default_rng(7)
    A1 = rng.standard_normal((4, 4)) * 0.3
    h = mimo_ar2ir(A1[None], 25)
    for n in range(25):
        np.testing.assert_allclose(h[n], np.linalg.matrix_power(A1, n), atol=1e-10)


def test_block_toeplitz_solve():
    # independent dense solve of the block lower-triangular Toeplitz system T_A h = E1
    rng = np.random.default_rng(3)
    d, p, m = 3, 3, 20
    A = rng.standard_normal((p, d, d)) * 0.2
    h = mimo_ar2ir(A, m)
    T = np.zeros((m * d, m * d))
    for i in range(m):
        T[i * d:(i + 1) * d, i * d:(i + 1) * d] = np.eye(d)
        for k in range(1, min(i, p) + 1):
            T[i * d:(i + 1) * d, (i - k) * d:(i - k + 1) * d] = -A[k - 1]
    E1 = np.zeros((m * d, d)); E1[:d] = np.eye(d)
    h_solve = np.linalg.solve(T, E1).reshape(m, d, d)
    np.testing.assert_allclose(h, h_solve, atol=1e-10)


def test_roundtrip_ir_to_ar():
    # recover the AR coefficients from the IR; nothing should appear beyond lag p
    rng = np.random.default_rng(5)
    d, p, m = 3, 2, 30
    A = rng.standard_normal((p, d, d)) * 0.25
    h = mimo_ar2ir(A, m)
    cal = np.zeros((m, d, d)); cal[0] = np.eye(d)
    for n in range(1, m):
        cal[n] = -sum(cal[k] @ h[n - k] for k in range(0, n))
    np.testing.assert_allclose(-cal[1:p + 1], A, atol=1e-10)
    assert np.max(np.abs(cal[p + 1:])) < 1e-10


def test_left_equals_right_recursion():
    # h_n = sum A_k h_{n-k} (left, used by the code) == sum h_{n-k} A_k (right)
    rng = np.random.default_rng(1)
    d, p, m = 3, 3, 40
    A = rng.standard_normal((p, d, d)) * 0.25
    h_left = mimo_ar2ir(A, m)
    h_right = np.zeros((m, d, d)); h_right[0] = np.eye(d)
    for n in range(1, m):
        h_right[n] = sum(h_right[n - k] @ A[k - 1] for k in range(1, min(n, p) + 1))
    np.testing.assert_allclose(h_left, h_right, atol=1e-10)


def test_scalar_reduces_to_ar1():
    # d=1 reduces to the scalar ar2ir: AR(1) -> a**n
    a = 0.7
    h = mimo_ar2ir(np.array([[[a]]]), 8).ravel()
    np.testing.assert_allclose(h, a ** np.arange(8), atol=1e-12)


def test_complex_coefficients():
    rng = np.random.default_rng(9)
    d, p, m = 2, 2, 30
    A = (rng.standard_normal((p, d, d)) + 1j * rng.standard_normal((p, d, d))) * 0.15
    h = mimo_ar2ir(A, m)
    assert _conv_resid(A, h, p, d) < 1e-10


def test_rejects_bad_shape():
    with pytest.raises(ValueError):
        mimo_ar2ir(np.ones((3, 2)), 10)        # not (p, d, d)
    with pytest.raises(ValueError):
        mimo_ar2ir(np.ones((2, 3, 4)), 10)     # non-square blocks
