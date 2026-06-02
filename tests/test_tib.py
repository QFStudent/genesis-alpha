"""
Tests for ga.filters.tib — TIB state-space model functions.

Correctness is checked via internal mathematical properties (TIB identity,
pole recovery, Blaschke null/unitary properties) rather than against any
external reference implementation.
"""

import numpy as np
import pytest

from ga.filters.tib import (
    TIBStateSpace,
    poles_chebyshev_roots,
    mimo_standard_null_vectors,
    blaschke_factor,
    blaschke_potapov_factor,
    normalize,
    null_basis_realization,
    krylov_basis,
    mimo_poles_to_ir,
)


# --------------------------------------------------------------------------
# Fixtures / helpers
# --------------------------------------------------------------------------

@pytest.fixture
def real_poles():
    return np.array([0.5, 0.3, -0.2, 0.1], dtype=complex)


def _unit_null_vectors(m, n, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal((m, n)) + 1j * rng.standard_normal((m, n))
    return v / np.linalg.norm(v, axis=0)


# --------------------------------------------------------------------------
# TIBStateSpace
# --------------------------------------------------------------------------

class TestTIBStateSpace:
    def test_p_q_properties(self):
        A = np.eye(3, dtype=complex)
        B = np.zeros((3, 2), dtype=complex)
        sys = TIBStateSpace(A_=A, B_=B)
        assert sys.p == 3       # number of poles (states)
        assert sys.q == 2       # number of inputs

    def test_A_dense_vs_sparse(self):
        A = np.eye(3, dtype=complex)
        B = np.zeros((3, 2), dtype=complex)
        sys = TIBStateSpace(A_=A, B_=B)
        assert hasattr(sys.A(), "toarray")          # sparse by default
        assert isinstance(sys.A(dense=True), np.ndarray)
        np.testing.assert_allclose(sys.A(dense=True), A)

    def test_str_does_not_raise(self):
        sys = TIBStateSpace(A_=np.eye(2, dtype=complex), B_=np.zeros((2, 1), dtype=complex))
        assert "TIBStateSpace" in str(sys)


# --------------------------------------------------------------------------
# poles_chebyshev_roots
# --------------------------------------------------------------------------

class TestChebyshevRoots:
    def test_count(self):
        assert len(poles_chebyshev_roots(5)) == 5

    def test_within_interval(self):
        p = poles_chebyshev_roots(8, a=0.0, b=1.0)
        assert np.all(p > 0.0) and np.all(p < 1.0)

    def test_ascending(self):
        p = poles_chebyshev_roots(6)
        assert np.all(np.diff(p) > 0)               # reversed -> ascending

    def test_rescaled_interval(self):
        p = poles_chebyshev_roots(4, a=-1.0, b=1.0)
        assert np.all(p > -1.0) and np.all(p < 1.0)


# --------------------------------------------------------------------------
# mimo_standard_null_vectors
# --------------------------------------------------------------------------

class TestStandardNullVectors:
    def test_shape(self):
        U = mimo_standard_null_vectors(6, 2)
        assert U.shape == (6, 2)

    def test_siso_is_ones_column(self):
        U = mimo_standard_null_vectors(4, 1)
        np.testing.assert_allclose(U, np.ones((4, 1)))


# --------------------------------------------------------------------------
# blaschke_factor (scalar)
# --------------------------------------------------------------------------

class TestBlaschkeFactor:
    def test_zero_at_own_pole(self):
        w = 0.4 + 0.1j
        assert abs(blaschke_factor(w, w)) < 1e-14

    def test_unimodular_on_unit_circle(self):
        w = 0.3 - 0.2j
        for theta in np.linspace(0, 2 * np.pi, 11):
            z = np.exp(1j * theta)
            assert abs(abs(blaschke_factor(w, z)) - 1.0) < 1e-12

    def test_real_pole_real_axis(self):
        # B_w(z) = (z - w) / (1 - w z) for real w, z
        w, z = 0.5, 0.3
        assert abs(blaschke_factor(w, z) - (z - w) / (1 - w * z)) < 1e-14


# --------------------------------------------------------------------------
# blaschke_potapov_factor (matrix)
# --------------------------------------------------------------------------

class TestBlaschkePotapovFactor:
    def test_null_property(self):
        # β_{w,u}(w) u = 0  (eq. 249 in Kong dissertation)
        w = 0.4 + 0.1j
        u = normalize(np.array([[1.0], [2.0], [0.5]], dtype=complex))
        bp = blaschke_potapov_factor(w, w, u)
        assert np.linalg.norm(bp @ u) < 1e-13

    def test_unitary_on_unit_circle(self):
        # With J = I, the BP factor is unitary for |z| = 1
        w = 0.3 + 0.2j
        u = normalize(np.array([[1.0], [1j], [0.5]], dtype=complex))
        z = np.exp(1j * 0.7)
        bp = blaschke_potapov_factor(w, z, u)
        np.testing.assert_allclose(bp @ bp.conj().T, np.eye(3), atol=1e-12)

    def test_shape(self):
        u = normalize(np.array([[1.0], [0.0]], dtype=complex))
        assert blaschke_potapov_factor(0.2, 0.5, u).shape == (2, 2)

    def test_rejects_short_vector(self):
        u = np.array([[0.5], [0.0]], dtype=complex)   # norm 0.5 < 1
        with pytest.raises(ValueError):
            blaschke_potapov_factor(0.2, 0.5, u)

    def test_rejects_long_vector(self):
        u = np.array([[2.0], [0.0]], dtype=complex)   # norm 2 > 1
        with pytest.raises(ValueError):
            blaschke_potapov_factor(0.2, 0.5, u)


# --------------------------------------------------------------------------
# normalize
# --------------------------------------------------------------------------

class TestNormalize:
    def test_unit_norm(self):
        x = np.array([3.0, 4.0])
        assert abs(np.linalg.norm(normalize(x)) - 1.0) < 1e-14

    def test_raises_on_zero(self):
        with pytest.raises(ValueError):
            normalize(np.zeros(3))


# --------------------------------------------------------------------------
# null_basis_realization — the core function
# --------------------------------------------------------------------------

class TestNullBasisRealization:
    def test_tib_identity_unit_vectors(self, real_poles):
        v = _unit_null_vectors(2, len(real_poles))
        sys = null_basis_realization(real_poles, v)
        A, B = sys.A(dense=True), sys.B(dense=True)
        resid = A @ A.conj().T + B @ B.conj().T - np.eye(len(real_poles))
        assert np.max(np.abs(resid)) < 1e-10

    def test_tib_identity_nonunit_vectors(self, real_poles):
        # Regression for the line-165 bug: B must use the *normalized* y1,
        # so the TIB identity must hold even when v columns are NOT unit norm.
        v = _unit_null_vectors(2, len(real_poles))
        v = v * np.array([2.0, 0.5, 3.0, 1.0])      # break unit norm
        sys = null_basis_realization(real_poles, v)
        A, B = sys.A(dense=True), sys.B(dense=True)
        resid = A @ A.conj().T + B @ B.conj().T - np.eye(len(real_poles))
        assert np.max(np.abs(resid)) < 1e-10

    def test_pole_recovery(self, real_poles):
        v = _unit_null_vectors(2, len(real_poles))
        sys = null_basis_realization(real_poles, v)
        A = sys.A(dense=True)
        ev = np.sort_complex(np.linalg.eigvals(A))
        np.testing.assert_allclose(ev, np.sort_complex(real_poles), atol=1e-10)

    def test_shapes(self, real_poles):
        m = 2
        v = _unit_null_vectors(m, len(real_poles))
        sys = null_basis_realization(real_poles, v)
        assert sys.A(dense=True).shape == (len(real_poles), len(real_poles))
        assert sys.B(dense=True).shape == (len(real_poles), m)

    def test_siso(self):
        lambd = np.array([0.6, 0.2, -0.3], dtype=complex)
        v = np.ones((1, 3), dtype=complex)
        sys = null_basis_realization(lambd, v)
        A, B = sys.A(dense=True), sys.B(dense=True)
        resid = A @ A.conj().T + B @ B.conj().T - np.eye(3)
        assert np.max(np.abs(resid)) < 1e-10

    def test_complex_conjugate_poles(self):
        lambd = np.array([0.5 + 0.3j, 0.5 - 0.3j, 0.2], dtype=complex)
        v = _unit_null_vectors(2, 3, seed=3)
        sys = null_basis_realization(lambd, v)
        A, B = sys.A(dense=True), sys.B(dense=True)
        resid = A @ A.conj().T + B @ B.conj().T - np.eye(3)
        assert np.max(np.abs(resid)) < 1e-10

    def test_rejects_unstable_pole(self):
        lambd = np.array([0.5, 1.5], dtype=complex)     # outside unit circle
        v = _unit_null_vectors(2, 2)
        with pytest.raises(ValueError):
            null_basis_realization(lambd, v)

    def test_rejects_zero_null_vector(self, real_poles):
        v = _unit_null_vectors(2, len(real_poles))
        v[:, 1] = 0.0
        with pytest.raises(ValueError):
            null_basis_realization(real_poles, v)

    def test_rejects_wrong_v_shape(self, real_poles):
        v = _unit_null_vectors(2, len(real_poles) - 1)   # wrong n
        with pytest.raises(ValueError):
            null_basis_realization(real_poles, v)


# --------------------------------------------------------------------------
# krylov_basis
# --------------------------------------------------------------------------

class TestKrylovBasis:
    def test_first_block_is_B(self):
        A = np.array([[0.5, 0.0], [0.1, 0.3]])
        B = np.array([[1.0], [0.0]])
        K = krylov_basis(A, B, 5)        # (p, n) since q == 1
        np.testing.assert_allclose(K[:, 0], B[:, 0])

    def test_recurrence_siso(self):
        A = np.array([[0.5, 0.0], [0.1, 0.3]])
        B = np.array([[1.0], [0.0]])
        n = 6
        K = krylov_basis(A, B, n)
        expected = B[:, 0].copy()
        for k in range(n):
            np.testing.assert_allclose(K[:, k], expected, atol=1e-12)
            expected = A @ expected

    def test_mimo_shape(self):
        A = np.eye(4)
        B = np.ones((4, 2))
        K = krylov_basis(A, B, 7)
        assert K.shape == (4, 2, 7)


# --------------------------------------------------------------------------
# mimo_poles_to_ir
# --------------------------------------------------------------------------

class TestMimoPolesToIR:
    def test_shape(self, real_poles):
        m, p, L = 2, 3, 20
        v = _unit_null_vectors(m, len(real_poles))
        C = np.random.default_rng(0).standard_normal((p, len(real_poles)))
        ir = mimo_poles_to_ir(real_poles, v, C, L)
        assert ir.shape == (p, m, L)

    def test_h0_equals_CB(self, real_poles):
        # h_0 = C A^0 B = C B  (complex when v is complex)
        m, p, L = 2, 3, 10
        v = _unit_null_vectors(m, len(real_poles))
        C = np.random.default_rng(1).standard_normal((p, len(real_poles)))
        sys = null_basis_realization(real_poles, v)
        B = sys.B(dense=True)
        ir = mimo_poles_to_ir(real_poles, v, C, L)
        np.testing.assert_allclose(ir[:, :, 0], C @ B, atol=1e-10)

    def test_h1_equals_CAB(self, real_poles):
        m, p, L = 2, 3, 10
        v = _unit_null_vectors(m, len(real_poles))
        C = np.random.default_rng(2).standard_normal((p, len(real_poles)))
        sys = null_basis_realization(real_poles, v)
        A, B = sys.A(dense=True), sys.B(dense=True)
        ir = mimo_poles_to_ir(real_poles, v, C, L)
        np.testing.assert_allclose(ir[:, :, 1], C @ A @ B, atol=1e-10)
