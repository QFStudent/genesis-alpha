"""
Tests for tib_from_state_space — the corrected BT -> TIB extraction.

The regression these lock in: a target system that is exactly TIB-representable
must survive the round-trip  target -> balanced_truncation -> tib_from_state_space
both as a valid TIB realization (A A* + B B* = I) and as the same transfer
function (impulse response). This is the path that extract_poles_and_nullvecs_from_bt
gets wrong (see docs TIBForm; scripts/compare_null_vectors.py quantifies the gap).
"""

import numpy as np
import pytest

from ga.filters.tib import null_basis_realization, mimo_poles_to_ir, poles_chebyshev_roots
from ga.reducers.balanced_truncation import balanced_truncation, tib_from_state_space


def _impulse_response(A, B, C, lags):
    """IR[:, :, k] = C A^k B."""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    p, m = C.shape[0], B.shape[1]
    ir = np.zeros((p, m, lags))
    Ak = np.eye(A.shape[0])
    for k in range(lags):
        ir[:, :, k] = C @ Ak @ B
        Ak = Ak @ A
    return ir


def _random_unit_nullvecs(m, n, seed):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal((m, n))
    return v / np.linalg.norm(v, axis=0)


@pytest.fixture
def target():
    """Order-12, 6-input, 4-output system with cross-coupled null vectors."""
    P, Q, P_OUT, LAGS = 12, 6, 4, 40
    poles = poles_chebyshev_roots(P)
    v_true = _random_unit_nullvecs(Q, P, seed=7)
    C_true = np.random.default_rng(7).standard_normal((P_OUT, P))
    sys = null_basis_realization(poles, v_true)
    ir = mimo_poles_to_ir(poles, v_true, C_true, LAGS)
    return dict(poles=poles, A=sys.A(dense=True), B=sys.B(dense=True),
                C=C_true, ir=ir, lags=LAGS)


class TestTibFromStateSpace:
    def test_roundtrip_reproduces_impulse_response(self, target):
        result = balanced_truncation(target["A"], target["B"], target["C"],
                                     order=len(target["poles"]))
        sys, transform = tib_from_state_space(result.A_r, result.B_r)
        C_tib = result.C_r @ transform
        ir = _impulse_response(sys.A(dense=True), sys.B(dense=True), C_tib, target["lags"])
        rel_err = np.linalg.norm(ir - target["ir"]) / np.linalg.norm(target["ir"])
        assert rel_err < 1e-10

    def test_result_is_input_balanced(self, target):
        result = balanced_truncation(target["A"], target["B"], target["C"],
                                     order=len(target["poles"]))
        sys, _ = tib_from_state_space(result.A_r, result.B_r)
        A, B = sys.A(dense=True), sys.B(dense=True)
        resid = A @ A.conj().T + B @ B.conj().T - np.eye(A.shape[0])
        assert np.max(np.abs(resid)) < 1e-10

    def test_poles_on_diagonal(self, target):
        result = balanced_truncation(target["A"], target["B"], target["C"],
                                     order=len(target["poles"]))
        sys, _ = tib_from_state_space(result.A_r, result.B_r)
        diag = np.sort(np.real(np.diag(sys.A(dense=True))))
        np.testing.assert_allclose(diag, np.sort(target["poles"]), atol=1e-9)

    def test_lower_triangular(self, target):
        result = balanced_truncation(target["A"], target["B"], target["C"],
                                     order=len(target["poles"]))
        sys, _ = tib_from_state_space(result.A_r, result.B_r)
        A = sys.A(dense=True)
        assert np.max(np.abs(np.triu(A, k=1))) < 1e-9

def _rot(r, th):
    return r * np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])


@pytest.fixture
def complex_target():
    """Real system with complex-conjugate poles (+ one real pole), in a random
    basis so the block structure is not handed to the algorithm for free."""
    import scipy.linalg as sla
    blocks = [_rot(0.8, 0.5), _rot(0.6, 1.2), np.array([[0.3]]), _rot(0.7, 0.9)]
    A0 = sla.block_diag(*blocks)
    n = A0.shape[0]
    rng = np.random.default_rng(0)
    Tm = rng.standard_normal((n, n))
    A = Tm @ A0 @ sla.inv(Tm)
    B = rng.standard_normal((n, 3))
    C = rng.standard_normal((2, n))
    LAGS = 40
    ir = _impulse_response(A, B, C, LAGS)
    return dict(A=A, B=B, C=C, ir=ir, lags=LAGS, n=n)


class TestComplexPoles:
    """Complex-conjugate poles must survive as 2x2 real diagonal blocks."""

    def _run(self, t):
        result = balanced_truncation(t["A"], t["B"], t["C"], order=t["n"])
        sys, transform = tib_from_state_space(result.A_r, result.B_r)
        return sys, result.C_r @ transform

    def test_roundtrip_reproduces_impulse_response(self, complex_target):
        sys, C_tib = self._run(complex_target)
        ir = _impulse_response(sys.A(dense=True), sys.B(dense=True), C_tib,
                               complex_target["lags"])
        rel_err = np.linalg.norm(ir - complex_target["ir"]) / np.linalg.norm(complex_target["ir"])
        assert rel_err < 1e-10

    def test_result_is_input_balanced(self, complex_target):
        sys, _ = self._run(complex_target)
        A, B = sys.A(dense=True), sys.B(dense=True)
        resid = A @ A.conj().T + B @ B.conj().T - np.eye(A.shape[0])
        assert np.max(np.abs(resid)) < 1e-10

    def test_realization_is_real(self, complex_target):
        sys, _ = self._run(complex_target)
        assert np.max(np.abs(np.imag(sys.A(dense=True)))) < 1e-12
        assert np.max(np.abs(np.imag(sys.B(dense=True)))) < 1e-12

    def test_eigenvalues_recovered(self, complex_target):
        sys, _ = self._run(complex_target)
        got = np.sort_complex(np.linalg.eigvals(sys.A(dense=True)))
        want = np.sort_complex(np.linalg.eigvals(complex_target["A"]))
        np.testing.assert_allclose(got, want, atol=1e-9)

    def test_block_lower_triangular(self, complex_target):
        # block-lower-triangular = lower-Hessenberg with isolated 2x2 blocks:
        # nothing above the first superdiagonal, but the superdiagonal itself
        # is nonzero (the 2x2 complex blocks).
        sys, _ = self._run(complex_target)
        A = sys.A(dense=True)
        assert np.max(np.abs(np.triu(A, k=2))) < 1e-9
        assert np.max(np.abs(np.triu(A, k=1))) > 1e-3
