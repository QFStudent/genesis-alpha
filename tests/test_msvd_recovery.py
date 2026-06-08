"""
Tests for the corrected Hankel-SVD recovery:

- msvdreduce_fixed  -- finishes the null-vector deflation that msvdreduce computes
  then discards ("Bug 2"), so it returns the orthonormal null vectors y_k.
- build_tib_from_directions -- rebuilds a TIB pair from poles + already-orthonormal
  directions, WITHOUT the Blaschke re-deflation that null_basis_realization would
  wrongly apply to them.

Regression these lock in: a known MIMO TIB system, reduced to its own order from its
impulse response alone, must come back -- exact poles, exact transfer function,
input-balanced -- when rebuilt with build_tib_from_directions; and the shortcut of
feeding the recovered directions back into null_basis_realization must NOT reproduce
it (double-deflation; see docs/derivations/03-null-vector-recovery.md).
"""

import numpy as np
import pytest

from ga.filters.tib import (null_basis_realization, build_tib_from_directions,
                            mimo_poles_to_ir, krylov_basis, poles_chebyshev_roots)
from ga.reducers.hankel import msvdreduce_fixed


def _ir(A, B, C, lags):
    K = np.asarray(krylov_basis(A, B, lags))
    return np.real(np.tensordot(C, K, (-1, 0)))


def _fit_C(A, B, lags, target):
    K = np.asarray(krylov_basis(A, B, lags)).reshape(A.shape[0], -1)
    return np.real(target.reshape(target.shape[0], -1) @ np.linalg.pinv(K))


def _relerr(ir, ref):
    return np.linalg.norm(ir - ref) / np.linalg.norm(ref)


@pytest.fixture
def target():
    """Known order-6, 2-input, 3-output MIMO TIB system (real poles)."""
    n, q, p, lags = 6, 2, 3, 80
    poles = poles_chebyshev_roots(n, 0.0, 0.85)
    rng = np.random.default_rng(3)
    v = rng.standard_normal((q, n)); v /= np.linalg.norm(v, axis=0)
    C = rng.standard_normal((p, n))
    sys = null_basis_realization(poles, v)
    ir = np.real(mimo_poles_to_ir(poles, v, C, lags))
    return dict(n=n, q=q, p=p, lags=lags, poles=poles,
                A=np.real(sys.A(dense=True)), B=np.real(sys.B(dense=True)), C=C, ir=ir)


class TestMsvdRecovery:
    def test_poles_recovered_full_order(self, target):
        w, _ = msvdreduce_fixed(target["ir"], target["n"])
        np.testing.assert_allclose(np.sort(np.real(w)), np.sort(target["poles"]), atol=1e-9)

    def test_recovered_null_vectors_are_unit(self, target):
        _, u = msvdreduce_fixed(target["ir"], target["n"])
        np.testing.assert_allclose(np.linalg.norm(np.real(u), axis=0), 1.0, atol=1e-9)

    def test_roundtrip_reproduces_impulse_response(self, target):
        w, u = msvdreduce_fixed(target["ir"], target["n"])
        sys = build_tib_from_directions(np.real(w), np.real(u))
        A, B = sys.A(dense=True), sys.B(dense=True)
        C = _fit_C(A, B, target["lags"], target["ir"])
        assert _relerr(_ir(A, B, C, target["lags"]), target["ir"]) < 1e-7

    def test_rebuilt_is_input_balanced(self, target):
        w, u = msvdreduce_fixed(target["ir"], target["n"])
        sys = build_tib_from_directions(np.real(w), np.real(u))
        A, B = sys.A(dense=True), sys.B(dense=True)
        resid = A @ A.conj().T + B @ B.conj().T - np.eye(A.shape[0])
        assert np.max(np.abs(resid)) < 1e-9

    def test_do_not_redeflate(self, target):
        # build_tib_from_directions rebuilds correctly; feeding the SAME recovered
        # directions to null_basis_realization re-applies the Blaschke deflation
        # (double-deflation) and must NOT reproduce the system.
        w, u = msvdreduce_fixed(target["ir"], target["n"])
        good = build_tib_from_directions(np.real(w), np.real(u))
        bad = null_basis_realization(np.real(w), np.real(u))
        Ag, Bg = good.A(dense=True), good.B(dense=True)
        Ab, Bb = np.real(bad.A(dense=True)), np.real(bad.B(dense=True))
        Cg = _fit_C(Ag, Bg, target["lags"], target["ir"])
        Cb = _fit_C(Ab, Bb, target["lags"], target["ir"])
        assert _relerr(_ir(Ag, Bg, Cg, target["lags"]), target["ir"]) < 1e-7
        assert _relerr(_ir(Ab, Bb, Cb, target["lags"]), target["ir"]) > 1e-2


class TestBuildFromDirections:
    """build_tib_from_directions produces a valid TIB pair for any orthonormal y."""

    def test_valid_tib_for_arbitrary_directions(self):
        n, q = 5, 3
        poles = poles_chebyshev_roots(n, 0.0, 0.8)
        rng = np.random.default_rng(11)
        y = rng.standard_normal((q, n)); y /= np.linalg.norm(y, axis=0)
        sys = build_tib_from_directions(poles, y)
        A, B = sys.A(dense=True), sys.B(dense=True)
        assert np.max(np.abs(A @ A.conj().T + B @ B.conj().T - np.eye(n))) < 1e-9   # balanced
        assert np.max(np.abs(np.triu(A, k=1))) < 1e-9                               # lower-triangular
        np.testing.assert_allclose(np.sort(np.real(np.diag(A))), np.sort(poles), atol=1e-9)

    def test_rejects_shape_mismatch(self):
        poles = poles_chebyshev_roots(4)
        with pytest.raises(ValueError):
            build_tib_from_directions(poles, np.ones((2, 3)))   # y has n=3, poles has n=4
