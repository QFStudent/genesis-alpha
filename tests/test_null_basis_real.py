"""
Tests for null_basis_realization_real -- the real TIB realization for poles that
may include complex-conjugate pairs.

It builds null_basis_realization on the pole magnitudes (real, input-balanced) and
injects each pair's phase with a 2x2 rotation read off the realized block, so:
A, B come out REAL, eig(A) equals the (complex) poles, the system stays
input-balanced and block-lower-triangular, and the impulse response is real --
unlike null_basis_realization on the complex poles directly, which is complex.
"""

import numpy as np
import scipy.linalg as la
import pytest

from ga.filters.tib import (null_basis_realization, null_basis_realization_real,
                            build_tib_from_directions, build_tib_from_directions_real,
                            krylov_basis, poles_chebyshev_roots)


def _conj_pairs(mags, angs):
    return np.array([z for m, a in zip(mags, angs)
                     for z in (m * np.exp(1j * a), m * np.exp(-1j * a))])


def _ir(A, B, C, lags):
    K = np.asarray(krylov_basis(A.astype(complex), B.astype(complex), lags))
    return np.tensordot(C, K, (-1, 0))


class TestNullBasisRealizationReal:
    def test_complex_pairs_realization(self):
        rng = np.random.default_rng(0)
        poles = _conj_pairs([0.5, 0.7, 0.6], [0.4, 0.9, 1.1])      # 3 pairs, n=6
        n = len(poles)
        v = rng.standard_normal((4, n))
        sys = null_basis_realization_real(poles, v)
        A = np.asarray(sys.A(dense=True))
        B = np.asarray(sys.B(dense=True))

        assert not np.iscomplexobj(A)                                           # A real
        assert not np.iscomplexobj(B)                                           # B real
        np.testing.assert_allclose(np.sort_complex(la.eigvals(A)),
                                   np.sort_complex(poles), atol=1e-9)            # eig(A) = poles
        assert np.max(np.abs(np.imag(la.eigvals(A)))) > 0.1                      # genuinely complex
        assert np.max(np.abs(A @ A.T + B @ B.T - np.eye(n))) < 1e-9             # input-balanced
        assert np.max(np.abs(np.triu(A, 2))) < 1e-9                             # block-lower-triangular

        C = rng.standard_normal((3, n))
        ir = _ir(np.real(A), np.real(B), C, 80)
        assert np.max(np.abs(np.imag(ir))) < 1e-9                               # real impulse response

    def test_mixed_real_and_complex(self):
        rng = np.random.default_rng(1)
        poles = np.concatenate([[0.3, -0.5], _conj_pairs([0.6, 0.7], [0.5, 1.0])])  # 2 real + 2 pairs
        n = len(poles)
        v = rng.standard_normal((3, n))
        sys = null_basis_realization_real(poles, v)
        A = np.asarray(sys.A(dense=True))
        B = np.asarray(sys.B(dense=True))

        assert not np.iscomplexobj(A)
        np.testing.assert_allclose(np.sort_complex(la.eigvals(A)),
                                   np.sort_complex(poles), atol=1e-9)
        assert np.max(np.abs(A @ A.T + B @ B.T - np.eye(n))) < 1e-9

    def test_unordered_input_is_handled(self):
        # poles given in a scrambled order (pair members not adjacent) must still work
        rng = np.random.default_rng(7)
        z0 = 0.6 * np.exp(1j * 0.5); z1 = 0.7 * np.exp(1j * 1.0)
        poles = np.array([z0, 0.3, np.conj(z1), z1, -0.4, np.conj(z0)])
        v = rng.standard_normal((3, 6))
        A = null_basis_realization_real(poles, v).A(dense=True)
        np.testing.assert_allclose(np.sort_complex(la.eigvals(A)),
                                   np.sort_complex(poles), atol=1e-9)

    def test_all_real_reduces_to_null_basis(self):
        # with only real poles, Q is the identity -> identical to null_basis_realization
        poles = poles_chebyshev_roots(5, 0.0, 0.8)
        rng = np.random.default_rng(2)
        v = rng.standard_normal((3, 5))
        A_real = np.asarray(null_basis_realization_real(poles, v).A(dense=True))
        A_base = np.real(np.asarray(null_basis_realization(poles, v).A(dense=True)))
        np.testing.assert_allclose(A_real, A_base, atol=1e-12)

    def test_recoverable_by_msvdreduce_complex(self):
        # end-to-end: a real complex-pole system from this builder must be recovered
        from ga.reducers.hankel import msvdreduce_complex
        rng = np.random.default_rng(3)
        poles = _conj_pairs([0.5, 0.7], [0.6, 1.1])               # 2 pairs, n=4
        n = len(poles)
        v = rng.standard_normal((3, n)); C = rng.standard_normal((2, n))
        sys = null_basis_realization_real(poles, v)
        A, B = np.real(sys.A(dense=True)), np.real(sys.B(dense=True))
        ir = np.real(_ir(A, B, C, 200))

        rsys, _ = msvdreduce_complex(ir, n)
        rA = np.asarray(rsys.A(dense=True))
        np.testing.assert_allclose(np.sort_complex(la.eigvals(rA)),
                                   np.sort_complex(poles), atol=1e-7)

    def test_rejects_unpaired_complex(self):
        poles = np.array([0.5 + 0.3j, 0.4, 0.2])                  # lone complex pole
        v = np.random.default_rng(0).standard_normal((2, 3))
        with pytest.raises(ValueError):
            null_basis_realization_real(poles, v)
