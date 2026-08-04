"""
Tests for ga.vc — Schlicht/Ludsteck time-varying coefficients estimator.

Checks mirror VC.m structure and the penalized-LS / ML objective in
Schlicht (2020) IZA DP 12920.
"""

import numpy as np
import pytest

from ga.vc.schlicht import (
    VCEstimateResult,
    make_P,
    make_block_X,
    make_S,
    ml_criterion,
    solve_coefficients,
    vc_estimate,
)
from ga.vc.simulate import (
    fitted_values,
    random_walk_residuals,
    simulate_vc_model,
    stacked_state,
)


class TestMatrixConstruction:
    def test_block_X_shape(self):
        x = np.arange(12, dtype=float).reshape(4, 3)
        X = make_block_X(x)
        assert X.shape == (4, 12)
        assert X[2, 6] == x[2, 0]
        assert X[2, 8] == x[2, 2]
        assert X[2, 5] == 0.0

    def test_P_applies_random_walk(self):
        n, T = 2, 4
        coeff = np.arange(T * n, dtype=float).reshape(T, n)
        a = stacked_state(coeff)
        Pa = make_P(n, T) @ a
        expected = (coeff[1:] - coeff[:-1]).reshape(-1)
        np.testing.assert_allclose(Pa, expected)


class TestPenalizedLeastSquares:
    @pytest.fixture
    def small_problem(self):
        rng = np.random.default_rng(0)
        T, n = 20, 2
        x = rng.standard_normal((T, n))
        y = rng.standard_normal(T)
        r = np.array([0.5, 2.0])
        return x, y, r

    def test_normal_equations(self, small_problem):
        x, y, r = small_problem
        a, M, X, P, S = solve_coefficients(x, y, r)
        tXy = X.T @ y
        np.testing.assert_allclose(M @ a, tXy, atol=1e-10)

    def test_Q_decomposition(self, small_problem):
        x, y, r = small_problem
        a, M, X, P, S = solve_coefficients(x, y, r)
        Pa = P @ a
        u = y - X @ a
        Q = float(u @ u + Pa @ S @ Pa)
        assert Q > 0
        T, n = x.shape
        assert ml_criterion(r, x, y) == pytest.approx(
            np.linalg.slogdet(M.toarray())[1]
            + (T - n) * np.log(Q)
            + (T - 1) * np.sum(np.log(r))
        )


class TestVCEstimate:
    def test_returns_expected_fields(self):
        sim = simulate_vc_model(T=40, n=3, sigma_u=1.0, sigma_v=np.array([0.05, 0.1, 0.15]))
        res = vc_estimate(sim["x"], sim["y"], optimize_ratios=False, variance_ratios=1.0)
        assert isinstance(res, VCEstimateResult)
        assert res.coeff.shape == (40, 3)
        assert res.sdb.shape == (40, 3)
        assert res.variance_ratios.shape == (3,)
        assert res.sdu > 0
        assert np.all(res.sdi > 0)

    def test_in_sample_fit_beats_ols_for_rw_data(self):
        sim = simulate_vc_model(
            T=60,
            n=3,
            sigma_u=0.5,
            sigma_v=np.array([0.08, 0.12, 0.06]),
            seed=7,
        )
        res = vc_estimate(sim["x"], sim["y"], variance_ratios=1.0)
        yhat_vc = fitted_values(res.coeff, sim["x"])
        coef_ols, _, _, _ = np.linalg.lstsq(sim["x"], sim["y"], rcond=None)
        yhat_ols = sim["x"] @ coef_ols
        sse_vc = np.sum((sim["y"] - yhat_vc) ** 2)
        sse_ols = np.sum((sim["y"] - yhat_ols) ** 2)
        assert sse_vc < sse_ols

    def test_recovers_variance_ratios_on_simulation(self):
        true_r = np.array([0.04, 0.16, 0.09]) / 0.25
        sim = simulate_vc_model(
            T=120,
            n=3,
            sigma_u=0.5,
            sigma_v=np.sqrt(true_r * 0.25),
            seed=11,
        )
        res = vc_estimate(sim["x"], sim["y"], variance_ratios=1.0)
        # Order-of-magnitude recovery (ML is noisy on short samples).
        for est, truth in zip(res.variance_ratios, true_r):
            assert est == pytest.approx(truth, rel=0.8, abs=0.15)

    def test_large_ratios_approach_constant_ols(self):
        sim = simulate_vc_model(T=50, n=2, sigma_u=1.0, sigma_v=np.array([1e-6, 1e-6]), seed=3)
        res = vc_estimate(sim["x"], sim["y"], variance_ratios=1e-6, optimize_ratios=False)
        ols, _, _, _ = np.linalg.lstsq(sim["x"], sim["y"], rcond=None)
        mean_vc = res.coeff.mean(axis=0)
        np.testing.assert_allclose(mean_vc, ols, rtol=0.15)

    def test_coefficient_path_tracks_true_dgp(self):
        sim = simulate_vc_model(
            T=80,
            n=2,
            sigma_u=0.3,
            sigma_v=np.array([0.05, 0.05]),
            seed=19,
        )
        res = vc_estimate(sim["x"], sim["y"], variance_ratios=1.0)
        for j in range(sim["b"].shape[1]):
            corr = np.corrcoef(res.coeff[:, j], sim["b"][:, j])[0, 1]
            assert corr > 0.6

    def test_random_walk_residuals_match_Pa(self):
        sim = simulate_vc_model(T=30, n=2, sigma_u=1.0, sigma_v=np.array([0.1, 0.1]))
        res = vc_estimate(sim["x"], sim["y"], optimize_ratios=False)
        Pa = random_walk_residuals(res.coeff)
        a = stacked_state(res.coeff)
        Pa_direct = (make_P(2, 30) @ a).reshape(29, 2)
        np.testing.assert_allclose(Pa, Pa_direct)


class TestValidation:
    def test_T_must_exceed_n(self):
        x = np.ones((3, 4))
        y = np.ones(3)
        with pytest.raises(ValueError, match="T > n"):
            vc_estimate(x, y)

    def test_mismatched_lengths(self):
        x = np.ones((10, 2))
        y = np.ones(9)
        with pytest.raises(ValueError, match="same number of rows"):
            vc_estimate(x, y)
