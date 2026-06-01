---
type: paper
tags: [market-impact, execution, factor-model, multi-asset]
sources: [Kong2018-TIBInfoGeometry]
updated: 2026-05-31
---

# Information Geometry and Dimensionality Reduction

**Authors:** Kong, C.  **Venue / Year:** PhD Dissertation, Stony Brook University, 2018  
**Advisor:** Andrew Mullhaupt (Quantitative Finance, Stony Brook)

## Contribution

Foundational paper for the TIB framework used throughout genesis-alpha. Four main contributions:
1. Generalises information distance (KL divergence / Fisher metric) from SISO to **MIMO LTI systems** — when cepstrum coefficients are used as coordinates the manifold becomes Euclidean, making distance computation tractable.
2. Develops the **Triangular Input Balanced (TIB) form** — a permutation-invariant state-space parameterisation whose system matrix is A = M⁻¹NQ with bidiagonal M, N. Numerically superior to generic state-space forms.
3. Proposes **Null Basis Parametrisation** for MIMO TIB reduction — the null vectors U (rows unitary, U[:q,:] unitary) encode input directions; each state has an interpretable temporal shape and input direction.
4. Develops efficient FFT-based algorithms for **large block-Toeplitz matrix** operations (square root, logarithm), enabling MIMO cepstrum computation and reduction on statistical manifolds.

## Methodology

- **Universe:** MIMO discrete-time LTI systems, general p inputs / q outputs
- **Core representation:** State-space (A, B, C, D); A = M⁻¹NQ where M, N are bidiagonal sparse matrices, Q is a rotation matrix keeping entries real even for complex poles
- **Information distance:** Geodesic distance on the Fisher-information Riemannian manifold; equivalent to H₂-norm of the log transfer function cepstrum
- **Model reduction methods compared:** Balanced truncation, POD, tangential interpolation (IRKA), and the new null-basis SVD method
- **Key benchmark:** H₂ norm and Hankel norm error vs. number of poles retained

## Key Results

- TIB form is equivalent to balanced realisation but admits a sparse bidiagonal parameterisation (Table 1–2 in Ch. 7): e.g. power-law decay system reduced from 100 → 5 poles with H₂ error < 0.01 (pre-cost, exact numbers in tables).
- Null Basis Reduction outperforms Hankel SVD on non-rational systems (infinite-dimensional) in Hankel norm (Figs 5–6).
- Fast block-Toeplitz logarithm: O(n log²n) vs O(n³) for direct method; error < 10⁻⁶ for matrices up to n = 1000 (Tables 3–8, Figs 7–9).
- Cepstrum reduction competitive with balanced truncation in H₂ norm (Figs 12–17).

## Limitations / Caveats

- Stability requirement: poles strictly inside unit disk; non-stationary processes (random walk) are handled via finite information-length extension but require care.
- Null-basis reduction produces **approximate** null vectors when starting from balanced-truncation output (BT doesn't preserve TIB structure — see `extract_poles_and_nullvecs_from_bt` in `ga/reducers/balanced_truncation.py`).
- Numerical results are on synthetic systems; no live-trading capacity estimates.

## Connection to genesis-alpha

- `ga/filters/tib.py` is a direct implementation of TIB form: `TIBSystem`, `siso_system_matrices_real`, `mimo_system_matrices`, `krylov_basis`, `blkhankel`, `msvdreduce`, `mimo_null_basis`.
- `ga/reducers/hankel.py` implements Ch. 7 reduction: `FastHankelProduct` (FFT Hankel-vector multiply), `reduce_fft_truncate`, `reduce_svd_truncate`, `info_svd_reduce`.
- The interpretability of TIB states (each state = decay rate + input direction) is valuable for the execution simulator: states can be labelled as "long-memory market impact", "short-term reversion", etc.
- Connects to `backtest/costs/` via the impulse response representation of market impact.

## Related Pages

[[concepts/TIBForm]], [[concepts/ModelReduction]], [[concepts/InformationGeometry]], [[concepts/MarketImpact]], [[papers/Mu2026-ModelReductionNotes]], [[papers/MullhauptRiedel2003-TIBBandMatrix]]
