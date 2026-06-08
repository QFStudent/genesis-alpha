---
type: concept
tags: [execution, market-impact]
sources: [Kong2018-TIBInfoGeometry, Mu2026-ModelReductionNotes]
updated: 2026-06-07
---

# Model Reduction for LTI Systems

## Definition

Given a high-order LTI system (A∈ℝⁿˣⁿ, B, C), find a low-order system (Aᵣ∈ℝʳˣʳ, Bᵣ, Cᵣ) with r ≪ n that approximates the input-output behaviour. The key metric is the H∞ or Hankel norm of the approximation error.

Three main families ([[papers/Kong2018-TIBInfoGeometry]], Ch. 2) — POD and tangential interpolation (IRKA) also exist but are not used in genesis-alpha:

| Method | Input | Output | Error norm |
|--------|-------|--------|------------|
| **Balanced Truncation (BT)** | State-space (A,B,C) | Reduced state-space | H∞ with provable bound |
| **Hankel SVD / msvdreduce** | Impulse response | Poles + null vectors | Hankel norm |
| **Info SVD (cepstrum)** | Impulse response | Poles | Information distance |
| POD | Snapshots | Reduced basis | L² projection error — *not used* |
| IRKA | Transfer function samples | Poles | H₂ — *not used* |

> 📐 **Full derivation:** `docs/derivations/01-hankel-svd-reduction.md` — the `msvdreduce` (Hankel-SVD) algorithm from Yu §6.2, with a boxed **BT vs `msvdreduce` separation**.

## Why It Matters for genesis-alpha

Market impact is modelled as the impulse response of a causal LTI filter from order flow to price. Reduction is necessary because:
- High-order impulse responses (100+ lags) are expensive to evaluate in real-time
- Overfitting risk with too many poles
- Interpretability: a 3–8 pole TIB model gives actionable state labels

The pipeline in genesis-alpha:
```
Market data → Estimate IR → msvdreduce / BT → TIB(A_r, B_r, C_r) → execution/simulator/
```

## Empirical Findings

**Balanced Truncation:**
- H∞ error bound: ‖Gₙ - Gᵣ‖_{H∞} ≤ 2(σᵣ₊₁ + ⋯ + σₙ) — tight in practice ([[papers/Mu2026-ModelReductionNotes]])
- Square-root (SR) method is numerically robust — avoids direct Grammian inversion
- Requires stable system (all poles inside unit disk); ill-conditioned for near-unstable systems

**Hankel SVD (msvdreduce):**
- Data-driven: works directly from impulse response, no state-space required
- For power-law decay systems: 5 poles captures > 99% of **Hankel energy** (Fig 2) and achieves H₂ error < 0.01 (Fig 1) ([[papers/Kong2018-TIBInfoGeometry]], Table 1 — both metrics reported separately)
- FFT-based fast variant (`reduce_fft_truncate`) scales to large n

**Info SVD (cepstrum reduction):**
- Minimises information distance (geodesic on Fisher manifold)
- Better than BT for non-rational/infinite-dimensional systems
- Implemented in `info_svd_reduce` in `ga/reducers/hankel.py`

> ⚠️ **`info_svd_reduce` does not actually realize this (gap found 2026-06-07).** The "minimises information distance" description is *aspirational, not implemented*. As written the function cannot recover poles — it runs the pole-finder directly on the cepstrum `aₖ = (1/k)Σλᵏ`, which the `1/k` makes **not** a low-rank-Hankel sequence — and it silently drops its own `rho`-damping, and shares the `msvdreduce` transpose-vs-pinv shift bug. A corrected companion `info_svd_reduce_fixed` (multiply by `k` → power sums `Σλᵏ`, then shift) recovers poles **only for identification / full order**, *not* stable order *reduction* (the power-sum Hankel has no energy ordering). **Genuine information-distance reduction (Kong2018) is a different, more involved algorithm — not yet implemented.** See `docs/derivations/01-hankel-svd-reduction.md` (SISO implementation-note) and [[concepts/InformationGeometry]].

> ⚠️ BT and msvdreduce target different norms (H∞ vs Hankel). For execution cost modelling, Hankel norm is more relevant since it bounds the energy of the response rather than the worst-case frequency gain.

## Key Papers

- [[papers/Kong2018-TIBInfoGeometry]] — TIB form, null basis reduction, fast Toeplitz algorithms
- [[papers/Mu2026-ModelReductionNotes]] — BT square-root algorithm derivation and sanity checks

## Open Questions

- Which reduction method (BT vs msvdreduce vs info SVD) is best for market impact IRs estimated from noisy data?
- How does capacity (AUM limit) interact with model order? Fewer poles = coarser impact model = potential systematic bias in cost estimates.
- Online re-estimation: how frequently should the reduced model be re-fit as market microstructure evolves?

## Related Pages

[[concepts/TIBForm]], [[concepts/MarketImpact]]
