---
type: concept
tags: [factor-model, execution, market-impact]
sources: [Kong2018-TIBInfoGeometry]
updated: 2026-06-07
---

# Information Geometry of LTI Systems

## Definition

Information geometry applies the framework of Riemannian differential geometry to probability distributions. For LTI systems, the key idea is:

- Each stable LTI transfer function H(z) defines a probability distribution over output sequences (via its spectral density).
- The **Fisher information metric** on this space gives a natural notion of distance between systems.
- The **information distance** I(f, g) is the geodesic length between two systems f and g on this Riemannian manifold.

The key result ([[papers/Kong2018-TIBInfoGeometry]], Ch. 3) is that when **cepstrum coefficients** {aₖ} are used as coordinates:

```
log H(z) = Σ aₖ zᵏ
```

the Fisher information matrix becomes the **identity** — i.e. the manifold is **Euclidean** in cepstrum coordinates. This makes distance computation tractable:

```
I²(H, white noise) = Σ |aₖ|²  =  ‖log H(z)‖²_H₂
```

## Connection to LTI Distance Measures

Information distance is sandwiched between Hellinger and total variation distance (Lemma 3.1.1–3.1.2 in Kong2018):

```
8 H(f, g) ≤ I(f, g)    and    K·I(f, g) ≤ 8·H(f, g)   for K < 1
```

This means minimising information distance is equivalent (up to constants) to minimising the probability of confusing two systems in a hypothesis test — a statistically principled objective for system approximation.

For MIMO systems with n inputs/outputs, the per-sample information distance generalises to:

```
I²(Tm, I) → (1/4n) ∫ ‖log Φ(eⁱᶿ)‖²_F dθ
```

where Φ is the spectral density matrix (block Toeplitz limit via the generalised Szegő theorem).

## Why It Matters for genesis-alpha

1. **Reduction objective**: `info_svd_reduce` in `ga/reducers/hankel.py` minimises information distance rather than H∞ or H₂. For market impact IRs that are non-rational or near-unit-root, this is more appropriate than H∞.
2. **Non-stationary robustness**: Information geometry handles non-stationary processes (random walk) gracefully — finite information distance can exist even without a spectral density.
3. **Cepstrum coordinates**: The Euclidean structure means cepstrum-domain model selection (e.g. comparing two impact models) is as simple as computing a squared L₂ norm. No geodesic integration needed.
4. **TIB connection**: The TIB form is information-geometrically motivated — it is the parameterisation in which the Fisher information matrix has the simplest structure, making it optimal for prior estimation and dimensionality reduction on the statistical manifold. See [[concepts/TIBForm]].

> ⚠️ **Implementation gap (found 2026-06-07):** the claim in #1 that `info_svd_reduce` "minimises information distance rather than H∞ or H₂" is **not realized in the current code**. As written the function (a) runs the pole-finder directly on the cepstrum `aₖ`, but `aₖ = (1/k)Σλᵏ` is *not* a low-rank-Hankel sequence so Ho–Kalman cannot recover poles from it; (b) computes a `rho`-damped IR but then ignores it; and (c) shares the `msvdreduce` transpose-vs-pinv shift bug. A corrected companion `info_svd_reduce_fixed` (undo the `1/k` → power sums `Σλᵏ`, then shift) recovers poles correctly **but only for identification / full order** — the power-sum Hankel has no energy ordering, so it does **not** give numerically stable order *reduction*. So the figures below are Kong (2018)'s results *in principle*; the repo does not yet implement genuine information-distance reduction. See `docs/derivations/01-hankel-svd-reduction.md` (SISO implementation-note) and [[concepts/ModelReduction]].

## Empirical Findings

- Cepstrum reduction achieves H₂ norm competitive with balanced truncation on rational systems ([[papers/Kong2018-TIBInfoGeometry]], Figs 12–15).
- For infinite-dimensional systems (non-rational IRs), cepstrum/info-SVD reduction outperforms BT in Hankel norm (Figs 14–17).
- Fast block-Toeplitz logarithm (O(n log²n)) enables MIMO cepstrum computation at scale — crucial for multi-asset impact modelling.

## Open Questions

- Is the information distance objective the right one for market impact? The connection to hypothesis testing says "minimise the chance of confusing two impact models" — is that the same as minimising trading cost error?
- How does cepstrum reduction interact with the TIB null-basis structure? The `info_svd_reduce` function currently returns poles only, not null vectors.
- **Implement genuine information-distance order reduction.** The current `info_svd_reduce` path is at best a full-order *identifier* (see ⚠️ above), not a reducer. A real info-distance reducer (Kong2018) needs a different algorithm than Ho–Kalman-on-cepstrum.

## Related Pages

[[concepts/TIBForm]], [[concepts/ModelReduction]], [[papers/Kong2018-TIBInfoGeometry]]
