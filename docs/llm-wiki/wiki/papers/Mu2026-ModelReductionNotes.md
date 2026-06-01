---
type: paper
tags: [market-impact, execution]
sources: [Mu2026-ModelReductionNotes]
updated: 2026-05-31
---

# Notes for Model Reduction

**Authors:** Mu, Y.  **Venue / Year:** Internal Notes, Jan 2026

## Contribution

Working notes deriving the **Square-Root (SR) Balanced Truncation** algorithm from first principles, with full proofs of biorthogonality and balanced Grammians. Written as the direct theoretical basis for `ga/reducers/balanced_truncation.py`.

## Methodology

Starting from the Hankel singular value identity σᵢ(H)² = λᵢ(Wo Wc), the notes derive:

**Square-Root Algorithm (Algorithm 1):**
1. Cholesky factorise Grammians: P = Zₚ Zₚ*, Q = Z_Q Z_Q*
2. SVD of bridge matrix: Zₚ* Z_Q = U Σ V* (block-partitioned into kept/truncated)
3. Define transformations: W = Zₚ U₁ Σ₁^{-1/2}, V = Z_Q V₁ Σ₁^{-1/2}
4. Reduced system: Aᵣ = W*AV, Bᵣ = W*B, Cᵣ = CV

**Error bound:** ‖Gₙ - Gᵣ‖_{H∞} ≤ 2(σᵣ₊₁ + ⋯ + σₙ)

**Sanity checks (from notes):**
1. U₁*U₁ ≈ I
2. V₁*V₁ ≈ I (tolerance 10⁻¹⁰ to 10⁻⁶)
3. ‖Zₚ* Z_Q V₁ − U₁ Σ₁‖ (SVD block identity check)
4. W*V ≈ Iᵣ (biorthogonality)
5. W*QW ≈ Σ₁, V*PV ≈ Σ₁ (balanced Grammians)
6. Condition numbers of W, V

## Key Results

- Full proof that W*V = Iᵣ using SVD block identity Zₚ* Z_Q V₁ = U₁ Σ₁.
- Full proof that W*QW = Σ₁ and V*PV = Σ₁.
- Notes also include an interpretability section on TIB states: each state corresponds to a specific temporal decay rate and input direction — enables labelling states as "long-memory" vs. "fast transient".

## Connection to genesis-alpha

This is the direct specification document for `ga/reducers/balanced_truncation.py`:
- `balanced_truncation()` implements Algorithm 1 step by step
- `BalancedReductionResult.error_bound()` implements the H∞ error bound
- `sanity_check_bt()` implements all 6 sanity checks
- `compute_grammians()`, `cholesky_with_fallback()` implement steps 1–2

## Related Pages

[[concepts/ModelReduction]], [[concepts/MarketImpact]], [[papers/Kong2018-TIBInfoGeometry]]
