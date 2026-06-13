---
type: concept
tags: [futures, multi-asset, factor-model, data]
sources: []
updated: 2026-06-13
---

# Energy-Based Covariate Selection (MIMO sysID)

## Definition

Principled selection of inputs for the MIMO TIB / LTI prediction model, based on
**energy** (squared ℓ²-norm / variance) — but the *right* energy.

**The weak criterion (what to avoid).** Ranking inputs by in-sample impulse-response
channel energy `Σ_τ h_{ij}(τ)²` (the [[concepts/FuturesCovariates]] channel heatmap) is
the weakest measure: it (a) conflates predictable signal with **overfit noise**, (b)
ignores **redundancy / collinearity** — two correlated inputs each look high-energy but
together add little, and (c) measures raw *drive*, not *predictive value*.

**The principle.** Rank inputs by the **predictable, incremental, out-of-sample** output
energy (variance) they contribute — not by raw IR magnitude.

## Methods (most-principled-for-LTI first)

1. **CCA / subspace ID / reduced-rank regression — the natural MIMO energy method.**
   Compute the **canonical correlations between the past** (stacked lagged inputs *and*
   outputs) **and the future outputs**, after pre-whitening. The squared canonical
   correlations are a **predictable-energy spectrum** (fraction of future-output variance
   explainable from the past per direction). Then:
   - the number of *significant* canonical correlations = the **effective state order**
     (the data-domain version of the Hankel-σ decay, see [[concepts/ModesAndHankel]]);
   - the canonical **loadings** say which inputs carry predictable energy — an input that
     doesn't load on a significant direction is dropped.
   This is **subspace identification (CCA/CVA, N4SID)** — the redundancy-aware, predictive
   generalization of Hankel-SVD; it yields covariate selection *and* model order from one
   spectrum. Ties to [[concepts/PyMORDataDrivenID]]. **Reduced-rank regression** is the
   same idea in regression form (keep the top-`r` predictable directions).

2. **Hankel-σ contribution (controllability × observability), not raw drive.** An input
   can pour energy into a state that produces no output. Rank inputs by how they load on
   the **dominant Hankel directions** (high-σ right singular vectors), which combine
   reachability and observability. The per-input controllability Grammian
   `tr(P_j)`, `P_j = Σ_k A^k B_jB_j^⊤A^{k⊤}`, measures excitation energy but is only
   meaningful *weighted by observability* — the Hankel singular value does both.

3. **Group, out-of-sample, incremental variance explained.** Drop input `j` if including
   its *entire IR column* doesn't reduce OOS residual output energy — i.e. **group
   partial-R² under cross-validation**. Being *incremental* handles redundancy
   automatically; *group* (whole filter, not per-lag) respects that each input is one
   channel.

4. **Group-LASSO + reduced-rank (the practical workhorse).** Penalize each input's IR-block
   norm (group sparsity → drops whole irrelevant inputs) plus a **nuclear-norm /
   reduced-rank** penalty (few-mode structure → the low-energy-tail truncation). Tune by
   CV. Operationalizes 1–3: selects inputs *and* order by penalizing energy and validating
   out-of-sample.

5. **Information distance / transfer entropy (info-geometry weighting).** Replace L²
   energy with **information distance** (cepstrum, [[concepts/InformationGeometry]]) or
   mutual/transfer entropy between lagged input `j` and the output — weights *relative*
   predictive content, better for non-Gaussian / low-SNR data. More principled, more work.

## Recommended workflow

1. **CCA / subspace predictable-energy spectrum** (past inputs+outputs → future outputs):
   read off the effective order *and* the per-input loadings → a ranked covariate list.
2. **CV group partial-R²** to confirm and catch redundancy the spectrum hides.
3. **Group-LASSO + reduced-rank regression**, CV-tuned, as the production selector.

All steps **out-of-sample** — in-sample energy overstates importance (overfitting), exactly
the failure mode behind the ill-conditioned high-order fits (see [[concepts/ModesAndHankel]],
the energy floor / identifiability limit).

## Why It Matters for genesis-alpha

This is the principled upgrade to the candidate-input list in
[[concepts/FuturesCovariates]]: replace ad-hoc channel energy with **predictable energy**,
get the **model order and the covariate set together**, and be **redundancy- and
OOS-aware** — directly feeding the pooled MIMO construction of
[[concepts/OHLCVPooledPrediction]] (shared pole bank, per-channel residue).

## Open Questions

- Which canonical/predictable-energy directions survive **net of cost** intraday (the fast
  ones may not, per the turnover/√-breadth analysis)?
- Best whitening / vol-normalization for the CCA step so the LTI stationarity assumption
  holds (links the vol-normalize/de-seasonalize cautions in [[concepts/FuturesCovariates]])?
- L² predictable energy vs information-distance selection — when do they disagree, and which
  generalizes better at low SNR?

## Related Pages

[[concepts/FuturesCovariates]], [[concepts/ModesAndHankel]], [[concepts/PyMORDataDrivenID]], [[concepts/InformationGeometry]], [[concepts/OHLCVPooledPrediction]], [[concepts/SysIDReturnPrediction]]
