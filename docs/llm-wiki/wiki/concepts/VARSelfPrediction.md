---
type: concept
tags: [multi-asset, factor-model, futures]
sources: [Mu2026-TIBNote]
updated: 2026-06-14
---

# (V)AR Self-Prediction — the TIB state-space as an inverse / whitening filter

## Definition

The genesis-alpha TIB state-space `x_{t+1} = A x_t + B u_t`, `ŷ_t = C x_t` is used in
**self-prediction** mode: the input `u_t` **and** the output are the *same* data — the
`p`-vector of (lagged) returns `y_t`. That makes the model a **vector autoregression (VAR)**:

$$\hat y_t = \sum_{k\ge1} (C A^{k-1} B)\, y_{t-k}, \qquad A_k := C A^{k-1} B \in \mathbb{R}^{p\times p}.$$

So **`CAⁿB` are the (matrix) AR coefficients** — square `p×p` — not the forward impulse
response. The model realizes the **prediction** polynomial `Â(z) = Σₖ A_k z⁻ᵏ`
(data → prediction `ŷ`); the associated **whitening / prediction-error** filter is
`𝒜(z) = I − Â(z)` (data → residual `ε`), and the **forward / synthesis** filter is its
inverse `𝒜(z)⁻¹` (innovation → data):

| filter | input → output | transfer | IR coefficients |
|---|---|---|---|
| forward / synthesis | innovation `ε` → data `y` | `𝒜(z)⁻¹` | MA(∞) coeffs (the impulse response) |
| **prediction** *(this model)* | data `y` → prediction `ŷ` | `Â(z) = Σₖ Aₖ z⁻ᵏ` | `[A₁, A₂, …]` = **AR coeffs** = `CAⁿB` |
| whitening / prediction-error | data `y` → residual `ε = y − ŷ` | `𝒜(z) = I − Â(z)` | `[I, −A₁, −A₂, …]` |

These are **two different filters**, related by `I − (·)`, not one. Your state space realizes
the **prediction** filter `Â(z)` — output `ŷ` (≈ the data), IR `= CAⁿB =` the AR coefficients
`[A₁, A₂, …]`. The **whitening / prediction-error** filter is `𝒜(z) = I − Â(z)` — output the
residual `ε = y − ŷ`, IR `[I, −A₁, −A₂, …]`. `ar2ir` inverts the *whitening* polynomial: it
takes the AR coefficients `Aₖ` (= your `CAⁿB`), forms `𝒜(z) = I − Σₖ Aₖ z⁻ᵏ`, and returns the
forward IR `𝒜(z)⁻¹`.

## Why it matters for genesis-alpha

- It's the natural form for predicting returns from their own (and cross-asset) lags:
  `ŷ_t = Σ_k A_k y_{t-k}` with shared poles across the `p×p` blocks.
- The TIB parameterization (poles + null vectors) **is** the balanced MIMO lattice / Schur
  form of this inverse filter, so the whole reduction stack (Hankel-SVD, BT, `msvdreduce`)
  operates directly on the AR-coefficient sequence `CAⁿB` — see [[concepts/TIBForm]],
  [[concepts/ModesAndHankel]].
- To go from the fitted **inverse** filter to a **forward impulse response** (for simulation,
  IRF analysis, or signal construction) you expand `𝒜(z)⁻¹` — the `ar2ir` step below.

## `ar2ir`: AR coefficients → forward IR

`ar2ir` (SISO, `dynamical_system/filters/tib.py`) maps scalar AR coefficients to the IR of
`1/A(z)` by solving a lower-triangular Toeplitz system — equivalently the recursion
`h₀ = 1`, `hₙ = Σ_{k≥1} aₖ hₙ₋ₖ`. The **MIMO** analog is the VAR→VMA expansion of `𝒜(z)⁻¹`:

$$h_0 = I, \qquad h_n = \sum_{k=1}^{n} A_k\, h_{n-k}\quad (A_k = 0\text{ for }k>p),$$

each `hₙ` a `p×p` matrix — a **block lower-triangular Toeplitz** solve (or block forward
substitution). It reduces to the scalar `ar2ir` when `p = 1`. The scalar derivation
(`T_a h = e₁`, the AR↔IR Toeplitz duality `T_h T_a = I`) is in [[papers/Mu2026-TIBNote]],
which marks the MIMO case "pass"; this block recursion fills it (`ga/filters/tib.py:
mimo_ar2ir`). **Full multivariate derivation** (block-Toeplitz duality `T_H T_𝒜 = I`,
`H = 𝒜⁻¹`, left=right recursion, validation): `docs/derivations/07-mimo-ar-to-ir.md`. Pipeline:

```
TIB (A,B,C)  →  A_k = C A^{k-1} B  (square AR coeffs)  →  mimo_ar2ir  →  forward IR  h_n
```

## Reducing for the data poles — and predicting from them

> **Data poles (definition).** The poles of the *data-generating* process — the eigenvalues
> that govern the data's own dynamics: each mode's oscillation frequency (`arg λ`) and
> decay / persistence (`|λ|`). Formally, the poles of the forward filter `𝒜(z)⁻¹` (the roots
> of `det 𝒜(z) = 0`), recovered as `eig(A_r)` after `ar2ir` + reduction. They are the modes
> *of the data itself* — distinct from the prediction filter `Â(z)`'s poles, which sit at 0
> for a finite AR (those encode only how fast the AR coefficients fall off, not the data's
> dynamics). The oscillatory effects you want to model live in the **data poles**.

**Model reduction (ERA / Hankel-SVD / BT) extracts the poles of whatever IR you feed it**,
so the order matters. The AR coefficients `{A_k} = CAⁿB` are the IR of the *prediction*
filter `Â(z)`; feeding them to reduction recovers `Â`'s poles — which for a finite AR(p) sit
at **0** (`Â(z)` is a polynomial in `z⁻¹`). The data's oscillatory poles live in the
*forward* filter `𝒜(z)⁻¹`, so you must **`ar2ir` first, then reduce**:

```
fit → AR coeffs {A_k}=CAⁿB → ar2ir → forward IR {h_n} → reduce → data poles (= poles of 𝒜⁻¹)
```

Illustration — AR(1) `y_t = a·y_{t-1} + ε`: the AR coeff `{a}` is FIR ("pole" 0); the forward
IR `{1, a, a², …}` has pole `a`. Reduce the AR coefficients directly and you find a pole at 0,
not `a` — for a finite AR you'd get nothing but poles at the origin.

**Predicting from the extracted data poles.** Reduction returns the reduced state-space
`(A_r, B_r, C_r)` directly, with `eig(A_r)` = the data poles. You predict by **propagating
this reduced state-space recursively**, driving it with the observed returns:

$$x_{t+1} = A_r x_t + B_r u_t, \qquad \hat y_t = C_r x_t .$$

**No separate inversion step is needed** — running the reduced recursion *is* the inverse /
prediction filter, with the data poles carried in `A_r`. (Writing out an explicit inverse
transfer function would only re-express this same recursion.) You stay in the original
*data-in → prediction-out* form, now low-order, with the data's dominant dynamics in `A_r`.
Full round-trip:

```
fit → AR coeffs CAⁿB → ar2ir → forward IR {h_n} → reduce → data poles + reduced (A_r, B_r, C_r)
   → propagate the reduced state-space on the observed returns → ŷ
```

## Connections (is this a known thing?)

Yes — it is **linear prediction (LPC)**: `𝒜(z)` is the *prediction-error / whitening filter*,
`𝒜(z)⁻¹` the *synthesis filter*. The same object is the **lattice filter / Levinson–Durbin /
Schur** recursion; **TIB is its balanced MIMO parameterization** (tangential Schur = the
multivariate Levinson recursion). So "a state space driven by the data whose IR is the AR
coefficients" is the natural state-space form of the inverse filter the Schur/TIB machinery
was built for — not an ad-hoc construction.

## Caveats / open questions

- **"AR" is exact only in self-prediction** (input = the modeled series' own past). With
  *exogenous* features (`p×q`, e.g. ES→NQ) it becomes a transfer-function / ARX model;
  `CAⁿB` is then a cross-IR, not AR coefficients, `𝒜(z)` isn't square, and `ar2ir` doesn't
  apply (see [[concepts/SysIDReturnPrediction]] for the I/O design question).
- **Whiteness is a property of the *fit*, not the structure.** A partial whitener (residual
  not perfectly white) is still a valid predictor — "it works" without exact innovation
  recovery.
- **Minimum-phase requirement.** `𝒜(z)⁻¹` decays into a usable IR only if `𝒜(z)` is
  minimum-phase (forward poles inside the unit circle). TIB enforces this by construction;
  a drifted fit that goes non-minimum-phase makes the IR blow up — a useful diagnostic.

## Related Pages

[[concepts/TIBForm]], [[concepts/ModesAndHankel]], [[concepts/SysIDReturnPrediction]],
[[concepts/FuturesCovariates]]
