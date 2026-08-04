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

## Practical case: ARX / MISO prediction (own-lags + exogenous inputs)

The theory above is the **square self-prediction VAR** (`CAⁿB` = AR coefficients, inverted with
`mimo_ar2ir`). A typical *prediction* model is narrower and needs much less machinery: a
**univariate target** from **its own lagged returns/volume plus other instruments' lagged
returns**. That is an **ARX** model (autoregressive + exogenous):

$$y_t = \underbrace{\sum_k a_k\, y_{t-k}}_{\text{own AR (scalar)}}
      + \underbrace{\sum_k b_k\, u_{t-k}}_{\text{exogenous: volume, other instruments}}
      + \varepsilon_t,
  \qquad H(z) = \frac{B(z)}{\mathcal A(z)},\ \ \mathcal A(z) = 1 - \textstyle\sum_k a_k z^{-k}.$$

**For prediction — no inversion, no innovations.** The one-step predictor is the conditional
expectation, and the innovation **drops out** of it (it is uncorrelated with the past, so
`E[ε_t | past] = 0`):

$$\hat y_t = \mathbb{E}[y_t\mid\text{past}] = \sum_k a_k y_{t-k} + \sum_k b_k u_{t-k}
   + \underbrace{\mathbb{E}[\varepsilon_t\mid\text{past}]}_{=\,0}
   = \sum_k a_k y_{t-k} + \sum_k b_k u_{t-k}.$$

Every right-hand-side term is **observed** (past `y`, past `u`), so you predict by plugging in
observed values — the fitted ARX regression **is** the optimal linear predictor. The innovation
`ε_t = y_t − ŷ_t` is the prediction *error*: a byproduct known only after the fact, **never an
input**. Own-lags don't reintroduce innovations because they're *observed data* — autoregression
runs the **inverse direction** (data in → residual out), so you feed data, not innovations. And
"no inversion" because you use the AR/ARX coefficients *directly*; you never form `1/𝒜(z)`.

> **Two caveats — so "no inversion, no innovations" is precise, not over-broad:**
> 1. **This is pure ARX.** With an **MA error term** (ARMAX, `… + Σ cⱼ ε_{t-j} + ε_t`), the
>    predictor *does* use the **past innovations** `ε_{t-j}`, which you compute recursively as
>    residuals — `E[ε_t|past]=0` still holds, but the lagged `ε_{t-j}` are no longer zero in the
>    prediction. So "no innovations" holds for **ARX, not ARMAX**. Your model (own lags +
>    exogenous, no MA term) is ARX.
> 2. **"No inversion" is about *prediction*.** If you also want the **forward IR to reduce**, you
>    *do* invert — but only the **scalar** `𝒜(z)` (own-AR) via the scalar `ar2ir`, then convolve
>    with the exogenous numerator `bₖ`. That inversion is real and scalar, but **separate from
>    prediction** (and never `mimo_ar2ir`).

**The autoregression is scalar.** The only AR / invertible part is the target's *own* past:
- **Data poles** = roots of the *scalar* `𝒜(z)` (own-AR) — read them off directly (or `eig` of
  the scalar AR companion);
- **Exogenous terms** (volume, other instruments) are the **numerator** `B(z)` — they add
  **zeros**, observed inputs, *not poles*; no innovations, no inversion.

So **`mimo_ar2ir` is not what an ARX model needs** — its multivariate (square `d×d`) AR
structure isn't present. At most you touch the **scalar `ar2ir`**, and only to get the explicit
input→target IR for reduction:
`IR = bₖ * ar2ir(aₖ)` — scalar expansion of `1/𝒜(z)`, convolved with the exogenous numerator.

| you want… | innovations? | `ar2ir`? |
|---|---|---|
| **predict** the target | no | no — regress on observed lags |
| **data poles** | no | root-find the scalar `𝒜(z)` |
| explicit **forward IR** (to reduce) | no | **scalar** `ar2ir` on `aₖ`, convolve with `bₖ` |
| `mimo_ar2ir` / multivariate inversion | — | **never** for ARX (own-AR is scalar) |

**"But a forward system needs innovations as inputs?"** — only the *self*-MA representation
(innovation → data). A model predicting the target from *other observed series* is a forward
**input→output** map driven by **observed inputs**, not innovations — so you skip both the
innovation computation and `mimo_ar2ir`. Computing a univariate target's innovations would mean
whitening it (`ε_t = y_t − Σ aₖ y_{t-k}`), which *is* the inverse filter you're trying to avoid.
For returns (already nearly white) the inputs ≈ their own innovations anyway — feed raw lagged
returns. And model **reduction works on a MISO IR (`p=1`)** exactly as on MIMO (the Hankel/SVD
machinery is dimension-agnostic — verified).

**When you DO need `mimo_ar2ir` / the full `d×d` VAR:** only if (a) you model all instruments
jointly as a *square* self-prediction VAR (the inverse-filter framing of this page), or (b) you
want the **structural / feedback-aware IRF** (response to an isolated innovation shock, tracing
feedback). Plain prediction needs neither — the ARX regression already conditions on the
observed inputs, with feedback absorbed implicitly into the coefficients.

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
