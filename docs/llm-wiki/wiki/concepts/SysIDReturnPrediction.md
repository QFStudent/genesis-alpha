---
type: concept
tags: [factor-model, backtest-methodology, multi-asset]
sources: []
updated: 2026-06-02
---

# Applying System Identification (TIB / Kernel) to Return Prediction

## Scope

How to use LTI system-identification — [[TIBForm|TIB]] and kernel/Bayesian ID ([[PyMORDataDrivenID]]) — for financial **return prediction**, and how to fuse it with nonlinear ML (trees / deep learning). Two questions: (1) what should the **inputs and outputs** be? (2) how to **combine** the linear sysID with a nonlinear model?

The fact these methods encode: they fit `y(t) = (g * u)(t) + noise` — output = convolution of the input with an impulse response `g` (the memory kernel). **The method's whole job is to learn the lag structure**, so the input should be the *raw causal series*, not pre-filtered features.

## Inputs & Outputs

**Output = forward returns** (yes). Refinements:
- multi-horizon → MIMO vector output = a *term structure of predictability*;
- vol-normalise the target (LTI assumes ~homoskedastic noise; returns have vol clustering).

**Inputs = raw, causal, lagged driver series — NOT TA indicators.** Why:
1. **SysID *is* the feature engineering.** MA / EMA / MACD are *linear filters* of past returns — exactly what the IR `g` represents. The learned IR reproduces any linear indicator as a special case and finds the *optimal* one. Feeding linear TA indicators is redundant and imposes a fixed, likely-suboptimal shape.
2. **Two paradigms.** TA feature-stacking = cross-sectional ML (many features, shallow lag). SysID = few raw channels, rich learned dynamics. Pre-baking lagged/filtered features front-loads the temporal modelling sysID is built to do.
3. **Recommended raw inputs (MIMO channels):**
   - lagged returns of the target → IR = optimal momentum / short-term-reversal filter;
   - **order flow / OFI / signed volume** → usually better than past returns (causal, higher SNR);
   - related-asset / factor returns → lead-lag & cross-impact (the MIMO TIB null vectors are the input directions across instruments).

**The linear/nonlinear boundary (honest limit).** LTI captures linear, time-invariant predictability only; nonlinear TA value (regimes, thresholds, ratios, interactions) is *unreachable* by LTI.
- Linear TA indicators (MA/EMA/MACD) → don't use as inputs (subsumed by the IR).
- Nonlinear features → either add as input channels (linear *in* the feature) or move to a nonlinear model (next section).

**Interpretation caveat:** for prediction the IR is a regularised **Wiener / predictive filter**; poles = predictability *timescales*, not physical dynamics — correlational, not generative. Returns are extreme low-SNR → exactly where the kernel / [[PyMORDataDrivenID|Bayesian-TIB]] decaying-lag prior pays off.

## Combining with a nonlinear model

If nonlinear TA features go to a tree/DL model, **fuse — don't naively chain.** Division of labour: **linear sysID = dynamic linear alpha** (`ŷ_lin` + TIB latent states `x(t)`); **nonlinear model = regimes / interactions.**

| approach | what it does | verdict |
|---|---|---|
| **A. Residual / boosting-offset** — fit nonlinear on `r = y − ŷ_lin` (or GBM `base_margin`) | nonlinear mops up what linear missed | theoretically clean (residual = innovation), but **fragile at low SNR** (residual ≈ noise → overfit) and bakes in an ordering assumption. *Not the default.* |
| **B. Stacking** — regularised meta-learner on out-of-fold preds | learns the optimal blend | **safe robust default**; no ordering assumption; handles correlated learners |
| **C. Feature fusion** — feed `ŷ_lin` + TIB states `x(t)` into the nonlinear model | nonlinear model **gates** the linear signal by regime | **most powerful when nonlinear value is conditional**; TIB states are a compact, denoised summary of history — better tree/NN food than raw lags |
| **D. Meta-labeling** — linear = side; nonlinear = P(correct) → sizing | nonlinear *sizes*, doesn't re-predict | finance-native, robust; easier (higher-SNR) learning target |

**Recommendation:** keep TIB/kernel as a **stable linear backbone**; default to **stacking (B)**; prefer **feature-fusion (C)** or **meta-labeling (D)** when the nonlinear contribution is conditional / sizing; **avoid pure residual-chasing (A)**. Layering is fine (e.g. meta-labeling for sizing + a small stacked additive term, non-negative weights).

> **Don't default to "nonlinear predicts the residual."** At return-prediction SNR, regime-gating and regularised blending beat residual-chasing.

## Methodology guards (matter more than the architecture)

- **Leakage from overlapping forward returns:** any meta-learner / stack trains on out-of-fold predictions with an **embargo ≥ the forecast horizon** (purged / combinatorial-purged CV).
- **Correlation / double-counting:** TA features are functions of the same returns the linear model uses → correlated base preds → **regularise the combiner** (ridge / non-negative / convex weights).
- **Combine at the signal level, size once** — one blended forecast → one position.
- **Validate the *combination* OOS** after costs, not the components in-sample.

## Why It Matters for genesis-alpha

Defines how the [[TIBForm|TIB]] / [[PyMORDataDrivenID|kernel]] machinery plugs into alpha research: raw causal inputs → linear dynamic forecast (+ interpretable mode states) → fused with nonlinear ML for regimes. The TIB latent states double as a denoised feature representation for downstream models.

## Open Questions

- Which raw driver carries the linear signal at our horizon — lagged returns vs order flow vs cross-asset?
- Does feature-fusion (TIB states into a GBM) beat stacking and meta-labeling OOS, after costs?
- Best purged-CV / embargo scheme for our horizons to keep the hybrid honest?
- Are TIB latent states better tree/NN features than raw lags or standard TA indicators?

## Related Pages

[[TIBForm]], [[PyMORDataDrivenID]], [[papers/PillonettoDeNicolao2010-KernelSysID]], [[MarketIntradayMomentum]], [[OHLCVPooledPrediction]], [[ModelReduction]], [[MarketImpact]]
