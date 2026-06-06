---
type: concept
tags: [factor-model, multi-asset, backtest-methodology]
sources: []
updated: 2026-06-05
---

# Own-asset OHLCV inputs & pooled MIMO prediction

## Scope

Practical construction for predicting many instruments' (execution-price) returns from **each instrument's own OHLCV-derived features**, using a **pooled** MIMO so the estimate stays robust at low SNR. Two parts: (1) what to feed in, (2) how to build the pooled multivariate model.

## Inputs: own-asset OHLCV — use returns, range, log-volume

OHLCV carries more than the close-to-close return; the useful, *distinct* channels:
- `Close − Open` → intrabar return (momentum/reversal),
- `High − Low` (range) → a **volatility** proxy,
- close position within the range → pressure,
- `Volume` → activity / conviction.

Checklist:
- **Stationarity:** use **returns, ranges, log-volume** — *not* raw O/H/L/C price *levels* (convolving non-stationary levels gives spurious fits).
- **Look-ahead (the #1 OHLCV trap):** High/Low/Close are realized only at the *end* of the bar → inputs must be **strictly prior** to the predicted-return window. Predict bar `t+1` (or the last-half-hour) from bars `≤ t` only.
- **Collinearity:** O/H/L/C are highly correlated → near-collinear channels → the fitted weights are unstable → **regularize `C`**.
- **Linear only:** an LTI/TIB filter captures the *linear* own-asset structure (reversal/momentum, vol-timing, volume-weighting); nonlinear OHLCV patterns (candles, interactions, regimes) need a nonlinear model — combine via [[concepts/SysIDReturnPrediction]] (stacking / feature-fusion / meta-labeling).

This is the **own-asset (diagonal)** signal — the robust place to start before any cross-asset/lead-lag (see [[concepts/MarketIntradayMomentum]], *Own-asset vs cross-asset*).

## Pooled MIMO construction (tied diagonal blocks)

The robust multivariate model = **own-asset (block-diagonal) + pooled (tied blocks)**.

- **Block-diagonal `A`, `B`, `C`** → each instrument is a decoupled own-asset subsystem (transfer matrix `H(z)` block-diagonal: output `i` depends only on instrument `i`'s features).
- **Tied** → all blocks are the *same matrices* (same poles, same null vectors, same `C`) ⇒ one shared filter ⇒ **pooling**.
- **Shared pole bank, replicated — not circulated.** Every instrument uses the *same* `m` timescales; copy the full `m`-pole bank into each block. Round-robin **canonical** null vectors instead give each instrument a *different* pole subset → incommensurable bases → **breaks pooling**.
- **Don't feed repeated poles into one `null_basis_realization` call** (equal poles make the Blaschke deflation singular). Build one block and reuse it.

**Recipe (you never build the big block-diagonal matrix):**
1. Build one shared block once: `A_b, B_b = null_basis_realization(poles_m, V_b)` — `poles_m` = `m` *distinct* timescales, `V_b` = (features × `m`) null vectors.
2. **Filter each instrument's features through that one block** → its state sequence `x_i(t)`.
3. **Stack** `{(x_i(t), y_i(t))}` across instruments and fit **one** `C` by regularized least squares. That *is* the tied-diagonal-block model — and equals the **pooled MISO** (one shared block per instrument, `C` on the concatenated panel).

**Partial pooling** (shrink each instrument's `C` toward the shared one) is usually the sweet spot — borrow strength but allow heterogeneity.

## Pooled vs dense — why pool

- **Pooled** = the full MIMO restricted to {block-diagonal, tied}; estimates `m` shared weights from `(#instruments) × N` rows → strong variance reduction (≈ `√#instruments` ×).
- **Dense** = each output reads *all* instruments' states with its *own* full `C` (no tying, `#inst·m` weights/output) → overfits at return SNR.
- Synthetic demo `notebooks/pooled_vs_dense_mimo.ipynb`: with a diagonal+tied truth, pooled ≈ population R² while dense is **negative OOS** (catastrophic at small N); and **even when real cross-asset signal exists, dense still loses to pooled**. ⇒ harvest cross-asset signal with a *small regularized / low-rank factor block*, not a dense fit.

## What you regularize

With fixed poles + null vectors, `A`, `B` are fixed and **only `C` is fit** — so `C` is the regularization target. Two structured `C`-priors: **smooth over modes** ([[concepts/PyMORDataDrivenID|Bayesian-TIB]], within-instrument) and **tie across instruments** (pooling).

## Related Pages

[[concepts/SysIDReturnPrediction]], [[concepts/MarketIntradayMomentum]], [[concepts/TIBForm]], [[concepts/PyMORDataDrivenID]]
