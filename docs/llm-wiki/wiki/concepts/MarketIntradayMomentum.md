---
type: concept
tags: [momentum, short-term-reversal, multi-asset, market-impact]
sources: [Baltussen2020-IntradayMomentumHedging]
updated: 2026-06-03
---

# Market Intraday Momentum (Hedging-Demand Channel)

## Definition

Time-series momentum at the **market level within a trading day**: the **rest-of-day return `r_ROD`** (previous close → 30 min before today's close) positively predicts the **last-half-hour return `r_LH`**, and then **reverts over the next ~3 days**. Documented across 60+ futures (equities, bonds, commodities, FX), 1974–2020 ([[papers/Baltussen2020-IntradayMomentumHedging]]).

## Mechanism

**Gamma hedging demand.** Short-gamma agents (option market makers, leveraged/inverse ETFs, variance swaps, vol-target / risk-parity programs) must trade *in the direction of* price moves to stay delta-neutral; the bulk of this hedging is concentrated near the close (better liquidity, overnight-risk and capital/margin incentives). This predictable flow pushes `r_LH` in the direction of `r_ROD` (momentum), and the **transitory** price pressure reverts over subsequent days. The effect is stronger when negative gamma exposure (NGE) is larger, and its magnitude tracks leveraged-ETF hedging demand.

## Why It Matters / the system-identification connection

This is a **return-predictability** result, and a natural fit for **system identification as predictive filtering** ([[concepts/SysIDReturnPrediction]]): forecast `r_LH` from a history of observable intraday returns. The paper uses a coarse, hand-crafted predictor (equal-weighted `r_ROD`); a learned linear dynamic filter (TIB / kernel) can do better.

> **Not a propagator / impact model.** A market-impact *propagator* ([[concepts/MarketImpact]]) maps observed *order flow* → price (causal, mechanical). Here we forecast *returns from returns* (correlational), so the estimated kernel is a **predictive (Wiener) filter**, not an impact propagator. The two share only the LTI / impulse-response *mathematics* — the same TIB/kernel tooling fits both, but they are different subjects. The hedging story below is the *economic mechanism* that explains the predictability, not the model we estimate.

**How to build signals (suggestions):**

1. **Learn the intraday filter instead of equal-weighting `r_ROD`.** Inputs = the sequence of intraday half-hour returns (`r_ON, r_FH, …, r_SLH`); output = `r_LH`. Estimate the predictive impulse response `g` (weights over intraday lags) with **kernel/Bayesian regularisation** (stable-spline prior, [[papers/PillonettoDeNicolao2010-KernelSysID]]) — low SNR ⇒ shrinkage matters. `r_ROD` is the flat-`g` special case; test whether a learned smooth `g` beats it OOS. (Directly the [[concepts/SysIDReturnPrediction]] recipe: "raw lagged drivers, let the IR be the indicator.")
2. **Model momentum + reversal jointly as a multi-horizon predictive filter.** Output = forward returns over several horizons (last-half-hour through the next few days); input = the predictor history. A [[concepts/TIBForm|TIB]]/kernel filter whose poles encode the *timescales* of the intraday-momentum and the few-day-reversal predictability captures the paper's two findings in one model and yields a natural exit / holding rule. (This becomes an actual impact-*propagator* estimation only if you feed *observed hedging order flow* as the input — see [[concepts/MarketImpact]].)
3. **Regime-gate by gamma exposure (nonlinear combine).** Use the linear sysID intraday-momentum signal as the backbone; **gate/size it by the NGE regime** (momentum strong when NGE negative) via feature-fusion or meta-labeling ([[concepts/SysIDReturnPrediction]]). The mechanism supplies a *structural* regime variable — far less overfit-prone than a data-mined gate.
4. **Joint MIMO across asset classes.** The effect is "everywhere"; a MIMO [[concepts/TIBForm|TIB]] pools the estimate across the 60+ contracts (shared hedging dynamics, cross-asset spillovers), reducing variance at low SNR — the joint-fitting advantage over per-instrument fits.
5. **Trade the reversal leg / use it as risk control.** The propagator's reversal portion is itself a fade signal and flags that the momentum is transitory (size and hold accordingly).

**Why the mechanism helps the modelling:** a *causal* hedging-flow story (not pure data-mining) justifies both the predictive-filter prior (smooth, decaying weights that revert) and the conditioning variable (NGE) — materially reducing overfitting risk in an otherwise very low-SNR signal.

## Empirical Findings

- `r_ROD` OOS R² ≈ **2.88%** (equity index futures); strategy Sharpe **0.87–1.73** gross at asset-class level ([[papers/Baltussen2020-IntradayMomentumHedging]]).
- Reversal over ~3 days; momentum magnitude scales with |NGE| and LETF demand.
- Distinct from cross-sectional intraday seasonality (Heston et al. 2010).

## Open Questions / Caveats

- Net of close-auction execution costs, how much of the gross Sharpe survives?
- Does a learned intraday IR (kernel/Bayesian) beat equal-weight `r_ROD` OOS after costs?
- Data dependency: NGE / LETF gamma data (e.g. SqueezeMetrics) — do we have access?
- Best regime-gating: hard NGE split vs meta-labeling vs feature-fusion?

## Related Pages

[[papers/Baltussen2020-IntradayMomentumHedging]], [[concepts/MarketImpact]], [[concepts/TIBForm]], [[concepts/SysIDReturnPrediction]], [[concepts/PyMORDataDrivenID]]
