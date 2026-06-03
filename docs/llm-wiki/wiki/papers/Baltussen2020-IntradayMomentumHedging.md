---
type: paper
tags: [momentum, short-term-reversal, multi-asset, futures, market-impact]
sources: [Baltussen2020-IntradayMomentumHedging]
updated: 2026-06-03
---

# Hedging Demand and Market Intraday Momentum

**Authors:** Baltussen, G., Da, Z., Lammers, S., Martens, M.  **Venue / Year:** SSRN working paper, Aug 2020 (later Journal of Financial Economics, 2021). SSRN abstract 3760365.

## Contribution

Documents "market intraday momentum" **everywhere** — across 60+ futures (equities, bonds, commodities, currencies), 1974–2020: the **rest-of-day return predicts the last-half-hour return**. The novel contribution is linking this to **gamma hedging demand** (short-gamma dealers / leveraged ETFs hedging into the close), evidenced by a negative-gamma-exposure proxy and LETF rebalancing demand; the effect **reverts over the next days**, indicating transitory price pressure rather than information.

## Methodology

- **Trading-day partition:** ON (overnight), FH (first 30 min), M (middle), SLH (second-to-last 30 min), LH (last 30 min). `ONFH = ON+FH`; **`ROD = ON+FH+M+SLH`** (rest of day).
- **Predictive regressions** of `r_LH` on: (i) `r_ONFH` (the Gao et al. 2018 predictor); (ii) `{r_ONFH, r_M, r_SLH}`; (iii) **`r_ROD`** (eqs 5–7). OOS R² via expanding window (≥500 obs), Clark–West test; pooled per asset class and per contract.
- **Mechanism tests:** negative gamma exposure (NGE) proxy from S&P 500 options; LETF hedging demand (market cap × leverage); 3-day reversal; a "futures-after-4pm" test (informed trading would extend predictability past option/LETF settlement — it does not).

## Key Results

- `r_ROD` positively and significantly predicts `r_LH` across all four asset classes and is the **best predictor** (highest t-stat and OOS R²). Equity index futures: `r_ROD` OOS R² **2.88%** vs `r_ONFH` **−1.71%**.
- Intraday-momentum strategy Sharpe **0.87–1.73** at asset-class level (gross).
- Momentum is **stronger when NGE is more negative**; LETF hedging demand drives its cross-sectional and time-series magnitude.
- Predictability **reverts over ~3 days** (significant mean-reversion in equities, bonds, commodities) → transitory price pressure, supporting the hedging (not informed-trading) channel.
- Distinct from Heston et al. (2010) cross-sectional intraday seasonality; lagged `r_LH(t−1)` does **not** predict market-level `r_LH` (negative for equity futures).

## Limitations / Caveats

- Results are **gross / pre-cost**; intraday + close-auction execution-heavy → transaction costs and close-auction microstructure are first-order in practice.
- Per-instrument SNR is low; Sharpe ratios are quoted at the pooled asset-class level.
- Mechanism evidence (NGE/LETF) is strongest for the S&P 500 and is data-intensive (options gamma, LETF AUM/leverage).
- Currency futures are the weakest (U-shape intraday volume least pronounced).

## Connection to genesis-alpha

This is a **return-predictability** result, well suited to **system identification as predictive filtering** ([[concepts/SysIDReturnPrediction]]) — *not* a market-impact propagator (we forecast returns from returns, not price from observed order flow; the shared content is only the LTI / impulse-response mathematics). The *economic mechanism* is temporary price pressure from short-gamma hedging flow, which reverts. See [[concepts/MarketIntradayMomentum]] for how to (1) replace the equal-weighted `r_ROD` with a *learned* intraday predictive filter (kernel/Bayesian, [[papers/PillonettoDeNicolao2010-KernelSysID]]), (2) model momentum + reversal jointly across horizons, and (3) regime-gate by NGE ([[concepts/SysIDReturnPrediction]]).

## Related Pages

[[concepts/MarketIntradayMomentum]], [[concepts/MarketImpact]], [[concepts/SysIDReturnPrediction]], [[concepts/TIBForm]]
