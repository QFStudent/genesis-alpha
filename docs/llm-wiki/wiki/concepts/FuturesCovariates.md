---
type: concept
tags: [futures, multi-asset, factor-model]
sources: []
updated: 2026-06-13
---

# Covariates for Multi-Asset Futures Prediction (TIB MIMO)

## Scope

Input (covariate) design for predicting nine futures with the TIB MIMO state-space
model on an **intraday** horizon (several entries/exits). Universe, by block:

| Block | Instruments |
|---|---|
| Equity index | **ES** (S&P 500), **NQ** (Nasdaq 100), **YM** (Dow 30), **DM** (S&P MidCap 400 — mid-cap) |
| Rates | **TY** (10y Note), **FV** (5y Note) |
| Volatility | **VX** (VIX futures) |
| FX | **URO** (EUR/USD), **JY** (JPY) |

## Design principles (for an LTI / TIB filter)

- **Raw, stationary, causal** inputs — let the impulse response learn the weighting;
  don't pre-shape with TA indicators (the linear filter reproduces any linear indicator
  and finds the optimal one — see [[concepts/SysIDReturnPrediction]]).
- **Log-returns, vol-normalized** (÷ trailing realized vol) and **de-seasonalized**
  intraday (remove the time-of-day mean). LTI assumes stationarity; intraday vol is
  U-shaped and regime-switching, so unnormalized inputs make the poles fit the
  seasonality instead of the dynamics.
- **MIMO**: feed *one shared input vector* — all nine vol-normalized lagged returns plus
  each instrument's order flow — and let the model learn which inputs drive which target
  ([[concepts/OHLCVPooledPrediction]]). The lists below are "which blocks matter most,"
  not separate feature sets.
- **Order-flow imbalance / signed volume** is the strongest causal intraday covariate.

## How the model is structured: the IR is a `p×q` matrix of shared-pole channels

With `q` inputs and `p` outputs, the model's impulse response is **matrix-valued**: at each
lag it is a `p×q` matrix, i.e. `p·q` scalar **SISO channels** `h_{ij}(τ)` = the response of
**output `i`** to a unit impulse on **input `j`**. By superposition each output is the sum of
its input channels,

$$y_i(t) = \sum_{j=1}^{q} (h_{ij} * u_j)(t).$$

These channels are **not independent systems** — they share one internal state and the
**same poles**: `h_{ij}(τ) = C_{i,:}\,A^{τ-1}\,B_{:,j}`, so every channel has identical decay
rates (eigenvalues of `A`) and differs only in *residue* — the output row of `C` times the
input column of `B` (the null-vector projection).

**Why this matters here.** It is the basis of **pooling**: the poles (a shared "pole bank")
are estimated *jointly* across all `p·q` channels and all instruments, while only the
per-channel residue (`C`, null vectors) varies — far fewer parameters, more robust pole
estimates, and the cross-asset lead–lag / error-correction below falls out of the shared
dynamics. This is the tied-block construction of [[concepts/OHLCVPooledPrediction]]. See
[[concepts/ModesAndHankel]] for the modal view `h(τ)=Σ_k λ_k^{τ}(C r_k)(w_k^{⊤}B)` (each
mode hits *every* channel, weighted by its output direction `C r_k` and input direction
`w_k^{⊤}B`).

## Shared core (every target)

own lagged returns · own order-flow imbalance (tick-rule if only trades) · own realized
vol / range · own bid-ask spread & book imbalance (if L1/L2) · prior-session / overnight
gap · time-of-day encoding (as an exogenous input, or de-seasonalize).

## Per-instrument drivers

**Equity block — ES / NQ / YM / DM** (highly cointegrated; intra-block lead–lag is the
key structure — see next section):
- **ES** — NQ, YM, DM lead–lag; VX (inverse vol regime); TY/FV (risk-on/off via rates);
  URO/JY as USD/risk proxy. ES is the most liquid, so its *own* order flow is the
  richest single input and it **leads** the block.
- **NQ** — ES, YM lead–lag; NQ–ES ratio (tech vs broad). **TY/FV especially** — NQ is
  long-duration / rate-sensitive (yields up → NQ down) more than ES. VX (higher vol beta).
- **YM** — ES (dominant), NQ lead–lag; **ES–YM spread** (mean-reverts intraday); TY/FV, VX.
  Mostly ES-beta + a big-cap-industrial residual.
- **DM** — ES (dominant), NQ/YM lead–lag; **DM–ES spread = size factor** (mid vs large-cap
  rotation). More cyclical/domestic → more rate-sensitive (TY/FV); higher beta and **less
  liquid than YM → stronger ES→DM lead–lag** (spread predicts DM catch-up more strongly).

**Rates block — TY / FV** (the curve):
- **TY** — **FV and US/bond (the curve)**; TY–FV spread; ES/NQ (risk-on/off into duration);
  front-end rates (SOFR/SR3) if available; URO/JY (rates↔FX); **macro-event clock** (8:30 ET
  releases, auctions, FOMC).
- **FV** — **TY (most correlated), TU(2y)/US** (belly of the curve); FV–TY, FV–TU spreads;
  front-end / Fed-path rates (SOFR) — FV is very Fed-path-sensitive; ES/NQ risk; same macro clock.

**Volatility — VX** (the least LTI-friendly target):
- **ES returns, esp. downside, and ES realized vol** — the dominant drivers (leverage /
  asymmetry); **VX term-structure slope** (front vs 2nd month — contango/backwardation and
  roll-down, huge into settlement); NQ vol, VVIX if available; time-to-expiry / roll.
- ⚠️ VX is mean-reverting, bounded, asymmetric — a plain linear filter fits it poorly.
  Model **log(VX)** changes, treat downside ES separately, or use VX mainly as an *input*
  to the others rather than a target.

**FX block — URO / JY:**
- **URO** — USD complex (JY, broad DXY; URO is ~half of DXY); **rate differential** (US
  front-end FV/TY/SOFR vs EU rates — carry/diff); ES (mild risk-on EUR bid); **session
  encoding** (London vs NY open) matters more in FX.
- **JY** — **US rates (TY/FV) are the strongest covariate**: USD/JPY tracks US yields, so
  the JY future is *positively* correlated with TY price (yields down → yen up → JY up);
  **risk-off / safe-haven** (ES down + VX up → yen strength, carry unwinds); broad USD (URO);
  EUR/JPY cross; session encoding.

## Cross-asset structure: lead–lag & error-correction (why spreads predict)

Each block is **cointegrated** (shares a common factor: the large-cap market for equities,
the level/curve for rates, the dollar for FX). Within a block, a properly-formed spread —
e.g. `log(YM) − β·log(ES)`, or the vol-matched cumulative-return spread — is **stationary
and mean-reverting**, so the *relative* performance of the two legs is predictable.

**ES–YM as the example.** It predicts YM's **relative (residual) return**, not the
market-beta move. Decompose `r_YM = β·r_market + residual`; the spread predicts the
**residual** (mean-reverting), and is silent on `r_market`. Sign: spread high (YM rich vs
ES) ⇒ expect YM to **underperform** ES next (≈ `−(spread − mean)`). Two mechanisms:
1. **Liquidity lead–lag** — ES is far more liquid, so price discovery leads in ES and **YM
   lags**; the spread measures YM's not-yet-completed catch-up. (ES–DM is stronger still:
   DM is less liquid.)
2. **Relative-value arbitrage** — transient order-flow dislocations get arbitraged back by
   index / stat-arb desks.

**The TIB connection.** An error-correction term is just a linear combination of lagged
ES and YM returns — exactly what the **MIMO filter learns endogenously** when both legs are
inputs (it also finds the right `β` and decay). The hand-built spread is the special case;
feeding the raw legs is the general version. The same logic covers NQ–ES, DM–ES (size),
TY–FV (curve), and URO–JY (dollar).

> ⚠️ **Stale-price artifact.** Asynchronous / stale prints create *spurious* spread
> predictability — the "reversion" is just the lagging leg updating, not tradeable alpha.
> Use synchronized, trade-time, liquid prices. The real edge is fast, heavily arbitraged,
> and cost-sensitive.

## Practical cautions for the LTI setup

- **Regime breaks** (8:30 ET releases, FOMC/ECB, VX settlement) violate stationarity —
  flag them as exogenous event dummies or gate the model around them rather than letting
  the filter average across them.
- **Vol-normalize and de-seasonalize** every input — a single LTI filter can't represent
  multiplicative intraday vol seasonality.
- **Costs / capacity** — intra-block reversion is fast and heavily arbitraged; evaluate
  net of cost (see [[concepts/MarketImpact]]) and consider EWMA / cost-aware smoothing of
  the position.

## Why It Matters for genesis-alpha

This is the input layer for the multi-asset prediction model. The MIMO TIB form is the
natural home: it ingests the shared vol-normalized lagged-return vector and **learns the
cross-asset lead–lag and error-correction structure endogenously** (the pooled,
tied-block construction of [[concepts/OHLCVPooledPrediction]]), so most of the per-asset
"covariate engineering" above is really about *which raw series to include* and *keeping
them stationary*, not hand-crafting features.

## Open Questions

- Which lead–lag edges survive net of cost intraday, and at what horizon (the reversion is
  fast)?
- How much does explicit event-gating beat letting the filter average across macro prints?
- Optimal vol-normalization / de-seasonalization window for the LTI stationarity assumption?

## Related Pages

[[concepts/SysIDReturnPrediction]], [[concepts/OHLCVPooledPrediction]], [[concepts/MarketIntradayMomentum]], [[concepts/MarketImpact]], [[concepts/TIBForm]], [[concepts/EnergyBasedCovariateSelection]]
