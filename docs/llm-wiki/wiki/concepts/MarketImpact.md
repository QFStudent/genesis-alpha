---
type: concept
tags: [market-impact, liquidity, execution]
sources: [Kong2018-TIBInfoGeometry, Mu2026-ModelReductionNotes]
updated: 2026-05-31
---

# Market Impact

## Definition
The adverse price movement caused by executing a trade. Split into two components:

- **Temporary impact:** Intraday price pressure that reverts after the trade is complete. Driven by liquidity supply/demand imbalance. Proportional to participation rate.
- **Permanent impact:** Lasting price change reflecting information content of the trade. Does not revert. Related to adverse selection.

Standard empirical model (square-root law):
```
Impact ≈ σ · √(Q / ADV)
```
where σ is daily volatility, Q is order size, ADV is average daily volume. Coefficient varies by venue and asset class (~0.1–1.0 in practice).

## Why It Matters for genesis-alpha

Market impact is the dominant cost for any signal with turnover. It enters the backtest in `backtest/costs/` and constrains realistic AUM in `risk/sizing/`. Ignoring it produces overly optimistic Sharpe ratios. The execution simulator (`execution/simulator/`) should implement a calibrated impact model.

Key questions for us:
- What impact model do we use in the backtest cost layer?
- How does impact scale with our target AUM?
- Are our signals capacity-constrained by impact at realistic position sizes?

## Empirical Findings

The LTI/TIB framework models market impact as the impulse response of a causal filter from order flow to price. Key findings from the wiki so far:

- A **5-pole TIB model** captures > 99% of the Hankel energy for power-law decay impact systems ([[papers/Kong2018-TIBInfoGeometry]], Table 1).
- The **H∞ error bound** for balanced truncation is tight: ‖Gₙ - Gᵣ‖_{H∞} ≤ 2 × sum of neglected Hankel singular values ([[papers/Mu2026-ModelReductionNotes]]).
- Hankel norm is more relevant than H∞ for execution cost modelling — it bounds response energy rather than worst-case gain ([[concepts/ModelReduction]]).

## Key Papers

- [[papers/Kong2018-TIBInfoGeometry]] — TIB form derivation; MIMO model reduction; interpretable state-space for impact modelling
- [[papers/Mu2026-ModelReductionNotes]] — BT square-root algorithm; error bound proof; sanity checks

## Open Questions

- Exact coefficient in the square-root law for our specific universe (futures vs equities)
- How does impact differ across futures (high liquidity, large notional) vs small-cap equities?
- Interaction between impact and signal decay: does the alpha survive net of impact costs?

## Related Pages

[[concepts/TIBForm]], [[concepts/ModelReduction]]
