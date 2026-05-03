# Genesis Alpha — Codebase Architecture

Production-level research environment for systematic futures and equity trading.

---

## Directory Structure

```
genesis-alpha/
├── data/               # Data layer — the foundation everything else depends on
│   ├── ingestion/      # Raw connectors: exchange APIs, vendors (Refinitiv, BBG)
│   ├── processing/     # Cleaning, normalization, corp actions, roll adjustments
│   ├── store/          # Unified data access layer (abstraction over storage backend)
│   └── schemas/        # Data contracts — enforce column names, dtypes, freq
│
├── universe/           # What you trade and when
│   ├── futures/        # Contract specs, roll schedules, continuous series construction
│   └── equities/       # Corporate actions, delistings, index membership history
│
├── alpha/              # Signal and factor library (production-grade, tested)
│   ├── signals/        # Individual signal implementations (stateless, vectorized)
│   ├── factors/        # Factor definitions (momentum, carry, value, etc.)
│   └── combination/    # Weighting, orthogonalization, alpha blending
│
├── research/           # Exploratory — notebooks only, nothing imported from here
│   ├── signals/
│   ├── strategies/
│   └── tearsheets/
│
├── backtest/           # Vectorized backtesting engine
│   ├── engine/         # Event loop / vectorized portfolio sim
│   ├── costs/          # Transaction cost models (linear, square-root, futures-specific)
│   ├── portfolio/      # Construction: optimization, risk parity, rank-weighting
│   └── analytics/      # Sharpe, drawdown, factor attribution, turnover
│
├── risk/               # Risk layer (separate from backtest — used live too)
│   ├── metrics/        # Vol, VaR, beta, correlation, Greeks
│   ├── limits/         # Hard limits checked pre-order
│   └── sizing/         # Position sizing rules
│
├── execution/          # Sim and live share the same interface
│   ├── simulator/      # Realistic fill model: latency, partial fills, market impact
│   ├── oms/            # Order management: state machine, rejects, partial fills
│   └── brokers/        # Broker connectors behind a common interface
│
├── config/             # All params live here, not in code
├── tests/
├── docs/
└── sample_data/
```

---

## Key Principles

### 1. Research → Production Promotion Path
The hardest thing to get right. Research notebooks should import from `alpha/` and `backtest/`, never the reverse. When a signal graduates from research, it gets refactored into `alpha/signals/` with tests, then the notebook becomes a tearsheet.

### 2. Data Access Through One Interface
Everything reads data through `data/store/` — never directly from files. This lets you swap storage backends (Parquet → Arctic → TimescaleDB) without touching research or alpha code.

### 3. Futures vs Equities Are Different in Ways That Matter
- **Futures:** continuous contract construction (Panama, back-adjust, ratio), roll schedules, notional sizing via point value, margin not cash
- **Equities:** corporate action adjustments, borrow costs, short constraints, index rebalance effects

Treat them as separate universe modules but share the same downstream interfaces.

### 4. Execution Simulator Must Match Live
`execution/simulator/` and `execution/brokers/` should implement the same interface. Your backtest engine calls the simulator; your live system calls the broker. This is what "live trading level simulator" really means — not just fill price, but order state machine, latency, queue position.

### 5. Config Over Code
Strategy parameters, universe filters, cost assumptions, risk limits — all in `config/`. Backtests should be fully reproducible by committing a config snapshot alongside results.

---

## Stack Recommendations

| Layer | Recommendation |
|---|---|
| Storage | Parquet + DuckDB for research; consider Arctic3 for live timeseries |
| Computation | NumPy/Pandas for research; Polars or Numba for hot paths in backtest |
| Pipeline orchestration | Prefect or Airflow for data ingestion jobs |
| Config | Pydantic v2 models — validates types, serializes to JSON |
| Testing | pytest + hypothesis for signal math; integration tests hit real data fixtures |

---

## Common Pitfall

The biggest mistake in quant codebases is letting research code bleed into production. Starting with this separation from day one — even if `alpha/` and `research/` initially contain the same ideas — pays off enormously as the codebase grows.
