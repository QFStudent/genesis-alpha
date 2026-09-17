# Transcription notes

All `.py` / `.toml` files transcribed from photographs are **verbatim** — no
markers, scaffolding, fixes or reconstructions inside the files themselves.
Some are therefore not syntactically valid Python; that is intentional. All
coverage and gap information lives here instead.

## Source A — Centiva GitHub Enterprise, `SFT/research`, `BasisArb/MainStrat`

| file | source | lines transcribed | sloc |
|---|---|---|---|
| `ga/rv/backtest.py` | Backtester.py 201 | all | 180/180 |
| `ga/rv/book.py` | Book.py 530 | all | 504/504 |
| `ga/rv/future.py` | Future.py 234 | all | 220/220 |
| `ga/rv/tradeList.py` | TradeList.py 203 | all | 184/184 |
| `ga/rv/simulationStats.py` | SimulationStats.py 201 | all | 174/174 |
| `ga/rv/processIntraDataForBasisArb.py` | 343 | all | 321/321 |
| `ga/rv/strategy.py` | Strategy.py 585 | all but 403-405 | 551/553 |
| `ga/rv/basis/basisAnalysis.py` | BasisAnalysis.py 91 | all but 33-38 | — |

### Gaps
- **`strategy.py` 403-405** — fell in the seam between two photos. Sits inside
  `updateMarginMultiplier`, between the `curMarg *= (1 - (self.gamma1 ...))`
  line and `if self.curStratRf * curSign * fut.riskFactor >= 0:`.
- **`basisAnalysis.py` 33-38** — between `msg['Subject'] = ...` and `s.quit()`.
  `s` is referenced by `s.quit()` and bound nowhere else in the file.
- **`strategy.py` `closeOutDateList`** — one date between `'2019-06-19'` and
  `'2019-09-13'` is absent (photo seam); a 2020-12 block between `'2020-09-17'`
  and `'2021-03-12'` is likewise absent.
- **`book.py` lines 161, 162, 164, 165, 239** — these five lines run off the
  right edge of the photo and are transcribed only as far as they were legible.
  All five are commented-out code.

## Source B — `yu.mu@ny5sftd01:~/yuresearch`, branch `paper`, `lib/sft_ym/cds/backtesting`

| file | source | lines transcribed |
|---|---|---|
| `tmp/beta_estimation.py` | `data/beta_estimation.py` 69 | all (69/69) |
| `tmp/residualization.py` | `data/residualization.py` 66 | all (66/66) |
| `tmp/data.py` | `data/data.py` 148 | all (174 sloc exact) |
| `tmp/oos_event_10bins.toml` | 103 | all |
| `tmp/sampler.py` | `data/sampler.py` | **1-331 only** |
| `tmp/backtesting.py` | `backtesting.py` 716 | **404-716 only** |

### Gaps
- **`tmp/backtesting.py`** — source lines **1-402 are absent** (~56% of the
  file): module docstring, all imports, `class BackTesting`, `__init__`,
  `get_data`, `load_data`, `gen_predictions`, `static_backtest`,
  `add_hp_2_prediction_path`, `debug`. The file therefore begins mid-method,
  inside `pick_prediction`, at source line 404 — it will not parse.
  One trailing comment at source line 440 runs off the right edge.
- **`tmp/sampler.py`** — source continues past **331**; `PreprocessSampler.fit`
  is cut off mid-body. Contains 7 classes: `ClockTimeSampler`,
  `ClockTimeSampler_cc_data78`, `ClockTimeSampler_complete`,
  `ClockTimeSampler_24HOUR`, `ClockTimeSampler_Liquid`, `SpikeTimeSampler`,
  `PreprocessSampler`.

## Not transcribed — insufficient coverage

- **Second `.toml` config.** Shares lines 74-103 verbatim with
  `oos_event_10bins.toml` but differs in `[sampler.kwargs]`:
  `exit_time = '16:00:00'`, `exeGap = '10min'`, `enterVWAP = '5min'`,
  `exitVWAP = '5min'`. Only lines 67-103 photographed; filename unknown.
- **`readAndProcessParamFile` fragment.** snake_case Python, apparently a port
  of BasisArb `Book.setBookParams` (same `maxBookBeta` / `maxBookRiskFactor` /
  `maxBookGmv` / `maxBookPos` keys, plus `libcts.timestamp`). Only lines
  112-140 photographed, most truncated at the right edge; filename unknown.

## Shelved observations

Recorded here only; **nothing was changed in the files.** Live
`pdb.set_trace()` calls in `tmp/backtesting.py` (source 595-597) and
`tmp/data.py` (source 99, 121). `add_trade_day` used at `sampler.py` 284 but
never imported or defined in 1-331. `preprocess_overlapping_data` in `data.py`
is annotated `-> pd.DataFrame` and returns nothing. `np.NaN` (removed in NumPy
2.0) across the Source A files and `tmp/data.py`. `interpolation=` kwarg
(removed in pandas 2.0) in `data.py` `rolling_winsorization`. `tau1` branches in
`tmp/backtesting.py` do not cover exactly `"20200301"`.
