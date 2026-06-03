# Activity Log

Append-only. Format: `## [YYYY-MM-DD] operation | description`

Quick tail: `grep "^## \[" log.md | tail -10`

---

## [2026-05-14] init | Wiki initialized — CLAUDE.md schema, index.md, log.md, wiki/ directory structure created

## [2026-06-01] update | Created concepts/BlaschkeFactor.md
- Explains b_w(z) definition, zero/analytic-pole/all-pass properties
- Pole (interpolation point) vs evaluation point: asymmetry, conjugation cue, role-based disambiguation in null-basis recursion
- Links to the line-185 bug and Olivi/Hanzon sources; index 12 pages

## [2026-06-01] ingest | 2 sources: Hanzon-Olivi-Peeters 2010 + Olivi 2010 HDR (raw/)
- Created HanzonOliviPeeters2010-TangentialSchur.md — Blaschke factor eq(26), elementary J-inner factor eq(28)
- Created Olivi2010-LosslessParametrization.md — exact Blaschke-Potapov def B_{w,u}(z) eq(1.11)
- Both confirm (pole=w, eval=z) convention → validates the line-185 bug fix and the two-sided ‖u‖=1 guard
- index.md: 11 pages, 5 sources

## [2026-06-01] update | Completed Mu2026-NullBasisProofs.md — full inductive proof derived
- Proved Lemma B: LₖLₖ* + PₖPₖ* = Iₘ (induction using JₖJₖ* = I − tₖ²yₖyₖ*)
- Proved Lemma A: A_{k−1}Lₖ* + B_{k−1}Pₖ* = 0 (induction using Lemma B and Jₖ* + ωₖI = (1+ωₖ)(I − yₖyₖ*))
- Completed inductive step: off-diagonal = 0 via Lemma A; diagonal = 1 via Lemma B
- Key insight: TIB condition requires only ‖yₖ‖ = 1, not the full Blaschke-Potapov construction

## [2026-06-01] ingest | null_basis_rep_proofs.pdf (raw/) — Mu 2026 proof of Theorem 6.5.1
- Created wiki/papers/Mu2026-NullBasisProofs.md
- Updated TIBForm.md: clarified why extract_poles_and_nullvecs_from_bt produces approximations
- Updated index.md: 9 pages, 3 sources ingested

## [2026-06-01] update | Fixed lint issues: bugs (#1-3) and missing pages (#4-5)
- Fixed MarketImpact.md sources frontmatter and updated date
- Resolved H₂/Hankel norm inconsistency in TIBForm and ModelReduction
- Fixed TIBForm→ModelReduction cross-reference
- Added POD/IRKA "not used" annotations in ModelReduction table
- Created wiki/concepts/InformationGeometry.md (Fisher metric, cepstrum, info_svd_reduce foundation)
- Created wiki/papers/MullhauptRiedel2003-TIBBandMatrix.md (stub — paper not yet ingested)
- Updated index.md: 8 pages, linked new entries

## [2026-06-01] lint | Full wiki health check — 10 issues found
- Bugs: MarketImpact sources frontmatter empty, stale updated date, H₂/Hankel norm inconsistency in TIBForm vs ModelReduction
- Missing pages: InformationGeometry concept, MullhauptRiedel2003 paper
- Missing cross-refs: TIBForm→ModelReduction, ModelReduction POD/IRKA dead links
- Data gaps: no empirical market impact papers (Almgren-Chriss, Gatheral), no noisy-IR estimation literature
- Suggested ingests: Almgren-Chriss 2000, Gatheral 2010, Laub-Heath-Paige-Ward

## [2026-05-31] ingest | Ingested 2 sources: Kong2018 dissertation (raw/) and Mu2026 model reduction notes (references/)
- Created wiki/papers/Kong2018-TIBInfoGeometry.md — TIB form, MIMO info geometry, null basis reduction, fast Toeplitz
- Created wiki/papers/Mu2026-ModelReductionNotes.md — BT square-root algorithm, proofs, sanity checks
- Created wiki/concepts/TIBForm.md — new concept page for TIB parameterisation
- Created wiki/concepts/ModelReduction.md — new concept page comparing BT, Hankel SVD, info SVD
- Updated wiki/concepts/MarketImpact.md — added empirical findings and key papers sections, linked new concepts

## [2026-06-02] update | Quantified the BT→TIB null-vector approximation gap in TIBForm.md
- Query: do random null vectors improve a TIB model over canonical? → built scripts/compare_null_vectors.py
- Part 1: canonical = decoupled axis-aligned input coupling (0% cross-channel); random = dense (84%); both valid TIB. Null-vector choice picks WHICH system, not its quality.
- Part 2: on an order-12/6-in/4-out cross-coupled target — oracle null vectors fit to ~1e-15 and bt_impulse_response to ~1e-14, but extract_poles_and_nullvecs_from_bt → null_basis_realization gives rel err 0.53, no better than canonical (0.35)/random (0.51)
- Root cause: extractor returns orthonormalized yₖ, fed in as raw vₖ → wrong transfer function (yₖ-vs-vₖ gap from Mu2026-NullBasisProofs)
- Updated TIBForm.md: ⚠️ block in "Why It Matters" #3 + new Open Question to fix the extraction; bumped updated date

## [2026-06-02] update | Resolved the BT→TIB extraction gap — added tib_from_state_space + tests
- Added ga/reducers/balanced_truncation.py::tib_from_state_space(A, B): input-balance (controllability Grammian → I via Cholesky-factor similarity) + orthogonal Schur → lower-triangular TIB form; returns (TIBStateSpace, C-transform); real-poles-only (ValueError on 2×2 Schur blocks)
- Reproduces target to ~1e-14; tests/test_bt_to_tib.py (5 tests: round-trip IR, input-balance, poles-on-diagonal, lower-triangular, complex-pole guard); full suite 39 passed
- null_basis_realization intentionally NOT changed (the fix needs no round-trip through it)
- No production consumers of extract_poles_and_nullvecs_from_bt; left in place, marked deprecated-for-fitting
- Updated TIBForm.md: ⚠️→✅ resolution in #3, Open Question struck through (remaining: complex poles via Q rotation); scripts/compare_null_vectors.py gains a 'corrected' reference row (~1e-14)

## [2026-06-02] update | Extended tib_from_state_space to complex poles
- Removed the real-pole guard: input-balance + real-Schur already produces a valid TIB realization for complex poles — complex-conjugate pairs come out as 2×2 real diagonal blocks (block-lower-triangular / lower-Hessenberg), realization stays real and input-balanced, IR reproduced to ~3e-15
- The orthogonal real-Schur step delivers exactly the block structure TIBForm's Q rotation describes, with no extra work (verified empirically before coding)
- tests/test_bt_to_tib.py: dropped the reject-complex test, added TestComplexPoles (5 tests: round-trip IR, input-balance, realness, eigenvalue recovery, block-lower-triangular); full suite 43 passed
- TIBForm.md: ✅ block now states complex poles handled; Open Question complex-pole caveat removed

## [2026-06-02] query | Created concepts/PyMORDataDrivenID.md — noise-robust derivation of poles/null vectors
- Query: is pyMOR (pymor.org) better, or are there alternatives, for deriving poles/null vectors from data at low SNR?
- Key framing filed: MOR != sysID. pyMOR (ERA/Loewner/AAA/BT/IRKA) is mature reduction + a validation reference, but does NOT address measurement noise
- The real lever is upstream system identification: N4SID/MOESP/CVA, N2SID (nuclear-norm, short batches), Bayesian/kernel-based ID, Hankel denoising (optimal SVHT Gavish-Donoho, Cadzow/SLRA), TLS-ERA, freq-domain VF/AAA/Loewner — then feed (A,B,C) into the exact tib_from_state_space
- Python tooling noted: SIPPY, control+slycot, PyDMD
- index.md: 13 pages; linked from TIBForm/ModelReduction/InformationGeometry/MarketImpact

## [2026-06-02] update | Added Bayesian-TIB fusion section to PyMORDataDrivenID.md
- Reframed (per user's view): TIB = linear basis expansion of the IR (fixed poles+null vectors define P basis functions; estimate C by LS, = fit_output_matrix)
- Unifying claim: TIB and kernel/Bayesian ID are both linear IR estimators / basis methods; TIB = hard-truncated explicit basis, kernel = soft-weighted implicit basis; same "stable decaying modes" prior
- Corrected earlier framing: they are NOT opposite ends of bias-variance; difference is regularisation mechanism (subset selection vs shrinkage = best-subset vs ridge)
- Fusion (Bayesian-TIB): estimate the same C but with a kernel prior over an overcomplete pole grid, marginal-likelihood tuned -> decouples basis richness from effective complexity; new Open Question to implement/benchmark it vs plain TIB-LS

## [2026-06-02] query | Created concepts/SysIDReturnPrediction.md — I/O design + linear/nonlinear fusion
- Q: for sysID (TIB/kernel) on financial data, what are inputs/outputs? lagged returns vs TA-lib features? how to combine with nonlinear (tree/DL) models?
- Filed: output = forward returns (multi-horizon = MIMO, vol-normalised); inputs = raw causal lagged drivers (lagged returns, order flow/OFI, cross-asset) NOT linear TA indicators (subsumed by the IR); sysID = the feature engineering; nonlinear TA value is unreachable by LTI
- Fusion architectures: A residual/boosting-offset (fragile at low SNR), B stacking (safe default), C feature-fusion (TIB states into nonlinear model -> regime gating), D meta-labeling (size, don't re-predict); avoid pure residual-chasing; strict purged/embargoed CV
- index.md: 14 pages; linked TIBForm/PyMORDataDrivenID/ModelReduction/MarketImpact

## [2026-06-03] ingest | PilloDeNico-Automatica-2010 (raw/) — Pillonetto & De Nicolao kernel-based system ID
- Created wiki/papers/PillonettoDeNicolao2010-KernelSysID.md — stable spline kernel K(s,t)=W(e^-βs,e^-βt) (integrated Wiener on exp-warped time axis); encodes smoothness AND BIBO-stability; nonparametric RKHS/GP estimate; 3 hyperparameters (λ,β,σ) by marginal likelihood; two-step "estimate then project" for reduced-order models
- Key results: beats ETFE/cubic-spline/Gaussian-kernel/PEM+AIC on 100-sample benchmarks, near PEM+oracle; order-30 random systems K=0.23 vs PEM+AIC 0.35/BIC 0.32/oracle 0.21; 95% CI coverage 93.7%; works on reduced/nonuniform grids
- Updated PyMORDataDrivenID.md (sources frontmatter + links + references) and SysIDReturnPrediction.md (related pages); the TC kernel in notebooks/kernel_bayesian_id.ipynb is the discrete simplification of this stable spline kernel
- index.md: 15 pages, 6 sources ingested

## [2026-06-03] ingest | heding demand and market intraday momentum.pdf (raw/) — Baltussen, Da, Lammers & Martens 2020
- Created wiki/papers/Baltussen2020-IntradayMomentumHedging.md — market intraday momentum everywhere (60+ futures, 1974-2020); rest-of-day return r_ROD predicts last-half-hour r_LH; gamma-hedging-demand mechanism (NGE proxy, LETF rebalancing); reverts over ~3 days (transitory price pressure); equity r_ROD OOS R^2 2.88%, Sharpe 0.87-1.73 gross
- Created wiki/concepts/MarketIntradayMomentum.md — phenomenon + mechanism + the sysID connection: it's a predictable-flow market-impact propagator. Suggestions to build signals: (1) learn the intraday IR (kernel/Bayesian) instead of equal-weight r_ROD; (2) model momentum+reversal as one TIB propagator; (3) regime-gate by NGE (feature-fusion/meta-labeling); (4) joint MIMO pooling across assets; (5) trade the reversal leg
- Updated MarketImpact.md (key papers + empirical: predictable hedging flow = temporary impact at market level) and SysIDReturnPrediction.md (related pages)
- index.md: 17 pages, 7 sources ingested

## [2026-06-03] update | Added "Data sources for hedging pressure" to MarketIntradayMomentum.md
- Tiered the data dependency: Tier 0 price-only (no alt data; r_ROD / learned IR), Tier 1 mechanical flow (LETF AUM+leverage rebalancing demand, NYSE/Nasdaq MOC imbalance feeds), Tier 2 dealer gamma (NGE/GEX from options OI + greeks + dealer-sign assumption; SqueezeMetrics/SpotGamma or OptionMetrics/OPRA; vanna/charm, 0DTE)
- Practical notes: sequence Tier 0->1->2; point-in-time discipline (look-ahead); gamma/LETF evidence is equity-index-specific
