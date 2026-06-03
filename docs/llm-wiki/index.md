# Wiki Index

*The LLM updates this on every ingest. Read this first before answering any query — use it to find relevant pages, then drill in.*

Last updated: 2026-06-03 | Pages: 17 | Sources ingested: 7

---

## Papers

- [[wiki/papers/Kong2018-TIBInfoGeometry|Kong 2018 — TIB & Info Geometry]] — PhD dissertation; TIB form derivation, MIMO information geometry, null-basis reduction, fast Toeplitz algorithms
- [[wiki/papers/Mu2026-ModelReductionNotes|Mu 2026 — Model Reduction Notes]] — Internal notes; BT square-root algorithm, biorthogonality/balanced-Grammian proofs, sanity checks
- [[wiki/papers/MullhauptRiedel2003-TIBBandMatrix|Mullhaupt & Riedel 2003 — TIB Band Matrix]] — ⚠️ stub; originating reference for bidiagonal M, N construction in tib.py; paper not yet ingested
- [[wiki/papers/Mu2026-NullBasisProofs|Mu 2026 — Null Basis Realization Proofs]] — inductive proof that null basis algorithm generates TIB matrices; confirms correctness of mimo_null_basis() and explains BT approximation gap
- [[wiki/papers/HanzonOliviPeeters2010-TangentialSchur|Hanzon, Olivi & Peeters 2010 — Tangential Schur]] — primary source for Blaschke factor (eq 26) and elementary J-inner factor (eq 28); confirms (pole, eval) convention
- [[wiki/papers/Olivi2010-LosslessParametrization|Olivi 2010 — Lossless Parametrization (HDR)]] — exact definition of Blaschke-Potapov factor B_{w,u}(z) (eq 1.11) matching the code; unit-vector-in-kernel requirement
- [[wiki/papers/PillonettoDeNicolao2010-KernelSysID|Pillonetto & De Nicolao 2010 — Kernel-Based System ID]] — stable spline kernel; nonparametric Bayesian IR estimation with marginal-likelihood hyperparameters; foundation of the kernel/Bayesian method and Bayesian-TIB fusion
- [[wiki/papers/Baltussen2020-IntradayMomentumHedging|Baltussen et al. 2020 — Intraday Momentum & Hedging]] — market intraday momentum everywhere (60+ futures); gamma-hedging mechanism; r_ROD predicts r_LH then reverts over days

---

## Concepts

- [[wiki/concepts/MarketImpact|Market Impact]] — cost of trading: permanent vs temporary impact, square-root law, venue effects; now includes LTI/TIB empirical findings
- [[wiki/concepts/TIBForm|TIB Form]] — Triangular Input Balanced state-space; interpretable states, sparse bidiagonal parameterisation, genesis-alpha codebase entry point
- [[wiki/concepts/ModelReduction|Model Reduction]] — BT vs Hankel SVD vs info SVD; error bounds; pipeline from IR estimation to execution simulator
- [[wiki/concepts/InformationGeometry|Information Geometry]] — Fisher metric on LTI systems; cepstrum as Euclidean coordinates; foundation for info_svd_reduce
- [[wiki/concepts/BlaschkeFactor|Blaschke Factor]] — b_w(z) and Blaschke-Potapov factor; pole vs evaluation point convention; why the args aren't symmetric (the line-185 bug)
- [[wiki/concepts/PyMORDataDrivenID|pyMOR & Data-Driven System ID]] — deriving poles/null vectors from noisy data; MOR vs sysID; why pyMOR helps reduction not SNR; N4SID/N2SID/Bayesian/Hankel-denoising upstream of tib_from_state_space
- [[wiki/concepts/SysIDReturnPrediction|SysID for Return Prediction]] — I/O design (raw lagged drivers vs TA features), linear/nonlinear boundary, and combining TIB/kernel with tree/DL models (residual vs stacking vs feature-fusion vs meta-labeling)
- [[wiki/concepts/MarketIntradayMomentum|Market Intraday Momentum]] — hedging-demand channel (r_ROD predicts r_LH + multi-day reversal); how to turn it into sysID/TIB signals: learned intraday IR, momentum+reversal propagator, NGE regime-gating, MIMO pooling

---

## Strategies

*(none yet)*

---

## Entities

*(none yet)*

---

## Syntheses

*(none yet — filed when a query produces a substantive multi-paper answer)*
