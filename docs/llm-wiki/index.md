# Wiki Index

*The LLM updates this on every ingest. Read this first before answering any query — use it to find relevant pages, then drill in.*

Last updated: 2026-08-03 | Pages: 24 | Derivations: 8 | Sources ingested: 8

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
- [[wiki/papers/Mu2026-TIBNote|Mu 2026 — TIB Note]] — TIB theory notes: orthonormal Blaschke bases, band-fraction filter `A=M⁻¹N`, H₂/H∞ criteria, and the **functions-of-Toeplitz** AR↔IR derivation (`T_h T_a = I`, `ar2ir` = `T_a h = e₁`); source for [[wiki/concepts/VARSelfPrediction]], MIMO case ("pass") filled by `mimo_ar2ir`

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
- [[wiki/concepts/OHLCVPooledPrediction|OHLCV inputs & pooled MIMO prediction]] — own-asset OHLCV features (returns/range/log-volume + look-ahead/stationarity checklist) and the pooled tied-diagonal-block construction (shared pole bank, build-one-block-pool-C) vs dense; demo notebook
- [[wiki/concepts/ModesAndHankel|Modes, poles, null vectors & the Hankel matrix]] — precise definitions: mode = (pole, null vector); state-space eigenvector vs null vector; modal IR decomposition; poles = shift-operator eigenvalues (not Hankel eigenvalues), HSVs = mode energies
- [[wiki/concepts/FuturesCovariates|Futures Prediction Covariates]] — per-instrument covariates for ES/NQ/YM/DM/TY/FV/VX/URO/JY in the TIB MIMO model; shared raw/vol-normalized inputs, cross-asset blocks, ES–YM (and intra-block) lead–lag / error-correction
- [[wiki/concepts/EnergyBasedCovariateSelection|Energy-Based Covariate Selection]] — principled input selection by *predictable* (not raw IR) energy: CCA/subspace predictable-energy spectrum, Hankel-σ contribution, CV group partial-R², group-LASSO + reduced-rank
- [[wiki/concepts/VARSelfPrediction|(V)AR Self-Prediction]] — the TIB state-space as an inverse/whitening filter when input=output=data: `CAⁿB` are square (V)AR coefficients (not the forward IR); LPC/lattice/Schur connection; `ar2ir` (and its MIMO block recursion) maps AR coeffs → forward IR; minimum-phase caveat

---

## Derivations

*Long-form math write-ups in `docs/derivations/` (outside the wiki vault) — source for a LaTeX write-up. Grounded in Yu (2014) dissertation §6 (`raw/xiaoyu_dissertation.pdf`) + [[wiki/papers/Kong2018-TIBInfoGeometry]] + [[wiki/papers/Mu2026-ModelReductionNotes]]. Cross-linked from the TIBForm / ModelReduction / ModesAndHankel concept pages.*

- [00 — Overview](../derivations/00-overview.md) — shared notation, source map, the unifying "Schur-triangularize → recover null structure" thread
- [01 — Hankel-SVD reduction (`msvdreduce`)](../derivations/01-hankel-svd-reduction.md) — block Hankel as the past→future operator, shift realization, Schur poles, `B = QᵀB`, why `(A₁,QᵀB)` isn't yet TIB; BT-vs-`msvdreduce` separation; SISO / `info_svd_reduce` bug notes
- [02 — `null_basis_realization` & Blaschke = Gram–Schmidt](../derivations/02-null-basis-realization.md) — forward-map review; why the Blaschke–Potapov deflation is Gram–Schmidt in the H² (reproducing-kernel) metric
- [03 — Null-vector recovery](../derivations/03-null-vector-recovery.md) — `msvdreduce` (tangential-Schur) vs `tib_from_state_space`; rebuild directly from `(poles, y)` — do **not** re-deflate
- [04 — System norms & H₂ optimality](../derivations/04-system-norms-and-h2-optimality.md) — ℓ²/H²/H∞ definitions and relationships; why H₂ is a sub-optimal reduction objective (information distance / cepstrum)
- [05 — Fast Hankel matvec (FFT)](../derivations/05-fast-hankel-matvec.md) — `O(k log k)` Hankel matvec (reversal → convolution → FFT), block extension, Lanczos partial SVD; `reduce_fft_truncate` float32 / triangular-Hankel bugs
- [06 — Real realization for complex poles](../derivations/06-real-realization-complex-poles.md) — magnitude+rotation construction (`null_basis_realization_real`); why only `A` is rotated (orthogonal `Q` cancels in `AA*+BB*=I`, `B` untouched); the Givens connection (one-sided to place a pair vs. two-sided similarity to preserve)
- [07 — MIMO AR ↔ IR (`mimo_ar2ir`)](../derivations/07-mimo-ar-to-ir.md) — full multivariate `ar2ir`: block-Toeplitz duality `T_H T_𝒜 = I`, `H(z)=𝒜(z)⁻¹`, recursion `h₀=I, hₙ=Σ Aₖhₙ₋ₖ` (left=right verified), block solve `T_𝒜 h = E₁`; fills the "MIMO: pass" in [[wiki/papers/Mu2026-TIBNote]]; validation table

---

## Strategies

*(none yet)*

---

## Entities

*(none yet)*

---

## Syntheses

- [[wiki/syntheses/SchlichtVC-ProdSimReplication|Schlicht VC vs Prod TOD-MISO “VC” — Replication Notes]] — two distinct “VC” objects (Schlicht `VC.m` vs prod `calc_tv_coefficients` TOD stack); `ga/vc/` vs `ga/sim/` mapping; replication confidence and open gaps (`calc_tv_coefficients` stub)
