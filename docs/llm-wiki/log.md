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

## [2026-06-03] update | Added own-asset-vs-cross-asset / MIMO section to MarketIntradayMomentum.md + Tier-0 notebook
- Created notebooks/intraday_ir_vs_rod.ipynb (synthetic): learned intraday IR (OLS + kernel/Bayesian) vs equal-weight r_ROD; kernel IR robust + best (OOS R^2 ~0.049 vs r_ROD 0.037), OLS IR overfits at small N (negative OOS R^2); sample-size sweep
- Added "Own-asset (diagonal) vs cross-asset (off-diagonal)" section + reworded suggestion #4. Clarifications: "diagonal" = transfer-function matrix H(z), not A; canonical null vectors realize the diagonal input->state map (block-diag C => q independent own-asset filters; full C sneaks in pole-constrained cross-asset); POOLING = tying C across assets (parameter sharing/prior), NOT a consequence of MIMO (untied MIMO = independent fits); clean pooling needs shared pole bank per asset; only C is regularized in fixed-basis fit (pooling and Bayesian-TIB are both C-priors)

## [2026-06-05] update | Added "Computing dealer gamma from OptionMetrics IvyDB" note (Tier 2) to MarketIntradayMomentum.md
- Build-your-own GEX/NGE pipeline: opprcd (gamma + OI) x secprd (spot) join; dgamma = gamma·OI·100·S²·0.01; GEX = calls − puts under the dealer-sign heuristic; gamma-flip spot; point-in-time (t−1 OI); limitations (EOD-only => no 0DTE/intraday; no dealer sign => CBOE/ISE Open-Close; sanity-check vs vendor)

## [2026-06-05] update | Created concepts/OHLCVPooledPrediction.md
- Focused page from the OHLCV->return + pooled-MIMO discussion: own-asset OHLCV inputs (returns/range/log-volume; look-ahead trap from H/L/C bar-end realization; stationarity; collinearity => regularize C; LTI = linear part only)
- Pooled tied-diagonal-block construction: block-diagonal+tied = pooling; shared pole bank REPLICATED not circulated; build-one-block -> states per instrument -> pooled C; vs dense (overfits, loses even with real cross-asset signal -> use small regularized/low-rank factor block). Refs notebooks/pooled_vs_dense_mimo.ipynb
- index.md: 18 pages; linked from SysIDReturnPrediction

## [2026-06-06] update | Created concepts/ModesAndHankel.md
## [2026-06-07] update | Added "What the Hankel singular values mean — energy" section to ModesAndHankel.md (signal-energy / reach-observe Grammian interpretation of σₖ); added docs/derivations/ pointer links to ModelReduction, TIBForm, ModesAndHankel
## [2026-06-07] update | Flagged ⚠️ implementation gap in ModelReduction.md and InformationGeometry.md: info_svd_reduce does not realize info-distance reduction (cepstrum not low-rank-Hankel; dropped rho-damping; transpose shift bug). Corrected companions reduce_svd_truncate_fixed / info_svd_reduce_fixed do identification only, not stable reduction
## [2026-06-07] update | Added Derivations section to index.md listing docs/derivations/ 00–05 (TIB model-reduction write-up + fast Hankel matvec); header now tracks "Derivations: 6"
## [2026-06-11] query | Created concepts/FuturesCovariates.md — input design for predicting ES/NQ/YM/DM/TY/FV/VX/URO/JY with the TIB MIMO model (raw vol-normalized lagged returns + order flow, cross-asset blocks, ES–YM lead-lag/error-correction, LTI cautions)
## [2026-06-13] update | Added "the IR is a p×q matrix of shared-pole channels" section to FuturesCovariates.md (matrix-valued MIMO IR = p·q SISO entries h_ij = C_i A^{τ-1} B_j sharing poles; basis of pooling; links ModesAndHankel / OHLCVPooledPrediction)
## [2026-06-13] query | Created concepts/EnergyBasedCovariateSelection.md — principled input selection by predictable (not raw-IR) energy: CCA/subspace predictable-energy spectrum (→ order + loadings), Hankel-σ contribution, CV group partial-R², group-LASSO + reduced-rank; back-linked from FuturesCovariates
- Precise vocabulary from the modes/eigenvectors discussion: mode = (pole, null vector) input-side descriptor; three distinct objects (pole=scalar/eig of A; state-space eigenvector in C^n; null vector in C^q = input direction); state-space eigenvector A v = λ v decouples modal coords; modal IR h(t)=Σ λ^t (C v_k)(w_k^T B) rank-1 per mode
- Hankel relationship: #modes=rank(H); HSVs = mode energies NOT poles; poles = eigenvalues of shift operator on Hankel range (not eig of H); null vectors from B-block; ties to msvdreduce. TIB A triangular => states coupled => states != pure modes (modes along eigenvectors)
- index.md: 19 pages; linked from TIBForm
- Expanded the state-space-eigenvector section: full dynamical meaning (modal-coordinate decoupling), the diagonalization, the rank-1 modal IR decomposition, and the right/left-eigenvector → output/input-direction table (null vector = input-space image of left eigenvector through B)
- Fixed notation clash (r_k = right eigenvector, w_k = left, v_k reserved for null vector); added the decoupling **proof** (z=R⁻¹x ⇒ z(t+1)=Λz(t); diagonalizable-only, else Jordan/t·λ^t) and the **stability** consequence (|z_k|=|λ_k|^t|z_k(0)| → 0 iff |λ|<1) under the modal-coordinate statement
- Added a **driven-vs-autonomous** note: x(t)→0 is the unforced coast-down; the driven model x(t+1)=Ax(t)+Bu(t) has a persisting forced term, so the forecast y=(h*u) doesn't vanish — pole decay = fading memory (EMA analogy), |λ|≥1 = infinite memory/blow-up
## [2026-06-13] update | Created docs/derivations/06-real-realization-complex-poles.md — keeping the null-basis realization real for complex poles (magnitude+rotation; null_basis_realization_real); why only A is rotated (orthogonal Q cancels in AA*+BB*=I, B untouched); the Givens connection (one-sided to place a pair vs two-sided similarity to preserve). Wired into derivations 00 index + 02 Related; pointer + Q-note added to concepts/TIBForm
## [2026-06-14] query | Created concepts/VARSelfPrediction.md — the TIB state-space as a (V)AR self-prediction / inverse (whitening) filter when input=output=data: CAⁿB are square p×p AR coefficients (not the forward IR); forward↔inverse (synthesis 1/𝒜 vs whitening 𝒜) duality table; ŷ-vs-ε framing note; ar2ir (SISO, dynamical_system) → forward IR via lower-tri Toeplitz / hₙ=Σaₖhₙ₋ₖ, and the MIMO block recursion hₙ=ΣAₖhₙ₋ₖ (h₀=I); LPC / lattice / Levinson–Schur connection (TIB = balanced MIMO inverse filter); caveats (AR exact only for self-prediction not p×q exogenous; whiteness is a fit property; minimum-phase needed for decaying IR). Linked from index; back-links TIBForm / ModesAndHankel / SysIDReturnPrediction / FuturesCovariates
## [2026-06-14] update | VARSelfPrediction.md: added "Reducing for the data poles — and predicting from them" — reduction extracts the poles of whatever IR you feed, so AR coeffs (FIR, poles ~0) must be ar2ir'd to the forward IR before reducing to get the data poles (AR(1) illustration); then predict by building the forward/generative TIB from (data poles, null vectors), inverting it (𝒜=(A_f−B_fC_f, B_f, −C_f, I); Âₖ=C_f(A_f−B_fC_f)^{k-1}B_f), and propagating the inverse on y; pole-distinction callout (predictor poles eig(A_f−B_fC_f) ≠ data poles eig(A_f); Kalman A_f−KC_f equivalent)
## [2026-06-14] update | VARSelfPrediction.md: corrected the prediction workflow per Yu — after ar2ir→reduce you propagate the reduced state-space (A_r,B_r,C_r) directly (x_{t+1}=A_r x_t+B_r u_t, ŷ=C_r x_t); the recursion IS the inverse/prediction filter (data poles carried in A_r), NO separate algebraic inversion step (removed the (A_f−B_fC_f) inverse-realization formula and the pole-distinction callout). Added a "Data poles (definition)" callout: poles of the data-generating process 1/𝒜(z) = roots of det 𝒜(z)=0 = eig(A_r) after ar2ir+reduction (oscillation arg λ / persistence |λ|), distinct from Â(z)'s poles (≈0 for finite AR)
## [2026-06-14] ingest | raw/tib_note.pdf (Yu Mu, TIB Theory, 7pp) → papers/Mu2026-TIBNote.md. Covers orthonormal Blaschke bases (general multi-pole Bᵢ(z) vs Laguerre/Kautz), band-fraction TIB filter A=M⁻¹N (bidiagonal M,N; ρₖ=√(1−|λₖ|²), μₖ=ρ_{k+1}/ρₖ, γₖ=λ̄ₖμₖ; "MIMO: pass"), H₂/H∞ criteria, realization theory, and the **functions-of-Toeplitz AR↔IR**: MA Y=T_hX, AR X=T_aY (T_a first col [1,−a₁,−a₂,…]), duality T_h·T_a=I (H=1/A), giving ar2ir = T_a h = e₁ ⟺ h₀=1, hₙ=Σaₖhₙ₋ₖ. CONFIRMED this is exactly the math behind ga's mimo_ar2ir (block recursion h₀=I, hₙ=ΣAₖhₙ₋ₖ) — the MIMO generalization the note marks "pass". Linked from/source for concepts/VARSelfPrediction (sources + ar2ir-section pointer); index Pages 22→23, Sources 7→8
## [2026-06-14] update | Created docs/derivations/07-mimo-ar-to-ir.md — full multivariate ar2ir derivation filling the "MIMO: pass" of Mu2026-TIBNote. MA Y=T_H X / AR X=T_𝒜 Y (T_𝒜 first block-col [I,−A₁,…,−A_p]); duality T_H T_𝒜=I ⟺ H(z)=𝒜(z)⁻¹; recursion h₀=I, h_m=Σ_{k} A_k h_{m-k} by z⁻ᵐ coefficient matching of 𝒜(z)H(z)=I; left vs right recursion equivalence (H=𝒜⁻¹ unique; verified 2.8e-17); block-Toeplitz solve T_𝒜 h = E₁ (MIMO analog of scalar T_a h=e₁); d=1 reduces to the note; validation table (𝒜H=I & H𝒜=I 5.6e-17, left=right 2.8e-17, scalar match, companion-IR/data-poles) + re-runnable snippet. Wired into derivations 00 index + llm-wiki index (Derivations 7→8) + VARSelfPrediction ar2ir pointer + Mu2026-TIBNote
## [2026-06-15] update | derivations/06-real-realization-complex-poles.md: added a "radial × angular" intuition box to §2 (magnitude→decay on the real axis, rotation→oscillation: spins the {r,r} block until eigenvalues split into {r·e^{±iθ}}; everything stays real, complex poles live only as eig(A) in 2×2 blocks). The rotation trick + the Givens comparison (§5: one-sided to *place* eigenvalues vs textbook two-sided similarity to *preserve*) were already documented here
## [2026-06-15] update | derivations/06-real-realization-complex-poles.md §3: added "Where does the imaginary part go? A rotation IS real complex multiplication" — ℂ≅ℝ², z=a+ib=r·e^{iθ} ↔ real matrix [[a,−b],[b,a]]=r·rotation(θ), i ↔ 90° turn; the imaginary part b=r·sinθ is stored as antisymmetric off-diagonal (geometry, not the symbol i), so there's nothing complex to "remove"; the rotation injects the angle/off-diagonal antisymmetry into the magnitude block (eig {r,r}→{r·e^{±iθ}}); why complex eigenvalues need a conjugate pair / 2×2 block; contrast with null_basis_realization's complex-diagonal entry
## [2026-06-19] update | VARSelfPrediction.md: added "Practical case: ARX / MISO prediction (own-lags + exogenous inputs)" — univariate target from own lagged returns/volume + other instruments = ARX, H(z)=B(z)/𝒜(z) with 𝒜(z)=1−Σaₖz⁻ᵏ. For prediction: all RHS lags observed → regress + predict, NO innovations / ar2ir / inversion. The AR part is SCALAR (target's own past) → data poles = roots of scalar 𝒜(z); exogenous = numerator B(z) (zeros, observed inputs, not poles). mimo_ar2ir never needed (no square multivariate AR); at most SCALAR ar2ir for the explicit input→target IR (IR = bₖ * ar2ir(aₖ)). Decision table; "forward needs innovations" only for self-MA (not input→output w/ observed inputs); reduction works on MISO p=1; full d×d VAR/mimo_ar2ir only for square self-prediction or structural IRF
## [2026-06-19] update | VARSelfPrediction.md ARX section sharpened: added the rigorous justification for "no inversion, no innovations" (one-step predictor = conditional expectation; innovation drops out since E[ε_t|past]=0 → ŷ_t=Σaₖy_{t-k}+Σbₖu_{t-k}, all observed; innovation is the prediction error/output, never an input; autoregression runs the inverse direction data→residual). Added a Two-Caveats callout: (1) ARMAX exception — with an MA error term the predictor DOES use past innovations ε_{t-j} (computed as residuals), so "no innovations" holds for ARX not ARMAX; (2) "no inversion" is prediction-only — extracting the forward IR to reduce needs the SCALAR ar2ir on 𝒜(z) (separate from prediction, never mimo_ar2ir)
## [2026-06-20] update | derivations/01-hankel-svd-reduction.md §1: added "Feedthrough caveat — the Hankel starts at CB, not D". A system with feedthrough D (y_t=Cx_t+Du_t, e.g. y_t=Cx_t+u_t ⇒ D=I) has a lag-0 IR term h_0=D, so ir=[D, CB, CAB, …]; the Hankel is built from h_1,h_2,… only (h_0=D is static, recovered separately as D=ir[0]). Strictly proper (D=0): ir[0]=CB, feed from ir[0]; with feedthrough (D≠0): ir[0]=D, ir[1]=CB → drop ir[0], feed ir[:,:,1:]. msvdreduce/blkhankel assume the array starts at CB (treat ir[:,:,0] as first Hankel block), so passing the full D≠0 IR makes the reducer mistake D for CB
