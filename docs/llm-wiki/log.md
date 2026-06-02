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
