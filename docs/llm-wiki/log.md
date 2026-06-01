# Activity Log

Append-only. Format: `## [YYYY-MM-DD] operation | description`

Quick tail: `grep "^## \[" log.md | tail -10`

---

## [2026-05-14] init | Wiki initialized — CLAUDE.md schema, index.md, log.md, wiki/ directory structure created

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
