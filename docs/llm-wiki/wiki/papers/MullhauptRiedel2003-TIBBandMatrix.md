---
type: paper
tags: [execution, market-impact]
sources: [MullhauptRiedel2003-TIBBandMatrix]
updated: 2026-06-01
---

# Band Matrix Representation of Triangular Input Balanced Form

**Authors:** Mullhaupt, A., Riedel, K.  **Venue / Year:** Technical Report, August 2003

> ⚠️ Stub — paper not yet ingested. Entry created from citations in `ga/filters/tib.py` and [[papers/Kong2018-TIBInfoGeometry]]. Ingest the actual report to fill in Key Results and Limitations.

## Contribution

The originating reference for the bidiagonal sparse matrix parameterisation of TIB state-space systems. Establishes that the system matrix A = M⁻¹NQ can be expressed using bidiagonal band matrices M, N (derived from the poles λ) and a rotation matrix Q, yielding a numerically superior and permutation-invariant state-space form.

## Known from Citations

From `ga/filters/tib.py` (`realtib` docstring) and [[papers/Kong2018-TIBInfoGeometry]]:
- M and N are bidiagonal sparse matrices constructed from poles λ via: c = 1/√(1 − |λ|²), s = λc
- Q is a block-diagonal rotation matrix handling complex conjugate pole pairs (2×2 Givens blocks)
- The construction satisfies A = M⁻¹NQ and is equivalent to a balanced state-space realisation but with O(n) sparse structure instead of O(n²) dense

## Connection to genesis-alpha

- `sptib()` and `siso_system_matrices()` in `ga/filters/tib.py` are direct implementations of the band matrix construction from this paper.
- `siso_system_matrices_real()` and `realtib()` implement the real-valued variant (with Q rotation for complex poles).

## Related Pages

[[concepts/TIBForm]], [[papers/Kong2018-TIBInfoGeometry]]
