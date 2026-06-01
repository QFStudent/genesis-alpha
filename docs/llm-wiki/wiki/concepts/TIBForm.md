---
type: concept
tags: [execution, market-impact, factor-model]
sources: [Kong2018-TIBInfoGeometry]
updated: 2026-05-31
---

# Triangular Input Balanced (TIB) Form

## Definition

A permutation-invariant state-space parameterisation of a MIMO LTI system where the system matrix takes the form:

```
A = M⁻¹ N Q
```

- **M, N**: sparse bidiagonal matrices derived from the poles λ via:
  - c = 1/√(1 - |λ|²), s = λc
  - M = diag(1/c) · bidiag(s*, c) · signature S
  - N = diag(1/c) · bidiag(c, s) · signature S
- **Q**: block-diagonal rotation matrix that keeps all entries real even for complex conjugate pole pairs (2×2 Givens blocks)
- **U**: null-vector matrix encoding input directions (rows of U[:q,:] are unitary)

The name "Triangular Input Balanced" refers to the triangular structure of the input null vectors and the balanced (information-geometric) parameterisation.

## Why It Matters for genesis-alpha

TIB is the core representation in `ga/filters/tib.py` and drives the entire model reduction pipeline:

1. **Interpretability**: Each TIB state corresponds to a specific *temporal decay rate* (ωₖ) and *input direction* (yₖ). You can say "this state models long-memory market impact driven by input channel 2" — impossible in a generic realization.
2. **Numerical stability**: The bidiagonal sparse structure of M, N avoids ill-conditioning from direct pole multiplication.
3. **Reduction compatibility**: The TIB form is the target representation for `msvdreduce` and `extract_poles_and_nullvecs_from_bt` — reduced systems feed back into TIB for online filtering.
4. **Execution simulator**: The `ir()` function computes the impulse response via the Krylov basis, directly usable in `execution/simulator/` for market impact simulation.

## Empirical Findings

- TIB-based Hankel SVD reduction achieves **H₂ error < 0.01** reducing 100 → 5 poles on power-law decay systems ([[papers/Kong2018-TIBInfoGeometry]], Table 1 / Fig 1); the same experiment shows > 99% Hankel energy captured (Fig 2). Both metrics are reported separately — see [[concepts/ModelReduction]] for the H₂ vs Hankel norm distinction.
- For rational systems, Null Basis Reduction in TIB outperforms standard Hankel SVD in Hankel norm, especially for non-rational (infinite-dimensional) systems.
- The bidiagonal sparse solve (via `scipy.sparse.linalg.spsolve`) is O(n) vs O(n³) for dense inversion.

## Key Papers

- [[papers/Kong2018-TIBInfoGeometry]] — original derivation of TIB form, null basis parameterisation, MIMO reduction algorithms

## Open Questions

- How does TIB reduction quality degrade with noisy impulse response estimates (relevant for live market data)?
- Optimal pole initialisation for market impact: Chebyshev roots (`poles_chebyshev_roots`) vs. data-driven initialisation?
- MIMO TIB for multi-asset impact: how many poles per asset pair before capacity degrades?

## Related Pages

[[concepts/ModelReduction]], [[concepts/InformationGeometry]], [[concepts/MarketImpact]], [[papers/MullhauptRiedel2003-TIBBandMatrix]]
