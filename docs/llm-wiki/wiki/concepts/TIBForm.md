---
type: concept
tags: [execution, market-impact, factor-model]
sources: [Kong2018-TIBInfoGeometry]
updated: 2026-06-02
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

## Intuition: how `null_basis_realization` builds the system

`null_basis_realization(lambd, v)` takes a pile of **(pole, null vector)** pairs and assembles them into a state-space `(A, B)`. Think of each pair as **one decaying mode**:

- **`λₖ`** = the mode's *decay rate*. |λ|<1; close to 1 → slow decay / long memory, close to 0 → fast transient. (Market impact: how long a trade's footprint lingers.)
- **`vₖ`** (null vector) = the mode's *input direction* — which combination of inputs excites that mode. Called "null" because it is a unit vector in the kernel of the inner function at `λₖ` (the Blaschke-Potapov factor annihilates it: `β_{λₖ,vₖ}(λₖ)·vₖ = 0`; see [[concepts/BlaschkeFactor]]).

**The core idea: it is recursive Gram-Schmidt for decaying modes.** The raw directions `vₖ` are not mutually orthogonal; using them directly would double-count energy. So the algorithm orthonormalises them as it adds them — that is what the `ys` are:

> **`yₖ` is the orthonormalised version of `vₖ`** — what is left after removing its overlap with all previously placed modes.

In flat Gram-Schmidt you subtract projections. Here the geometry depends on the *pole locations*, so instead you apply the product of Blaschke-Potapov factors as the deflation operator:

```
M = ∏_{i<k} β_{λᵢ, yᵢ}(λₖ*)      # curved-space "subtract previous projections"
yₖ = normalize(M @ vₖ)            # orthonormalised direction
```

**The cast of intermediate variables:**

| Symbol | Code | Role | Flat-space analogy |
|---|---|---|---|
| `vₖ` | `v[:, k]` | raw input direction (null vector) | the vector being added |
| `yₖ` | `ys[:, k]` | **orthonormalised** direction, ‖yₖ‖=1 | Gram-Schmidt output |
| `tₖ` | `t` | `√(1−\|λₖ\|²)` — energy/innovation gain | scale giving the state *unit* energy |
| `Jₖ` | `Js[k]` | `I − (1+λ̄ₖ)yₖyₖ*` — rank-1 deflation op | "remove yₖ going forward" |
| `Pₖ` | `P` | `J_{k−1}⋯J₁` cumulative deflation | all previous deflations stacked → builds B row |
| `Lₖ` | `L` | accumulates A's coupling columns | how new mode couples to old → builds A row |

**Why `tₖ = √(1−|λₖ|²)`:** energy normaliser. A mode `xₜ = λxₜ₋₁ + tₖeₜ` has unit stationary variance exactly when `tₖ = √(1−|λₖ|²)`, so each state stores precisely one unit of energy — this is what makes the realisation *balanced*. Slow modes (λ near 1) get a *small* tₖ because they accumulate energy over time.

**Why A is lower-triangular:** each new mode is driven by the *earlier* modes (the `a_row = tₖ yₖ* Lₖ` coupling) plus its own pole `λₖ` on the diagonal — but not by later ones. The triangular shape is a fingerprint of the recursive build order.

**The invariant:** every step preserves `AₖAₖ* + BₖBₖ* = I` (input-balanced). Crucially this holds for *any* unit `yₖ`, so the balance identity cannot detect whether the Blaschke orthogonalisation was done correctly — the identity guarantees *balance*, the Blaschke product guarantees the *right system*. (This is exactly why the line-185 argument-swap bug passed the sanity check while realising the wrong transfer function — see [[papers/Mu2026-NullBasisProofs]] and [[concepts/BlaschkeFactor]].)

One-line model: **recursive Gram-Schmidt that adds decaying modes one at a time — `vₖ` is the raw input direction, `yₖ` its orthonormalised version, `tₖ` gives it unit energy, and `Jₖ`/`Pₖ`/`Lₖ` carry forward "what's already placed" so the next mode knows what to be orthogonal to.**

## Why It Matters for genesis-alpha

TIB is the core representation in `ga/filters/tib.py` and drives the entire model reduction pipeline:

1. **Interpretability**: Each TIB state corresponds to a specific *temporal decay rate* (ωₖ) and *input direction* (yₖ). You can say "this state models long-memory market impact driven by input channel 2" — impossible in a generic realization.
2. **Numerical stability**: The bidiagonal sparse structure of M, N avoids ill-conditioning from direct pole multiplication.
3. **Reduction compatibility**: The TIB form is the target representation for `msvdreduce` and `extract_poles_and_nullvecs_from_bt` — reduced systems feed back into TIB for online filtering. `extract_poles_and_nullvecs_from_bt` produces **approximate** null vectors because BT output doesn't provide exact null basis pairs; the TIB condition AₖAₖ* + BₖBₖ* = I is only guaranteed when exact pairs are used (see [[papers/Mu2026-NullBasisProofs]]).
   > ⚠️ **The approximation gap is severe, not cosmetic** (quantified 2026-06-02, `scripts/compare_null_vectors.py`). On an order-12, 6-input, 4-output target built from known cross-coupled null vectors: the target is exactly TIB-representable (oracle null vectors → IR rel err `~1e-15`) and balanced truncation reproduces it (`bt_impulse_response` → `~1e-14`), but routing BT output through `extract_poles_and_nullvecs_from_bt` → `null_basis_realization` gives **rel err 0.53 — no better than canonical (0.35) or random (0.51)**. Root cause: the extractor returns the *orthonormalized* yₖ (normalised rows of the Schur-transformed B_r), but `null_basis_realization` treats its input as *raw* vₖ and re-applies the Blaschke-Potapov deflation — so the recovered system has the right poles and satisfies TIB balance but realises the **wrong transfer function**. Pole recovery itself is fine (`~4e-15`), and the extractor only handles **real** poles (real Schur + diagonal read-off).
   >
   > ✅ **Fixed** by `tib_from_state_space(A, B)` (in `ga/reducers/balanced_truncation.py`, 2026-06-02): instead of guessing null vectors, it re-coordinates the (reduced) system into TIB form directly — **input-balance** (similarity by the Cholesky factor of the controllability Grammian → controllability Grammian = I, i.e. AA*+BB*=I) followed by an **orthogonal real-Schur** rotation to (block-)lower-triangular. Reproduces the target to `~1e-14` (round-trip tests in `tests/test_bt_to_tib.py`; demo row in `scripts/compare_null_vectors.py`). **Handles complex poles**: real poles are 1×1 diagonal entries, complex-conjugate pairs come out as **2×2 real diagonal blocks** (block-lower-triangular / lower-Hessenberg) — the orthogonal real-Schur step delivers exactly the block structure that TIBForm's `Q` rotation describes, with no extra work, and the realization stays real and input-balanced.
4. **Execution simulator**: The `ir()` function computes the impulse response via the Krylov basis, directly usable in `execution/simulator/` for market impact simulation.

## Empirical Findings

- TIB-based Hankel SVD reduction achieves **H₂ error < 0.01** reducing 100 → 5 poles on power-law decay systems ([[papers/Kong2018-TIBInfoGeometry]], Table 1 / Fig 1); the same experiment shows > 99% Hankel energy captured (Fig 2). Both metrics are reported separately — see [[concepts/ModelReduction]] for the H₂ vs Hankel norm distinction.
- For rational systems, Null Basis Reduction in TIB outperforms standard Hankel SVD in Hankel norm, especially for non-rational (infinite-dimensional) systems.
- The bidiagonal sparse solve (via `scipy.sparse.linalg.spsolve`) is O(n) vs O(n³) for dense inversion.

## Key Papers

- [[papers/Kong2018-TIBInfoGeometry]] — original derivation of TIB form, null basis parameterisation, MIMO reduction algorithms

## Open Questions

- ~~Fix the BT → TIB null-vector extraction so a fitted model reproduces its target.~~ **Resolved 2026-06-02** by `tib_from_state_space` (input-balance + real-Schur; see ✅ in *Why It Matters* #3). Handles both real and complex poles (complex pairs → 2×2 real diagonal blocks). `extract_poles_and_nullvecs_from_bt` is left in place but should be considered deprecated for fitting.
- How does TIB reduction quality degrade with noisy impulse response estimates (relevant for live market data)?
- Optimal pole initialisation for market impact: Chebyshev roots (`poles_chebyshev_roots`) vs. data-driven initialisation?
- MIMO TIB for multi-asset impact: how many poles per asset pair before capacity degrades?

## Related Pages

[[concepts/ModelReduction]], [[concepts/InformationGeometry]], [[concepts/BlaschkeFactor]], [[concepts/ModesAndHankel]], [[concepts/MarketImpact]], [[papers/MullhauptRiedel2003-TIBBandMatrix]]
