---
type: concept
tags: [execution, factor-model]
sources: []
updated: 2026-06-06
---

# Modes, poles, null vectors & the Hankel matrix

Precise vocabulary for the TIB / realization picture — what a "mode" is, how the pole, the (state-space) eigenvector, and the null vector differ, and how they relate to the Hankel matrix. These are easy to conflate; this page pins them down.

## A mode

A **mode** is one eigen-component of the dynamics; the number of modes = system order = `n`. In the MIMO/TIB setting a mode carries *directions*, not just a rate. Mode `k` is:
- **pole `λₖ`** — the decay rate; an **eigenvalue of `A`** (the diagonal of the lower-triangular TIB `A`);
- **null vector `vₖ ∈ ℂ^q`** — the **input-space direction** that excites it;
- (+ an output direction, via `C`).

So the **`(pole, null vector)` pair is the input-side descriptor of a mode**. For SISO (`q=1`) the null vector is a trivial unit scalar → a mode ≈ its pole.

## Three different objects — don't conflate

| object | space | meaning |
|---|---|---|
| pole `λₖ` | scalar ℂ | decay rate = eigenvalue of `A` |
| state-space eigenvector | state ℂⁿ | `A vₖ = λₖ vₖ`; direction along which the state evolves as pure `λₖᵗ` |
| **null vector `vₖ`** | input ℂ^q | which input combination drives the mode |

The null vector is **not** the eigenvector of `A` — they live in different spaces (`ℂ^q` vs `ℂⁿ`).

## State-space eigenvector

`A vₖ = λₖ vₖ`, `vₖ ∈ ℂⁿ`. Diagonalize `A = VΛV⁻¹` and write `x(t) = Σ zₖ(t) vₖ`; the modal coordinates **decouple**, `zₖ(t) = λₖᵗ zₖ(0)`. So `vₖ` is the state-space *shape* that evolves coherently at rate `λₖ`. The **right** eigenvector `vₖ` gives the **output** direction `C vₖ`; the **left** eigenvector `wₖ` (rows of `V⁻¹`) gives the **input** direction `wₖᵀ B` ≈ the null vector.

**Modal impulse response** (everything in one place):

$$h(t) = C A^t B = \sum_k \lambda_k^{\,t}\,(C v_k)(w_k^\top B)$$

— a sum of **rank-1** terms: rate × output-direction × input-direction.

## Relationship to the Hankel matrix

The Hankel `H` is built from the Markov parameters `hₖ = C Aᵏ B` and factors as `H = O·C` (observability × controllability).

- **# modes = `rank(H)`** = number of nonzero singular values.
- **Hankel singular values ≠ poles.** They are the mode **energies / importances** — what balanced truncation and Hankel-norm approximation rank and truncate (see [[concepts/ModelReduction]]). (For a symmetric SISO Hankel, its eigenvalues are `±σₖ` = energies, still not poles.)
- **Poles = eigenvalues of the *shift* operator on the Hankel range**, *not* eigenvalues of `H`. Ho–Kalman / ERA: `H = UΣVᵀ`, build the one-block-shifted Hankel `H↑`, `A = Σ^{-1/2}Uᵀ H↑ V Σ^{-1/2}`, poles `= eig(A)`.
- **Singular vectors `U, V`** span the observability / controllability subspaces (the reduced state coordinates); the **null vectors** come from the *input block* of the controllability factor.

This is exactly what `msvdreduce` does: `S` = Hankel singular values (energy / truncation), `Vh` = controllability subspace, the shift of `Vh` → **poles**, the first `q`-block → **null vectors**. Poles and null vectors come from the **shift + B-block**, *not* from any "Hankel eigen-pair."

## In the TIB realization

`A` is **lower-triangular** (real poles) / block-lower-triangular (complex pairs), so:
- its eigenvalues = the poles (the diagonal), realization-invariant;
- but the states are **coupled** (triangular) → a single state coordinate is *not* a pure mode. The decoupled modes sit along the **eigenvectors** of `A` (triangular combinations of states). A *diagonal* realization would make states = modes; TIB trades that for the balanced triangular form (see [[concepts/TIBForm]]).

## Related Pages

[[concepts/TIBForm]], [[concepts/ModelReduction]], [[concepts/BlaschkeFactor]], [[concepts/InformationGeometry]]
