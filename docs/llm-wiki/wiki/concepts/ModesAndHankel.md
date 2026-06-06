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

## State-space eigenvector — the full picture

**Definition.** The state-space eigenvectors are the (right) eigenvectors of `A`:

$$A\,v_k = \lambda_k\, v_k, \qquad v_k \in \mathbb{C}^n .$$

They live in the **state space** `ℂⁿ` (dimension `n` = order = number of states) — *not* in the input space `ℂ^q` (where the null vectors live), and not in the output space.

**What it means dynamically.** `vₖ` is the direction in state space that the dynamics act on as **pure scaling by `λₖ`** — it does not mix into other directions. Diagonalize `A = V Λ V⁻¹` (`V = [v₁ … vₙ]`, `Λ = diag(λₖ)`) and write the state in that basis,

$$x(t) = \sum_k z_k(t)\, v_k .$$

Under the autonomous dynamics `x(t+1) = A x(t)` each **modal coordinate decouples**:

$$z_k(t) = \lambda_k^{\,t}\, z_k(0) .$$

So if you start the state *exactly* along `vₖ`, it **stays along `vₖ`** and merely decays/oscillates as `λₖᵗ`. That is the precise sense of "mode": a **shape (`vₖ`) in state space that evolves coherently at its own rate (`λₖ`)**, independent of the other modes.

**The modal impulse-response decomposition** ties all the objects together. With `Aᵗ = V Λᵗ V⁻¹ = Σₖ λₖᵗ vₖ wₖᵀ`, where `wₖᵀ` are the rows of `V⁻¹` — the **left** eigenvectors (`wₖᵀ A = λₖ wₖᵀ`) — sandwiching with `C … B` gives

$$h(t) = C A^t B = \sum_{k=1}^{n} \lambda_k^{\,t}\,\underbrace{(C v_k)}_{\text{output direction}}\,\underbrace{(w_k^\top B)}_{\text{input direction}} .$$

Each mode contributes a **rank-1** term: **rate × output-direction × input-direction**. This is where every object slots in:

| object | space | role in the mode |
|---|---|---|
| pole `λₖ` | scalar ℂ | the **rate** (decay / oscillation); `= eig(A)` |
| **right** eigenvector `vₖ` | state ℂⁿ | the mode's state-space **shape**; gives the **output direction** `C vₖ ∈ ℂ^p` |
| **left** eigenvector `wₖ` | state ℂⁿ | the dual direction; gives the **input direction** `wₖᵀ B ∈ ℂ^q` |
| **null vector** (input) | input ℂ^q | the input combination that drives the mode = `wₖᵀ B` (normalized) |

**So the null vector is the input-space *image* of the left eigenvector through `B`** — related to, but living in a different space than, the state-space eigenvector. That is the precise reason **"null vector ≠ eigenvector of `A`"**: the null vector is in input space `ℂ^q`, the eigenvector in state space `ℂⁿ`. Symmetrically, the **output direction `C vₖ`** is the output-space image of the *right* eigenvector through `C`.

**Two readings of one mode:**
- *State view* — a direction `vₖ` in `ℂⁿ` that decays as `λₖᵗ`.
- *Input/output view* — a rank-1 channel `λₖᵗ · (output direction)(input direction)` from inputs to outputs.

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
