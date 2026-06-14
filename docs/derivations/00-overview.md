# TIB Model Reduction — Derivation Notes (Overview)

Source material for a LaTeX write-up on the TIB model-reduction machinery. These
notes derive, in one consistent notation, the three pieces and how they connect:

| File | Topic | Answers |
|---|---|---|
| [`01-hankel-svd-reduction.md`](01-hankel-svd-reduction.md) | The `msvdreduce` technique (Hankel-SVD / ERA), cleanly separated from Balanced Truncation | "How does `msvdreduce` work, and how is it different from BT?" |
| [`02-null-basis-realization.md`](02-null-basis-realization.md) | `null_basis_realization` review; why its Blaschke deflation **is** Gram–Schmidt; its connection to `msvdreduce` | "What does the forward map do, and why is the deflation an orthogonalization?" |
| [`03-null-vector-recovery.md`](03-null-vector-recovery.md) | Recovering the null structure from a reduced model: `msvdreduce` vs `tib_from_state_space` | "Two routes to TIB from a reduced realization, and how they relate." |
| [`04-system-norms-and-h2-optimality.md`](04-system-norms-and-h2-optimality.md) | Background: $\ell^2$ / $H^2$ / $H^\infty$ norms, and why $H_2$ is a sub-optimal reduction objective | "What do these norms mean, and why prefer information distance / Hankel norm over $H_2$?" |
| [`05-fast-hankel-matvec.md`](05-fast-hankel-matvec.md) | FFT block-Hankel matvec and how it makes the partial SVD fast | "How is the Hankel SVD computed without forming $H$?" |
| [`06-real-realization-complex-poles.md`](06-real-realization-complex-poles.md) | Keeping the null-basis realization **real** for complex poles (magnitude + rotation); why only $A$ is rotated; the Givens connection | "How do complex poles give a real `(A,B)`, and why isn't $B$ rotated?" |

## Primary sources

- **Yu (2014)** — Xiao Yu, *Identification and Model Reduction of MIMO systems in
  Triangular Input Balanced Form*, PhD dissertation, Stony Brook (advisor
  Mullhaupt). `docs/llm-wiki/raw/xiaoyu_dissertation.pdf`. The governing reference
  for everything here: **§6.1** (reduction review), **§6.2** (Fast Partial Block
  Hankel SVD = `msvdreduce`), **Ch 5** (band-fraction / TIB construction, Thm 38),
  **Ch 3–4** (SISO TIB, lossless parametrization, tangential Schur). Equation
  numbers below refer to this document.
- **Kong (2018)** — `Chikong_dissertation.pdf`; parallel derivation of the TIB form
  and the null-basis construction (Thm 6.5.1).
- **Mu (2026)** — *Notes for Model Reduction* (`Mu2026-ModelReductionNotes`),
  Balanced-Truncation only — the anchor for the BT/`msvdreduce` separation; and
  *Some Proofs: Null Basis Realization* (`Mu2026-NullBasisProofs`), the proof that
  the forward map is TIB-balanced.

Cross-links into the wiki: `concepts/ModelReduction`, `concepts/TIBForm`,
`concepts/ModesAndHankel`, `concepts/BlaschkeFactor`.

## Notation (shared by all three files)

| Symbol | Meaning |
|---|---|
| $q,\;p,\;n$ | number of inputs, outputs, and states (= poles = system order) |
| $(A,B,C,D)$ | state space: $x_{t+1}=Ax_t+Bu_t,\;y_t=Cx_t+Du_t$ |
| $h_k = CA^{k-1}B$ | Markov parameters / impulse response (block, $p\times q$), $k\ge 1$ |
| $H$ | block Hankel matrix of the $h_k$ |
| $\sigma_k$ | Hankel singular values (mode **energies**), $\sigma_k=\sqrt{\lambda_k(PQ)}$ |
| $P,\,Q$ | reachability / observability Grammians (BT) |
| $\lambda_k$ | poles of the original system; $\omega_k$ poles of the **reduced** system |
| $v_k\in\mathbb{C}^q$ | **raw** null vector (input direction) of mode $k$ |
| $y_k\in\mathbb{C}^q$ | **orthonormalized** null vector, $\lVert y_k\rVert=1$ (Blaschke-deflated $v_k$) |
| $t_k=\sqrt{1-\lvert\lambda_k\rvert^2}$ | per-mode energy/innovation gain |
| $b_w(z)=\frac{z-w}{1-\bar w z}$ | scalar Blaschke factor |
| $\beta_{w,u}(z)=I+(b_w(z)-1)\,uu^{*}$ | Blaschke–Potapov factor, $\lVert u\rVert=1$ |
| $J_k=I-(1+\bar\lambda_k)\,y_ky_k^{*}$ | forward rank-1 deflation operator |

Sign/convention note: the code's `blaschke_factor(w, z)` puts the **system pole** at
`w` (a zero of $b_w$ as a function; the analytic pole $1/\bar w$ sits outside the
disk). See `concepts/BlaschkeFactor` for the "which slot is the pole" caveat.

## The unifying thread

Every method here ends in the same two-step move:

1. **Schur-triangularize $A$.** Whatever realization you start from (a BT-reduced
   $(A_r,B_r)$, or the Hankel-SVD shift realization), a real Schur factorization
   $A=QA_1Q^{\mathsf T}$ turns it into an equivalent realization $(A_1,Q^{\mathsf T}B)$
   with $A_1$ (block-)lower-triangular, so the **poles are read off the diagonal**.
2. **Recover the null structure.** $(A_1,Q^{\mathsf T}B)$ is *not yet* TIB (Yu, eq 451).
   You then either (a) peel the orthonormal input directions $y_k$ off the rows of
   the input matrix by **un-deflating** ($\,$`msvdreduce`$\,$, Yu eqs 452–453), or
   (b) re-coordinate into the band-fraction TIB form directly via input-balance +
   orthogonal Schur ($\,$`tib_from_state_space`$\,$, Yu Ch 5 eqs 394–401).

So `msvdreduce` and the BT path differ only in **how they obtain $(A,B)$** (a data
Hankel SVD vs. Grammian balancing) and **how they recover the null structure**
(explicit deflation vs. direct re-coordination). `null_basis_realization` is the
**forward** map that both are inverting.

> **Implementation-note boxes.** Where the current code in `ga/` diverges from the
> theory, it is flagged in a clearly-marked box (Bug 1, Bug 2). These are kept
> separate from the derivation so they can be dropped for a clean theoretical
> reference. No code is changed by these notes.
