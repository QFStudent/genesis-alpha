# Real realization for complex poles — the rotation trick

How to make the **forward** map produce a *real* `(A,B)` when some poles are
complex, staying entirely within the null-basis construction. This is what
`null_basis_realization_real` (`ga/filters/tib.py`) implements, and it answers two
questions that recur in practice:

- *Why does feeding complex poles to `null_basis_realization` give complex
  predictions, and how do we avoid it?*
- *Why does the fix rotate only $A$, never $B$ — and is that a Givens rotation?*

Notation is shared with [`00-overview.md`](00-overview.md): poles $\lambda_k$,
orthonormalized null vectors $y_k$, $t_k=\sqrt{1-|\lambda_k|^2}$, and the TIB
input-balance invariant $AA^{*}+BB^{*}=I$.

## 1. The problem: complex poles ⇒ complex realization

`null_basis_realization(\lambda, v)` places each pole **literally on a complex
diagonal** ($A_{kk}=\lambda_k$) and deflates one pole at a time. For a
complex-conjugate pair this gives a genuinely **complex** $(A,B)$: the Blaschke
orthonormalization is order-dependent, so the two rows of $B$ belonging to a pair
do **not** come out as conjugates of each other (empirically off by $0.2$–$0.5$).
That breaks the only structure that could make the impulse response real:

$$h_k = C A^{k} B \in\mathbb{R} \iff (A,B,C)\ \text{is conjugate-consistent.}$$

With non-conjugate $B$-rows, **no** $C$ — real or complex — makes $C A^{k}B$ real.
So a real complex-pole system cannot be obtained by putting complex poles on a
diagonal; the conjugate pair must be kept together as a **$2\times2$ real block**.

## 2. The construction: magnitude system + rotation

Split each complex pole $\lambda = r\,e^{i\theta}$ into its **magnitude** $r=|\lambda|$
(a *real* number) and its **phase** $\theta=\arg\lambda$, and handle them in two stages.

1. **Radial stage — build on the magnitudes.** Form the radial pole list (signed
   value for a real pole, magnitude $r$ — *repeated*, for each member of a pair) and
   order it so each conjugate pair is adjacent. Run `null_basis_realization` on it.
   Because every radial pole is real, the result $(A_{\text{mag}},B_{\text{mag}})$ is
   **real, input-balanced, and lower-triangular**, and each conjugate pair sits as a
   $2\times2$ diagonal block
   $$A_{\text{mag}}[j{:}j{+}2,\,j{:}j{+}2]=\begin{bmatrix} r & 0 \\ x & r\end{bmatrix},
     \qquad x = A_{\text{mag}}[j{+}1,\,j].$$
   The two diagonal entries are equal ($=r$) because $|\lambda|=|\bar\lambda|$; the
   lower entry $x$ is the within-pair coupling the deflation produced.

2. **Angular stage — inject the phase by an orthogonal rotation.** Build a
   block-diagonal $Q$ ($1$ for each real pole, a $2\times2$ rotation $G(\psi)$ per
   pair) and set
   $$\boxed{\,A = A_{\text{mag}}\,Q,\qquad B = B_{\text{mag}}\,.}$$

The eigenvalues of $A_{\text{mag}}Q$ become the conjugate pairs, $B$ is unchanged,
and (§4) the input balance is preserved. Verified: pole error $\sim10^{-15}$,
$\lVert AA^{\mathsf T}+BB^{\mathsf T}-I\rVert\sim10^{-16}$, $\operatorname{Im}h_k=0$.

## 3. The rotation angle (read off the realized block)

$A_{\text{mag}}Q$ is block-lower-triangular ($A_{\text{mag}}$ lower-triangular, $Q$
block-diagonal), so its spectrum is the union of the diagonal $2\times2$ products
$A_{ii}\,G(\psi)$. For one pair,

$$\begin{bmatrix} r & 0 \\ x & r\end{bmatrix}
  \begin{bmatrix}\cos\psi & -\sin\psi\\ \sin\psi & \cos\psi\end{bmatrix}
  \quad\Longrightarrow\quad
  \det = r^{2},\qquad \operatorname{tr} = 2r\cos\psi - x\sin\psi .$$

The eigenvalues satisfy $\mu^{2}-(\operatorname{tr})\mu+r^{2}=0$. We want
$\mu=r\,e^{\pm i\theta}$, i.e. $\operatorname{tr}=2r\cos\theta$ (the determinant is
already $r^{2}=|\lambda|^{2}$). Solving $2r\cos\psi - x\sin\psi = 2r\cos\theta$:

$$\boxed{\,\psi = -\operatorname{atan2}(x,\,2r)\;+\;\arccos\!\frac{2r\cos\theta}{\sqrt{4r^{2}+x^{2}}}\,}$$

Two points:

- **The angle is read off the realized block, not from a pole-only formula**, because
  $x$ is created by the deflation of *preceding* poles (and the null vectors).
  - If a pair has **no preceding poles**, $x=0$, the block is $r\,I$, and the formula
    collapses to $\psi=\theta$ — *just rotate by the pole angle.*
  - Once other poles are deflated first, $x\ne0$ and $\psi\ne\theta$; you must read $x$.
- **Contrast with the SISO `siso_system_matrices_real`.** That routine uses a closed
  form in the pole alone, $\psi=\arccos\frac{2r}{1+r^{2}}-\arccos\frac{2r\cos\theta}{1+r^{2}}$,
  because its *bidiagonal* construction has a fixed coupling ($x=-0.51$ in one test).
  The SISO closed form and the $x$-based formula are the two branches of the same
  equation $2r\cos\psi-x\sin\psi=2r\cos\theta$; both land on the same conjugate pair.

## 4. Why only $A$ is rotated, never $B$

Three ways to see it, strongest last.

**(i) The phase lives in $A$.** Poles are $\operatorname{eig}(A)$; $B$ does not enter
the eigenvalues. Moving a pole from $r$ to $r\,e^{\pm i\theta}$ is purely a
pole-location change, so it acts on the matrix whose spectrum we are editing. $B$ is
the input/reachability map (it came from the null vectors) — there is no frequency to
put there.

**(ii) It is a construction, not a similarity.** A similarity $T^{-1}AT,\ T^{-1}B$
*preserves* eigenvalues and the transfer function. Here $A=A_{\text{mag}}Q$ is a
**one-sided** product that deliberately *changes* the eigenvalues, so $(A,B)$ is a new
system, not a re-coordinatization. (If it were a similarity $Q^{\mathsf T}A_{\text{mag}}Q$,
the poles would stay at the magnitudes — wrong.)

**(iii) Input balance forces $B$ to stay put.** The TIB invariant is
$AA^{*}+BB^{*}=I$, and $A_{\text{mag}}A_{\text{mag}}^{*}+B_{\text{mag}}B_{\text{mag}}^{*}=I$
holds from the radial stage. With $Q$ **orthogonal**,

$$A A^{*} = (A_{\text{mag}}Q)(A_{\text{mag}}Q)^{*}
          = A_{\text{mag}}\,\underbrace{QQ^{*}}_{=I}\,A_{\text{mag}}^{*}
          = A_{\text{mag}}A_{\text{mag}}^{*},$$

so $AA^{*}+BB^{*}=A_{\text{mag}}A_{\text{mag}}^{*}+B_{\text{mag}}B_{\text{mag}}^{*}=I$
**with $B$ untouched.** The $QQ^{*}=I$ cancellation inside $AA^{*}$ is exactly why $B$
does not move; rotating $B$ as well would break the identity (or demand a
compensating change) and yield a non-balanced — i.e. non-TIB — system.

## 5. Is it a Givens rotation?

**Structurally, yes — that is the building block.** $Q$ is block-diagonal with
$2\times2$ planar rotations, one per conjugate pair acting on that pair's two
coordinates. A $2\times2$ planar rotation on a coordinate pair *is* a Givens
rotation, so $Q$ is a product of Givens rotations — one per oscillatory mode.

**In purpose, it is the opposite of the textbook use.**

| | textbook Givens (QR, Hessenberg, real Schur) | the TIB $Q$ |
|---|---|---|
| applied | **two-sided**, $G^{\mathsf T}AG$ (similarity) | **one-sided**, $A_{\text{mag}}\,Q$ |
| effect on spectrum | **preserves** eigenvalues (introduces zeros) | **changes** eigenvalues (injects phase) |
| goal | triangularize / standardize a block | place a real block's eigenvalues at $r\,e^{\pm i\theta}$ |

The only departure from a standard Givens step: $\psi$ is not chosen to annihilate an
entry, but solved from $2r\cos\psi-x\sin\psi=2r\cos\theta$ so the conjugate pair lands
at the target angle. Same primitive (a planar rotation), used to *place* a pair rather
than to *zero out* a coordinate.

## Implementation note — `null_basis_realization_real`

`ga/filters/tib.py: null_basis_realization_real(lambd, v)` (companion to
`null_basis_realization`; the original is left intact). It orders real poles first
then conjugate pairs adjacent, runs the radial stage, reads $x=A_{\text{mag}}[j{+}1,j]$
per pair, forms $Q$ with the $\psi$ above, and returns
`TIBStateSpace(A = A_mag @ Q, B = B_mag)`. Output: real $(A,B)$,
$\operatorname{eig}(A)=$ the poles, input-balanced, block-lower-triangular, real IR;
recoverable by `msvdreduce_complex`. Tests: `tests/test_null_basis_real.py`.

This is the **forward** analog, in the null-basis representation, of the rotation
$Q$ already present in the band-fraction `TIBSystem` (`A = M^{-1}NQ`) and in
`siso_system_matrices_real`. It keeps complex poles real where
`null_basis_realization` (complex diagonal) cannot.

## Related

- [`02-null-basis-realization.md`](02-null-basis-realization.md) — the forward map this extends (real-pole case).
- [`03-null-vector-recovery.md`](03-null-vector-recovery.md) — the real-Schur route (`tib_from_state_space`), the *recovery*-side analog of keeping complex pairs as $2\times2$ blocks.
- `concepts/TIBForm` (the $Q$ in $A=M^{-1}NQ$), `concepts/BlaschkeFactor`.
