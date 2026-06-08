# `null_basis_realization` — Review, and Why Blaschke Deflation **Is** Gram–Schmidt

> Reviews the forward map `null_basis_realization` in `ga/filters/tib.py`, places it
> as the map that `msvdreduce` inverts, and gives a precise explanation of why its
> Blaschke–Potapov deflation step is a Gram–Schmidt orthogonalization — just in the
> curved $H^2$ metric rather than the flat $\mathbb{C}^q$ one. Grounding: Yu Ch 3–4
> (SISO TIB, orthogonal-function POV, lossless parametrization), Kong (2018) Thm
> 6.5.1, and `Mu2026-NullBasisProofs`. Notation: [`00-overview.md`](00-overview.md).

## 1. What the forward map does

`null_basis_realization(lambd, v)` takes $n$ **(pole, raw null vector)** pairs
$\{(\lambda_k, v_k)\}$, $v_k\in\mathbb{C}^q$, and assembles a state-space pair
$(A,B)$ in **Triangular Input Balanced** form: $A\in\mathbb{C}^{n\times n}$
(block-)lower-triangular with the poles on the diagonal, $B\in\mathbb{C}^{n\times q}$,
satisfying the input-balance identity

$$A_kA_k^{*} + B_kB_k^{*} = I_k \quad\text{at every step } k.$$

Each pair is one decaying mode: $\lambda_k$ is the decay rate ($|\lambda_k|<1$),
$v_k$ the input direction that excites it. The recursion (Kong Thm 6.5.1;
`Mu2026-NullBasisProofs` Algorithm 1) is, for $k=1,\dots,n$:

$$
\begin{aligned}
&y_k = \operatorname{normalize}\!\Big(M_k\,v_k\Big),\qquad
M_k = \prod_{i=1}^{k-1}\beta_{\lambda_i,\,y_i}(\lambda_k^{*})^{*}, \quad (M_1=I)\\[2pt]
&t_k = \sqrt{1-|\lambda_k|^2},\qquad
J_k = I - (1+\bar\lambda_k)\,y_ky_k^{*},\\[2pt]
&L_k = \big[\,J_{k-1}L_{k-1}\ \big|\ t_{k-1}y_{k-1}\,\big],\qquad
P_k = J_{k-1}P_{k-1}\quad (P_1=I),\\[2pt]
&\text{append to } A:\ \text{row } (t_k\,y_k^{*}L_k,\ \lambda_k);\qquad
\text{append to } B:\ \text{row } b_k = t_k\,y_k^{*}P_k .
\end{aligned}
$$

The cast of variables (cross-referenced to `concepts/TIBForm`):

| Symbol | Code | Role |
|---|---|---|
| $v_k$ | `v[:, k]` | raw input direction (null vector) being added |
| $y_k$ | `ys[:, k]` | **orthonormalized** direction, $\lVert y_k\rVert=1$ — the Gram–Schmidt output |
| $t_k$ | `t` | $\sqrt{1-\lvert\lambda_k\rvert^2}$ — unit-energy normalizer |
| $J_k$ | `Js[k]` | $I-(1+\bar\lambda_k)y_ky_k^{*}$ — rank-1 deflation carried into the **state recursion** |
| $P_k$ | `P` | $J_{k-1}\cdots J_1$ — cumulative deflation, builds the $B$ row |
| $L_k$ | `L` | accumulates $A$'s sub-diagonal coupling |

Two facts to keep separate, because they use *different* deflations:

- The **input-direction orthogonalization** uses the Blaschke–Potapov product $M_k$
  (§3 below). This is what makes $y_k$ the right direction (the transfer function).
- The **state recursion** uses $J_k$ and its cumulative product $P_k$. This is what
  makes the realization *balanced*. The balance identity follows from $\lVert
  y_k\rVert=1$ alone (`Mu2026-NullBasisProofs`, Lemmas A & B; the key identity is
  $J_kJ_k^{*}=I-t_k^2 y_ky_k^{*}$). Balance does **not** test whether the $y_k$
  directions are right — which is exactly why a wrong Blaschke step passes the
  balance check while realizing the wrong transfer function (the historical
  argument-swap bug; see `concepts/BlaschkeFactor`).

## 2. The forward / inverse picture (connection to `msvdreduce`)

`null_basis_realization` is the **forward** map

$$\{(\lambda_k, v_k)\}\ \xrightarrow{\ \text{Blaschke deflation}\ }\ (A,B).$$

`msvdreduce` (and the BT path) need the **inverse**: from a reduced realization
$(A,B)$, recover the poles and the input directions. The poles are immediate (Schur
diagonal). The directions are encoded in the $B$ rows as

$$b_k = t_k\,y_k^{*}P_k,\qquad P_k = J_{k-1}\cdots J_1,$$

so the *first* row gives $y_1$ directly ($P_1=I$), but every later row has been
pre-multiplied by the cumulative deflation $P_k$. To read $y_k$ off row $k$ you must
**undo $P_k$**, peeling one $J_i^{-1}$ at a time:

$$J_i^{-1} = I - \Big(1+\tfrac{1}{\bar\lambda_i}\Big)y_iy_i^{*}.$$

That incremental un-deflation is precisely the `msvdreduce` `u`-loop (Yu eqs
452–453). The forward map deflates with coefficient $(1+\bar\lambda_i)$; the inverse
un-deflates with $(1+1/\bar\lambda_i)$. The mechanics and the code bug live in
[`03-null-vector-recovery.md`](03-null-vector-recovery.md); the point here is the
**duality**: `null_basis_realization` and the `u`-loop are inverse Gram–Schmidt
passes through the same Blaschke geometry.

## 3. Why the Blaschke deflation is Gram–Schmidt

### 3.1 Flat Gram–Schmidt, written as a projector

Classical Gram–Schmidt orthonormalizes $v_1,v_2,\dots$ by removing, at step $k$, the
overlap with the already-placed unit vectors:

$$y_k = \operatorname{normalize}\big(\Pi_{<k}\,v_k\big),\qquad
\Pi_{<k} = \prod_{i<k}\big(I - y_iy_i^{*}\big) = I - \sum_{i<k} y_iy_i^{*}.$$

Each factor $I-y_iy_i^{*}$ is the **orthogonal projector** that kills the $y_i$
component. Gram–Schmidt = "apply the product of the previous projectors, then
normalize."

### 3.2 The modes do not live in $\mathbb{C}^q$ — they live in $H^2$

A mode is not a bare vector; it is a **rational vector function**. Mode $k$ sits at
pole $\lambda_k$ with input direction $v_k$, and the object being orthogonalized is
its representative in the vector Hardy space $H^2(\mathbb{C}^q)$ — the
reproducing-kernel Hilbert space of the disk, with Szegő kernel
$\kappa_w(z)=\tfrac{1}{1-\bar w z}$ (Yu §2.3/§3.2; `concepts/InformationGeometry`,
`Olivi2010-LosslessParametrization`). In that space, two modes at *different* poles
are **not** orthogonal even when their direction vectors $v$ are orthogonal in
$\mathbb{C}^q$: the inner product couples direction and pole location through
$\kappa$. So flat $\mathbb{C}^q$ projection is the wrong orthogonalization — you need
the projector of the $H^2$ metric.

### 3.3 The Blaschke–Potapov factor **is** that projector

The elementary inner factor (Olivi 2010 eq 1.11)

$$\beta_{w,u}(z) = I + \big(b_w(z)-1\big)\,uu^{*},\qquad b_w(z)=\frac{z-w}{1-\bar w z},\ \lVert u\rVert=1,$$

is exactly the curved-metric deflation along direction $u$ at interpolation point
$w$. The decisive identity: **evaluate it at its own pole** $z=w$. Since
$b_w(w)=0$,

$$\boxed{\;\beta_{w,u}(w) = I + (0-1)\,uu^{*} = I - uu^{*}\;}$$

— the *flat* Gram–Schmidt projector. And it has the **null property**
$\beta_{w,u}(w)\,u = b_w(w)\,u = 0$ (Kong eq 249; `concepts/BlaschkeFactor`): each
factor annihilates its own direction at its own pole, the curved analog of "$y_i$ is
removed from the span." Away from $w$ ($b_w(z)\ne 0$), the factor is still a
rank-1 deflation along $u$, but **weighted by $b_w(z)$** — and that weight is the
pole-geometry correction that the $H^2$ metric demands. Flat Gram–Schmidt is the
special case "weight $=0$, evaluate at the pole."

### 3.4 The recursion is Gram–Schmidt with this projector

Line up the forward step against flat Gram–Schmidt:

$$
\underbrace{y_k = \operatorname{normalize}\Big(\,\textstyle\prod_{i<k}\beta_{\lambda_i,y_i}(\lambda_k^{*})^{*}\ v_k\Big)}_{\text{Blaschke deflation (this code)}}
\qquad\longleftrightarrow\qquad
\underbrace{y_k = \operatorname{normalize}\Big(\,\textstyle\prod_{i<k}(I-y_iy_i^{*})\ v_k\Big)}_{\text{flat Gram–Schmidt}} .$$

It is term-for-term the same algorithm: *replace the flat projector
$I-y_iy_i^{*}$ by the Blaschke–Potapov factor $\beta_{\lambda_i,y_i}$, and probe the
product at the new mode's (reflected) location $\lambda_k^{*}$ before applying it to
$v_k$.* Evaluating at $\lambda_k^{*}$ is the RKHS way of asking "how much does the
already-placed deflation operator see the new mode": the conjugate $\lambda_k^{*}$
(not $\lambda_k$) appears because interpolation conditions for inner functions live
at points reflected across the circle (`concepts/BlaschkeFactor`, §"Why this
evaluate-at-other-poles structure exists"). The output $y_k$ is the new direction
with its overlap against all earlier modes removed — orthonormalized in the metric
that the dynamics actually impose.

**One-line statement.** `null_basis_realization` is recursive Gram–Schmidt for
decaying modes; the only change from the textbook version is that the orthogonal
projector $I-y_iy_i^{*}$ is upgraded to the Blaschke–Potapov factor
$\beta_{\lambda_i,y_i}$, which *equals* that projector at its own pole and carries the
correct pole-dependent weighting elsewhere.

## Related

- [`01-hankel-svd-reduction.md`](01-hankel-svd-reduction.md) — `msvdreduce`, the reduction this feeds.
- [`03-null-vector-recovery.md`](03-null-vector-recovery.md) — the inverse pass (the `u`-loop) and its bug.
- `concepts/TIBForm`, `concepts/BlaschkeFactor`, `papers/Mu2026-NullBasisProofs`.
