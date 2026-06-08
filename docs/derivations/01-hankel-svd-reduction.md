# Hankel-SVD Reduction (`msvdreduce`) — Derivation

> Derives the `msvdreduce` algorithm in `ga/reducers/hankel.py` from Yu §6.2 ("Fast
> Partial Block Hankel SVD"), with the reduction-method context of §6.1. A boxed
> section at the end separates this method cleanly from **Balanced Truncation**,
> which is the subject of the user's own notes (`Mu2026-ModelReductionNotes`).
> Notation: see [`00-overview.md`](00-overview.md).

## 1. The problem

We are given the impulse response (Markov parameters) of a stable MIMO LTI system,

$$h_k = CA^{k-1}B \in \mathbb{C}^{p\times q}, \qquad k = 1, 2, \dots, $$

and want a **reduced** realization $(\hat A,\hat B,\hat C)$ of order $r\ll n$ whose
impulse response $\hat h_k=\hat C\hat A^{k-1}\hat B$ is close to $h_k$. The natural
criterion is $H_2$:

$$\min_{\deg \hat H = r}\ \sum_{k\ge 1}\lVert h_k-\hat h_k\rVert_F^2 .$$

Crucially, **only the impulse response is given** — no state-space $(A,B,C)$. This
is the defining difference from Balanced Truncation (§5). `msvdreduce` is a
*realize-then-reduce* method: it builds a state space directly from the data Hankel
matrix.

## 2. The block Hankel matrix and its partial SVD

Stack the Markov parameters into the block Hankel matrix (Yu eq 440):

$$H = \begin{bmatrix} h_1 & h_2 & \cdots & h_n\\ h_2 & h_3 & & h_{n+1}\\ \vdots & & & \vdots\\ h_n & h_{n+1} & \cdots & h_{2n-1}\end{bmatrix},\qquad h_i\in\mathbb{C}^{p\times q}.$$

The classical fact behind every Hankel realization method: $H$ factors as
observability $\times$ controllability,

$$H = \mathcal{O}\,\mathcal{C},\qquad \mathcal{O}=\begin{bmatrix}C\\ CA\\ \vdots\end{bmatrix},\quad \mathcal{C}=\begin{bmatrix}B & AB & \cdots\end{bmatrix},$$

so $\operatorname{rank}H = $ McMillan degree, and the **range of $H$** is the
observability subspace. Reduction = take the dominant part of $H$.

**Fast partial SVD (Yu §6.2, eqs 435–442).** For long impulse responses we do not
form $H$ densely. $H$ multiplies a vector in $O(n\log n)$ via the FFT (eqs 436–439;
implemented by `FastHankelProduct` in `hankel.py`). The block case is a sum of $q$
scalar Hankel products (eq 441). Since $H$ is generally not symmetric but $H^{*}H$ is
Hermitian PSD, its **singular values equal the eigenvalues of $H^{*}H$**, so a
**Lanczos** partial eigendecomposition of $H^{*}H$ (and $HH^{*}$) yields the rank-$\gamma$
partial SVD (eq 442)

$$H \approx \widetilde H = \widetilde U\,\widetilde S\,\widetilde V^{*},\qquad \widetilde S = \operatorname{diag}(\sigma_1\ge\cdots\ge\sigma_\gamma)>0 .$$

The $\sigma_k$ are the **Hankel singular values** — mode energies; their tail bounds
the reduction error. (The reference implementation `msvdreduce` calls a dense
`la.svd`; the "fast/partial" machinery is the production path and does not change the
algebra below.)

## 3. The shift realization

Yu's key identification (eq 443): in the truncated factorization, the right factor
$\widetilde V^{*}$ plays the role of the **controllability factor**,

$$\widetilde V^{*} \approx \begin{bmatrix} B & AB & \cdots & A^{k-1}B\end{bmatrix} =: \mathcal{C}.$$

This identification is the engine of the whole method — the shift step below only
recovers $A$ because of it — so it is worth deriving rather than asserting.

**The Hankel factorization, oriented.** Take the factorization $H=\mathcal{O}\mathcal{C}$
from §2 and write the block index explicitly:

$$H = \underbrace{\begin{bmatrix} C\\ CA\\ CA^2\\ \vdots\end{bmatrix}}_{\mathcal{O}\ \text{(observability)}}\ \underbrace{\begin{bmatrix} B & AB & A^2B & \cdots\end{bmatrix}}_{\mathcal{C}\ \text{(controllability)}},\qquad
H_{ij} = (CA^{\,i-1})(A^{\,j-1}B) = CA^{\,i+j-2}B = h_{i+j-1}.$$

The orientation is the whole point: $\mathcal{O}$ is the **left** (tall) factor, indexed
by output-time; $\mathcal{C}$ is the **right** (wide) factor, indexed by input-time.
Reading off the two subspaces of $H$,

$$\operatorname{range}(H)=\operatorname{col}(\mathcal{O})=\text{observability subspace},
\qquad
\operatorname{rowspace}(H)=\operatorname{row}(\mathcal{C})=\text{controllability subspace}.$$

So the row space — which the **right** singular vectors span — is the controllability
side. ($§2$ noted the range/observability half; this is the row-space/controllability
half.)

**What the Hankel matrix actually is: past inputs → future outputs.** The factorization
above is not a formal trick — $H$ is the matrix of a genuine, physical operator, and the
Markov parameters $h_k$ are its entries. Run the experiment: drive the system with inputs
only in the *past* ($u_s$, $s\le-1$), switch the input **off at $t=0$**, and watch the
output in the *future* ($t\ge0$). With $y_t=\sum_{k\ge1}h_k u_{t-k}$ (causal, $D=0$) and
$u_s=0$ for $s\ge0$, only the terms reaching before $0$ survive ($k\ge t+1$):

$$y_t \;=\; \sum_{k\ge t+1} h_k\,u_{t-k} \;\overset{j=k-t}{=}\; \sum_{j\ge1} h_{t+j}\,u_{-j},\qquad t\ge0 .$$

The coefficient linking past input $u_{-j}$ to future output $y_t$ is $h_{t+j}$ —
depending only on $t+j$, i.e. **constant along anti-diagonals**: exactly the Hankel
matrix $H_{ij}=h_{i+j-1}$. So $H$ *is* the map (past inputs) $\mapsto$ (future outputs).

The **state is the bottleneck**, which is *why* it factors as $\mathcal{O}\mathcal{C}$:
the past charges the state,
$x_0=\sum_{k\ge1}A^{k-1}Bu_{-k}=\mathcal{C}\,u_{\text{past}}$, and that state then rings
out into the future, $y_t=CA^{t}x_0$, i.e. $y_{\text{future}}=\mathcal{O}\,x_0$;
composing gives $y_{\text{future}}=\mathcal{O}\mathcal{C}\,u_{\text{past}}=H\,u_{\text{past}}$:

$$u_{\text{past}}\ \xrightarrow{\ \mathcal{C}\ }\ x_0\ \xrightarrow{\ \mathcal{O}\ }\ y_{\text{future}},\qquad \Gamma=\mathcal{O}\,\mathcal{C}.$$

The state is the only channel between past and future, which is why $\operatorname{rank}H
= $ state dimension.

**MIMO shapes — the factors are block matrices.** With $q$ inputs, $p$ outputs, $n$
states, the pieces carry block structure but compose cleanly. Here $B$ is $n\times q$,
$C$ is $p\times n$, each Markov parameter $h_k=CA^{k-1}B$ is a $p\times q$ **block**, and
$N$ counts the block rows/columns:

| object | shape |
|---|---|
| $u_{\text{past}}=(u_{-1},u_{-2},\dots)$, each $u_{-j}\in\mathbb{R}^{q}$ | $qN\times 1$ |
| $\mathcal{C}=[\,B\ \ AB\ \ A^{2}B\ \cdots]$ — block **columns**, width $q$ | $n\times qN$ |
| $x_0=\mathcal{C}\,u_{\text{past}}$ | $n\times 1$ |
| $\mathcal{O}=[\,C;\ CA;\ CA^{2};\ \cdots]$ — block **rows**, height $p$ | $pN\times n$ |
| $y_{\text{future}}=\mathcal{O}\,x_0$ | $pN\times 1$ |
| $H=\mathcal{O}\mathcal{C}$ — blocks $h_{i+j-1}$, each $p\times q$ | $pN\times qN$ |

The inner dimensions cancel:
$\mathcal{O}_{(pN\times n)}\,\mathcal{C}_{(n\times qN)}=H_{(pN\times qN)}$, and within each
term $x_0=\sum_j (A^{j-1}B)_{(n\times q)}\,(u_{-j})_{(q\times1)}$ the $q$ contracts away.
The decisive observation: the blocked, **wide** dimensions ($qN$, $pN$) sit on the
*outside*, while the dimension the product contracts over — the bottleneck — is the
**un-blocked scalar $n$**. Hence

$$\operatorname{rank}H = \operatorname{rank}(\mathcal{O}\mathcal{C}) \le n$$

no matter how many inputs/outputs you stack: $p$ and $q$ only widen the input/output
*ports*, while everything still funnels through the same $n$-dimensional state. (So the
block-Hankel SVD has at most $n$ nonzero singular values regardless of $p,q$ — which is
what makes the order-$n$ truncation meaningful.)

> **This is *not* the transfer function.** The full input–output map — convolution over
> *all* time, $y_t=\sum_k h_k u_{t-k}$ — is the **Toeplitz** operator $T_{ts}=h_{t-s}$,
> constant along *diagonals*, whose symbol on the circle is the transfer function
> $H(z)=\sum_k h_k z^k$. The **Hankel** operator is the **past→future corner** of that
> same Toeplitz operator (cut time at $0$; keep the block mapping strictly-past inputs
> to non-negative-time outputs). Causality is what turns its constant *diagonals* into
> constant *anti-diagonals*. The Toeplitz operator is the whole behaviour; the Hankel
> corner isolates **memory** — how the past determines the future — the finite-rank
> ($=$ state) part that model reduction must preserve.

**Why the right factor is the controllability side (operator view).** Since $H$ is the
matrix of that past-inputs-to-future-outputs operator, the SVD orientation does the rest.
For **any** operator's SVD $H=U\Sigma V^{*}$, the columns of $V$ are the singular
directions in the **input space** and the columns of $U$ those in the **output space**.
The Hankel operator's input space is "past inputs" (the controllability side) and its
output space is "future outputs" (the observability side), hence

$$V \longleftrightarrow \text{controllability (input directions)},\qquad
U \longleftrightarrow \text{observability (output directions)}.$$

This is the same fact that put the **output** directions in $U$ and the **input**
directions in $V$ in the energy picture (`concepts/ModesAndHankel`): right singular
vectors are always input directions, and the Hankel operator's inputs *are* the
system's past inputs.

> **Why right singular vectors are input directions (and left, output).** This is the
> meaning of the SVD, not a convention. View $H=U\Sigma V^{*}$ as a map
> $\mathbb{C}^{n}\!\to\!\mathbb{C}^{m}$ (input space $\to$ output space). The SVD says
> $H$ acts **diagonally** between two orthonormal bases — one per space:
>
> $$H\,v_k = \sigma_k\,u_k \qquad(\text{feed in } v_k\ \Rightarrow\ \text{out comes } \sigma_k u_k).$$
>
> That is, $H$ sends each input axis $v_k$ straight onto one output axis $u_k$, scaled
> by $\sigma_k$. It is immediate from $V^{*}v_k=e_k$ (orthonormal columns of $V$):
>
> $$H v_k = U\Sigma V^{*}v_k = U\Sigma e_k = U(\sigma_k e_k) = \sigma_k\,u_k .$$
>
> Three equivalent ways to read off which factor is which:
> - **Order of operations** ($H=U\Sigma V^{*}$, right-to-left): $V^{*}$ reads the
>   **input** coordinates $\langle v_k,x\rangle$, $\Sigma$ scales by $\sigma_k$, $U$
>   assembles the **output** along $u_k$. So $V$ is input-side, $U$ output-side.
> - **Which space each lives on**: the $v_k$ are the eigenvectors of $H^{*}H$ — an
>   $n\times n$ operator on the **input** space — while the $u_k$ are the eigenvectors
>   of $HH^{*}$, an $m\times m$ operator on the **output** space
>   ($H^{*}H=V\Sigma^{*}\Sigma V^{*}$, $HH^{*}=U\Sigma\Sigma^{*}U^{*}$).
> - **Shapes**: $V$ is $n\times n$, so its columns are vectors in the input space
>   $\mathbb{C}^{n}$; $U$ is $m\times m$, columns in the output space $\mathbb{C}^{m}$.
>
> ($\sigma_k=\lVert Hv_k\rVert$ is then literally the gain $H$ applies to input
> direction $v_k$ — the same energy/gain reading as in `concepts/ModesAndHankel`.)

**From "spans" to "equals" — why we may *set* $\mathcal{C}=\widetilde V^{*}$.** The
operator view gives only "$V$ spans the controllability *subspace*." Taking
$\widetilde V^{*}$ to literally *be* a controllability matrix uses the freedom in the
factorization:

1. For a true order-$n$ system, $H$ has rank $n$, so the rows of $\widetilde V^{*}$ (top
   $n$ right singular vectors) form an **orthonormal basis of the row space of $H$** =
   the controllability subspace. Hence $\widetilde V^{*}=T\,\mathcal{C}$ for some
   invertible $n\times n$ change of basis $T$.
2. But $T\mathcal{C}=\begin{bmatrix}TB & TA\,B & \cdots\end{bmatrix}
   =\begin{bmatrix}\tilde B & \tilde A\tilde B & \cdots\end{bmatrix}$ with
   $\tilde A=TAT^{-1},\ \tilde B=TB$ — *itself* the controllability matrix of the
   **similarity-transformed realization** $(\tilde A,\tilde B,\tilde C)$. A similarity
   transform leaves the transfer function unchanged, so choosing $\mathcal{C}=\widetilde
   V^{*}$ is just **choosing coordinates**: the realization whose controllability matrix
   is $\widetilde V^{*}$. We then *define* $A,B$ from it (the shift below, and $B=$ first
   block), and by construction $\widetilde V^{*}=[B\ AB\ \cdots]$.

**That realization is input-balanced.** Because $\widetilde V$ has orthonormal columns,
$\widetilde V^{*}\widetilde V=I_n$, so the (finite-horizon) reachability Grammian is

$$\mathcal{C}\mathcal{C}^{*} = \widetilde V^{*}\widetilde V
= BB^{*}+ABB^{*}A^{*}+\cdots+A^{k-1}BB^{*}(A^{*})^{k-1} = I .$$

Reachability Grammian $=I$ **is** the definition of input-balanced — this is the precise
content of "the rows of $\widetilde V^{*}$ are orthonormal $\Rightarrow$ the realization
is input-balanced." And it is exactly the fact that makes the shift recover $A$ in the
next step (the cross-terms in eqs 444–447 collapse to $A\cdot I=A$).

> **Coordinate choice vs. reduction — what is and isn't invariant.** Two distinct moves
> happen inside this SVD, with *opposite* effects on the poles:
> - **Choosing coordinates** (the input-balanced form here, vs internally balanced, vs
>   TIB, …) is a *similarity* $A\mapsto TAT^{-1}$. It leaves the characteristic
>   polynomial — hence $\operatorname{eig}(A)$ — **unchanged**: the eigenvalues of $A$
>   *are* the poles of the transfer function, an intrinsic, coordinate-free property.
>   Whichever factorization we pick, the recovered poles are identical.
> - **Reduction** (truncating the SVD to rank $r<n$) is a *projection*, not a similarity,
>   and it **does** move the poles — it keeps $r$ *approximate* ones. That truncation,
>   not the coordinate choice, is what sets the reduced poles.
>
> So the poles read off in §4 are invariant to the input-balanced choice made here, yet
> are determined by the order-$r$ truncation upstream. (It is also why the orthogonal
> Schur *similarity* in §4 reads poles off $\operatorname{diag}A_1$ without changing them
> — `concepts/ModesAndHankel`, "realization-invariant".)

**Recovering $A$ by shift-invariance (eqs 444–447).** Write $\widetilde V^{*}$ in block
columns of width $q$. Let $\widetilde V^{*}_{\rightarrow}$ be $\widetilde V^{*}$ with the
**first** block-column dropped ($[AB\ A^2B\ \cdots]$) and $\widetilde V^{*}_{\leftarrow}$
with the **last** dropped ($[B\ AB\ \cdots]$). Then

$$\widetilde V^{*}_{\rightarrow} = A\,\widetilde V^{*}_{\leftarrow}\quad\Longrightarrow\quad A \approx \widetilde V^{*}_{\rightarrow}\,\bigl(\widetilde V^{*}_{\leftarrow}\bigr)^{+},$$

where $(\cdot)^{+}$ is the Moore–Penrose pseudo-inverse. Yu's eqs 444–447 derive this
explicitly: $\widetilde V^{*}_{\rightarrow}\widetilde V_{\leftarrow}
= A\bigl(BB^{*}+ABB^{*}A^{*}+\cdots\bigr)\approx A$, the approximation holding **exactly
because the reachability Grammian $\approx I$**. So the right inverse of
$\widetilde V^{*}_{\leftarrow}$ is its conjugate-transpose *only in the limit of a perfect
input-balanced factor*; for the truncated factor it is the pseudo-inverse.

> **Implementation-note — Bug 1 (shift uses transpose instead of pseudo-inverse).**
> The committed `msvdreduce` (HEAD) computes the shift as
> ```python
> A = Vh[:, di:].dot(Vh.T[:-di, :])          # uses the transpose
> ```
> i.e. $A=\widetilde V^{*}_{\rightarrow}\,(\widetilde V^{*}_{\leftarrow})^{\mathsf T}$.
> The truncated $\widetilde V^{*}_{\leftarrow}$ has **orthonormal rows, not orthonormal
> columns**, so its right inverse is the pseudo-inverse, not the transpose. With the
> transpose the recovered $A$ is wrong: on a 16-pole MIMO test the pole error is
> $8.75\times10^{-2}$ (with spurious complex eigenvalues, $\max|\operatorname{Im}|=0.16$).
> The working tree already corrects this to
> ```python
> A = Vh[:, di:].dot(np.linalg.pinv(Vh[:, :-di]))   # pseudo-inverse
> ```
> which drops the pole error to $2.37\times10^{-11}$. *(Status: this fix is present in
> the working tree but uncommitted; not changed by these notes.)*

**Recovering $B$ and $C$.** With $\widetilde V^{*}\approx[B\ AB\ \cdots]$, the input matrix
is the first block,

$$B = \widetilde V_{\,:,\,1:q},$$

(no extra $\Sigma^{1/2}$ scaling — Yu folds the singular values so that
$\widetilde V^{*}$ *is* the controllability factor). In the code, `B = Q.T.dot(Vh[:, :di])`
after the Schur rotation $Q$ of the next step. Finally $C$ is obtained by an $H_2$
least-squares fit in the reduced Krylov basis (eqs 454–455):

$$K = \begin{bmatrix}B & AB & \cdots & A^{k-1}B\end{bmatrix},\qquad C = \begin{bmatrix}h_1 & h_2 & \cdots & h_k\end{bmatrix}K^{*}.$$

(`msvdreduce` itself returns only the poles and null vectors; $C$ is fit downstream,
e.g. via `krylov_basis` + least squares.)

## 4. Schur step → poles

The shift realization $A$ is a full matrix. A **real Schur** factorization makes it
(block-)lower-triangular (Yu eq 451):

$$A = Q\,A_1\,Q^{\mathsf T},\qquad A_1\ \text{(block-)lower-triangular},$$

and $(A_1,\,Q^{\mathsf T}B)$ is an equivalent realization (impulse response is invariant
under the similarity $A\mapsto TAT^{-1},\,B\mapsto TB,\,C\mapsto CT^{-1}$, eqs 448–450,
here with the orthogonal $T=Q^{\mathsf T}$). The **reduced poles are the diagonal of
$A_1$**:

$$\omega_k = (A_1)_{kk}.$$

In `hankel.py` this is `lschur` (Schur + reversal to lower-triangular), and `w =
np.diag(A)`. Yu notes explicitly (eq 451) that **$(A_1,Q^{\mathsf T}B)$ is not in general
a TIB pair** — turning it into one (recovering the null vectors $y_k$) is the
tangential-Schur step, carried out by the `u`-loop. That recovery, and the second
bug in it, is the subject of [`03-null-vector-recovery.md`](03-null-vector-recovery.md).

## 5. The algorithm, mapped to the code

`ga/reducers/hankel.py: msvdreduce(ir, order)`:

| Step | Math | Code |
|---|---|---|
| Block Hankel | $H$ (eq 440) | `H = blkhankel(ir)` |
| (Partial) SVD | $H\approx\widetilde U\widetilde S\widetilde V^{*}$ (442) | `U, S, Vh = la.svd(H)`; `Vh = Vh[:order]` |
| Shift → $A$ | $A=\widetilde V^{*}_{\rightarrow}(\widetilde V^{*}_{\leftarrow})^{+}$ (444–447) | `A = Vh[:,di:].dot(np.linalg.pinv(Vh[:,:-di]))` — **Bug 1** |
| Schur → poles | $A=QA_1Q^{\mathsf T}$, $\omega_k=\operatorname{diag}A_1$ (451) | `A, Q = lschur(A)`; `w = np.diag(A)` |
| Input matrix | $B=\widetilde V_{:,1:q}$ (rotated) | `B = Q.T.dot(Vh[:, :di])` |
| Null vectors | tangential-Schur deflation (452–453) | the `u`-loop — **Bug 2**, see file 03 |

> **Implementation-note — the SISO siblings (`reduce_svd_truncate`, `info_svd_reduce`)
> share Bug 1, and add more.** `ga/reducers/hankel.py` also has SISO reducers that are
> an earlier, rougher draft of this same shift-realization idea. They carry **Bug 1** —
> the shift uses `v[:order,1:] @ v.T[:-1,:order]`, a transpose where the pseudo-inverse
> is required (on a proper Hankel, $\sim\!10^{-1}$ pole error vs $\sim\!10^{-9}$ with
> `pinv`) — plus two more faults:
>
> - **`reduce_svd_truncate`** builds `la.hankel(h)` with a *single* argument, which is a
>   **triangular** Hankel (lower-right triangle zeroed); it silently discards the late
>   impulse-response samples. Harmless for a fast-decaying IR, but $\approx 0.3$ pole
>   error when the tail carries energy. Needs a properly populated square Hankel
>   `la.hankel(h[:m], h[m-1:2m-1])` **and** the pseudo-inverse.
> - **`info_svd_reduce`** (the cepstrum / information-distance path,
>   `concepts/InformationGeometry`) additionally (i) computes a `rho`-damped IR `hr` but
>   then runs on the *un-damped* `h_` — the damping is silently dropped; and (ii) runs
>   the pole-finder **directly on the cepstrum** $a_k$. But for a rational system
>   $a_k=\tfrac1k\sum_i\lambda_i^k$, and the $1/k$ makes $a_k$ **not** a low-rank-Hankel
>   sequence, so Ho–Kalman cannot recover the poles from it. One must first undo the
>   $1/k$ to the **power sums** $g_k=k\,a_k=\sum_i\lambda_i^k$ — a sum of geometrics,
>   Prony-recoverable. With the damping applied and this $k$-multiplication, full-order
>   recovery returns to $\sim\!10^{-8}$.
> - **Caveat even after fixing:** the power-sum Hankel has **no energy ordering** (every
>   pole enters with residue 1), so SVD-truncating it for *aggressive* order reduction is
>   unstable (can return out-of-disk poles). This path is sound for **identification /
>   full-order recovery**, not deep reduction — for that, use the energy-ranked
>   impulse-response `msvdreduce`.
>
> Verified corrected companions live in `ga/reducers/hankel.py` (originals left intact):
> `reduce_svd_truncate_fixed` (full-order recovery $\sim\!9\times10^{-10}$) and
> `info_svd_reduce_fixed` ($\sim\!10^{-8}$). `reduce_fft_truncate` is a separate FFT/QR
> path, not audited here.

## 6. Separation from Balanced Truncation

Balanced Truncation (BT) and `msvdreduce` are both "model reduction" and both rank
modes by the **Hankel singular values**, but they are otherwise different methods
operating on different inputs. BT is the subject of `Mu2026-ModelReductionNotes` and
`ga/reducers/balanced_truncation.py`; `msvdreduce` is the subject of this file.

> **BT vs. `msvdreduce` (Hankel-SVD).**
>
> | | Balanced Truncation (Yu §6.1; Mu notes) | `msvdreduce` / Hankel-SVD (Yu §6.2) |
> |---|---|---|
> | **Input** | state space $(A,B,C)$ | impulse response $\{h_k\}$ only |
> | **Core object** | Grammians $P=APA^{*}+BB^{*}$, $Q=A^{*}QA+C^{*}C$ (eqs 402–403) | data block Hankel $H=\mathcal{O}\mathcal{C}$ (eq 440) |
> | **What is ranked** | Hankel singular values $\sigma_k=\sqrt{\lambda_k(PQ)}$ (eq 404) | Hankel singular values $\sigma_k$ of $H$ (eq 442) — *same quantity* |
> | **Mechanism** | balance (square-root: $Z_P^{*}Z_Q=U\Sigma V^{*}$, $W=Z_PU_1\Sigma_1^{-1/2}$, $V=Z_QV_1\Sigma_1^{-1/2}$), then **truncate**: $A_r=W^{*}AV$ | **realize** from $H$ (shift), then Schur-triangularize |
> | **Error control** | $H_\infty$ bound $\lVert G-G_r\rVert_\infty\le 2\sum_{i>r}\sigma_i$ (eq 405) | $H_2$ / tangential interpolation (Walsh's theorem, Thm 40; MIMO Thm 41 / Lemma 42) |
> | **Reduced realization** | projection $W_1^{\mathsf T}V_1=I$, $A_r=W_1^{\mathsf T}AV_1$ (eqs 411–414); balanced but not TIB | input-balanced shift realization; not yet TIB |
> | **Code** | `balanced_truncation()` | `msvdreduce()` |
>
> The shared object is the **Hankel singular value spectrum** — both methods learn it
> (BT from the Grammian product $PQ$, `msvdreduce` from the data Hankel $H$) and both
> truncate on it. They differ in what they *start from*: BT needs a state-space model
> already in hand; `msvdreduce` is data-driven and builds the state space from the
> impulse response. For market-impact work where the IR is estimated directly from
> data, `msvdreduce` is the natural entry point; BT applies once a model exists and
> carries the provable $H_\infty$ bound.

**Error-norm caveat.** BT targets $H_\infty$ (worst-case frequency gain) with a hard
$2\sum\sigma_i$ bound; `msvdreduce`/Hankel-SVD targets $H_2$ / Hankel norm via the
interpolation conditions of §6.1 (Walsh's theorem: the $H_2$-optimal reduced model
interpolates $H$ at the reflected poles $1/\bar z_k$, eqs 417–419; the MIMO
*tangential* version is Thm 41 / Lemma 42). They generally produce different reduced
models even at the same order. See `concepts/ModelReduction` for the wiki summary.

## 7. The Schur step in detail: `B = Qᵀ B`, and why it isn't yet TIB

This section proves the two facts §4 states without derivation — that $B$ transforms as
$Q^{\mathsf T}B$ under the Schur step, and that the resulting $(A_1,\,Q^{\mathsf T}B)$ is
not yet a TIB pair — and closes with the recovery (full detail in
[`03-null-vector-recovery.md`](03-null-vector-recovery.md)).

**The shift $A$ is dense; Schur triangularizes it.** The shift estimate
$A=\widetilde V^{*}_{\rightarrow}(\widetilde V^{*}_{\leftarrow})^{+}$ is just *some*
representative of the input-balanced coordinate class: it has the right eigenvalues (the
poles) but no imposed sparsity — a full, dense $n\times n$ matrix. TIB needs $A$
(block-)lower-triangular with the poles on the diagonal, so we apply a **real Schur**
factorization $A=QA_1Q^{\mathsf T}$. Real Schur stays real: real poles become $1\times1$
diagonal entries, complex-conjugate pairs become $2\times2$ real diagonal blocks
(block-lower-triangular).

**Why $B$ transforms as $Q^{\mathsf T}B$.** The Schur step relabels the *state* by the
orthogonal $Q$, and $B$ follows from that change of coordinates. With the realization
$x_{t+1}=Ax_t+Bu_t,\ y_t=Cx_t$, define new coordinates $\tilde x_t=Q^{\mathsf T}x_t$ (so
$x_t=Q\tilde x_t$, using $Q^{-1}=Q^{\mathsf T}$) and substitute:

$$Q\,\tilde x_{t+1}=AQ\,\tilde x_t+Bu_t \;\Longrightarrow\; \tilde x_{t+1}=\underbrace{Q^{\mathsf T}AQ}_{A_1}\,\tilde x_t+\underbrace{Q^{\mathsf T}B}_{\tilde B}\,u_t,\qquad y_t=\underbrace{CQ}_{\tilde C}\,\tilde x_t .$$

So $B\mapsto Q^{\mathsf T}B$ because that is how the input enters the *new* state
coordinates — the general rule "$\tilde x=Tx \Rightarrow (A,B,C)\mapsto(TAT^{-1},TB,CT^{-1})$"
with $T=Q^{\mathsf T}$. It is the **unique** transform that preserves the impulse
response, which telescopes ($Q^{\mathsf T}Q=I$):

$$\tilde C\tilde A^{\,k-1}\tilde B=(CQ)(Q^{\mathsf T}AQ)^{k-1}(Q^{\mathsf T}B)=CQ\,(Q^{\mathsf T}A^{k-1}Q)\,Q^{\mathsf T}B=CA^{k-1}B=h_k .$$

Because $Q$ is **orthogonal**, this similarity also preserves the reachability Grammian,
$\tilde P=Q^{\mathsf T}PQ=Q^{\mathsf T}IQ=I$, so input-balance survives the Schur step.

**Why $(A_1,Q^{\mathsf T}B)$ is not yet a TIB pair.** After Schur, $A_1$ is lower-
triangular with the poles on its diagonal and input-balance is preserved — so what is
missing? TIB is **not** merely "lower-triangular + input-balanced"; it is the canonical
form parameterized by (poles, *unit null vectors* $\{y_k\}$), and that imposes a coupling
between $A$ and $B$ which Schur does not produce. In the null-basis construction
([`02-null-basis-realization.md`](02-null-basis-realization.md)),

$$\text{row } k \text{ of } A:\ \ a_k=t_k\,y_k^{*}L_k,\qquad \text{row } k \text{ of } B:\ \ b_k=t_k\,y_k^{*}P_k,$$

the **same** $y_k$ appears in both. A generic $Q^{\mathsf T}B$ violates this. For mode 2,
TIB forces $A_{21}=t_1t_2\,(y_2^{*}y_1)$ and $b_2=t_2\,y_2^{*}J_1$ to share one unit
$y_2$: solving $y_2^{*}=b_2J_1^{-1}/t_2$ from the $B$-row generically gives
$\lVert y_2\rVert\ne1$ and an $A_{21}$ that disagrees. Schur triangularized the **state**
($A$) but did nothing to express the **input directions** through the null-vector
recursion, nor did it list the $y_k$. So $(A_1,Q^{\mathsf T}B)$ hands you the **poles**
(the diagonal) but neither the null-vector parameters nor the $A$–$B$ coupling that "TIB
pair" means. (And the input-balance inherited from a *truncated* SVD is only approximate,
so even the balance needs restoring.) Hence Yu's "not in general a TIB pair."

**Recovering the TIB pair (tangential-Schur).** The fix reads the null vectors off $B$
and rebuilds the canonical pair. Since $b_k=t_k\,y_k^{*}P_k$ with
$P_k=J_{k-1}\cdots J_1$, the first row gives $y_1$ directly and each later row is the next
$y_k$ pre-multiplied by the cumulative deflation; peel one $J_i$ at a time (Yu eqs
452–453):

$$u_i=\operatorname{normalize}(b_1^{*}),\qquad B\leftarrow B_{2:}\Big(I-\big(1+\tfrac{1}{\bar\lambda_i}\big)u_iu_i^{*}\Big),\quad\text{drop row 1,}$$

where $I-(1+1/\bar\lambda_i)u_iu_i^{*}=J_i^{-1}$ is the Sherman–Morrison inverse of the
forward deflation $J_i$. This recovers all $\{y_k\}$; the strict TIB pair is then built
**directly** from $(\text{poles},\{y_k\})$ by the forward recursion *with the $y_k$ as the
final orthonormal directions* (Yu eqs 454–455) — **not** by re-running
`null_basis_realization` on them, which would re-apply the Blaschke deflation to vectors
that are *already* deflated (double-deflation) and realize the **wrong** transfer
function. Built directly, the pair is exactly input-balanced and reproduces the system;
$C$ is fit by least squares. (Verified at full order on a known MIMO TIB system: poles to
$2\times10^{-15}$, the direct rebuild to $\sim\!3\times10^{-9}$ in IR — whereas the naive
re-deflation gives $\sim\!0.13$.) The conditioning ($J_i^{-1}$ blows up as
$\lambda_i\to0$), the code's **Bug 2** (the deflation is computed but discarded, with the
wrong coefficient $1+\bar\lambda$ instead of $1+1/\bar\lambda$), the **do-not-re-deflate**
rebuild, and the robust alternative `tib_from_state_space` (input-balance + orthogonal
Schur + band factorization, no explicit null vectors) are covered in
[`03-null-vector-recovery.md`](03-null-vector-recovery.md).

## Related

- [`02-null-basis-realization.md`](02-null-basis-realization.md) — the forward map this inverts.
- [`03-null-vector-recovery.md`](03-null-vector-recovery.md) — the `u`-loop (Bug 2) and the BT-side counterpart.
- `concepts/ModesAndHankel`, `concepts/ModelReduction`.
