# Recovering the Null Structure from a Reduced Model: `msvdreduce` vs `tib_from_state_space`

> Two routes from a reduced realization $(A,B)$ to TIB form. One (`msvdreduce`, Yu
> eqs 451–453) recovers the orthonormal input directions $y_k$ **explicitly** by
> un-deflating the $B$ rows; the other (`tib_from_state_space`, Yu Ch 5 eqs 394–401)
> re-coordinates into the band-fraction TIB form **implicitly** and never names the
> null vectors. Both are the same first move — Schur-triangularize $A$ — followed by
> different second moves. Notation: [`00-overview.md`](00-overview.md). The forward
> map being inverted is [`02-null-basis-realization.md`](02-null-basis-realization.md).

## 1. The shared inverse problem

A reduced realization $(A,B)$ (from `msvdreduce`'s shift step, or from
`balanced_truncation`) is balanced but **not** TIB. Yu states the situation exactly
(eq 451): after the real Schur factorization

$$A = Q\,A_1\,Q^{\mathsf T},\qquad A_1\ \text{(block-)lower-triangular},$$

the pair $(A_1,\,Q^{\mathsf T}B)$ is an equivalent realization with the **poles on the
diagonal of $A_1$**, but it "is not in general a TIB pair." Both methods below begin
here; they differ in how they finish the conversion to TIB.

## 2. Route A — `msvdreduce`: explicit null-vector recovery by un-deflation

Recall (file 02) that a TIB $B$ has rows

$$b_k = t_k\,y_k^{*}P_k,\qquad P_k = J_{k-1}\cdots J_1,\qquad J_i = I-(1+\bar\lambda_i)y_iy_i^{*}.$$

Row 1 gives $y_1$ directly; later rows carry the cumulative deflation $P_k$. To peel
$y_k$ off row $k$ we undo one $J_i$ at a time. The inverse is a Sherman–Morrison
rank-1 update (using $y_i^{*}y_i=1$):

$$J_i^{-1} = \big(I-(1+\bar\lambda_i)y_iy_i^{*}\big)^{-1}
= I + \frac{1+\bar\lambda_i}{1-(1+\bar\lambda_i)}\,y_iy_i^{*}
= I - \Big(1+\tfrac{1}{\bar\lambda_i}\Big)\,y_iy_i^{*}.$$

So the recovery loop is (Yu eqs 452–453): take the current first row as the next
direction, then right-multiply the remaining rows by $J_i^{-1}$ to strip its factor:

$$
u_i = \operatorname{normalize}\big(b_1^{*}\big),\qquad
B \leftarrow B_{2:}\,\Big(I - \big(1+\tfrac{1}{\bar\lambda_i}\big)u_iu_i^{*}\Big),\quad
\text{drop row 1.}
$$

Then $b_k J_i^{-1} = t_k y_k^{*}J_{k-1}\cdots J_{i+1}$ — one factor gone — so after the
$i$-th pass the new first row is $t_{i+1}y_{i+1}^{*}$ and $u_{i+1}=y_{i+1}$ exactly.
This is the inverse Gram–Schmidt pass dual to the forward Blaschke deflation: the
forward map deflates with $(1+\bar\lambda_i)$, the inverse un-deflates with
$(1+1/\bar\lambda_i)$.

**Verification.** On a TIB system built by `null_basis_realization` (order 6, $q=3$),
recovering the directions from $B$ and comparing to the forward $y_k$ gives
$\max_k\big(1-|\langle y_k,u_k\rangle|\big)$:

| recovery | result |
|---|---|
| no deflation (current code behaviour) | $6.1\times10^{-1}$ |
| deflation with the code's coefficient $(1+\bar\lambda)$ | $6.4\times10^{-1}$ |
| deflation with Yu's coefficient $(1+1/\bar\lambda)$ | $0.0$ (exact) |

(The first direction matches in all cases; only $k\ge2$ depend on the deflation.)

**Rebuilding the TIB pair — build directly, do not re-deflate.** The recovered $u_i$ are
the *already-orthonormalized* directions $y_k$. To rebuild the canonical TIB $(A,B)$ you
run the forward recursion of [`02-null-basis-realization.md`](02-null-basis-realization.md)
**with the $y_k$ used as the final orthonormal directions** — i.e. skip the Blaschke
orthogonalization step (Yu eqs 454–455). You must **not** feed the recovered $y_k$ into
`null_basis_realization`, because that function treats its input as *raw* $v_k$ and
re-applies the Blaschke–Potapov deflation: deflating an already-deflated vector
(double-deflation) yields a pair that is still input-balanced but realizes the **wrong**
transfer function. This is the same "orthonormalized $y$ vs. raw $v$" trap as the
historical BT-extraction gap (see [[bt-to-tib-extraction-gap]] / `concepts/TIBForm`).
Verified at full order on a known MIMO TIB system ($q=2,p=3,n=6$):

| rebuild from recovered $(\text{poles}, y_k)$ | IR rel err |
|---|---|
| **directly** (forward recursion, $y_k$ as final directions) | $3\times10^{-9}$ ✓ |
| via `null_basis_realization` (re-deflates) | $0.13$ ✗ |
| (for comparison) bare equivalent realization $(A_1, Q^{\mathsf T}B)$ + least-squares $C$ | $2\times10^{-15}$ |

The last row is the reminder that if you only need *a* realization of the transfer
function — not the canonical TIB form — the Schur pair $(A_1,Q^{\mathsf T}B)$ with a fitted
$C$ already reproduces it to machine precision; the null-vector recovery is only needed to
land in TIB coordinates.

> **Implementation-note — Bug 2 (deflation computed, then discarded; wrong coefficient).**
> The `u`-loop in `ga/reducers/hankel.py: msvdreduce` is:
> ```python
> for i in range(order):
>     ui = (B[0] / np.linalg.norm(B[0], ord=2)).reshape(-1, 1)
>     u[:, i] = ui
>     tmp = np.eye(di) - (1 + np.conj(w[i])) * ui.dot(ui.T)   # computed...
>     B = B[1:, :]                                            # ...never used; raw rows taken
> ```
> Two faults: **(a)** the deflation `tmp` is computed but never applied — `B = B[1:]`
> just drops the row and takes the *raw* remaining rows, so every returned $u_i$ for
> $i\ge2$ is an un-deflated row, not a mode direction (error $0.61$ above). **(b)**
> even the computed `tmp` uses $(1+\bar w_i)$, which is the *forward* $J_i$ coefficient,
> not the inverse $J_i^{-1}$ coefficient $(1+1/\bar w_i)$ (error $0.64$). The correct
> loop applies $J_i^{-1}$:
> ```python
> for i in range(order):
>     ui = (B[0] / np.linalg.norm(B[0], ord=2)).reshape(-1, 1)
>     u[:, i] = ui
>     B = B[1:, :] @ (np.eye(di) - (1 + 1/np.conj(w[i])) * ui.dot(ui.conj().T))
> ```
> which recovers the directions exactly (error $0.0$). *(Not changed by these notes;
> Bug 2 is present in both HEAD and the working tree.)*
>
> **Conditioning caveat.** $J_i^{-1}$ contains $1/\bar\lambda_i$, so the explicit
> un-deflation blows up as a pole approaches the origin ($\lambda_i\to0$). This is an
> intrinsic ill-conditioning of *explicit* null-vector recovery — and a concrete
> reason to prefer Route B when you only need a TIB realization, not the vectors.

## 3. Route B — `tib_from_state_space`: implicit re-coordination into band-fraction TIB

Yu's Ch 5 gives the structural target (Thm 38): every TIB pair has a **band-fraction**
representation $A=M^{-1}N$, $B=M^{-1}\!\begin{bmatrix}U\\0\end{bmatrix}$ with $M,N$
bidiagonal and $U$ unitary. For a real system with complex-conjugate poles, the
real construction (eqs 394–401) is: take an input-balanced real pair $(A,B)$, apply
the **real Schur** $A=QTQ^{*}$ ($T$ upper quasi-triangular, $2\times2$ blocks for
conjugate pairs), factor $T=LG$ ($L$ lower-triangular, $G$ block-unitary via the LQ
of each $2\times2$ block), and read off the **real TIB pair $(L,\,Q^{*}B)$**
(post-multiplying by the unitary $G$). The balance $AA^{*}+BB^{*}=I$ is preserved
because $Q$ and $G$ are orthogonal.

`ga/reducers/balanced_truncation.py: tib_from_state_space(A, B)` implements exactly
this, getting the input balance first (the reduced model isn't input-balanced) via a
Grammian similarity:

| Step | Math | Code |
|---|---|---|
| Reachability Grammian | $P = APA^{*}+BB^{*}$, $P=ZZ^{*}$ | `P = solve_lyapunov_discrete(A, B@B.T)`; `Z = cholesky(P)` |
| Input-balance | $A_1=Z^{-1}AZ,\ B_1=Z^{-1}B$ ⟹ $A_1A_1^{*}+B_1B_1^{*}=I$ | `A1 = la.solve(Z, A@Z)`; `B1 = la.solve(Z, B)` |
| Orthogonal real Schur | $A_1=UTU^{*}$, flip to lower: $A_{\text{tib}}=U_f^{*}A_1U_f$ | `_, U = la.schur(A1, 'real')`; `Uf = U[:, ::-1]` |
| TIB pair | $(A_{\text{tib}},\,U_f^{*}B_1)$ | `A_tib = Uf.T@A1@Uf`; `B_tib = Uf.T@B1` |
| Output map | $C_{\text{tib}} = C\,(ZU_f)$ | `transform = Z @ Uf` |

The result is exactly TIB, real, and reproduces the transfer function to $\sim10^{-14}$
(round-trip tests in `tests/test_bt_to_tib.py`; demo in
`scripts/compare_null_vectors.py`). Real poles appear as $1\times1$ diagonal entries,
complex-conjugate pairs as $2\times2$ real diagonal blocks (block-lower-triangular) —
the orthogonal real-Schur step delivers precisely the block structure that TIBForm's
$Q$ rotation describes. **No null vectors are ever extracted**; the input directions
stay implicit in the band-fraction $B$.

## 4. The connection

Both routes are the same two moves, differing only in the second:

```
reduced (A,B)  ──real Schur──▶  (A₁, QᵀB),  A₁ (block-)lower-tri,  poles = diag      [Yu eq 451]
                                        │  "not in general a TIB pair"
              ┌─────────────────────────┴─────────────────────────┐
   Route A (msvdreduce)                                  Route B (tib_from_state_space)
   peel explicit yₖ by un-deflation                      input-balance + orthogonal Schur
   Jᵢ⁻¹ = I − (1+1/λ̄ᵢ)yᵢyᵢ*  (Yu 452–453)               → band-fraction TIB (Yu Ch5, 394–401)
   → explicit (ωₖ, uₖ); ill-conditioned for small λ      → implicit TIB (A,B); robust (~1e-14)
```

The parallel extends to the *failure mode*. The buggy `msvdreduce` `u`-loop (take raw
$B$ rows, no un-deflation) is structurally the same mistake as
`extract_poles_and_nullvecs_from_bt` on the BT side — both return the **normalized
rows of the Schur-rotated $B$**, i.e. the un-deflated rows, which carry the wrong
input directions and so realize the wrong transfer function (the $0.53$-rel-err
failure documented in `concepts/TIBForm`). There are then two correct fixes:

- **Fix the explicit recovery** — apply $J_i^{-1}$ (Route A done right). Recovers the
  vectors, but inherits the $1/\bar\lambda_i$ conditioning.
- **Sidestep it** — re-coordinate directly into TIB (Route B, `tib_from_state_space`).
  You lose the explicit $(\omega_k,u_k)$ list but gain a numerically robust TIB pair
  that exactly reproduces the transfer function.

For genesis-alpha's fitting pipeline, Route B is the recommended path
(`concepts/TIBForm`, *Why It Matters* #3); Route A is the right tool only when the
explicit per-mode input directions are themselves the deliverable.

## Related

- [`01-hankel-svd-reduction.md`](01-hankel-svd-reduction.md) — where the `msvdreduce` $(A,B)$ comes from.
- [`02-null-basis-realization.md`](02-null-basis-realization.md) — the forward Blaschke deflation these invert.
- `concepts/TIBForm` (*Why It Matters* #3), `tests/test_bt_to_tib.py`, `scripts/compare_null_vectors.py`, `notebooks/model_reduction_validation.ipynb`.
