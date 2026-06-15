# MIMO AR ↔ IR (`ar2ir` for the multivariate case)

The full multivariate derivation of `ar2ir` — converting matrix **AR coefficients** to the
**impulse response** — filling the **"MIMO: pass"** placeholder in Yu Mu's TIB note
(`raw/tib_note.pdf`, "Functions of Toeplitz Matrices"; [[wiki/papers/Mu2026-TIBNote]]). It is
the math behind `ga/filters/tib.py: mimo_ar2ir`, and the bridge in
[[concepts/VARSelfPrediction]] from `CAⁿB` (AR coefficients) to the forward IR you feed to
model reduction. Everything reduces to the note's scalar derivation when `d = 1`.

## Notation

| symbol | meaning |
|---|---|
| `d` | number of series (inputs = outputs in the self-prediction VAR) |
| `y_t, x_t ∈ ℝ^d` | data (output) and innovation (input) at time `t` |
| `h_k ∈ ℝ^{d×d}` | matrix impulse response / MA (Markov) coefficients, `h_0 = I` |
| `A_k ∈ ℝ^{d×d}` | matrix AR coefficients (`k = 1..p`), e.g. `A_k = CA^{k-1}B` |
| `H(z) = Σ_{k≥0} h_k z^{-k}` | MA / IR transfer function |
| `𝒜(z) = I − Σ_{k≥1} A_k z^{-k}` | AR / whitening polynomial (matrix-valued) |

## 1. The two processes (block-Toeplitz form)

**MA / forward (IR) process** — output is a causal convolution of the input with the IR:

$$y_t = \sum_{k\ge 0} h_k\, x_{t-k},\qquad\text{i.e.}\qquad Y = T_H X,$$

with `T_H` the **block lower-triangular Toeplitz** matrix `(T_H)_{ij} = h_{i-j}` (`= 0` for `i<j`):

$$T_H=\begin{bmatrix} h_0 & & & \\ h_1 & h_0 & & \\ h_2 & h_1 & h_0 & \\ \vdots & & & \ddots\end{bmatrix}.$$

**AR / inverse (whitening) process** — innovation is the data minus its prediction:

$$y_t = \sum_{k=1}^{p} A_k\, y_{t-k} + x_t \;\Longleftrightarrow\; x_t = y_t - \sum_{k=1}^{p} A_k\, y_{t-k},\qquad X = T_{\mathcal A} Y,$$

with `T_𝒜` block lower-triangular Toeplitz, first block-column `(I, −A_1, …, −A_p, 0, …)`, i.e.
`(T_𝒜)_{ij} = 𝒜_{i-j}` with `𝒜_0 = I`, `𝒜_k = −A_k`:

$$T_{\mathcal A}=\begin{bmatrix} I & & & \\ -A_1 & I & & \\ -A_2 & -A_1 & I & \\ \vdots & & & \ddots\end{bmatrix}.$$

This is the block ("MIMO") version of the scalar `T_h`, `T_a` in the note; `A_k` multiplies
`y_{t-k}` **on the left**, the only place the matrix order matters.

## 2. The AR ↔ IR duality `T_H T_𝒜 = I`

Substituting `X = T_𝒜 Y` into `Y = T_H X` gives `Y = T_H T_𝒜 Y` for every `Y`, hence

$$\boxed{\,T_H\,T_{\mathcal A} = I\,}$$

the MIMO analog of the note's `T_h·T_a = I`. Block-Toeplitz multiplication corresponds to
transfer-function multiplication, so equivalently

$$H(z)\,\mathcal A(z) = I \quad\Longleftrightarrow\quad H(z) = \mathcal A(z)^{-1}.$$

**The forward IR is the inverse of the AR/whitening polynomial.** (Both `𝒜H = I` and `H𝒜 = I`
hold because `𝒜(z)` is square and invertible — used in §3.)

## 3. The recursion (`mimo_ar2ir`)

Equate coefficients of `z^{-m}` in `𝒜(z) H(z) = I`:

$$\Big(I - \sum_{k\ge1} A_k z^{-k}\Big)\Big(\sum_{n\ge0} h_n z^{-n}\Big) = I
\;\Longrightarrow\; h_m - \sum_{k=1}^{m} A_k\, h_{m-k} = \delta_{m0}\,I,$$

giving the **block forward-substitution recursion**

$$\boxed{\,h_0 = I,\qquad h_m = \sum_{k=1}^{m} A_k\, h_{m-k}\quad(A_k=0\text{ for }k>p)\,}$$

each `h_m` a `d×d` matrix. This is exactly `mimo_ar2ir`.

**Left vs. right.** Equating coefficients in the *other* product `H(z)𝒜(z) = I` gives the
**right** recursion `h_m = Σ_k h_{m-k} A_k`. Because `H = 𝒜^{-1}` is unique, the two produce the
**same** `{h_m}` — verified numerically to `2.8×10⁻¹⁷` (§5). The code uses the left form; either
is correct. (For `d = 1` they are literally identical since scalars commute.)

## 4. Block-Toeplitz solve form (the note's `T_a h = e₁`, MIMO)

Stacking `h = (h_0; h_1; …; h_{m-1})` and `E_1 = (I; 0; …; 0)`, the recursion is the block
forward-substitution of the **block lower-triangular Toeplitz system**

$$T_{\mathcal A}\, h = E_1,\qquad
T_{\mathcal A}=\begin{bmatrix} I \\ -A_1 & I \\ -A_2 & -A_1 & I \\ \vdots & & & \ddots\end{bmatrix},\quad
E_1=\begin{bmatrix} I\\0\\0\\\vdots\end{bmatrix}.$$

This is the block analog of the note's final scalar equation `T_a h = e₁`. The scalar version
uses `scipy.linalg.solve_toeplitz`; the block version is solved by the substitution recursion
of §3 (SciPy's Toeplitz solver is scalar-only). Same equation either way.

## 5. Validation (re-runnable)

With random `A_k` (`d=3`, `p=3`, scaled for stability), `m=40`:

| check | what it confirms | result |
|---|---|---|
| `𝒜(z)·H(z) = I` (conv of `[I,−A_1,…]` with `h`, left) | `h` inverts `𝒜` | `5.6e-17` |
| `H(z)·𝒜(z) = I` (right) | both one-sided inverses agree | `5.6e-17` |
| left recursion `== ` right recursion | order-of-multiplication is immaterial | `2.8e-17` |
| `d=1` vs scalar `h_0=1, h_n=Σ a_k h_{n-k}` | reduces to the note's `ar2ir` | exact |
| `mimo_ar2ir` `==` VAR companion IR `C_f comp^n B_f` | IR carries the **data poles** `eig(comp)` | exact |

The last row connects to [[concepts/VARSelfPrediction]]: the data poles are
`eig(companion) = ` roots of `det 𝒜(z) = 0 = ` poles of `H = 𝒜^{-1}`, so reducing `{h_n}`
recovers them. Reproduce with:

```python
import numpy as np, scipy.linalg as la
from ga.filters.tib import mimo_ar2ir
rng = np.random.default_rng(1); d, p, m = 3, 3, 40
A = rng.standard_normal((p, d, d)) * 0.25
h = mimo_ar2ir(A, m)
Acal = [np.eye(d)] + [-A[k] for k in range(p)]
err = max(np.max(np.abs(sum(Acal[j] @ h[n-j] for j in range(min(n, p)+1))
                         - (np.eye(d) if n == 0 else 0))) for n in range(m))
print("A(z)H(z)=I err:", err)            # ~1e-16
```

## 6. Pipeline & code

`ga/filters/tib.py: mimo_ar2ir(ar_coeffs, m)` — `ar_coeffs` shape `(p, d, d)`, returns
`h` shape `(m, d, d)`. In context:

```
CAⁿB (square AR coeffs) → mimo_ar2ir → forward IR {h_n} → reduce → data poles (= poles of 𝒜⁻¹)
```

## Related

- [[wiki/papers/Mu2026-TIBNote]] — the source note (scalar derivation; "MIMO: pass").
- [[concepts/VARSelfPrediction]] — why `CAⁿB` are AR coefficients and the forward IR is `𝒜⁻¹`.
- [`01-hankel-svd-reduction.md`](01-hankel-svd-reduction.md) — the reduction this feeds.
