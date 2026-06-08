# Fast Hankel Matrix–Vector Multiplication (FFT) — Derivation

> Derives the $O(k\log k)$ FFT-based Hankel matvec behind `FastHankelProduct` in
> `ga/reducers/hankel.py`, maps the code to the math, flags the two bugs in
> `reduce_fft_truncate`, and extends to the block (MIMO) case that the fast
> `msvdreduce` needs. This is the multiplication that lets Lanczos (`eigsh`/`svds`)
> compute a partial Hankel SVD without ever forming the dense matrix — Yu §6.2 "Fast
> Partial Block Hankel SVD" (eqs 435–442). Notation: [`00-overview.md`](00-overview.md).

## 1. The problem

A (scalar) Hankel matrix $H\in\mathbb{R}^{k\times k}$ has constant anti-diagonals,

$$H_{ij}=h_{i+j},\qquad i,j=0,\dots,k-1,$$

so it is determined by the $2k-1$ numbers $h_0,\dots,h_{2k-2}$. The matvec is

$$(Hx)_i=\sum_{j=0}^{k-1} h_{i+j}\,x_j .$$

Done densely this is $O(k^2)$ per product. Lanczos methods (`eigsh` on $H$ or
$H^{*}H$, `svds` on $H$) need many matvecs, so an $O(k\log k)$ matvec is the whole game
for large $k$.

## 2. Reversal turns the Hankel matvec into a convolution

A Hankel matvec is a Toeplitz/convolution in disguise. Reverse the input,
$x'_m=x_{k-1-m}$, and substitute $m=k-1-j$:

$$(Hx)_i=\sum_{j=0}^{k-1}h_{i+j}x_j=\sum_{m=0}^{k-1}h_{\,i+k-1-m}\,x'_m=(h * x')_{\,i+k-1},$$

where $*$ is the ordinary **linear convolution**. So the Hankel matvec is a single
convolution of the data $h$ (length $2k-1$) with the reversed input $x'$ (length $k$),
read off at indices $i+k-1$ for $i=0,\dots,k-1$ — i.e. the middle slice
$[\,k-1:2k-1\,]$ of the convolution. (A Toeplitz matvec is the same with $x$ *not*
reversed; reversing the rows of $H$ swaps the two — Hankel $\leftrightarrow$ Toeplitz.)

## 3. FFT: convolution in $O(k\log k)$

Linear convolution is pointwise multiplication in the Fourier domain. Zero-pad both
sequences to a common length $N\ge \operatorname{len}(h)+\operatorname{len}(x')-1=3k-2$
(in practice the next power of two), then

$$h * x' = \mathrm{IFFT}\big(\mathrm{FFT}(h)\odot\mathrm{FFT}(x')\big),$$

with $\odot$ the elementwise product. Cost: three FFTs $+$ one multiply $=O(k\log k)$.
Take the slice $[k-1:2k-1]$ and you have $Hx$.

**Circulant view (equivalent).** Embedding $H$ into a circulant $C$ of size $N$ and
using $C=F^{-1}\operatorname{diag}(Fc)F$ (a circulant is diagonalized by the DFT) gives
the same thing: $Cv=\mathrm{IFFT}(\mathrm{FFT}(c)\odot\mathrm{FFT}(v))$. The
first-column vector $c$ of the circulant encodes the Hankel entries; padding $v$ with
zeros isolates the linear (non-circular) part of the convolution.

## 4. Mapping to `FastHankelProduct`

The code (`ga/reducers/hankel.py`) is the circulant form, with $N=2k-1$:

```python
chat = np.hstack([self.h[-1], np.zeros(self.n-1), self.h[:-1]])   # length 2k-1
Pi   = np.eye(self.n)[::-1]            # reversal
xhat = np.hstack([Pi.dot(x), np.zeros(self.n-1)])                 # [reverse(x), 0..0]
y    = ifft(fft(chat) * fft(xhat))
return y[:self.n].real
```

Working out the circular convolution $y_k=\sum_m \texttt{chat}_m\,\texttt{xhat}_{(k-m)\bmod N}$
with this `chat` gives

$$y_i=(Hx)_i=\sum_{j} \widehat H_{ij}\,x_j,\qquad \widehat H_{ij}=\begin{cases}h_{i+j}, & i+j\le k-1,\\[2pt] 0,& i+j\ge k.\end{cases}$$

> **Note — `FastHankelProduct` implements the *triangular* Hankel.** Because `chat`
> only carries $h_0,\dots,h_{k-1}$ (the last $k-1$ slots are zero), the operator equals
> `scipy.linalg.hankel(h)` — the **single-argument**, lower-right-triangle-zeroed
> Hankel — *not* the fully populated one. Verified: `FastHankelProduct(h) @ x` matches
> `la.hankel(h) @ x` to $\sim\!10^{-16}$. This is the **triangular-Hankel defect**: the
> late samples $h_k,\dots,h_{2k-2}$ are discarded, which corrupts pole recovery when the
> impulse-response tail still carries energy (harmless for fast decay). A *fully
> populated* fast Hankel matvec instead pads `chat` with all $2k-1$ entries (or uses the
> linear-convolution form of §2–§3 with $h_0,\dots,h_{2k-2}$); that is what the fast
> `msvdreduce` below uses.

## 5. Block (MIMO) extension

For the block Hankel $H\in\mathbb{R}^{pk\times qk}$ with $p\times q$ blocks
$H_{ij}=h_{i+j}$ (Markov parameters), partition $x\in\mathbb{R}^{qk}$ into $k$ blocks
$x_j\in\mathbb{R}^{q}$. Then

$$(Hx)_i=\sum_{j=0}^{k-1} h_{i+j}\,x_j
\;\Longrightarrow\;
(Hx)_i^{(a)}=\sum_{b=1}^{q}\ \underbrace{\sum_{j} h_{i+j}^{(a,b)}\,x_j^{(b)}}_{\text{scalar Hankel matvec}}$$

for each output channel $a=1,\dots,p$. So a block Hankel matvec is **$pq$ scalar
Hankel matvecs** — each an FFT convolution of the scalar sequence $\{h_t^{(a,b)}\}_t$
against $\{x_j^{(b)}\}_j$ — summed over the input channels $b$ (Yu eq 441). Cost
$O(pq\,k\log k)$. Because each scalar Hankel is symmetric, the adjoint $H^{*}$ (needed by
`svds`) is the same routine with the channel roles swapped, $a\leftrightarrow b$.

## 6. Why it matters: fast partial Hankel SVD

$H^{*}H$ is Hermitian PSD, so its eigenvalues are $\sigma_k^2$ and a few iterations of
**Lanczos** (`scipy.sparse.linalg.eigsh` / `svds`) recover the top-$r$ singular triples
using only fast matvecs — the rank-$r$ partial SVD $H\approx\widetilde U\widetilde S
\widetilde V^{*}$ of [`01-hankel-svd-reduction.md`](01-hankel-svd-reduction.md) §2,
**without forming $H$**. That is Yu's "Fast Partial Block Hankel SVD" (§6.2), and it is
what makes Hankel-SVD reduction scale to the long ($k\sim10^3$–$10^4$) impulse responses
of market-impact data.

> **Implementation-note — `reduce_fft_truncate` bugs.** The SISO
> `reduce_fft_truncate(h, p)` builds `FastHankelProduct(h)` and runs `eigsh`:
> - **`dtype='float32'` (default in `FastHankelProduct.__init__`)** — the FFT actually
>   computes in float64, but `eigsh` trusts the declared dtype and iterates in single
>   precision, costing ~5–6 digits: pole error $6.4\times10^{-4}$ vs the dense
>   `reduce_svd_truncate`'s $6.0\times10^{-9}$. Declaring float64 recovers
>   $5.4\times10^{-9}$. **This is the primary bug.**
> - **Triangular Hankel** (§4 note) — inherited; corrupts slow-decay / tail-heavy IRs.
>
> The eigenvector shift `A = v.T[:p,1:] @ v[:-1,:p]` itself is fine for the symmetric
> Hankel (its eigenvectors are the singular vectors up to sign). The corrected fast
> path — float64 operator, fully-populated block Hankel, Lanczos `svds` — is implemented
> as `msvdreduce_fast` in `ga/reducers/hankel.py` (companion to `msvdreduce`).

## Related

- [`01-hankel-svd-reduction.md`](01-hankel-svd-reduction.md) — the (dense) Hankel-SVD reduction this accelerates.
- `concepts/ModelReduction` — `reduce_fft_truncate` is listed there as the FFT-based fast variant.
