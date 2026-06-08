# System Norms ($\ell^2$, $H^2$, $H^\infty$) and Why $H_2$ Is a Sub-optimal Reduction Objective

> Background for the reduction write-up. §1 defines and relates the three norms that
> show up across these notes; §2 explains the sense(s) in which choosing $H_2$ as the
> *objective* of model reduction is sub-optimal — the point that motivates the
> information-distance / cepstrum reduction (`info_svd_reduce`). Notation:
> [`00-overview.md`](00-overview.md).

## 1. Three norms: $\ell^2$, $H^2$, $H^\infty$

The first thing to keep straight is **what each norm measures**. $\ell^2$ norms a
**signal** (a time-domain sequence); $H^2$ and $H^\infty$ norm a **system** (a
transfer function, i.e. a function analytic in the unit disk). They are tightly
linked, but they live on different objects.

### 1.1 $\ell^2$ — signal energy (time domain)

For a sequence (signal) $x=(x_t)_{t}$,

$$\lVert x\rVert_{\ell^2}^2 = \sum_t \lVert x_t\rVert^2 .$$

This is **signal energy** — the squared $\ell^2$-norm, the discrete analog of
$\int\lVert x(t)\rVert^2\,dt$ (see `concepts/ModesAndHankel`, *What the Hankel singular
values mean*). Nothing thermodynamic; just the total size of a signal.

### 1.2 $L^2(\mathbb{T})$ and $H^2$ — energy of a function on the circle

Move to the frequency domain via the $z$-transform on the unit circle
$z=e^{i\theta}$. For a matrix-valued function $F$ on the circle,

$$\lVert F\rVert_{L^2}^2 = \frac{1}{2\pi}\int_{-\pi}^{\pi}\lVert F(e^{i\theta})\rVert_F^2\,d\theta
= \sum_{k\in\mathbb{Z}}\lVert c_k\rVert_F^2 \quad(\text{Parseval, two-sided coeffs}).$$

The **Hardy space** $H^2$ is the subspace of $L^2(\mathbb{T})$ of functions whose
**negative Fourier coefficients vanish** — equivalently, boundary values of functions
**analytic inside the disk**:

$$G(z)=\sum_{k\ge 0} h_k\, z^k,\qquad
\lVert G\rVert_{H^2}^2 = \frac{1}{2\pi}\int_{-\pi}^{\pi}\lVert G(e^{i\theta})\rVert_F^2\,d\theta
= \sum_{k\ge 0}\lVert h_k\rVert_F^2 .$$

For a system, $G$ is the **transfer function** and the $h_k=CA^{k-1}B$ are the Markov
parameters; analyticity in the disk $=$ **causality** of the impulse response. Two
structures distinguish $H^2$ from a bare sequence space:

- **Causality** is built in (no negative-lag part). $\ell^2(\mathbb{N})$ (causal
  signals) $\leftrightarrow H^2$; $\ell^2(\mathbb{Z})$ (two-sided) $\leftrightarrow
  L^2(\mathbb{T})$.
- **$H^2$ is a reproducing-kernel Hilbert space**: point evaluation $G\mapsto G(w)$ is
  bounded, with Szegő kernel $\kappa_w(z)=\tfrac{1}{1-\bar w z}$. This is the geometry
  in which the Blaschke–Potapov deflation is a Gram–Schmidt step
  ([`02-null-basis-realization.md`](02-null-basis-realization.md) §3).

### 1.3 $H^\infty$ — peak gain (worst case)

$$\lVert G\rVert_{H^\infty} = \sup_{|z|=1}\ \sigma_{\max}\!\big(G(e^{i\theta})\big)
= \sup_{u\ne 0}\ \frac{\lVert y\rVert_{\ell^2}}{\lVert u\rVert_{\ell^2}},\qquad y = G\!*\!u .$$

$H^\infty$ is the **worst-case energy gain** of the system — the $\ell^2\!\to\!\ell^2$
**induced (operator) norm** of the convolution. Unlike $H^2$ it is therefore
**submultiplicative**, $\lVert GH\rVert_\infty\le\lVert G\rVert_\infty\lVert
H\rVert_\infty$, which is why it is the natural norm for robustness/feedback bounds.

### 1.4 How they relate

**Parseval ties $\ell^2$ and $H^2$ — they are the *same number*.** For a causal IR
$\{h_k\}_{k\ge0}$ with transfer function $G$,

$$\boxed{\;\lVert h\rVert_{\ell^2}^2 = \sum_{k\ge0}\lVert h_k\rVert_F^2 = \lVert G\rVert_{H^2}^2\;}$$

So the $H_2$ reduction criterion $\sum_k\lVert h_k-\hat h_k\rVert^2$ is simultaneously
the $\ell^2$ norm of the impulse-response error (time domain) and the $H^2$ norm of the
transfer-function error (analytic domain).

**$H^2$ vs $H^\infty$ — average vs peak.** $H_2$ is an **energy/average** over the
circle; $H_\infty$ is the **peak** over the circle. $H_2$ is the IR energy (= output
power for unit-variance white-noise input); $H_\infty$ is the worst-case sinusoidal
amplification. They are different norms and generally rank reduced models
differently.

| | $\ell^2$ | $H^2$ | $H^\infty$ |
|---|---|---|---|
| **object** | signal (sequence) | transfer function (analytic in disk) | transfer function (analytic in disk) |
| **definition** | $\sum_t\lVert x_t\rVert^2$ | $\frac{1}{2\pi}\!\int\lVert G\rVert_F^2\,d\theta=\sum_k\lVert h_k\rVert_F^2$ | $\sup_\theta\sigma_{\max}(G(e^{i\theta}))$ |
| **domain** | time | frequency / analytic | frequency / analytic |
| **interpretation** | signal energy | IR energy / avg power / white-noise output variance | peak gain / worst-case energy amplification |
| **operator norm?** | — | no (not induced) | yes ($\ell^2\!\to\!\ell^2$), submultiplicative |
| **extra structure** | — | causality; RKHS (Szegő kernel, point eval) | — |
| **role in reduction** | the IR-error energy | `msvdreduce` / IRKA objective; Walsh interpolation | BT error bound $2\sum_{i>r}\sigma_i$ |

(See `concepts/ModelReduction` for the method-by-norm summary;
[`01-hankel-svd-reduction.md`](01-hankel-svd-reduction.md) §6 for the BT-vs-`msvdreduce`
split that turns on $H_\infty$ vs $H_2$.)

## 2. Why $H_2$ is a sub-optimal reduction objective

"Sub-optimal" is said in three different senses. The one that matters most for
*data-estimated* systems (market-impact IRs) is the first.

### 2.1 The metric sense — $H_2$ is the wrong distance on the space of systems

A linear system is really a **probability model**: it defines a spectral density / a
distribution over output processes. The right question for approximation is not "are
the impulse responses close as vectors?" but "**are the two systems hard to tell apart
statistically?**" Those differ.

- **$H_2$** measures Euclidean distance on the *raw* transfer function,
  $\lVert H-\hat H\rVert_{H_2}^2=\sum_k\lVert h_k-\hat h_k\rVert^2$ — a geometric
  distance in coefficient space.
- **Information distance** (Fisher–Rao / $\approx$ KL) is the geodesic distance on the
  **statistical manifold** of systems. Kong (2018) (`concepts/InformationGeometry`,
  Ch 3) shows that in **cepstrum coordinates** $\log H(z)=\sum_k a_k z^k$ the Fisher
  information matrix is the **identity**, and

$$I^2(H,\text{white}) = \sum_k |a_k|^2 = \big\lVert \log H\big\rVert_{H_2}^2 .$$

The statistically natural distance is the $H_2$ norm of **$\log H$ (the cepstrum)**,
*not* of $H$ itself. So plain $H_2$ optimizes Euclidean distance on $H$ when it should
optimize Euclidean distance on $\log H$.

**Why the $\log$ matters.** It makes the metric **relative/multiplicative** instead of
**absolute**:

- $H_2$ weights errors by the gain, so it spends its budget on the **loud bands** and
  barely notices a factor-of-2 error in a **quiet band**.
- Information distance treats a relative error the same everywhere — a $2\times$ misfit
  in a low-gain band is as statistically confusable as in a high-gain band, and the
  cepstrum norm sees it. The Hellinger / total-variation sandwich ($8H\le I$,
  $K\!\cdot\!I\le 8H$; Kong Lemma 3.1.1–3.1.2) ties information distance directly to
  the **probability of confusing two systems in a hypothesis test**.

**Where it bites for genesis-alpha.** Noisy, non-stationary, near-unit-root, or
non-rational IRs — exactly market impact (power-law decay, random-walk-like).
Information distance stays **finite even without a proper spectral density**;
$H_2$ can blow up or be dominated by the near-unit-root pole. Mullhaupt–Choi's result
(Yu §6.2) that *prediction information lives in the unstable/outer part* is the same
idea: the log/cepstrum view isolates the predictively relevant part that raw $H_2$
smears together. This is exactly why `info_svd_reduce` (cepstrum reduction) exists in
`ga/reducers/hankel.py` alongside `msvdreduce`.

### 2.2 The computational sense — $H_2$-optimal reduction is non-convex

Even if you *want* the $H_2$ optimum, it is hard to reach. The $H_2$-optimality
conditions (Walsh's theorem / Meier–Luenberger interpolation at the mirrored poles
$1/\bar z_k$; Yu §6.1 Thm 40–41, eqs 417–419 and 423–434) characterize only
**stationary points**. Iterative solvers — IRKA, and the MIRIAm algorithm Yu cites —
converge to *local* optima with **no global guarantee** (Yu: "a proof of the
convergence is not provided"). Yu's own Hankel-SVD method "guarantees a small $H_2$
norm, **but not necessarily minimal** $H_2$ norm for the given order" (§6.2). So
choosing $H_2$ buys a non-convex problem you typically don't solve exactly.

### 2.3 The well-posedness sense — Hankel/$H^\infty$ has an exact optimum, $H_2$ doesn't

Adamjan–Arov–Krein (AAK) theory gives a **closed-form** optimal order-$r$ Hankel-norm
approximation — the error equals the $(r{+}1)$-th Hankel singular value, constructed
directly from the SVD. $H_2$ has no such closed-form optimum. Yu's hybrid scheme
(§6.3) exploits this: find the optimal **inner part in the $H^\infty$ sense (AAK)**,
then fit the rest with the nice TIB geometry. Picking $H_2$ as the global objective
gives up that well-posed, computable optimum.

### 2.4 Summary

| Sense | "$H_2$ sub-optimal" means | Better choice |
|---|---|---|
| **Metric (primary)** | $H_2$ is Euclidean on $H$, not the statistical (Fisher) distance | information distance $=\lVert\log H\rVert_{H_2}$ (cepstrum; `info_svd_reduce`) |
| Computational | $H_2$-optimal reduction is non-convex; only local optima | accept a stationary point, or change objective |
| Well-posedness | $H_2$ has no closed-form optimum | Hankel norm / $H^\infty$ via AAK (exact, computable) |

For systems *estimated from data*, the principled objective is to minimize statistical
confusability (information distance); raw $H_2$ is a convenient but sub-optimal stand-in.

## Related

- [`01-hankel-svd-reduction.md`](01-hankel-svd-reduction.md) — `msvdreduce` ($H_2$/Hankel) and the BT ($H^\infty$) separation.
- `concepts/InformationGeometry`, `papers/Kong2018-TIBInfoGeometry` — the information-distance / cepstrum story.
- `concepts/ModelReduction` — methods ranked by target norm; `concepts/ModesAndHankel` — the energy meaning of the Hankel singular values.
