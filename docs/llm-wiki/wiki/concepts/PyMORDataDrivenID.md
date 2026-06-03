---
type: concept
tags: [market-impact, execution, data]
sources: [PillonettoDeNicolao2010-KernelSysID]
updated: 2026-06-03
---

# pyMOR & Data-Driven / Noise-Robust System Identification

## Definition

The problem: **derive a reduced state-space realisation — poles `{λₖ}` and null vectors `{vₖ}` — from *measured* data** (an estimated impulse response, or noisy input/output time series), especially in a **low signal-to-noise** environment such as market-impact estimation from order-flow data.

Key distinction that organises everything on this page:

> **Model order reduction (MOR) ≠ system identification (sysID).**
> MOR compresses a *known, clean* model (or clean measurement data). sysID *estimates* a model from *noisy* data. The genesis-alpha pipeline currently does MOR on a (assumed clean) IR; the hard part at low SNR is the sysID it skips.

## Why It Matters for genesis-alpha

The current paths to `(poles, null vectors)` are both **deterministic realisation theory**:

- `ga/reducers/hankel.py` → `msvdreduce` / `reduce_svd_truncate`: Hankel-SVD = the **Eigensystem Realization Algorithm (ERA)** / Ho–Kalman.
- `ga/reducers/balanced_truncation.py` → balanced truncation → `tib_from_state_space` (the corrected BT→TIB extraction; see [[TIBForm]]).

Both assume a **clean** impulse response. At low SNR the failure is structural: the Hankel singular-value spectrum has **no clean gap**, truncation absorbs **noise modes as if they were dynamics**, and the estimated poles are biased. This is precisely the [[TIBForm]] open question *"how does TIB reduction quality degrade with noisy IR estimates?"*. `info_svd_reduce` (cepstrum / [[InformationGeometry]]) is a different projection but still SVD-on-noisy-data.

**The reframing:** poles and null vectors are only as good as the `(A,B,C)` estimate they come from. `tib_from_state_space` is *exact* — so the work is entirely **upstream**, in how `(A,B,C)` is estimated.

```
noisy I/O data ──[robust sysID]──▶ (A,B,C) estimate
                                       │
                           tib_from_state_space   (exact, real+complex poles)
                                       ▼
                       poles + data-driven null vectors
```

At low SNR the **joint MIMO** estimate (pooling across instruments, shared modes) is where subspace/Bayesian methods buy variance reduction that per-instrument fitting cannot — the representational power the *canonical* null vectors throw away (each input gets a disjoint pole subset; the MIMO fit decouples into independent per-input fits with no pooling benefit).

## What pyMOR ([pymor.org](https://pymor.org/)) provides — and its limit

pyMOR is a mature MOR library (parametric PDEs + LTI systems). System-theoretic reductors: **BT, IRKA / SOR-IRKA, Loewner, AAA, ERA, Hankel-norm-type, spectral-factor, passivity-preserving (pH-IRKA)**.

| pyMOR method | Data-driven? | Relevance |
|---|---|---|
| **ERA** | from impulse-response (Markov params) | mature, MIMO-correct version of our `msvdreduce` — **same noise sensitivity** |
| **BT / IRKA** | need a clean state-space | better-tested than our `balanced_truncation.py`; not estimators |
| **Loewner / AAA** | from frequency-response samples | useful *if* working in frequency domain; assume relatively clean samples (AAA regularisable) |

> ⚠️ **pyMOR will not solve the low-SNR problem.** Its data-driven methods (ERA, Loewner, AAA) presume the data samples are good; none model a noise process. Use pyMOR to (a) get better-engineered reduction and (b) cross-check our `msvdreduce` / `tib_from_state_space`, **not** as a noise remedy.

## The real lever: noise-aware methods (ordered by value here)

1. **Subspace identification — N4SID / MOESP / CVA.** Estimate `(A,B,C)` from noisy I/O via QR + truncated SVD of block-Hankel data, with explicit process/measurement-noise handling. The standard low-SNR MIMO workhorse; direct drop-in upstream of `tib_from_state_space`. Python: **SIPPY**, or `control` + `slycot`.
2. **Nuclear-norm / regularised subspace (N2SID)** for **short, noisy batches** — convex, far more stable than plain SVD truncation when the SV gap is washed out.
3. **Bayesian / kernel-based ID** (stable-spline kernels — [[papers/PillonettoDeNicolao2010-KernelSysID|Pillonetto & De Nicolao 2010]]) — best bias/variance trade-off for short, low-SNR records; encodes stability + smoothness as priors. **Directly fusible with TIB — see *Relationship to TIB* below.**
4. **Denoise the Hankel before SVD** — cheapest, highest-leverage upgrade to the *existing* pipeline: **optimal singular-value hard thresholding (Gavish–Donoho)** for data-driven order selection (replace the fixed `order`), and/or **Cadzow / structured low-rank approximation (SLRA)** to reproject onto the rank-`r` Hankel manifold.
5. **Errors-in-variables / TLS-ERA** if staying with realisation-style methods.
6. **Frequency domain: Vector Fitting / AAA / Loewner with regularisation** — if a (smoothed) transfer-function estimate is available; pyMOR covers AAA/Loewner.

## Relationship to TIB: the Bayesian-TIB fusion

**Reframing (the genesis-alpha view of TIB).** For fixed poles `{λₖ}` and null vectors `{vₖ}`, TIB fixes `(A,B)` and the state Krylov sequences `ψₖ(j,t) = (Aᵗ B)[k,j]` are a set of **`P` basis functions** for the impulse response; the empirical IR is `ĝ(o,j,t) = Σₖ C[o,k] ψₖ(j,t)`, with `C` estimated by least squares. So **TIB is a linear basis expansion of the IR** ([[TIBForm]]; this is what `fit_output_matrix` does).

Under this view TIB and kernel/Bayesian ID are the **same kind of object** — both are *linear estimators of the same empirical IR*, i.e. both are basis methods. They differ only in how complexity is regularised:

| | basis | complexity control |
|---|---|---|
| **TIB** | explicit, finite: the `P` decaying modes you choose (poles + null vectors) | **hard truncation** (which / how many poles) + plain LS on `C` |
| **Kernel/Bayesian** | implicit, large/∞: the kernel's Mercer eigenbasis | **soft shrinkage** (stability/decay prior), strength set by marginal likelihood |

> **TIB = hard-truncated explicit basis; kernel = soft-weighted implicit basis.** Both encode the *same* prior — "the IR is a sum of stable, decaying modes." They are **not** opposite ends of a bias–variance spectrum (each spans it via its own knob — `P` vs shrinkage `λ`); the real difference is the *regularisation mechanism*: **subset selection (TIB) vs shrinkage (kernel)** — the classic best-subset-vs-ridge dichotomy.

**The fusion (Bayesian-TIB).** Keep the TIB basis; replace plain least-squares for `C` with a kernel-regularised (Bayesian) fit — estimating the *same* `C`, just penalised:

$$\widehat C = \arg\min_C \ \|g_{\text{emp}} - \Phi\,C\|^2 \;+\; \sigma^2\,\operatorname{tr}\!\big(C\,P^{-1}C^\top\big)$$

where `Φ` is the TIB basis (from poles + null vectors), `P` a prior covariance on the **TIB mode coefficients** (the kernel, now on coefficients rather than raw IR lags), and `σ²` / prior-scale tuned by **marginal likelihood**. Plain TIB is the no-penalty limit (`P → ∞`).

**Why it helps:** it decouples *basis richness* from *effective complexity*. Today TIB's only knob is `P` (and the pole choice): lowering bias means enlarging `P`, but plain LS on many coefficients overfits noise — fatal at low SNR. With the prior you use a **large / overcomplete pole grid** (low bias) and let evidence-tuned shrinkage control variance, while staying inside the interpretable TIB basis. In code it is a small change to the `C`-fit: ridge / Bayesian regression in the TIB feature space instead of `np.linalg.lstsq`.

➡️ **Bayesian-TIB = TIB's basis + the kernel's penalty.** Same estimation target (`C` / the empirical IR), self-tuning regularisation, interpretable modes retained — so kernel ID and TIB are not rivals but two implementations of one idea (impose stable decaying modes: *pick* them, or *shrink* over them).

## Recommendation

- **Do not** swap our reducers for pyMOR for SNR reasons — it won't help there; do use pyMOR's **ERA/BT** as a tested reference to validate against.
- **For low SNR, add a sysID front-end:** start with **N4SID (SIPPY/slycot) → `tib_from_state_space`**, plus **optimal SVHT order selection + a Cadzow/SLRA Hankel-denoise** as a quick win on the current pipeline. For short/very noisy data, go **Bayesian kernel-based**.

## Python tooling

- Reduction / data-driven MOR: **pyMOR** (ERA, Loewner, AAA, BT, IRKA), `control` + `slycot`.
- Subspace sysID: **SIPPY**, `slycot` (`n4sid`-style).
- DMD variants (modal extraction, noise-robust): **PyDMD** (optDMD, TLS-DMD, BOP-DMD).

## Open Questions

- Benchmark `noisy IR → N4SID → tib_from_state_space` vs `msvdreduce` for pole/IR recovery across SNR levels — quantify the crossover where sysID beats plain reduction.
- **Implement & benchmark Bayesian-TIB** (kernel prior on `C` over an overcomplete pole grid, marginal-likelihood-tuned) vs plain TIB-LS (`fit_output_matrix`) across SNR — does evidence-tuned shrinkage on a rich basis beat hard pole truncation at low SNR? Which prior on the mode coefficients (ridge vs decay-weighted vs stable-spline-on-modes)?
- Best order-selection rule for market-impact IR: optimal SVHT vs Hankel-energy threshold vs information criteria.
- Does joint-MIMO pooling across instruments measurably reduce pole-estimate variance at realistic order-flow SNR? (ties to the canonical-vs-data-driven null-vector advantage)
- Frequency- vs time-domain estimation for impact: is a smoothed transfer-function + Loewner/AAA more robust than time-domain subspace ID?

## References / External Sources

- **Pillonetto & De Nicolao (2010)**, *A new kernel-based approach for linear system identification*, Automatica 46:81–93 — the stable-spline kernel; foundation of the kernel/Bayesian method and the Bayesian-TIB fusion here. **Ingested:** [[papers/PillonettoDeNicolao2010-KernelSysID]].
- pyMOR reductors: https://docs.pymor.org/latest/autoapi/pymor/reductors/index.html
- N2SID — Nuclear Norm Subspace Identification for short data batches: https://arxiv.org/pdf/1401.4273
- Hankel matrix denoising for subspace state-space ID (modal analysis): https://www.researchgate.net/publication/400328237
- Performance analysis of N4SID: https://www.researchgate.net/publication/3758793
- Gavish & Donoho (2014), *The Optimal Hard Threshold for Singular Values is 4/√3* (SVHT)

## Related Pages

[[TIBForm]], [[ModelReduction]], [[InformationGeometry]], [[MarketImpact]]
