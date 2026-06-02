---
type: concept
tags: [execution, market-impact]
sources: [Olivi2010-LosslessParametrization, HanzonOliviPeeters2010-TangentialSchur, Kong2018-TIBInfoGeometry]
updated: 2026-06-01
---

# Blaschke Factor & Blaschke-Potapov Factor

## Definition

**Scalar Blaschke factor** (`blaschke_factor(w, z)` in `ga/filters/tib.py`):
```
b_w(z) = (z − w) / (1 − w̄·z)
```
A Möbius (fractional-linear) function of the variable `z`, parameterised by `w` (|w| < 1):

- **Zero** at `z = w` (numerator vanishes).
- **Analytic pole** at `z = 1/w̄`, which lies *outside* the unit disk (since |w| < 1) — so b_w is holomorphic inside the disk.
- **All-pass / unimodular**: on the unit circle |z| = 1, |b_w(z)| = 1.

**Matrix Blaschke-Potapov factor** (`blaschke_potapov_factor(w, z, u)`), Olivi 2010 eq. (1.11):
```
B_{w,u}(z) = I + (b_w(z) − 1) u u*,   u a unit vector
```
With J = I and ‖u‖ = 1 this is the elementary inner / Potapov factor. The general J-inner form is Hanzon-Olivi-Peeters 2010 eq. (28): `I + (b_w(z) − 1)·(x x* J)/(x* J x)`.

## ⚠️ Naming caveat: "pole" vs analytic pole

The code/docstrings call `w` the **"pole"**, meaning the **system pole** (an eigenvalue of A you prescribed, |w|<1). As a *function*, b_w actually has its **zero** at `w` and its analytic pole at `1/w̄` (outside the disk). The cleaner term (used by Hanzon-Olivi-Peeters) is **interpolation point**. Do not read "pole" as "the analytic pole of b_w".

## The two arguments are NOT symmetric

```
b_w(z) = (z − w) / (1 − w̄·z)
b_z(w) = (w − z) / (1 − z̄·w)      # different function!
```
Numerators differ by a sign; denominators differ unless w, z are real. So **b_w(z) ≠ b_z(w)** — swapping arguments yields a genuinely different value (e.g. real poles 0.5, 0.3 → +0.235 one way, −0.235 the other).

**Tell-tale cue:** the *pole / interpolation point* is the argument that appears **conjugated in the denominator** (`w̄`). The *evaluation point* appears un-conjugated in the denominator (and bare in the numerator). If unsure which slot is which, look for the conjugate.

## How to tell pole from evaluation point in the null-basis recursion

In `null_basis_realization`, **both** arguments come from the prescribed pole set {λ₁,…,λₙ}, so you cannot distinguish them by origin — only by **role**. At step k:
```
yₖ ∝ [ ∏_{i=1}^{k−1} β_{λᵢ, yᵢ}(λₖ*) ] · vₖ
            └ pole = λᵢ ┘  └ eval = λₖ* ┘
```

| Question | Pole (arg 1) | Evaluation point (arg 2) |
|---|---|---|
| Which loop index? | `i` (inner loop, *old* poles) | `k` (current *new* pole, conjugated) |
| How many in the product? | k−1, one per existing pole | exactly 1, shared by all factors |
| Role | *identity* of each factor | *where you probe* the product |
| In the formula | conjugated (`w̄`) | not conjugated in denominator |

Code: `blaschke_potapov_factor(lambd[i], z, ys[:, i])` with `z = conj(lambd[k])` → `lambd[i]` = old pole (parameter), `z` = new point (evaluation).

## Why this evaluate-at-other-poles structure exists

Orthogonalisation. Each factor satisfies the **null property** `β_{λᵢ,yᵢ}(λᵢ)·yᵢ = 0` (Kong eq. 249): a factor annihilates its own direction at its own pole. Forming the product of all previous factors and evaluating at the new point λₖ*, then applying to vₖ, deflates the new direction against the already-placed ones — a recursive Gram-Schmidt in the reproducing-kernel space. The conjugate λₖ* (not λₖ) comes from how interpolation conditions for inner functions live at reflected points in the disk.

## Why It Matters for genesis-alpha

- Getting the (pole, eval) order right is essential: swapping them (the line-185 bug, now fixed) still produces unit vectors and still passes the TIB sanity check (Lemma B), but silently realises the **wrong transfer function**. See [[papers/Mu2026-NullBasisProofs]].
- The unit-vector requirement on `u` (Olivi eq. 1.11, "u in ker Q(w)") justifies the two-sided ‖u‖ = 1 guard in `blaschke_potapov_factor`.
- The null property `β_{w,u}(w)·u = 0` is tested by `test_null_property` in `tests/test_tib.py`; unimodularity on |z|=1 by `test_unitary_on_unit_circle`.

## Related Pages

[[concepts/TIBForm]], [[papers/Olivi2010-LosslessParametrization]], [[papers/HanzonOliviPeeters2010-TangentialSchur]], [[papers/Kong2018-TIBInfoGeometry]], [[papers/Mu2026-NullBasisProofs]]
