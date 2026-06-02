---
type: paper
tags: [execution, market-impact]
sources: [Olivi2010-LosslessParametrization]
updated: 2026-06-01
---

# Parametrization of Rational Lossless Matrices with Applications to Linear Systems Theory

**Author:** Olivi, M.  **Venue / Year:** Habilitation à Diriger des Recherches (HDR), Université de Nice Sophia Antipolis, 25 Oct 2010

## Contribution

Olivi's habilitation thesis — a comprehensive treatment of rational lossless (all-pass / inner) matrix parameterisation and its use in linear systems theory. Gives the cleanest statement of the **elementary inner (Potapov) factor** that the genesis-alpha TIB code is built on, in the J = I form that matches the implementation directly.

## Key Definition (exact match to the code)

**Elementary inner / Potapov factor** (eq. 1.11):
```
B_{w,u}(z) = I + (b_w(z) − 1) u u*,   b_w(z) = (z − w)/(1 − w̄ z)
```
where **w** is a zero of the inner function Q(z) inside the open unit disk, and **u** is a **unit vector in ker Q(w)**.

This is *identically* `blaschke_potapov_factor(w, z, u)` in `ga/filters/tib.py`:
- arg1 `w` = pole / interpolation point
- arg2 `z` = evaluation point (inside `b_w(z)`)
- arg3 `u` = unit vector

**Potapov factorization** (Prop. 1.1.1, eq. 1.12):
```
Q(z) = B_{wₙ,uₙ}(z) · B_{wₙ₋₁,uₙ₋₁}(z) ⋯ B_{w₁,u₁}(z) · Q₀
```
with wᵢ in the open unit disk and uᵢ unit vectors. The product structure (and the fact each factor needs a **unit** uᵢ) is exactly what the null-basis recursion exploits.

## Connection to genesis-alpha

- **Definitive source** for the `blaschke_potapov_factor` definition and the unit-vector requirement on `u` — directly justifies the two-sided `‖u‖ = 1` guard we added (rejecting both too-short and too-long vectors).
- The unit-vector-in-kernel condition (`B_{w,u}(w)·u = 0`, Kong eq. 249) is the property tested by `test_null_property` in `tests/test_tib.py`.
- The Potapov product order (eq. 1.12) underpins the `M = M @ β†` accumulation in `null_basis_realization`.

## Related Pages

[[concepts/BlaschkeFactor]], [[papers/HanzonOliviPeeters2010-TangentialSchur]], [[papers/Kong2018-TIBInfoGeometry]], [[papers/Mu2026-NullBasisProofs]], [[concepts/TIBForm]]
