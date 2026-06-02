---
type: paper
tags: [execution, market-impact]
sources: [HanzonOliviPeeters2010-TangentialSchur]
updated: 2026-06-01
---

# Balanced Realizations of Discrete-Time Stable All-Pass Systems and the Tangential Schur Algorithm

**Authors:** Hanzon, B., Olivi, M., Peeters, R.L.M.  **Venue / Year:** arXiv:1012.3272, Dec 2010

## Contribution

Connects two parameterisations of multivariable stable all-pass (lossless) systems: (1) the **tangential Schur algorithm** (linear fractional transformations, reproducing kernel Hilbert spaces, Schur parameters + interpolation points) and (2) **state-space balanced realisations**. The recursive balanced canonical form is generalised from scalar to multivariable, with the realisation expressed directly in terms of tangential Schur parameters. This is the **primary upstream reference** for the Blaschke-Potapov machinery used in genesis-alpha (cited as "[Olivi 2010]" in Kong's dissertation eq. 243).

## Key Definitions (the ones we depend on)

**Blaschke factor** (eq. 26):
```
b_w(z) := (z − w) / (1 − w̄ z),   w ∉ T   (T = unit circle)
```
This is exactly `blaschke_factor(w, z)` in `ga/filters/tib.py` — **first argument = pole w, argument = evaluation point z**.

**Elementary J-inner factor** (eq. 28), the general (J ≠ I) form:
```
Φ(z) = X_x(b_w(z)) = I_{2p} + (b_w(z) − 1) · (x x* J)/(x* J x)
```
With **J = I** and a **unit vector** x (so x*Jx = ‖x‖² = 1) this reduces to the Blaschke-Potapov factor
`β_{w,u}(z) = I + (b_w(z) − 1) u u*` used throughout the TIB code.

By Potapov's theorem, every J-inner matrix function of McMillan degree n factors into n such elementary factors.

## Connection to genesis-alpha

- Confirms the **(pole, eval-point) argument convention** of `blaschke_factor(w, z)` and `blaschke_potapov_factor(w, z, u)`. This is exactly the convention that the line-185 bug in `null_basis_realization` violated (it passed `(eval, pole)`); the fix restores `β_{ωᵢ, yᵢ}(ωₖ*)` = `(pole=λᵢ, eval=z)`.
- The tangential Schur algorithm is the "reversed-direction" sibling of the null-basis recursion (Kong §6.5, Lemma 6.5.2) — both build balanced (A,B) but update opposite rows.
- Source authority for [[concepts/TIBForm]] and [[papers/Mu2026-NullBasisProofs]].

## Related Pages

[[concepts/BlaschkeFactor]], [[papers/Olivi2010-LosslessParametrization]], [[papers/Kong2018-TIBInfoGeometry]], [[papers/Mu2026-NullBasisProofs]], [[concepts/TIBForm]]
