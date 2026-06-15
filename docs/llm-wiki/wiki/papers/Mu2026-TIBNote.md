---
type: paper
tags: [factor-model, multi-asset]
sources: [Mu2026-TIBNote]
updated: 2026-06-14
---

# TIB Note — TIB Theory (orthonormal bases, band-fraction filter, Toeplitz AR↔IR)

**Author:** Yu Mu.  **Source:** `docs/llm-wiki/raw/tib_note.pdf` (7 pp).

## Contribution

Working notes that assemble the TIB framework: (1) the orthonormal rational-basis view of
the transfer function, (2) the band-fraction TIB filter `A = M⁻¹N`, (3) the optimality
criteria, (4) realization theory, and (5) the **functions-of-Toeplitz-matrices** view of
AR/MA processes — which is where `ar2ir` (AR coefficients → impulse response) is derived.

## Contents

- **Orthonormal bases.** `H(z) = Σᵢ θᵢ Bᵢ(z)`; fixed-pole, Laguerre, Kautz bases each
  incorporate only one pole, whereas the general Blaschke basis
  `Bᵢ(z) = (√(1−|ξᵢ|²)/(z−ξᵢ)) ∏ₖ (1−ξ̄ₖz)/(z−ξₖ)` lets a *variety* of poles in — the TIB
  construction. (Cf. [[concepts/TIBForm]], [[concepts/BlaschkeFactor]].)
- **TIB filter (band fraction).** `A = M⁻¹N`, state advanced by `M z_{t+1} = N z_t + M B ε_t`
  (solve the banded system rather than invert `M`). Single-input `M, N` are bidiagonal with
  `ρₖ=√(1−|λₖ|²), μₖ=ρ_{k+1}/ρₖ, γₖ=λ̄ₖμₖ` — the `siso_system_matrices` parameterization.
  The **MIMO** case is marked "pass" (left open).
- **Optimality.** `H₂ = Σ |hₖ−ĥₖ|²`; `H∞` = the operator norm of the Hankel difference. (Cf.
  `docs/derivations/04-system-norms-and-h2-optimality.md`.)
- **Realization theory.** Any proper `p×q` rational `F(λ) = D + C(λI−A)⁻¹B`.
- **Functions of Toeplitz matrices (AR ↔ IR).** MA process `y_t = Σ hₖ x_{t-k}` ↔ block
  Toeplitz `Y = T_h X`; AR process `y_t = Σ aₖ y_{t-k} + x_t` ↔ `X = T_a Y` with `T_a` having
  first column `[1, −a₁, −a₂, …]`. The duality `T_h · T_a = I` (i.e. `H(z)A(z)=1`,
  `H=1/A`) yields **`ar2ir`** as the lower-triangular Toeplitz solve `T_a h = e₁`, equivalently
  `h₀=1, hₙ=Σ_{k≥1} aₖ hₙ₋ₖ`.

## Connection to genesis-alpha

This note is the **source for the AR↔IR / `ar2ir`** math used in
[[concepts/VARSelfPrediction]]: `CAⁿB` are AR coefficients, and the forward impulse response
is `1/𝒜(z)` recovered by `ar2ir`. The note's scalar `ar2ir` (`T_a h = e₁`) is implemented in
`dynamical_system/filters/tib.py: ar2ir`; the note's **"MIMO: pass"** placeholder is filled by
`ga/filters/tib.py: mimo_ar2ir` — the block recursion `h₀=I, hₙ=Σ Aₖhₙ₋ₖ` (verified:
`A(z)H(z)=I` to ~1e-16, matches scalar `ar2ir` at `d=1`, and the forward IR carries the VAR
companion's data poles). The **full multivariate derivation** of the note's "MIMO: pass" case
is `docs/derivations/07-mimo-ar-to-ir.md`.

## Related Pages

[[concepts/VARSelfPrediction]], [[concepts/TIBForm]], [[concepts/ModesAndHankel]],
[[concepts/BlaschkeFactor]], [[papers/Kong2018-TIBInfoGeometry]]
