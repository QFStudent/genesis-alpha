---
type: paper
tags: [execution, market-impact]
sources: [Mu2026-NullBasisProofs]
updated: 2026-06-01
---

# Some Proofs: Null Basis Realization

**Authors:** Mu, Y.  **Venue / Year:** Internal Notes, Jan 2026  
**Completed by:** wiki derivation, 2026-06-01

## Contribution

Complete proof by induction that the null basis recursive algorithm (Algorithm 1) generates matrices (Aₖ, Bₖ) that are **triangular input balanced** — i.e. satisfy AₖAₖ* + BₖBₖ* = Iₖ at every step k. This is Theorem 6.5.1 from [[papers/Kong2018-TIBInfoGeometry]] Ch. 6.5.

---

## The Algorithm

**Dimensions:** Aₖ ∈ ℂᵏˣᵏ, Bₖ ∈ ℂᵏˣᵐ, Lₖ ∈ ℂᵐˣ⁽ᵏ⁻¹⁾, Pₖ ∈ ℂᵐˣᵐ, yₖ ∈ ℂᵐ with ‖yₖ‖ = 1.

Given null basis pairs {ωₖ, vₖ}, vₖ ∈ ℂᵐ, for k = 1…n:

**Initialization (k = 1):**
```
y₁ = v₁ / ‖v₁‖,   t₁ = √(1 − |ω₁|²),   A₁ = [ω₁],   B₁ = [t₁ y₁*]
L₁ = []  (m×0 empty),   P₁ = Iₘ,   J₁ = Iₘ − (1 + ω₁*) y₁ y₁*
```

**Recursive step (k → k+1):**
1. `Lₖ = [J_{k−1} L_{k−1} | t_{k−1} y_{k−1}]`  — m×(k−1) matrix
2. `Pₖ = J_{k−1} P_{k−1}`  — m×m (product of all previous Jᵢ)
3. `tₖ = √(1 − |ωₖ|²)`
4. `yₖ` = normalised image of vₖ under Blaschke-Potapov product (ensures ‖yₖ‖ = 1)
5. `Jₖ = Iₘ − (1 + ωₖ*) yₖ yₖ*`
6. Let `a = tₖ yₖ* Lₖ`  (1×(k−1)),  `b = tₖ yₖ* Pₖ`  (1×m)

```
Aₖ = [[A_{k−1},  0 ],    Bₖ = [B_{k−1}]
      [a,        ωₖ]]          [b      ]
```

---

## Proof of Theorem 6.5.1

**Goal:** AₖAₖ* + BₖBₖ* = Iₖ for all k ≥ 1.

**Proof by induction.**

### Base case k = 1

A₁A₁* + B₁B₁* = |ω₁|² + t₁² y₁* y₁ = |ω₁|² + (1 − |ω₁|²)·1 = 1 ✓

### Inductive step

Assume A_{k−1}A_{k−1}* + B_{k−1}B_{k−1}* = I_{k−1}. Expanding AₖAₖ* + BₖBₖ*:

```
AₖAₖ* + BₖBₖ* = [[A_{k−1}A_{k−1}* + B_{k−1}B_{k−1}*,   A_{k−1}a* + B_{k−1}b*],
                   [aA_{k−1}* + bB_{k−1}*,                 aa* + bb* + |ωₖ|²    ]]
```

By hypothesis the top-left block equals I_{k−1}. It remains to prove:

- **(1) Off-diagonal = 0:** A_{k−1}a* + B_{k−1}b* = 0
- **(2) Diagonal = 1:** aa* + bb* + |ωₖ|² = 1

These reduce to two structural lemmas about Lₖ and Pₖ.

---

### Lemma B: LₖLₖ* + PₖPₖ* = Iₘ

**Proof by induction on k.**

**Base (k=1):** L₁L₁* + P₁P₁* = 0 + Iₘ = Iₘ ✓

**Step (k → k+1):** Assume LₖLₖ* + PₖPₖ* = Iₘ. Then:

```
L_{k+1}L_{k+1}* = [JₖLₖ | tₖyₖ][JₖLₖ | tₖyₖ]*
                 = JₖLₖLₖ*Jₖ* + tₖ² yₖyₖ*

P_{k+1}P_{k+1}* = JₖPₖPₖ*Jₖ*
```

Adding: `L_{k+1}L_{k+1}* + P_{k+1}P_{k+1}* = Jₖ(LₖLₖ* + PₖPₖ*)Jₖ* + tₖ² yₖyₖ* = JₖJₖ* + tₖ² yₖyₖ*`

**Key computation of JₖJₖ*:** With Jₖ = Iₘ − (1 + ωₖ*)yₖyₖ* and Jₖ* = Iₘ − (1 + ωₖ)yₖyₖ*:

```
JₖJₖ* = (I − (1+ωₖ*)yₖyₖ*)(I − (1+ωₖ)yₖyₖ*)
       = I − (1+ωₖ)yₖyₖ* − (1+ωₖ*)yₖyₖ* + (1+ωₖ*)(1+ωₖ)yₖyₖ*yₖyₖ*
```

Since ‖yₖ‖ = 1: yₖyₖ*yₖyₖ* = yₖ(yₖ*yₖ)yₖ* = yₖyₖ*. Collecting coefficients of yₖyₖ*:

```
coefficient = −(1+ωₖ) − (1+ωₖ*) + (1+ωₖ*)(1+ωₖ)
            = −1−ωₖ−1−ωₖ* + 1+ωₖ+ωₖ*+|ωₖ|²
            = −1 + |ωₖ|² = −tₖ²
```

Therefore **JₖJₖ* = Iₘ − tₖ² yₖyₖ*** and:

```
L_{k+1}L_{k+1}* + P_{k+1}P_{k+1}* = (Iₘ − tₖ²yₖyₖ*) + tₖ²yₖyₖ* = Iₘ  ✓
```

---

### Lemma A: A_{k−1}Lₖ* + B_{k−1}Pₖ* = 0  (for k ≥ 2)

**Proof by induction on k.**

**Base (k=2):** L₂ = [t₁y₁] (m×1), P₂ = J₁.

```
A₁L₂* + B₁P₂* = ω₁(t₁y₁)* + (t₁y₁*)J₁*
               = t₁[ω₁y₁* + y₁*J₁*]
               = t₁ y₁*[ω₁I + J₁*]
```

With J₁* = I − (1+ω₁)y₁y₁*:

```
ω₁I + J₁* = ω₁I + I − (1+ω₁)y₁y₁* = (1+ω₁)(I − y₁y₁*)
```

So: `t₁ y₁*(1+ω₁)(I − y₁y₁*) = t₁(1+ω₁)(y₁* − y₁*y₁y₁*) = t₁(1+ω₁)(y₁* − y₁*) = 0  ✓`

(used ‖y₁‖ = 1 → y₁*y₁ = 1)

**Step (k → k+1):** Assume A_{k−1}Lₖ* + B_{k−1}Pₖ* = 0. Show AₖL_{k+1}* + BₖP_{k+1}* = 0.

Partition L_{k+1}* (k×m) and P_{k+1}* (m×m):

```
L_{k+1}* = [[Lₖ*Jₖ*],    P_{k+1}* = Pₖ*Jₖ*
             [tₖyₖ* ]]
```

Then:

```
AₖL_{k+1}* = [[A_{k−1},  0 ]] · [[Lₖ*Jₖ*]] = [[A_{k−1}Lₖ*Jₖ*              ]]
              [[a,        ωₖ]]   [[tₖyₖ* ]]    [[aLₖ*Jₖ* + ωₖtₖyₖ*          ]]

BₖP_{k+1}* = [[B_{k−1}]] · Pₖ*Jₖ* = [[B_{k−1}Pₖ*Jₖ*]]
              [[b      ]]               [[bPₖ*Jₖ*      ]]
```

**Top block:** `(A_{k−1}Lₖ* + B_{k−1}Pₖ*)Jₖ* = 0 · Jₖ* = 0` ✓  (induction hypothesis)

**Bottom block:**

```
aLₖ*Jₖ* + ωₖtₖyₖ* + bPₖ*Jₖ*
= tₖyₖ*LₖLₖ*Jₖ* + ωₖtₖyₖ* + tₖyₖ*PₖPₖ*Jₖ*
= tₖyₖ*(LₖLₖ* + PₖPₖ*)Jₖ* + ωₖtₖyₖ*
= tₖyₖ*Jₖ* + ωₖtₖyₖ*          (by Lemma B)
= tₖyₖ*(Jₖ* + ωₖI)
= tₖyₖ*(1+ωₖ)(I − yₖyₖ*)       (same algebra as base case)
= tₖ(1+ωₖ)(yₖ* − yₖ*yₖyₖ*)
= tₖ(1+ωₖ)(yₖ* − yₖ*) = 0  ✓
```

---

### Completing the inductive step

**(1) Off-diagonal = 0:**

```
A_{k−1}a* + B_{k−1}b* = tₖ(A_{k−1}Lₖ* + B_{k−1}Pₖ*)yₖ = tₖ · 0 · yₖ = 0  ✓
```

(by Lemma A)

**(2) Diagonal = 1:**

```
aa* + bb* = tₖ²yₖ*LₖLₖ*yₖ + tₖ²yₖ*PₖPₖ*yₖ
          = tₖ²yₖ*(LₖLₖ* + PₖPₖ*)yₖ
          = tₖ²yₖ*Iₘyₖ          (by Lemma B)
          = tₖ²‖yₖ‖² = tₖ²
```

Therefore: `aa* + bb* + |ωₖ|² = tₖ² + |ωₖ|² = (1 − |ωₖ|²) + |ωₖ|² = 1  ✓`

**QED.**

---

## Summary of Proof Structure

The entire proof rests on two lemmas that are proven by induction simultaneously:

| Lemma | Statement | Key identity used |
|-------|-----------|-------------------|
| **B** | LₖLₖ* + PₖPₖ* = Iₘ | JₖJₖ* = Iₘ − tₖ² yₖyₖ* |
| **A** | A_{k−1}Lₖ* + B_{k−1}Pₖ* = 0 | Jₖ* + ωₖI = (1+ωₖ)(I − yₖyₖ*) → yₖ*(I − yₖyₖ*) = 0 |

Both lemmas only require **‖yₖ‖ = 1**. More precisely:

- **TIB condition** (AₖAₖ* + BₖBₖ* = Iₖ): follows from ‖yₖ‖ = 1 alone.
- **Pole locations** {ω₁, ..., ωₖ}: these sit on the diagonal of the block lower-triangular Aₖ regardless of how yₖ is chosen.
- **Correct transfer function**: this is what the Blaschke-Potapov construction actually guarantees. The null basis pairs {ωₖ, vₖ} describe a target MIMO H(z). The BP factors "subtract out" the contribution of the already-placed poles ω₁…ω_{k−1}, mapping vₖ → yₖ such that the running (A, B) realises exactly H(z). Using any other normalised yₖ still gives a valid TIB system with the right poles, but the B matrix would encode different input directions — a different transfer function.

---

## Connection to genesis-alpha

- `mimo_null_basis()` in `ga/filters/tib.py` implements Algorithm 2: `cumJ` = Pₖ, `y[i]` = yₖ, `J[i]` = Jₖ, `L` = Lₖ.
- This proof confirms `mimo_null_basis()` produces an exactly TIB (A, B) when given exact normalised null basis vectors vₖ.
- `extract_poles_and_nullvecs_from_bt()` approximates vₖ from BT output — since ‖vₖ‖ = 1 is maintained by the normalisation step, Lemma A and B still hold structurally, but the resulting (A, B) represents a different system than the original BT reduction.

## Related Pages

[[concepts/TIBForm]], [[concepts/InformationGeometry]], [[papers/Kong2018-TIBInfoGeometry]], [[papers/Mu2026-ModelReductionNotes]]
