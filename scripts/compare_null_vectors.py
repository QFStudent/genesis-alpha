"""
Canonical vs random null vectors in a TIB realization.

Demonstrates what the null-vector choice actually controls in
`null_basis_realization`: the *input coupling structure* of the MIMO system,
not its validity. Both choices yield a valid, stable, balanced TIB system
(A A* + B B* = I); they differ in *which* transfer function gets realized.

Canonical null vectors (`mimo_standard_null_vectors`) assign each decay mode
to a single input channel in round-robin order (mode k <- channel k % q), so
input coupling is block-sparse and channels are decoupled at the input.
Random unit null vectors mix every channel into every mode -> dense coupling.

Run:  python scripts/compare_null_vectors.py
Plots are written next to this script if matplotlib is available; the
numerical summary always prints to stdout.
"""

import numpy as np

from ga.filters.tib import (
    poles_chebyshev_roots,
    mimo_standard_null_vectors,
    null_basis_realization,
    mimo_poles_to_ir,
    krylov_basis,
)
from ga.reducers.balanced_truncation import (
    balanced_truncation,
    extract_poles_and_nullvecs_from_bt,
    bt_impulse_response,
    tib_from_state_space,
)

P = 12          # number of poles / states
Q = 6           # number of input channels  (q > 5: the MIMO regime)
P_OUT = 4       # number of outputs (for the fitted case)
LAGS = 40       # impulse-response length
SEED = 0


def canonical_null_vectors(p: int, q: int) -> np.ndarray:
    """(q, p) null-vector matrix: mode k driven by channel k % q (round-robin)."""
    return mimo_standard_null_vectors(p, q).T          # (p, q) -> (q, p)


def random_null_vectors(p: int, q: int, seed: int = SEED) -> np.ndarray:
    """(q, p) matrix of real unit null vectors, one per pole."""
    rng = np.random.default_rng(seed)
    v = rng.standard_normal((q, p))
    return v / np.linalg.norm(v, axis=0)


def tib_residual(sys) -> float:
    """max |A A* + B B* - I| — should be ~0 for any unit null vectors."""
    A, B = sys.A(dense=True), sys.B(dense=True)
    resid = A @ A.conj().T + B @ B.conj().T - np.eye(A.shape[0])
    return float(np.max(np.abs(resid)))


def cross_channel_energy_fraction(ir: np.ndarray, q: int) -> float:
    """
    Fraction of impulse-response energy that lands *off* a mode's assigned
    channel. With C = I, output k is state k; canonical assigns state k to
    channel k % q, so any energy in ir[k, j != k % q, :] is cross-channel.
    """
    p = ir.shape[0]
    assigned = np.arange(p) % q
    total = float(np.sum(np.abs(ir) ** 2))
    on = sum(float(np.sum(np.abs(ir[k, assigned[k], :]) ** 2)) for k in range(p))
    return (total - on) / total


def b_nonzeros(sys, tol: float = 1e-10) -> int:
    return int(np.sum(np.abs(sys.B(dense=True)) > tol))


def fit_output_matrix(A: np.ndarray, B: np.ndarray,
                      target_ir: np.ndarray, lags: int) -> np.ndarray:
    """
    Best-case IR for a given (A, B): least-squares output matrix C.

    The impulse response is linear in C  (ir[o,j,t] = sum_s C[o,s] (A^t B)[s,j]),
    so for fixed (A, B) the optimal C is a linear least-squares solve. Fitting C
    this way isolates the effect of the *null vectors*: it asks "with the best
    possible output map, how well can these input directions reproduce the
    target?" — so any residual error is due to (A, B), i.e. the null vectors.
    """
    K = np.asarray(krylov_basis(A, B, lags), dtype=float)   # (P, Q, lags)
    p = A.shape[0]
    M = K.reshape(p, -1)                                    # (P, Q*lags)
    Y = target_ir.reshape(target_ir.shape[0], -1)           # (P_OUT, Q*lags)
    C = np.linalg.lstsq(M.T, Y.T, rcond=None)[0].T          # (P_OUT, P)
    return (C @ M).reshape(target_ir.shape)


def rel_error(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-15))


def main() -> None:
    poles = poles_chebyshev_roots(P)           # real, decaying modes in (0, 1)
    C = np.eye(P)                              # outputs = states, so we can read
                                               # per-mode input coupling directly

    cases = {
        "canonical": canonical_null_vectors(P, Q),
        "random": random_null_vectors(P, Q),
    }

    systems, irs = {}, {}
    print(f"TIB system: P={P} poles, Q={Q} input channels, {LAGS} lags\n")
    print(f"{'case':<10} {'TIB resid':>12} {'B nnz':>8} {'B nnz max':>10} "
          f"{'cross-chan E':>14}")
    print("-" * 58)
    for name, v in cases.items():
        sys = null_basis_realization(poles, v)
        ir = mimo_poles_to_ir(poles, v, C, LAGS)       # (P, Q, LAGS)
        systems[name], irs[name] = sys, ir
        print(f"{name:<10} {tib_residual(sys):>12.2e} {b_nonzeros(sys):>8d} "
              f"{P * Q:>10d} {cross_channel_energy_fraction(ir, Q):>13.1%}")

    print(
        "\nReading:\n"
        "  * TIB resid ~0 for BOTH -> both are valid balanced realizations.\n"
        f"  * Canonical B has {P} nonzeros (1 per mode); random fills all {P * Q}.\n"
        "  * Canonical cross-channel energy ~0 (each mode lives on one channel);\n"
        "    random spreads each mode across all channels.\n"
        "  -> The null-vector choice sets WHICH system you realize, not its quality."
    )

    fit_irs = fit_comparison()
    _maybe_plot(systems, irs, fit_irs)


def fit_comparison() -> dict:
    """
    PART 2 — when there IS a target: which null vectors reproduce it?

    Ground truth: an order-P TIB system with *non-axis-aligned* (cross-coupled)
    null vectors and a random output matrix. Its impulse response is the "data".

    Two reference points establish what is achievable:
        oracle    -> the TRUE null vectors  (target is exactly TIB-representable;
                     also validates the least-squares C fit)
        BT direct -> balanced_truncation's own reduced IR (the data is fittable)

    Then, holding the poles fixed, we vary ONLY the null vectors and fit C
    optimally each time:
        canonical -> axis-aligned input structure
        random    -> arbitrary input structure
        fitted    -> null vectors recovered from the data via
                     `extract_poles_and_nullvecs_from_bt`
    """
    poles = poles_chebyshev_roots(P)               # real poles (BT extraction
                                                   # only handles real poles)
    rng = np.random.default_rng(7)
    v_true = random_null_vectors(P, Q, seed=7)     # true cross-coupled directions
    C_true = rng.standard_normal((P_OUT, P))
    target_ir = mimo_poles_to_ir(poles, v_true, C_true, LAGS)   # (P_OUT, Q, LAGS)

    # Recover (poles, null vectors) from the data via balanced truncation of the
    # (here, known) full-order system -- the genesis-alpha BT -> TIB pipeline.
    true_sys = null_basis_realization(poles, v_true)
    result = balanced_truncation(
        true_sys.A(dense=True), true_sys.B(dense=True), C_true, order=P,
    )
    poles_fit, nullvecs_fit = extract_poles_and_nullvecs_from_bt(result)
    pole_recovery = rel_error(np.sort(poles_fit), np.sort(poles))

    def tib_fit(v):
        sys = null_basis_realization(poles, v)     # poles fixed: only v varies
        return fit_output_matrix(sys.A(dense=True), sys.B(dense=True),
                                 target_ir, LAGS)

    fit_irs = {"target": target_ir}
    print("\n\n=== PART 2: which null vectors reproduce a target IR? ===")
    print(f"target: order-{P} system, {P_OUT} outputs, {Q} inputs, "
          f"cross-coupled true null vectors")
    print(f"BT pole recovery (sorted, rel err): {pole_recovery:.2e}\n")

    # Reference points: what is achievable.
    fit_irs["oracle"] = tib_fit(v_true)
    bt_ir = bt_impulse_response(result, LAGS)
    # The fix: re-coordinate the BT system into TIB form directly.
    corr_sys, transform = tib_from_state_space(result.A_r, result.B_r)
    corr_K = np.asarray(krylov_basis(corr_sys.A(dense=True),
                                     corr_sys.B(dense=True), LAGS))
    corr_ir = np.tensordot(result.C_r @ transform, corr_K, (-1, 0))
    fit_irs["corrected"] = corr_ir
    print(f"{'reference':<26} {'IR rel err':>12}")
    print("-" * 40)
    print(f"{'oracle (true v)':<26} {rel_error(fit_irs['oracle'], target_ir):>12.2e}")
    print(f"{'BT reduced (direct)':<26} {rel_error(bt_ir, target_ir):>12.2e}")
    print(f"{'corrected (tib_from_ss)':<26} {rel_error(corr_ir, target_ir):>12.2e}")

    # TIB realization from each null-vector choice (best-C fit).
    cases = {
        "canonical": canonical_null_vectors(P, Q),
        "random": random_null_vectors(P, Q, seed=99),
        "fitted (BT extract)": nullvecs_fit.T,     # (P, Q) -> (Q, P)
    }
    print(f"\n{'TIB null vectors':<22} {'IR fit rel err (best C)':>24}")
    print("-" * 48)
    for name, v in cases.items():
        fit_irs[name.split()[0]] = tib_fit(v)      # key: canonical/random/fitted
        print(f"{name:<22} {rel_error(fit_irs[name.split()[0]], target_ir):>24.2e}")

    print(
        "\nReading:\n"
        "  * oracle ~0 and BT-direct ~0: the target IS exactly representable, and\n"
        "    balanced truncation reproduces it -- so the data is perfectly fittable.\n"
        "  * corrected (tib_from_state_space) ~0: re-coordinating the BT system into\n"
        "    TIB form (input-balance + Schur) recovers the target exactly. THIS is\n"
        "    the right BT -> TIB path.\n"
        "  * BUT the three null-vector choices below are far off, and 'fitted' is NO\n"
        "    BETTER than canonical/random: extract_poles_and_nullvecs_from_bt returns\n"
        "    the orthonormalized y_k of a Schur (not input-balanced) realization, so\n"
        "    feeding them to null_basis_realization gives the WRONG transfer function\n"
        "    (the approximation gap in Mu2026-NullBasisProofs, quantified).\n"
        "  -> Use tib_from_state_space, not extract_poles_and_nullvecs_from_bt."
    )
    return fit_irs


def _maybe_plot(systems: dict, irs: dict, fit_irs: dict) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed - skipping plots)")
        return

    import os
    out_dir = os.path.dirname(os.path.abspath(__file__))

    # Figure 1: |B| coupling pattern, canonical vs random.
    fig, axes = plt.subplots(1, 2, figsize=(8, 5))
    for ax, name in zip(axes, systems):
        B = np.abs(systems[name].B(dense=True))
        im = ax.imshow(B, aspect="auto", cmap="viridis")
        ax.set_title(f"|B|  ({name})")
        ax.set_xlabel("input channel")
        ax.set_ylabel("mode (state)")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    f1 = os.path.join(out_dir, "compare_B_coupling.png")
    fig.savefig(f1, dpi=120)

    # Figure 2: IR of one mode across all channels, canonical vs random.
    mode = Q + 1                                   # a mode assigned to channel 1
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    for ax, name in zip(axes, irs):
        ir = irs[name]
        for j in range(Q):
            ax.plot(ir[mode, j, :], label=f"ch {j}")
        ax.set_title(f"IR of mode {mode}  ({name})")
        ax.set_xlabel("lag")
    axes[0].set_ylabel("response")
    axes[1].legend(fontsize=8, ncol=2)
    fig.tight_layout()
    f2 = os.path.join(out_dir, "compare_ir.png")
    fig.savefig(f2, dpi=120)

    # Figure 3: fitted vs target IR for one (output, input) pair.
    out, inp = 0, 0
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(fit_irs["target"][out, inp, :], "k-", lw=2.5, label="target")
    for name in ("canonical", "random", "fitted"):
        ax.plot(fit_irs[name][out, inp, :], "--", label=name)
    ax.set_title(f"best-C IR fit to target  (output {out}, input {inp})")
    ax.set_xlabel("lag")
    ax.set_ylabel("response")
    ax.legend()
    fig.tight_layout()
    f3 = os.path.join(out_dir, "compare_fit.png")
    fig.savefig(f3, dpi=120)

    print(f"\nPlots written:\n  {f1}\n  {f2}\n  {f3}")


if __name__ == "__main__":
    main()
