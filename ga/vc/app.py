#!/usr/bin/env python3
"""
Standalone CLI for Schlicht/Ludsteck VC estimation.

Example:
    python -m ga.vc.app --y data/y.csv --x data/x.csv --output results/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ga.vc.schlicht import vc_estimate


def _load_matrix(path: Path) -> np.ndarray:
    df = pd.read_csv(path)
    return df.to_numpy(dtype=float)


def _load_vector(path: Path) -> np.ndarray:
    df = pd.read_csv(path)
    if df.shape[1] == 1:
        return df.iloc[:, 0].to_numpy(dtype=float)
    raise ValueError(f"{path}: y file must have exactly one column")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Schlicht (1989) time-varying coefficients (VC) estimation",
    )
    parser.add_argument("--x", type=Path, required=True, help="Regressors CSV (T rows, n cols)")
    parser.add_argument("--y", type=Path, required=True, help="Dependent variable CSV (T rows, 1 col)")
    parser.add_argument("--output", type=Path, required=True, help="Output directory")
    parser.add_argument(
        "--variance-ratios",
        type=float,
        default=1.0,
        help="Initial variance ratio r = Var(v)/Var(u) (scalar or use --variance-ratios-file)",
    )
    parser.add_argument("--variance-ratios-file", type=Path, default=None)
    parser.add_argument(
        "--fixed-ratios",
        action="store_true",
        help="Skip ML optimization; use variance ratios as fixed",
    )
    args = parser.parse_args(argv)

    x = _load_matrix(args.x)
    y = _load_vector(args.y)
    if args.variance_ratios_file is not None:
        vr = _load_vector(args.variance_ratios_file)
    else:
        vr = args.variance_ratios

    result = vc_estimate(
        x,
        y,
        variance_ratios=vr,
        optimize_ratios=not args.fixed_ratios,
    )

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(result.coeff, columns=[f"b{i}" for i in range(result.coeff.shape[1])]).to_csv(
        out / "coefficients.csv", index=False
    )
    pd.DataFrame(result.sdb, columns=[f"se_b{i}" for i in range(result.sdb.shape[1])]).to_csv(
        out / "coefficient_std.csv", index=False
    )
    pd.DataFrame({"y_hat": (x * result.coeff).sum(axis=1), "y": y}).to_csv(
        out / "fit.csv", index=False
    )

    summary = {
        "sdu": result.sdu,
        "sdi": result.sdi.tolist(),
        "variance_ratios": result.variance_ratios.tolist(),
        "residual_var": result.residual_var,
        "criterion": result.criterion,
        "T": int(x.shape[0]),
        "n": int(x.shape[1]),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
