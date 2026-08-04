"""Tests for ga.vc.app CLI."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from ga.vc.simulate import simulate_vc_model


def test_cli_end_to_end(tmp_path: Path):
    sim = simulate_vc_model(T=35, n=2, sigma_u=1.0, sigma_v=np.array([0.1, 0.2]), seed=5)
    x_path = tmp_path / "x.csv"
    y_path = tmp_path / "y.csv"
    out_path = tmp_path / "out"
    pd.DataFrame(sim["x"], columns=["x0", "x1"]).to_csv(x_path, index=False)
    pd.DataFrame({"y": sim["y"]}).to_csv(y_path, index=False)

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "ga.vc.app",
            "--x",
            str(x_path),
            "--y",
            str(y_path),
            "--output",
            str(out_path),
            "--fixed-ratios",
            "--variance-ratios",
            "1.0",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    summary = json.loads(proc.stdout)
    assert summary["T"] == 35
    assert summary["n"] == 2
    assert (out_path / "coefficients.csv").exists()
    assert (out_path / "summary.json").exists()

    coeff = pd.read_csv(out_path / "coefficients.csv").to_numpy()
    assert coeff.shape == (35, 2)
