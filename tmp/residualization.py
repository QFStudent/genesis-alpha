import os
import sys
from pathlib import Path

import fire
import logbook
import numpy as np
import pandas as pd
from logbook import Logger, StreamHandler
from tqdm import tqdm

from lib.sft_ym.cds.backtesting.config import load_config
from lib.sft_ym.cds.backtesting.data.data import calc_ret

_instruments = ["ES", "YM", "NQ", "MI"]

_logger = Logger("Residualization")
log_handler = StreamHandler(sys.stdout)
log_handler.push_application()


def residualization(config: str, symbol: str, field: str):
    # load config
    cfg = load_config(config)

    # load minute bars
    df = pd.read_parquet(os.path.join(cfg.input_prefix, symbol + "_adjusted_price.parq"))
    ret = calc_ret(df[field])

    # load index return
    index_path = os.path.join(cfg.input_prefix, "_".join(_instruments) + field + "_minute_ret" + ".parq")
    _logger.info(f"Loading Minute Bar Index Returns from {index_path}.")
    idx_ret = pd.read_parquet(index_path)

    # load beta estimation
    beta_path = os.path.join(cfg.input_prefix, "_".join(_instruments), symbol + "_beta" + ".parq")
    _logger.info(f"Loading Beta Estimation from {beta_path}.")
    beta = pd.read_parquet(beta_path)

    residuals = []
    for row in tqdm(beta.iterrows()):
        if not np.isnan(row[1].values[0]):
            date = row[0].strftime("%Y-%m-%d")

            # align return and index
            daily_df = pd.concat([ret.loc[date], idx_ret.loc[date, "IDX"].squeeze() * row[1].values[0]], axis=1)
            daily_df.fillna(0.0, inplace=True)

            daily_residuals = daily_df[field] - daily_df["IDX"]

            residuals.append(daily_residuals)
        else:
            continue

    residuals_df = pd.concat(residuals, axis=0)

    Path(cfg.output_prefix).mkdir(parents=True, exist_ok=True)  # create directory if it doesn't exist
    residual_path = os.path.join(cfg.output_prefix, symbol + field + "_beta_residuals" + ".parq")

    _logger.info(f"Save residuals for {symbol} to {residual_path}.")
    pd.DataFrame(residuals_df).to_parquet(residual_path)


if __name__ == "__main__":
    with logbook.StreamHandler(sys.stdout).applicationbound():
        fire.Fire(residualization)
