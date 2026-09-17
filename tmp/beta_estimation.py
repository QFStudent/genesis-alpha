import os
import sys
from pathlib import Path

import fire
import logbook
import pandas as pd

from lib.sft_ym.cds.backtesting.config import load_config
from lib.sft_ym.cds.backtesting.data.data import calc_ret

_instruments = ["ES", "YM", "NQ", "MI"]


def gen_synthetic_index(config: str, field: str):
    """
    Generate synthetic index from the adjusted minute bar data.
    """
    # load config
    cfg = load_config(config)

    # load adjusted price from config and calculate returns from it
    dfs = []
    for inst in _instruments:
        df = pd.read_parquet(os.path.join(cfg.input_prefix, inst + "_adjusted_price.parq"))
        ret = calc_ret(df[field])
        dfs.append(ret)

    # calculate minute bar returns and minute bar index
    minute_rets = pd.concat(dfs, axis=1)
    minute_rets.fillna(0.0, inplace=True)

    index = minute_rets.mean(axis=1)
    minute_rets["IDX"] = index
    minute_rets.columns = _instruments + ["IDX"]

    # Calculate daily return from minute return
    daily_rets = minute_rets.groupby(minute_rets.index.date).sum()
    daily_rets["IDX"] = daily_rets.mean(axis=1)
    daily_rets.columns = _instruments + ["IDX"]
    daily_rets.index = pd.to_datetime(daily_rets.index)  # index becomes index after groupby

    minute_rets.to_parquet(os.path.join(cfg.input_prefix, "_".join(_instruments) + field + "_minute_ret.parq"))
    daily_rets.to_parquet(os.path.join(cfg.input_prefix, "_".join(_instruments) + field + "_daily_ret.parq"))


def get_beta(config: str, instrument: str, field: str) -> None:
    cfg = load_config(config)

    daily_ret = pd.read_parquet(os.path.join(cfg.input_prefix, "_".join(_instruments) + field + "_daily_ret.parq"))

    # calculate rolling window 1 year beta
    cov_df = daily_ret["IDX"].rolling(window=252).cov(daily_ret[instrument]).shift(1)  # shift one step to remove beta lookahead
    var_df = daily_ret["IDX"].rolling(window=252).var().shift(1)
    Beta = 0.3 + 0.7 * (cov_df / var_df)

    upload_path = os.path.join(cfg.input_prefix, "_".join(_instruments))
    Path(upload_path).mkdir(parents=True, exist_ok=True)  # create directory if it doesn't exist
    pd.DataFrame(Beta).to_parquet(os.path.join(upload_path, instrument + field + "_beta" + ".parq"))


if __name__ == "__main__":
    with logbook.StreamHandler(sys.stdout).applicationbound():
        fire.Fire(
            {
                "get_beta": get_beta,
                "gen_synthetic_index": gen_synthetic_index,
            }
        )
