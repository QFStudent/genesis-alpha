import logging
from datetime import date, datetime, time
from typing import List

import numpy as np
import pandas as pd
from tqdm import tqdm

from lib.sft_ym.cds.backtesting.backtesting import BackTesting, create_multi_index_data
from lib.sft_ym.cds.backtesting.data.data import calc_ret

_format = "%Y-%m-%d %H:%M:%S%z"

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s:%(name)s:%(message)s")
file_handler = logging.FileHandler("/home/yu.mu/yuresearch/lib/sft_ym/cds/log/sampler.log")
file_handler.setFormatter(formatter)
_logger.addHandler(file_handler)


class ClockTimeSampler:
    """
    Sampling data based on the clock time
    """

    def __init__(self, freq: str, offset: int = 0):
        self.freq = freq
        self.offset = offset

    def fit(self, df: pd.DataFrame):
        if not isinstance(df.index, pd.DatetimeIndex):
            raise ValueError("Data index needs to be DatetimeIndex.")

        if isinstance(df, pd.Series):
            df = pd.DataFrame(df)

        if not "DATE" in df.columns:
            df["DATE"] = df.index.date

        # subsample using groupby
        grouper = df.groupby(
            [
                pd.Grouper(key="DATE"),
                pd.Grouper(
                    level="DATETIME", freq=self.freq, origin="start", closed="right", label="right", offset=pd.Timedelta(self.offset, "min")
                ),
            ]
        )

        # here we can also use last to calc ohlc to reduce high freq effects.
        sub_df = pd.concat(
            [
                grouper["OPEN"].first().rename("OPEN"),
                grouper["HIGH"].max().rename("HIGH"),
                grouper["LOW"].min().rename("LOW"),
                grouper["LAST"].last().rename("CLOSE"),
                grouper["VOLUME"].sum(),
            ],
            axis=1,
        )
        sub_df = sub_df.reset_index().set_index("DATETIME")[["OPEN", "HIGH", "LOW", "CLOSE", "VOLUME"]]
        sub_df["RET"] = calc_ret(sub_df["CLOSE"])
        return sub_df


class ClockTimeSampler_cc_data78:
    """
    Sampling data based on the clock time
    """

    def __init__(self, freq: str, offset: int = 0):
        self.freq = freq
        self.offset = offset

    def fit(self, df: pd.DataFrame):
        if not isinstance(df.index, pd.DatetimeIndex):
            raise ValueError("Data index needs to be DatetimeIndex.")

        if not "DATE" in df.columns:
            df["DATE"] = df.index.date

        # subsample using groupby
        grouper = df.groupby(
            [
                pd.Grouper(key="DATE"),
                pd.Grouper(freq=self.freq, origin="start", closed="right", label="right", offset=pd.Timedelta(self.offset, "min")),
            ]
        )

        # here we can also use last to calc ohlc to reduce high freq effects.
        sub_df = pd.concat(
            [
                grouper["LAST"].last().rename("LAST"),
                grouper["MID"].last().rename("MID"),
                grouper["WEIGHTED"].last().rename("WEIGHTED"),
                grouper["HS"].last().rename("HS"),
            ],
            axis=1,
        )

        sub_df = sub_df.reset_index().set_index("DATETIME")[["LAST", "MID", "WEIGHTED", "HS"]]
        return sub_df


class ClockTimeSampler_complete:
    """
    Sampling data based on the clock time
    """

    def __init__(self, freq: str, offset: int = 0):
        self.freq = freq
        self.offset = offset

    def fit(self, df: pd.DataFrame):
        if not isinstance(df.index, pd.DatetimeIndex):
            raise ValueError("Data index needs to be DatetimeIndex.")

        if not "DATE" in df.columns:
            df["DATE"] = df.index.date

        # subsample using groupby
        grouper = df.groupby(
            [
                pd.Grouper(key="DATE"),
                pd.Grouper(freq=self.freq, origin="start", closed="right", label="right", offset=pd.Timedelta(self.offset, "min")),
            ]
        )

        # here we can also use last to calc ohlc to reduce high freq effects.
        sub_df = pd.concat(
            [
                grouper["LAST"].last().rename("LAST"),
                grouper["MID"].last().rename("MID"),
                grouper["WEIGHTED"].last().rename("WEIGHTED"),
                grouper["HS"].last().rename("HS"),
                grouper["VOLUME"].sum().rename("VOLUME"),
                grouper["HIGH"].max().rename("HIGH"),
                grouper["LOW"].min().rename("LOW"),
                grouper["OPEN"].first().rename("OPEN"),
                # grouper[["WEIGHTED", "VOLUME"]].apply(lambda x: (x["WEIGHTED"] * x["VOLUME"]).sum() / x["VOLUME"].sum()).rename("VWAP"),
                grouper["VWAP"].last().rename("VWAP1S"),   # one second vwap
                grouper["VWAP1"].last().rename("VWAP1M"),   # one min vwap
            ],
            axis=1,
        )

        sub_df = sub_df.reset_index().set_index("DATETIME")[
            ["OPEN", "HIGH", "LOW", "LAST", "VOLUME", "MID", "WEIGHTED", "HS", "VWAP1S", "VWAP1M"]
        ]
        return sub_df


class ClockTimeSampler_24HOUR:
    """
    Sampling data based on the clock time
    """

    def __init__(self, freq: str, offset: int = 0):
        self.freq = freq
        self.offset = offset

    def fit(self, df: pd.DataFrame):
        if not isinstance(df.index, pd.DatetimeIndex):
            raise ValueError("Data index needs to be DatetimeIndex.")

        if not "DATE" in df.columns:
            df["DATE"] = df.index.date

        # subsample using groupby
        grouper = df.groupby(
            [
                pd.Grouper(key="DATE"),
                pd.Grouper(freq=self.freq, origin="start", closed="right", label="right", offset=pd.Timedelta(self.offset, "min")),
            ]
        )

        # here we can also use last to calc ohlc to reduce high freq effects.
        sub_df = pd.concat(
            [
                grouper["LAST"].last().rename("LAST"),
                grouper["MID"].last().rename("MID"),
                grouper["WEIGHTED"].last().rename("WEIGHTED"),
                grouper["HS"].last().rename("HS"),
                grouper["VOLUME"].sum().rename("VOLUME"),
                grouper["HIGH"].max().rename("HIGH"),
                grouper["LOW"].min().rename("LOW"),
                grouper["OPEN"].first().rename("OPEN"),
                # grouper[["WEIGHTED", "VOLUME"]].apply(lambda x: (x["WEIGHTED"] * x["VOLUME"]).sum() / x["VOLUME"].sum()).rename("VWAP"),
                grouper["VWAP"].last().rename("VWAP1S"),   # one second vwap
                grouper["VWAP1"].last().rename("VWAP1M"),   # one min vwap
            ],
            axis=1,
        )

        sub_df = sub_df.reset_index().set_index("DATETIME")[
            ["OPEN", "HIGH", "LOW", "LAST", "VOLUME", "MID", "WEIGHTED", "HS", "VWAP1S", "VWAP1M"]
        ]
        return sub_df


class ClockTimeSampler_Liquid:
    """
    Sampling data based on the clock time
    """

    def __init__(self, freq: str, offset: int = 0):
        self.freq = freq
        self.offset = offset

    def fit(self, df: pd.DataFrame):
        if not isinstance(df.index, pd.DatetimeIndex):
            raise ValueError("Data index needs to be DatetimeIndex.")

        if not "DATE" in df.columns:
            df["DATE"] = df.index.date

        df["TIME"] = df.index.time

        # remove iliquid periods
        df = df.groupby(pd.Grouper(key="DATE")).apply(lambda x: x[(x["TIME"] <= time(15, 0)) & (x["TIME"] >= time(11, 0))])
        df = df.drop(columns="DATE").reset_index().set_index("DATETIME")

        # subsample using groupby
        grouper = df.groupby(
            [
                pd.Grouper(key="DATE"),
                pd.Grouper(freq=self.freq, origin="start", closed="right", label="right", offset=pd.Timedelta(self.offset, "min")),
            ]
        )

        # here we can also use last to calc ohlc to reduce high freq effects.
        sub_df = pd.concat(
            [
                grouper["LAST"].last().rename("LAST"),
                grouper["MID"].last().rename("MID"),
                grouper["WEIGHTED"].last().rename("WEIGHTED"),
                grouper["HS"].last().rename("HS"),
                grouper["VOLUME"].sum().rename("VOLUME"),
                grouper["HIGH"].max().rename("HIGH"),
                grouper["LOW"].min().rename("LOW"),
                grouper["OPEN"].first().rename("OPEN"),
                # grouper[["WEIGHTED", "VOLUME"]].apply(lambda x: (x["WEIGHTED"] * x["VOLUME"]).sum() / x["VOLUME"].sum()).rename("VWAP"),
                grouper["VWAP"].last().rename("VWAP1S"),   # one second vwap
                grouper["VWAP1"].last().rename("VWAP1M"),   # one min vwap
            ],
            axis=1,
        )

        sub_df = sub_df.reset_index().set_index("DATETIME")[
            ["OPEN", "HIGH", "LOW", "LAST", "VOLUME", "MID", "WEIGHTED", "HS", "VWAP1S", "VWAP1M"]
        ]
        return sub_df


class SpikeTimeSampler:
    """
    Sampling data based on the spike time
    """

    def __init__(self, spike_time: List[str], offset: int):
        self.spike_time = [datetime.strptime(x, "%H:%M:%S").time() for x in spike_time]
        self.offset = offset

    def fit(self, df: pd.DataFrame):
        if isinstance(df, pd.Series):
            df = pd.DataFrame(df)

        if not isinstance(df.index, pd.DatetimeIndex):
            raise ValueError("Data index needs to be DatetimeIndex.")

        if not "DATE" in df.columns:
            df["DATE"] = df.index.date

        df["TIME"] = df.index.time
        df["INTERVAL"] = 0

        # remove iliquid periods
        df = df.groupby(pd.Grouper(key="DATE")).apply(
            lambda x: pd.concat([x[x["TIME"] <= time(17, 0)], x[x["TIME"] >= time(18, 0)]], axis=0)
        )

        df = df.drop(columns="DATE").reset_index().set_index("DATETIME")
        df = add_trade_day(df)

        def assign_interval():
            n = 1
            for i, st in enumerate(self.spike_time[:-1]):
                st_filt = (df["TIME"] > self.spike_time[i]) & (df["TIME"] <= self.spike_time[i + 1])
                df.loc[st_filt, "INTERVAL"] = n
                n += 1

        assign_interval()

        df = df.reset_index()
        grouper = df.groupby([pd.Grouper(key="TRADEDAY"), pd.Grouper(key="INTERVAL")])

        sub_df = pd.concat(
            [
                grouper["DATETIME"].last(),
                grouper["LAST"].last(),
                grouper["MID"].last(),
                grouper["WEIGHTED"].last(),
                grouper["HS"].last(),
                grouper["VOLUME"].sum().rename("VOLUME"),
                grouper["HIGH"].max().rename("HIGH"),
                grouper["LOW"].min().rename("LOW"),
                grouper["OPEN"].first(),
                grouper[["WEIGHTED", "VOLUME"]].apply(lambda x: (x["WEIGHTED"] * x["VOLUME"]).sum() / x["VOLUME"].sum()).rename("VWAP"),
                grouper["VWAP"].last().rename("VWAP1S"),   # one second vwap
                grouper["VWAP1"].last().rename("VWAP1M"),   # one min vwap
            ],
            axis=1,
        )
        sub_df = sub_df.set_index("DATETIME")[["OPEN", "HIGH", "LOW", "LAST", "VOLUME", "MID", "WEIGHTED", "HS", "VWAP1S", "VWAP1M"]]

        return sub_df


class PreprocessSampler:
    """
    Sampling data after preprocessing
    """

    def __init__(self, freq: str, offset: int = 0) -> None:
        self.freq = freq
        self.offset = offset

    def fit(self, df: pd.DataFrame):
        if not isinstance(df.index, pd.DatetimeIndex):
            raise ValueError("Data index needs to be DatetimeIndex.")

