import pickle
from typing import List

import numpy as np
import pandas as pd
import logging

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('/home/yu.mu/yuresearch/lib/sft_ym/cds/log/data.log')
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
_logger.addHandler(file_handler)
_logger.addHandler(stream_handler)


def rolling_winsorization(df: pd.DataFrame,
                          lookback: int,
                          limits: List,
                          shift_size: int) -> pd.DataFrame:
    """
    Rolling window winsorization for the data

    Parameters
    ----------
    df: univariate DataFrame, shape = (n, )
    lookback: lookback window
    """
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    df.fillna(0.0, inplace=True)
    lb = df.rolling(lookback).quantile(limits[0], interpolation="nearest").shift(shift_size)
    ub = df.rolling(lookback).quantile(limits[1], interpolation="nearest").shift(shift_size)

    large_ol_mask = (df > ub)[df.columns[0]]
    small_ol_mask = (df < lb)[df.columns[0]]

    df.loc[large_ol_mask] = ub[large_ol_mask]
    df.loc[small_ol_mask] = lb[small_ol_mask]
    df = df.iloc[(lookback+shift_size):]
    return df


def calc_ret(price: pd.Series) -> pd.Series:
    """
    Calculate return process from price
    """
    ret = np.log(price).diff()
    ret.fillna(0.0, inplace=True)
    return ret


def preprocess(y: pd.DataFrame, std_lookback: int, winsor_lookback: int, winsor_limits: List[float], preprocess_shift) -> pd.DataFrame:
    """
    Preprocess data using rolling window winsorization
    and rolling window standardization to improve
    stationarity of our data.
    """
    for col in y.columns:
        if y[col].isna().sum() > 0:
            _logger.info(f"There are {y[col].isna().sum()} Nans.")
            y[col].fillna(0.0, inplace=True)

        std = y[col].rolling(std_lookback).std().shift(preprocess_shift)
        _logger.info(f"Rolling standardizing with shift {preprocess_shift} for column {col}.")
        y[col] = y[col] / std
        y[col + "_std"] = std
        y[col].fillna(0.0, inplace=True)
    return y


def write_pickle(data: dict, path: str) -> None:
    with open(path, "wb") as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)


def read_pickle(path: str) -> dict:
    dbfile = open(path, "rb")
    return pickle.load(dbfile)


def roll_preprocess(y: pd.Series, lookback: int, limits: List,
                    field: str='WEIGHTED', shift_size: int=1):
    """
    preprocess data in rolling manner by
    using both winsorization and standardization
    """
    y_ = rolling_winsorization(y, lookback, limits, shift_size)

    # if y_.squeeze().isna().any():
    #     raise ValueError('Data has Nan Values, rolling std will fail.')

    # std = y_.rolling(lookback).std().shift(shift_size)     # first lookback + shift_size will be NaNs
    std = y_.rolling(lookback).std()      # first lookback + shift_size will be NaNs

    if (std.squeeze() == 0).any():
        import pdb; pdb.set_trace()

        raise ValueError('Rolling std is zero.')

    y_ = y_ / std
    y_.fillna(0.0, inplace=True)
    y_ = y_.iloc[lookback:]

    return y_[field], std


def vol_curve_preprocess(y: pd.Series, lookback: int, limits: List,
                         field: str="WEIGHTED", shift_size: int=1):
    """
    Preprocess by winsorization and normalizing by vol curve
    """
    name = y.name
    y_ = rolling_winsorization(y, lookback, limits, shift_size)
    y_['TIME'] = y_.index.time
    y_new = y_.reset_index()
    y_table = y_new.pivot(index="DATETIME", values=name, columns="TIME")
    vol_curve = y_table.std(axis=0)
    import pdb; pdb.set_trace()

    print("done")


def expanding_window_preprocess(train: pd.Series, lookback: int, limits: List,
                                field: str='WEIGHTED', shift_size: int=1):
    """
    This is for expanding window CV preprocess, within each window, we also
    want to do rolling preprocess so that this will be consistent with the
    alpha evaluation step
    """
    pass


def preprocess_overlapping_data(dp: pd.DataFrame, lookback: int, limits: List) -> pd.DataFrame:
    ys = []
    for col in dp.columns:
        y_ = rolling_winsorization(dp[col], lookback, limits, 1)
        ys.append(y_)

    ys_df = pd.concat(ys, axis=1)
