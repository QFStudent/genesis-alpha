import base64
import io
import logging
import os
from dataclasses import dataclass
from datetime import datetime
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lib.sft_ym.dp.dp import calc_previous_bday
from lib.sft_ym.cds.backtesting.config import load_config as bt_load_config, Config
from pathlib import Path


_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s:%(name)s:%(message)s")


@dataclass
class ModelMatrices:
    coef_dict: dict
    z_dict: dict
    M_dict: dict
    N_dict: dict
    U_dict: dict
    Q_dict: dict
    scale_dict: dict


def load_model_matrices(date: datetime.date, cfg: Config, universe: List):
    """Load model matrices for both prod and research run."""
    previous_date = calc_previous_bday(date, 'ES')
    _logger.info(f"Load model matrices for {str(previous_date)}")
    strats_name = f"strategy={cfg.expr_name}"
    path = os.path.join(
        cfg.output_prefix,
        cfg.system,
        strats_name,
        'model_selection',
        '{sym}',
        str(previous_date.year),
        str(previous_date.month).zfill(2),
        str(previous_date.day).zfill(2),
    )

    _logger.info(f'Constructed prediction path is: {path}.')
    coef_dict = {}
    z_dict = {}
    M_dict = {}
    N_dict = {}
    U_dict = {}
    Q_dict = {}
    scale_dict = {}

    pre_date_prefix = previous_date.strftime("%Y%m%d")
    for sym in universe:
        pred_path_ = path.format(sym=sym)
        coef_dict[sym] = pd.read_csv(
            os.path.join(pred_path_, f"{pre_date_prefix}.coef.csv"), index_col=0
        ).squeeze().values
        z_dict[sym] = pd.read_csv(
            os.path.join(pred_path_, f"{pre_date_prefix}.z0.csv"), index_col=0
        ).squeeze().values
        M_dict[sym] = pd.read_csv(
            os.path.join(pred_path_, f"{pre_date_prefix}.M.csv"), index_col=0
        ).values
        N_dict[sym] = pd.read_csv(
            os.path.join(pred_path_, f"{pre_date_prefix}.N.csv"), index_col=0
        ).values
        U_dict[sym] = pd.read_csv(
            os.path.join(pred_path_, f"{pre_date_prefix}.U.csv"), index_col=0
        ).values
        Q_dict[sym] = pd.read_csv(
            os.path.join(pred_path_, f"{pre_date_prefix}.Q.csv"), index_col=0
        ).values
        scale_dict[sym] = pd.read_csv(
            os.path.join(pred_path_, f"{pre_date_prefix}.scale.csv"), index_col=0
        )

    return ModelMatrices(
        coef_dict=coef_dict,
        z_dict=z_dict,
        M_dict=M_dict,
        N_dict=N_dict,
        U_dict=U_dict,
        Q_dict=Q_dict,
        scale_dict=scale_dict,
    )


def check_coef_consistency(prod_sys: ModelMatrices, res_sys: ModelMatrices, universe: List):
    coef_err = {}
    z0_err = {}
    scale_err = {}
    for sym in universe:
        coef_err[sym] = np.abs(prod_sys.coef_dict[sym] - res_sys.coef_dict[sym]).sum()
        z0_err[sym] = np.abs(prod_sys.z_dict[sym] - res_sys.z_dict[sym]).sum()
        scale_err[sym] = np.abs(prod_sys.scale_dict[sym] - res_sys.scale_dict[sym]).sum()
    coef_err_df = pd.DataFrame.from_dict(coef_err, orient="index", columns=["coef"])
    z0_err_df = pd.DataFrame.from_dict(z0_err, orient="index", columns=["z0"])
    scale_err_df = pd.DataFrame.from_dict(scale_err, orient="index")
    return pd.concat([coef_err_df, z0_err_df, scale_err_df], axis=1)


def check_coefficients_and_pred_matches(date: datetime.date, rcfg: Config, pcfg: Config):
    """Reconcile model coefficients and predictions between research and prod."""
    if rcfg.universe != pcfg.universe:
        raise ValueError('Research universe is not the same as production universe')

    res_sys = load_model_matrices(date, rcfg, rcfg.universe)
    prod_sys = load_model_matrices(date, pcfg, pcfg.universe)
    err_df = check_coef_consistency(prod_sys, res_sys, rcfg.universe)
    _logger.info(f'Coefficient consistency check:\n{err_df}')

    if not np.allclose(err_df.values, 0.0):
        raise ValueError('research coefficients differentiate from production coefficients.')

    _logger.info(f'Research coefficients match prod for date {str(date)}.')


def check_adjusted_price(date: datetime.date, rcfg: Config, pcfg: Config):
    """Reconcile adjusted price for research and prod are the same"""
    if rcfg.universe != pcfg.universe:
        raise ValueError('Research universe is not the same as production universe')

    date_format = date.strftime("%Y%m%d")
    for sym in rcfg.universe:
        research_path = os.path.join(
            rcfg.input_prefix,
            sym,
            date_format,
            f"{sym}_{date_format}_adjusted_price.csv",
        )
        prod_path = os.path.join(
            pcfg.input_prefix,
            sym,
            date_format,
            f"{sym}_{date_format}_adjusted_price.csv",
        )
        r_df = pd.read_csv(research_path, index_col=0)
        p_df = pd.read_csv(prod_path, index_col=0)
        r_df.index = pd.to_datetime(r_df.index)
        p_df.index = pd.to_datetime(p_df.index)
        if not np.allclose(r_df['WEIGHTED'], p_df['WEIGHTED']):
            raise ValueError('research long bar differentiate from production long bar.')

    _logger.info(f'Research data matches prod data for date {str(date)}.')


def main():
    rec_list = ["yu.mu@centivacapital.com"]
    date = datetime.strptime('20241012', '%Y%m%d').date()

    # set up log files
    log_path = f'/ccny5/dat/strat/sft/yu/cds/log/directional/{date.strftime("%Y%m%d")}'
    Path(log_path).mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(os.path.join(log_path, "research_consistency.log"))
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    _logger.addHandler(file_handler)
    _logger.addHandler(stream_handler)

    previous_date = calc_previous_bday(date, 'ES')

    # rconfig = '/home/yu.mu/yuresearch/lib/sft_ym/cds/prod/config/sim_research_consistency.toml'
    rconfig = '/home/yu.mu/yuresearch/lib/sft_ym/cds/prod/config/direct_research_consistency.toml'
    pconfig = '/home/yu.mu/yuresearch/lib/sft_ym/cds/prod/config/direct_prod_dev1.toml'

    rcfg = bt_load_config(rconfig)
    pcfg = bt_load_config(pconfig)

    # check_adjusted_price(previous_date, rcfg, pcfg)

    check_coefficients_and_pred_matches(previous_date, rcfg, pcfg)


if __name__ == "__main__":
    main()
