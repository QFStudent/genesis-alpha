import os 
import datetime
import io 
import pandas as pd 
import numpy as np 
from dataclasses import dataclass


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
    """
    Load model matrices for both prod and research run 
    """
    previous_date = calc_previous_bday(date, 'ES')
    _logger.info(f"Load model matrices for {str(previous_date)}")
    strats_name = f"strategy={cfg.expr_name}"
    path = os.path.join(cfg.output_prefix, 
                        cfg.system, 
                        strats_name, 
                        'model_selection', 
                        '{sym}', 
                        str(previous_date.year), 
                        str(previous_date.month).zfill(2), 
                        str(previous_date.day).zfill(2), 
                        )
    
    _logger.info(f'Constructed prediction path is: {path}.')
    price_path = os.path.join(cfg.output_prefix, cfg.system, strats_name, 'risk', 
                              str(previous_date.year), str(previous_date.month).zfill(2))
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
        coef_dict[sym] = pd.read_csv(os.path.join(pred_path_, f"{pre_date_prefix}.coef.csv"), index_col=0).squeeze().values
        z_dict[sym] = pd.read_csv(os.path.join(pred_path_, f"{pre_date_prefix}.z0.csv"), index_col=0).squeeze().values 
        M_dict[sym] = pd.read_csv(os.path.join(pred_path_, f"{pre_date_prefix}.M.csv"), index_col=0).values 
        N_dict[sym] = pd.read_csv(os.path.join(pred_path_, f"{pre_date_prefix}.N.csv"), index_col=0).values 
        U_dict[sym] = pd.read_csv(os.path.join(pred_path_, f"{pre_date_prefix}.U.csv"), index_col=0).values 
        Q_dict[sym] = pd.read_csv(os.path.join(pred_path_, f"{pre_date_prefix}.Q.csv"), index_col=0).values 
        scale_dict[sym] = pd.read_csv(os.path.join(pred_path_, f"{pre_date_prefix}.scale.csv"), index_col=0)
    
    return ModelMatrices(coef_dict=coef_dict, 
                         z_dict=z_dict, 
                         M_dict=M_dict, 
                         N_dict=N_dict, 
                         U_dict=U_dict, 
                         Q_dict=Q_dict, 
                         scale_dict=scale_dict)

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


def load_prod_price_date(cfg: Config, previous_date: datetime.date, 
                         cur_date: datetime.date, 
                         market: str): 
    """
    Load previous trading day electron price data 
    and predictions for each tod 
    """
    prod_data_path = "/cc/dat/strat/sft/elektron/"
    dateformat = previous_date.strftime("%Y%m%d")
    cur_dateformat = cur_date.strftime("%Y%m%d")
    
    if market == 'US': 
        tod_list = TOD_ 
    elif market == "FX": 
        tod_list = FX_TOD_ 
    elif market == "EU": 
        tod_list = EU_TOD_ 

    electronPrice = {} 
    for symbol in cfg.universe: 
        electronPrice[symbol] = {} 
        for tod in tod_list: 
            if tod >= '18:00:00' and tod <= '23:59:00': 
                pathDateFormat = cur_dateformat
            else: 
                pathDateFormat = dateformat 
            
            try: 
                df = pd.read_csv(os.path.join(prod_data_path, f"Futures_{pathDateFormat}_{tod[:2] + tod[3:5]}.csv.gz"))
                eprice = df.loc[df.index.str.contains(symbol)].iloc[0]['Tomorrow']
            except: 
                eprice = np.nan 
            
            tod_dt = datetime.strptime(tod, "%H:%M:%S").time() 
            date_tod_dt = datetime.combine(previous_date, tod_dt)
            electronPrice[symbol][date_tod_dt] = eprice
    
    electronPrice_df = pd.DataFrame.from_dict(electronPrice, orient="columns")
    return electronPrice_df

def load_prod_prediction(cfg: Config, previous_date: datetime.date): 
    pass 


def load_research_price_data(cfg: Config, previous_date: datetime.date, market: str): 
    if market == 'US': 
        research_path = f"rth_cmefut_prod/{sym}"
    elif market == 'FX': 
        research_path = f'rth_fx_1min_200_starttime/{sym}'
    elif market == 'EU': 
        research_path = f"rth_eufut_prod/{sym}"

    dateformat = previous_date.strftime("%Y%m%d")
    price = {} 

    for symbol in cfg.universe: 
        path = os.path.join(research_path.format(sym=symbol), dateformat, f"{symbol}_{dateformat}_adjusted_price.csv")
        df = pd.read_csv(path, index_col=0)
        df.index = pd.to_datetime(df.index)
        tod_price = extract_tod_price(df.loc[str(previous_date)], market)
        price[symbol] = tod_price 
    return pd.DataFrame.from_dict(price, orient="columns")


