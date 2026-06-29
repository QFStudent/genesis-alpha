import os
import sys
import schedule
import time
import pandas as pd
import numpy as np
from datetime import datetime
from lib.sft_ym.cds.backtesting.config import load_config
from lib.sft_ym.dp.dp import calc_previous_bday, calc_next_bday
import logging
from pathlib import Path


_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s:%(name)s:%(message)s")


US_TOD = ['09:40:00', '11:40:00', '13:40:00', '15:40:00', '16:00:00']
FX_TOD = ['03:00:00', '06:00:00', '09:00:00', '12:00:00', '15:00:00', '18:00:00', '21:00:00']
EU_TOD = ['03:10:00', '05:10:00', '07:10:00', '09:10:00', '11:10:00', '11:30:00']  # ETC


class Trader:
    def __init__(self, config: str, date: str, market: str) -> None:
        self.cfg = load_config(config)
        self.curDate = date
        self.nextDate = (
            calc_next_bday(self.curDate, 'MP')
            if market == 'FX'
            else calc_next_bday(self.curDate, 'ES')
        )
        self.market = market

        if self.market == 'US':
            self.tod_idx = US_TOD
        elif self.market == 'EU':
            self.tod_idx = EU_TOD
        elif self.market == 'FX':
            self.tod_idx = FX_TOD

        self.curDateFormat = self.curDate.strftime('%Y%m%d')
        self.nextDateFormat = self.nextDate.strftime('%Y%m%d')
        self.preDate = calc_previous_bday(self.curDate, 'ES')
        self.preDateFormat = self.preDate.strftime('%Y%m%d')
        self.model_selection_path = os.path.join(
            self.cfg.output_prefix,
            self.cfg.system,
            f"strategy={self.cfg.expr_name}",
            "model_selection",
            '{sym}',
            str(self.preDate)[0:4],
            str(self.preDate)[5:7],
            str(self.preDate)[8:10],
        )

        self.live_data_path = '/cc/dat/strat/sft/elektron'
        # self.live_data_path = '/ccny5/dat/strat/sft/yu/dp/sim_rth/'
        # self.pred_path = '/ccny5/dat/strat/sft/yu/dp/prod_predictions/'
        self.pred_path = f'/cc/dat/strat/sft/yu_alpha/{self.market}'
        self.universe = self.cfg.universe
        self.system = self.cfg.system
        self.initializeModelMatrices()
        self.setup_log()

    def setup_log(self):
        self.log_path = (
            f'/ccny5/dat/strat/sft/yu/cds/log/directional/{self.market}/'
            f'{self.curDate.strftime("%Y%m%d")}'
        )
        Path(self.log_path).mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(os.path.join(self.log_path, "daily_trade.log"))
        file_handler.setFormatter(formatter)
        _logger.addHandler(file_handler)

    def load_last_price(self):
        dfs = []
        for u in self.cfg.universe:
            date_format = datetime.strptime(str(self.cfg.preDate), "%Y-%m-%d").strftime("%Y%m%d")
            path_ = os.path.join(
                self.cfg.input_prefix,
                u,
                date_format,
                f"{u}_{date_format}_adjusted_price.csv",
            )
            df = pd.read_csv(path_, index_col=0)
            df.index = pd.to_datetime(df.index)
            dfs.append(df['WEIGHTED'])
        minute_prices = pd.concat(dfs, axis=1)
        minute_prices.columns = self.cfg.universe
        last_price = minute_prices.iloc[-1]

        if last_price.isna().any():
            raise ValueError('Last price from yesterday is Nan, something is wrong.')
        return last_price

    @staticmethod
    def _predict(z, C, invM, N, Q, U, y):
        zii = N.dot(Q.dot(z))
        zii += U.dot(y.reshape(-1, 1)).ravel()
        z = invM.dot(zii)
        yhat = C.dot(z)
        return yhat, z

    def initializeModelMatrices(self):
        self.trade_price = {}
        self.trade_ret = []
        self.scaled_ret = []
        self.coefDict = {}
        self.zDict = {}
        self.invMDict = {}
        self.NDict = {}
        self.UDict = {}
        self.QDict = {}
        self.scaleDict = {}

        for sym in self.universe:
            modelPath = self.model_selection_path.format(sym=sym)

            coef_df = pd.read_csv(
                os.path.join(modelPath, f"{self.preDateFormat}.coef.csv"), index_col=0
            )
            self.coefDict[sym] = coef_df.values

            z0_df = pd.read_csv(
                os.path.join(modelPath, f"{self.preDateFormat}.z0.csv"), index_col=0
            )
            self.zDict[sym] = z0_df.values.ravel()

            M = pd.read_csv(os.path.join(modelPath, f"{self.preDateFormat}.M.csv"), index_col=0)
            self.invMDict[sym] = np.linalg.inv(M.values)

            N = pd.read_csv(os.path.join(modelPath, f"{self.preDateFormat}.N.csv"), index_col=0)
            self.NDict[sym] = N.values

            U = pd.read_csv(os.path.join(modelPath, f"{self.preDateFormat}.U.csv"), index_col=0)
            self.UDict[sym] = U.values

            Q = pd.read_csv(os.path.join(modelPath, f"{self.preDateFormat}.Q.csv"), index_col=0)
            self.QDict[sym] = Q.values

            scale = pd.read_csv(
                os.path.join(modelPath, f"{self.preDateFormat}.scale.csv"), index_col=0
            )
            if not (scale.columns == ['std', 'lb', 'ub']).all():
                raise ValueError('Scale format is incorrect.')
            self.scaleDict[sym] = scale

        self.last_price = self.load_last_price()
        for sym in self.cfg.universe:
            self.trade_price[sym] = [self.last_price.loc[sym]]

    def make_prediction(self):
        predictions = {}
        for i, fut in enumerate(self.universe):
            pred, self.zDict[fut] = Trader._predict(
                self.zDict[fut],
                self.coefDict[fut],
                self.invMDict[fut],
                self.NDict[fut],
                self.QDict[fut],
                self.UDict[fut],
                self.scaled_ret[-1][i],
            )
            predictions[fut] = pred * self.scaleDict[fut]["std"].squeeze()
        return predictions

    def getLast(self, curtime: str):
        file_path = os.path.join(
            self.live_data_path,
            f"Futures_{self.pathDateFormat}_{curtime[:2] + curtime[3:5]}.csv.gz",
        )

        # read data in a while loop with waiting for 5 seconds
        c = 0
        while c < 10:
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, index_col=1)
                break
            else:
                time.sleep(1)
                c = c + 1
                continue

        # generate row of data for each future
        data_row = {}
        for sym in self.universe:
            if 'elektron' in self.live_data_path:
                data_row[sym + '_last'] = df.loc[df.index.str.contains(sym)].iloc[0]['Tommorrow']
            else:
                data_row[sym + '_last'] = df.loc[sym, 'LAST']
        return data_row

    def processBarData(self, data_row: dict):
        self.trade_ret.append([])
        self.scaled_ret.append([])

        for sym in self.universe:
            self.trade_price[sym].append(data_row[sym + "_last"])
            price = self.trade_price[sym]
            # forward fill nans
            if price[-1] != price[-1]:
                _logger.info(f'{sym} price is Nan, ffill.')
                price[-1] = price[-2]

            # calculate returns
            ret = np.log(price[-1]) - np.log(price[-2])
            self.trade_ret[-1].append(ret)

            # preprocess
            if self.system == 'SISO':
                lb = self.scaleDict[sym]["lb"].squeeze()
                ub = self.scaleDict[sym]["ub"].squeeze()
                std = self.scaleDict[sym]["std"].squeeze()
                ret = np.clip(ret, lb, ub)
                ret = ret / std
            self.scaled_ret[-1].append(ret)

        self.alphas = self.make_prediction()

    def fill_missing(self):
        data_row = {}
        for sym in self.universe:
            data_row[sym + '_last'] = self.trade_price[sym][-1]
        return data_row

    def job(self, curtimeformat: str = None):
        if curtimeformat is None:
            curtime = datetime.now().time()
            curtimeformat = curtime.strftime('%H:%M:%S')

        if curtimeformat >= '18:00:00' and curtimeformat <= '23:59:00':
            self.pathDateFormat = self.nextDateFormat
        else:
            self.pathDateFormat = self.curDateFormat

        outPath = os.path.join(
            self.pred_path,
            f'Futures_{self.pathDateFormat}_{curtimeformat[:2] + curtimeformat[3:5]}.csv.gz',
        )

        try:
            data_row = self.getLast(curtimeformat)
        except Exception:
            _logger.info(f'Missing data at {curtimeformat} and forward fill the price.')
            data_row = self.fill_missing()

        self.processBarData(data_row)
        pd.DataFrame.from_dict(self.alphas, orient='index', columns=['prediction']).to_csv(outPath)
        _logger.info(f'Scaled return at {curtimeformat} is: {self.scaled_ret}')

    def sim_onData(self):
        for tformat in self.tod_idx:
            print(tformat)
            # self.wait_job(tformat)
            self.job(tformat)

    def wait_job(self, timeformat: str):
        """
        Wait before live data being loaded
        """
        while True:
            try:
                self.job(timeformat)
                print("job processed without wait.")
                break
            except Exception:
                print("wait..")
                time.sleep(1)
                continue

    def serial_process(self):
        starttime = time.monotonic()
        while True:
            curtime = datetime.now().time()
            curtimeformat = curtime.strftime('%H:%M:%S')

            if curtimeformat in self.tod_idx:
                # re-run missed hour
                if self.tod_idx.index(curtimeformat) > 0 and len(self.trade_ret) == 0:
                    # if it's second or higher tod but trade ret hasn't been loaded yet, re-run the first few tods.
                    for i in range(self.tod_idx.index(curtimeformat) + 1):
                        self.wait_job(self.tod_idx[i])
                else:
                    self.wait_job(curtimeformat)
            else:
                print(curtimeformat)
                time.sleep(1 - (time.monotonic() - starttime) % 1)

    def simple_serial_process(self):
        while True:
            curtime = datetime.now().time()
            curtimeformat = curtime.strftime('%H:%M:%S')
            starttime = time.monotonic()

            if curtimeformat in self.tod_idx:
                # re-run missed hour
                if self.tod_idx.index(curtimeformat) > 0 and len(self.trade_ret) == 0:
                    # if it's second or higher tod but trade ret hasn't been loaded yet, re-run the first few tods.
                    for i in range(self.tod_idx.index(curtimeformat)):
                        self.job(self.tod_idx[i])
                self.job(curtimeformat)
                time.sleep(1 - (time.monotonic() - starttime) % 1)
            else:
                print(curtimeformat)
                time.sleep(1 - (time.monotonic() - starttime) % 1)

    def runData(self):
        for tformat in self.tod_idx:
            print(tformat)
            self.wait_job(tformat)
            self.job(tformat)


if __name__ == "__main__":
    market = sys.argv[1]
    # prod config
    config = (
        f"/home/yu.mu/yuresearch/lib/sft_ym/cds/prod/config/{market.lower()}/"
        f"direct_prod_dev1_{market.lower()}.toml"
    )
    curDate = datetime.now().date()

    # dates = ['20241104', '20241105', '20241106', '20241107', '20241108',
    #          '20241111', '20241112', '20241113', '20241114', '20241115',
    #          '20241118', '20241119', '20241120', '20241121', '20241122']
    # for date in dates:
    #     curDate = datetime.strptime(date, '%Y%m%d').date()
    #     inst = Trader(config, curDate, market)
    #     inst.sim_onData()
    #     inst.serial_process()
    inst = Trader(config, curDate, market)
    inst.simple_serial_process()

