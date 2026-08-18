'''
Created on Oct 7, 2019

@author: john.casale
'''
import pandas as pd
import numpy as np
import csv, yaml, datetime, os
from BasisArb.MainStrat import Book
import matplotlib.pyplot as plt

class Backtester(object):

    def __init__(self, sampleInterval, paramFilePath, region='US', **kwargs):
        fullAudit = kwargs.get('fullAudit',False)
        impactLimitOn = kwargs.get('impactLimit',True)
        execThresh = kwargs.get('execThresh',1.0)
        randSeed = kwargs.get('randSeed',0)
        closeOutExpiry = kwargs.get('closeOutExpiry', 0.0)
        maxBfList = kwargs.get('maxBfList',[])
        #
        self.extAlphaMult = kwargs.get('extAlphaMult', 0)
        self.extFileName = kwargs.get('extFileName', '')
        #self.extFileName = '_new'+self.extFileName if self.extAlphaMult > 0 else ''
        self.barDataFilePath = r'/datadrive/Futures/BasisArb/%s_futures_%ssec_bars%s.csv' % (region, str(sampleInterval), self.extFileName)
        self.startDate = kwargs.get('startDate','2016-01-01')
        self.endDate = kwargs.get('endDate','2099-12-31')
        self.region = region
        # load all book and strategy parameters
        with open(paramFilePath) as f:
            self.paramStr = f.read()
        with open(paramFilePath) as f:
            params = yaml.load_all(f,Loader=yaml.FullLoader)
            self.book = Book.Book(self.startDate,region,next(params),params,sampleInterval,fullAudit,impactLimitOn,execThresh,randSeed,closeOutExpiry, maxBfList, self.extAlphaMult)
        #
        print('Using external alpha multiplier:', self.extAlphaMult)
        print('Finished initializing all objects, starting backtest')

    def loadPositionFile(self, dateToLoad):
        bbgToGeneric = {'ES':'ES', 'NQ':'NQ', 'DM':'YM', 'FA':'DM', 'RTY':'TFS'}
        dateStr = str(dateToLoad).replace('-','')
        prodFile = r'/mnt/tradingtech/data/archive/logs/sftprod/%s/%s/%s.eod.posn.csv' % (dateStr[0:4], dateStr[0:8], dateStr[0:8])
        df = pd.read_csv(prodFile)
        #
        closePx = self.book.closePxUnadj.loc[str(dateToLoad)].to_dict()
        adjFactor = self.book.adjFactor.loc[str(dateToLoad)].to_dict()
        #
        df['fut'] = df['ticker'].apply(lambda x: bbgToGeneric[x[:-2]])
        df['close'] = df['fut'].apply(lambda x: closePx[x])
        df['adjFactor'] = df['fut'].apply(lambda x: adjFactor[x])
        df['key'] = df['book'] + '_' + df['fut']
        df = df.set_index('key')
        #
        for name, fut in self.book.bookFutDict.items():
            curDict = df.loc['main_' + name].to_dict()
            fut.currentPosition = curDict['posn']
            fut.todayClosePx = curDict['close']
            fut.todayAdjFactor = fut.yestAdjFactor = curDict['adjFactor']
        for strat in self.book.stratList:
            for fut in strat.futList:
                curDict = df.loc[strat.stratName + '_' + fut.name].to_dict()
                fut.currentPosition = curDict['posn']
                fut.todayClosePx = curDict['close']
                fut.todayAdjFactor = fut.yestAdjFactor = curDict['adjFactor']
        #
        print('Finished loading all positions:')
        print([k+':'+str(v.currentPosition) for k,v in self.book.bookFutDict.items()])

    def runBacktest(self):
        t1 = datetime.datetime.now()
        prevDate = None
        with open(self.barDataFilePath) as f:
            #data = csv.DictReader(f)
            headers = f.readline().rstrip().split(',')
            smallHeaders = headers[2:]
            for row in f.readlines():
                curDate = row[0:10] #row['Date']
                if curDate < self.startDate:
                    continue
                if curDate > self.endDate:
                    break
                # if we are on a new day, process data
                if curDate != prevDate:
                    #print(curDate, prevDate)
                    # process trading data from yesterday
                    if not prevDate is None:
                        self.book.progressForwardOneDay(prevDate)
                        print('Finished day:',prevDate)
                    # initialize today's data
                    self.book.initializeTodayValues(curDate)
                #
                #self.book.processBarData(dict(zip(headers,row.rstrip().split(','))))
                self.book.processBarData(row[0:10], row[11:19], dict(zip(smallHeaders, [np.NaN if not x else float(x) for x in row[20:].rstrip().split(',')])))
                self.book.checkAndProcessBookTrades()
                prevDate = curDate
            #
            self.book.progressForwardOneDay(prevDate)
            print('Finished day:',prevDate)
        #
        self.book.getOutput()
        t2 = datetime.datetime.now()
        print('Total runtime: %s seconds' % ('{:.1f}'.format((t2-t1).total_seconds())))

    def runBacktestMultipleParam(self, paramPath):
        paramDates = [x[0:8] for x in os.listdir(paramPath)]
        t1 = datetime.datetime.now()
        prevDate = None
        with open(self.barDataFilePath) as f:
            #data = csv.DictReader(f)
            headers = f.readline().rstrip().split(',')
            smallHeaders = headers[2:]
            for row in f.readlines():
                curDate = row[0:10] #row['Date']
                if curDate < self.startDate:
                    continue
                if curDate > self.endDate:
                    break
                # if we are on a new day, process data
                if curDate != prevDate:
                    # process trading data from yesterday
                    if not prevDate is None:
                        self.book.progressForwardOneDay(prevDate)
                        print('Finished day:',prevDate)
                    # check if we changed param file
                    if curDate.replace('-','') in paramDates:
                        self.book.updateParams(paramPath + curDate.replace('-','') + '.Jus.yaml')
                        print('Finished updating params on:',curDate)
                    # initialize today's data
                    self.book.initializeTodayValues(curDate)
                #
                #self.book.processBarData(dict(zip(headers,row.rstrip().split(','))))
                self.book.processBarData(row[0:10], row[11:19], dict(zip(smallHeaders, [float(x) for x in row[20:].rstrip().split(',')])))
                self.book.checkAndProcessBookTrades()
                prevDate = curDate
            #
            self.book.progressForwardOneDay(prevDate)
            print('Finished day:',prevDate)
        #
        self.book.getOutput()
        t2 = datetime.datetime.now()
        print('Total runtime: %s seconds' % ('{:.1f}'.format((t2-t1).total_seconds())))

    def getBookLevelTradeData(self):
        return self.getBookLevelPnlAll(), self.getBookLevelPnlIntra(), self.getBookLevelPnlExtra(), self.getBookLevelGmv(), self.getBookLevelT0()

    def getBookLevelPnlAll(self):
        return self.book.bookTradeList.pnlAll

    def getBookLevelPnlIntra(self):
        return self.book.bookTradeList.pnlIntra

    def getBookLevelPnlExtra(self):
        return self.book.bookTradeList.pnlExtra

    def getBookLevelGmv(self):
        return self.book.bookTradeList.gmv

    def getBookLevelT0(self):
        return self.book.bookTradeList.t0

    def getStats(self):
        return self.book.getOutput()

    def getStratLevelTradeLimits(self,stratName):
        for strat in self.book.stratList:
            if strat.stratName == stratName:
                return strat.stratTradeList.tradeLimitDict
        return {}

    def getStratLevelTradeData(self,stratName):
        for strat in self.book.stratList:
            if strat.stratName == stratName:
                return strat.stratTradeList.pnlAll, strat.stratTradeList.pnlIntra, strat.stratTradeList.pnlExtra, strat.stratTradeList.gmv, strat.stratTradeList.t0
        return [], [], [], [] ,[]

    def plotStrategyPnl(self,stratName):
        plt.figure()
        for strat in self.book.stratList:
            if strat.stratName == stratName:
                strat.stratTradeList.pnlAll['Total'].cumsum().plot(label=strat.stratName+'_Total',rot=45)
                strat.stratTradeList.pnlIntra['Total'].cumsum().plot(label=strat.stratName+'_Intra',rot=45)
                strat.stratTradeList.pnlExtra['Total'].cumsum().plot(label=strat.stratName+'_Extra',rot=45)
        plt.tight_layout()
        plt.legend()
        plt.show()

    def printStats(self):
        return self.book.getOutput()

    def getParams(self):
        return self.paramStr

    def outputOrdersFile(self):
        futToBbg = {'ES':'ES','NQ':'NQ','YM':'DM','DM':'FA','TFS':'RTY'}
        ord = pd.DataFrame(self.book.bookTradeList.rawTradeList, columns=['date','time','fut','size','avg_fill_px'])
        ord['side'] = ord['size'].apply(lambda x: 'BUY' if x > 0 else 'SELL')
        ord['fill_qty'] = ord['size'].abs()
        ord['symbol'] = ord['fut'].apply(lambda x: futToBbg[x] + 'M0')
        ord.to_csv(r'/datadrive/BasisArb/Logging/orders.csv')
