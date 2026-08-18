'''
Created on Oct 7, 2019

module for individual Book class

@author: john.casale
'''
import numpy as np
import pandas as pd
import datetime, math, yaml
from BasisArb.MainStrat import Strategy, Future, TradeList

sign = lambda x: -1 if x < 0 else (1 if x > 0 else (0 if x == 0 else np.NaN))

class Book:

    def __init__(self, startDate, region, params, strats, sampleInterval, fullAudit, impactLimitOn, execThresh=1.0, randSeed=0, closeOutExpiry=0.0, maxBfList=[], extAlphaMult=0):
        self.sampleInterval = sampleInterval
        self.fullAudit = fullAudit
        self.impactLimitOn = impactLimitOn
        self.startDate = startDate
        self.region = region
        self.randSeed = randSeed
        self.closeOutExpiry = closeOutExpiry
        #
        self.setBookParams(params, maxBfList)
        #
        self.volPosMult = 1.0
        self.vixDict = pd.read_csv(r'/mnt/nfs1/user/sft/sftdata/CurrentFiles/vix.csv').set_index('Date')['Close']
        self.vixDict.index = pd.to_datetime(self.vixDict.index)
        self.vixDict = self.vixDict.reindex(pd.date_range(self.vixDict.index[0], self.vixDict.index[-1]))
        self.vixDict.index = [x.strftime('%Y-%m-%d') for x in self.vixDict.index]
        self.vixDict = self.vixDict.shift(1).ffill().bfill().to_dict()
        #
        self.sizePosMult = 1.0
        self.sizeMaxMult = 1.0
        retName = 'SIZE'
        self.sizeDict = pd.read_csv(r'/datadrive/Equity/Barra/USMED/FactorReturn/BarraFactorReturn.csv').set_index('DataDate')[retName]
        #pd.read_csv(r'/home/fut_strat/media/disk/_cygdrive_C_JOHNCA1_temp/move.csv').set_index('Date')['MOVE']
        self.sizeDict = self.sizeDict[self.sizeDict.index >= '2016-01-01']
        #self.sizeDict = self.sizeDict / self.sizeDict.shift(1) - 1.0
        self.sizeDict = self.sizeDict.rolling(window=10,min_periods=2).std().shift(1).ffill().bfill().to_dict()
        #self.sizeDict = self.sizeDict.shift(1).ffill().bfill().to_dict()
        #
        # close prices
        try:
            self.closePxAdj = pd.read_csv(params['closePriceAdjFile']).set_index('Date').ffill()
        except:
            self.closePxAdj = pd.read_csv(r'/mnt/nfs1/user/sft/sftdata/CurrentFiles/%s_futures_close_adj.csv' % (region)).set_index('Date').ffill()
        try:
            self.closePxUnadj = pd.read_csv(params['closePriceUnadjFile']).set_index('Date').ffill()
        except:
            self.closePxUnadj = pd.read_csv(r'/mnt/nfs1/user/sft/sftdata/CurrentFiles/%s_futures_close_unadj.csv' % (region)).set_index('Date').ffill()
        self.adjFactor = self.closePxAdj / self.closePxUnadj
        # events
        self.eventMarginMult = 1.0
        self.readAndProcessEventList(params['eventMaxMargin'])
        self.fx = self.readFxToDataframe()
        # MANUAL BT PARAMS
        self.reduceRiskOnlyEod = False
        self.reduceRiskOnly = False
        self.reduceRiskOnlyStartTime = '16:00:00'
        # initialize risk positions
        self.curBookBeta = 0
        self.curBookGmv = 0
        self.curBookRf = 0
        self.curBookBf = [0] * len(self.barraFactorNameList)
        self.globalBookMargin = 1.0
        # load all strategy parameters
        self.stratList = []
        for strat in strats:
            self.stratList.append(Strategy.Strategy(self,strat,fullAudit,execThresh,randSeed,closeOutExpiry))
        self.numStrats = len(self.stratList)
        # save book level future list
        self.bookFutDict = {}
        for fut in self.allFutNames:
            self.bookFutDict[fut] = Future.Future(None,fut,{},{},self.startDate,self.fx)
        self.numFuts = len(self.allFutNames)
        # initiate book level trade list
        self.bookTradeList = TradeList.TradeList(self.allFutNames,self.closePxAdj,self.closePxUnadj)
        self.orderAggFeatures = ['BuyQty','SellQty','NetQty','IntQty','BidSize','AskSize','ExecQty','BuyPx','SellPx','IntPx','ExtPx']
        self.futRowDict = dict(zip(self.allFutNames,range(self.numFuts)))
        self.useOrderAggregator = False
        self.execThresh = execThresh
        np.random.seed(randSeed)
        #
        self.computeBetaMatrix(self.allFutNames,self.closePxAdj)
        self.computeBetaPnl = False
        #
        self.extAlphaMult = extAlphaMult

    def setBookParams(self, params, maxBfList):
        self.maxBookBeta = params['maxBookBeta']
        self.maxBookRiskFactor = params['maxBookRiskFactor']
        self.maxBookBarraFactor = params.get('maxBookBarraFactor', [])#
        self.maxBookGmv = params['maxBookGmv']
        self.timeBetweenTrades = params['timeBetweenTrades']
        self.runStratsIndependent = params['runStratsIndependent']
        self.i1EwmSecs = params['i1EwmSecs']
        self.maxBookPos = params['maxBookPos']
        self.allFutNames = list(self.maxBookPos.keys())
        self.riskFactors = params['riskFactor']
        self.capI2 = params['capI2']
        #
        if 'PCA' in params.get('barraFactor1Name',''):
            self.barraFactorNameList = [params['barraFactor1Name'], params['barraFactor2Name']] if 'barraFactor2Name' in params.keys() else [params['barraFactor1Name']]
            self.maxBookBarraFactor = [params['maxBookBarraFactor1'], params['maxBookBarraFactor2']] if 'barraFactor2Name' in params.keys() else [params['maxBookBarraFactor1']]
            self.barraFactorNameToPos = dict(zip(self.barraFactorNameList, range(len(self.barraFactorNameList))))
            static = pd.read_csv(r'/datadrive/Futures/BasisArb/AP_PCA_static.csv').set_index('Future')
            static.columns = ['PCA1','PCA2','PCA3','PCA4'] if '_' not in params['barraFactor1Name'] else ['PCA_1','PCA_2','PCA_3','PCA_4']
            self.barraFactorDict, self.barraFactors = {}, {}
            for k in self.barraFactorNameList:
                if 'PCA' in k:
                    self.barraFactors[k] = static[k].to_dict() #{'JNM':0.0, 'JTI':0.0, 'YAP':0.0, 'HCEI':0.0, 'HSI':0.0, 'KS':0.0}#
                elif 'DELTA' in k:
                    self.barraFactors[k] = {'JNM':1.0, 'JTI':1.0, 'YAP':1.0, 'HCEI':1.0, 'HSI':1.0, 'KS':1.0}
        else:
            self.barraFactorNameList = params.get('barraFactorName', [])#[]#['GROWTH','MOMENTUM']#
            if 'barraFactor1Name' in params.keys():
                self.maxBookBarraFactor = [params['maxBookBarraFactor1'], params['maxBookBarraFactor2']]
                self.barraFactorNameList = [params['barraFactor1Name'], params['barraFactor2Name']]
            self.barraFactorNameToPos = dict(zip(self.barraFactorNameList, range(len(self.barraFactorNameList))))
            self.barraFactorDict = {}
            self.barraFactors = {}
            barraModel = 'USMED' if self.region == 'US' else 'GEMLT'
            barraFactorsInput = pd.read_csv(r'/datadrive/Equity/Barra/%s/Future_exposures.csv' % barraModel)
            temp = barraFactorsInput[barraFactorsInput['Factor'] == 'SIZE'].drop_duplicates(['DataDate', 'Future']).pivot('DataDate', 'Future', 'Exposure')
            temp = temp.shift(1).ffill().bfill()
            self.riskFactorsDateList = temp.T.to_dict()
            for factor in self.barraFactorNameList:
                temp = barraFactorsInput[barraFactorsInput['Factor'] == factor].drop_duplicates(['DataDate', 'Future']).pivot('DataDate', 'Future', 'Exposure')
                temp = temp.shift(1).ffill().bfill()
                self.barraFactorDict[factor] = temp.T.to_dict()
                self.barraFactors[factor] = self.barraFactorDict[factor][min(self.barraFactorDict[factor].keys())]
        #
        print(self.barraFactorNameList, self.maxBookBarraFactor)
        self.calcI4 = params.get('calcI4', False)
        self.ewmI4 = str(params.get('ewmI4', 300))
        self.ewmI4Alpha = 1.0 - np.exp(-2.0 * self.sampleInterval / float(self.ewmI4))
        self.futVolI4 = [0.0] * len(self.allFutNames)
        #
        self.startTime = params['startTime']
        self.endTime = params['endTime']
        self.timezone = params['timezone']
        #
        self.vixMinMult = params.get('vixMinMult',1.0)
        self.sizeMinMult = params.get('sizeMinMult',1.0)
        print('Vix min mult =',self.vixMinMult,'Size min mult =',self.sizeMinMult)
        # impact
        self.impactEwms = params['impact'].keys()
        self.impactLimits = params['impact'].values()
        self.impactEwmAlphas = [1.0 - np.exp(-2.0 * self.sampleInterval / impEwmSecs) for impEwmSecs in self.impactEwms]
        self.mktVolEwmDict = {fut+'_mkt_imp_'+str(x):0.01 for x in self.impactEwms for fut in self.allFutNames}
        self.bookVolEwmDict = {fut+'_mkt_imp_'+str(x):0 for x in self.impactEwms for fut in self.allFutNames}
        # dynamic params
        self.startMultTime = params['startMultTime'] #: "09:30:00"
        self.endMultTime = params['endMultTime'] #: "16:00:00"
        timeList = [x.strftime('%H:%M:%S') for x in pd.date_range(pd.to_datetime('2017-01-01 '+self.startMultTime), pd.to_datetime('2017-01-01 '+self.endMultTime), freq='1S')]
        timeListMorn = [x.strftime('%H:%M:%S') for x in pd.date_range(pd.to_datetime('2017-01-01 09:30:00'), pd.to_datetime('2017-01-01 '+self.startMultTime), freq='1S')]
        # triangle curve
        #self.timeBetaMultDict = dict(zip(timeList,np.append(np.linspace(params['startBetaMult'],1.5,np.sum([x <= '12:00:00' for x in timeList])), np.linspace(1.5, params['endBetaMult'],  <<<LINE TRUNCATED IN SOURCE PHOTO>>>
        #self.timeLambdaMultDict = dict(zip(timeList,np.append(np.linspace(params['startLambdaMult'],0.5,np.sum([x <= '12:00:00' for x in timeList])), np.linspace(0.5, params['endLambdaM  <<<LINE TRUNCATED IN SOURCE PHOTO>>>
        # step curve
        #self.timeBetaMultDict = dict(zip(timeList,np.append([1.]*1800 + [1.125]*1800 + [1.25]*1800 + [1.375]*1800 + [1.5]*1800, np.linspace(1.5, params['endBetaMult'], np.sum([x >= '12:  <<<LINE TRUNCATED IN SOURCE PHOTO>>>
        #self.timeLambdaMultDict = dict(zip(timeList,np.append([1.]*1800 + [0.875]*1800 + [0.75]*1800 + [0.625]*1800 + [0.5]*1800, np.linspace(0.5, params['endLambdaMult'], np.sum([x >=  <<<LINE TRUNCATED IN SOURCE PHOTO>>>
        # eod linear curve
        self.timeBetaMultDict = dict(zip(timeList,np.linspace(params['startBetaMult'],params['endBetaMult'],len(timeList))))
        if self.region == 'US' and params.get('morningBetaMult', 1.0) > 1.0:
            self.timeBetaMultDict.update(dict(zip(timeListMorn,np.linspace(params['morningBetaMult'],params['startBetaMult'],len(timeListMorn)))))
        self.timeLambdaMultDict = dict(zip(timeList,np.linspace(params['startLambdaMult'],params['endLambdaMult'],len(timeList))))
        self.timeBarra1MultDict = dict(zip(timeList,np.linspace(params['startBetaMult'],params['endBetaMult'],len(timeList))))
        self.timeBarra2MultDict = dict(zip(timeList,np.linspace(params['startBetaMult'],params['endBetaMult'],len(timeList))))
        self.useLambda1Mult = ('lambda1' in params['multNames'])
        self.useLambda2Mult = ('lambda2' in params['multNames'])
        self.useMaxBetaMult = ('beta' in params['multNames'])
        self.useMaxGmvMult = ('gmv' in params['multNames'])
        #
        self.t0MarginInc = params['t0MarginInc']
        self.t0CurveDict = pd.read_csv(params['t0CurveFile']).set_index('Time')['Percent'].to_dict()
        if self.runStratsIndependent:
            self.targetT0 = np.inf
        else:
            self.targetT0 = params['targetT0']

    def updateParams(self, paramFilePath):
        with open(paramFilePath) as f:
            newParams = yaml.load_all(f,Loader=yaml.FullLoader)
            params = next(newParams)
            strats = newParams
            #
            self.setBookParams(params)
            # load all strategy parameters
            curStratNameList = [x.stratName for x in self.stratList]
            for strat in strats:
                curStratName = strat['stratName']
                if curStratName in curStratNameList:
                    curStratIndex = curStratNameList.index(curStratName)
                    self.stratList[curStratIndex].updateParams(strat)
                else:
                    self.stratList.append(Strategy.Strategy(self,strat,self.fullAudit,self.execThresh,self.randSeed,self.closeOutExpiry))
            #
            self.numStrats = len(self.stratList)

    def readAndProcessEventList(self,eventMaxMargin):
        self.eventMultDict = {}
        eventList = pd.read_csv(r'/datadrive/BasisArb/BasisFiles/EventFile_US.csv')
        eventList = eventList[eventList['Date'] >= '2016-01-01']
        eventList['Datetime'] = pd.to_datetime(eventList['Date']+' '+eventList['Time'])
        for d in eventList['Datetime']:
            timeList = [x.strftime('%Y-%m-%d %H:%M:%S') for x in pd.date_range(d,d+datetime.timedelta(minutes=30),freq=str(self.sampleInterval)+'S')]
            self.eventMultDict.update(dict(zip(timeList,np.linspace(eventMaxMargin,1.0,len(timeList)))))

    def readFxToDataframe(self, fxPath=r'/datadrive/Futures/fx.csv'):
        fx = pd.read_csv(fxPath).set_index('Date')
        fxOut = pd.DataFrame()
        for col in fx.columns:
            p1 = col[0:3]
            p2 = col[3:6]
            fxOut[p1+p2] = fx[p1+p2]
            if p2+p1 not in fxOut.columns:
                fxOut[p2+p1] = 1.0 / fx[p1+p2]
            if p1+p1 not in fxOut.columns:
                fxOut[p1+p1] = 1.0
            if p2+p2 not in fxOut.columns:
                fxOut[p2+p2] = 1.0
        #
        fxOut.index = pd.to_datetime(fxOut.index)
        fxOut = fxOut.shift(1)
        fxOut = fxOut.reindex(pd.date_range(fxOut.index[0], fxOut.index[-1] + datetime.timedelta(days=90))).ffill()
        fxOut.index = [x.strftime('%Y-%m-%d') for x in fxOut.index]
        return fxOut

    def initializeTodayValues(self,curDate):
        self.todayT0 = 0
        self.bookBetaPnl = 0
        self.volPosMult = max(min(12.0 / self.vixDict[curDate], 1.0), self.vixMinMult)
        #self.sizePosMult = max(min(self.sizeDict.get(curDate, 1) / 0.002, self.sizeMaxMult), self.sizeMinMult)
        #max(min(0.0025 / self.sizeDict.get(curDate, 1), self.sizeMaxMult), self.sizeMinMult)
        #max(min(166.667*self.sizeDict.get(curDate, 0) + 1.0 - (1.0*(self.sizeDict.get(curDate, 0) > 0)*33.3*self.sizeDict.get(curDate, 0)), self.sizeMaxMult), self.sizeMinMul  <<<LINE TRUNCATED IN SOURCE PHOTO>>>
        #max(min(self.sizeDict.get(curDate, 1) / 0.002, self.sizeMaxMult), self.sizeMinMult)
        #
        if len(self.barraFactorNameList) > 0:
            try:
                if self.region == 'US':
                    self.riskFactors = self.riskFactorsDateList[curDate]
            except:
                asdf = 10
            #
            for factor in self.barraFactorNameList:
                try:
                    self.barraFactors[factor] = self.barraFactorDict[factor][curDate]#, self.barraFactorDict[factor])#{x: 0 for x in self.allFutNames})
                except:
                    asdf = 10
        #
        for strat in self.stratList:
            strat.initializeTodayValues(curDate, self.barraFactors)
        #
        setCloseOutPct = (self.closeOutExpiry > 0) and (curDate in self.stratList[0].closeOutDateList)
        for fut in self.bookFutDict.values():
            fut.initializeTodayValues(curDate, self.barraFactors, self.riskFactors[fut.name])
            if setCloseOutPct:
                fut.closeOutPositionPct = self.closeOutExpiry
            else:
                fut.closeOutPositionPct = 0.0

    def processBarData(self, curDate, curTime, row):
        # save date and time data
        self.curDate = curDate
        self.curTime = curTime
        self.curTimeSecs = 3600 * int(self.curTime[0:2]) + 60 * int(self.curTime[3:5]) + int(self.curTime[6:8])
        self.curTimeBetaMult = self.timeBetaMultDict.get(self.curTime, 1.0)
        self.curTimeLambdaMult = self.timeLambdaMultDict.get(self.curTime, 1.0)
        self.curTimeBarra1Mult = self.timeBarra1MultDict.get(self.curTime, 1.0)
        self.curTimeBarra2Mult = self.timeBarra2MultDict.get(self.curTime, 1.0)
        #
        self.reduceRiskOnly = self.reduceRiskOnlyEod and (self.curTime >= self.reduceRiskOnlyStartTime)
        #
        self.eventMarginMult = self.eventMultDict.get(self.curDate+' '+self.curTime, 1.0)
        self.maxBetaMult = self.curTimeBetaMult if self.useMaxBetaMult else 1.0
        self.maxGmvMult = max(self.curTimeBetaMult, 0.9) if self.useMaxGmvMult else 1.0
        self.maxRfMult = self.curTimeBetaMult if self.useMaxBetaMult else 1.0
        self.maxBarraMult = [self.curTimeBarra1Mult if self.useMaxBetaMult else 1.0, self.curTimeBarra2Mult if self.useMaxBetaMult else 1.0] #[1.0, self.curTimeBarra2Mult if self.useMaxBetaMult else 1.0]#[self.curTimeBarra1Mult if self.useMaxBetaMult else 1.0, 1.0]# [1.0,1.0]#
        # compute avg return
        if self.computeBetaPnl:
            self.avgRet = np.nanmean([row[futName+'_mid'] / (fut.lastValid['Mid'] * fut.yestAdjFactor / fut.todayAdjFactor) - 1.0 for futName, fut in self.bookFutDict.items()]) if (self.curTime == self.startTime) else np.nanmean([row[futName+'_mid'] / fut.lastValid['Mid'] - 1.0 for futName, fut in self.bookFutDict.items()])
            self.avgRet = self.avgRet if (self.avgRet == self.avgRet) else 0
        else:
            self.avgRet = 0
        # compute i4 if needed
        if self.calcI4:
            for fut, futPos in zip(self.allFutNames, range(len(self.allFutNames))):
                self.futVolI4[futPos] = self.ewmI4Alpha * row[fut+'_vol'] + (1 - self.ewmI4Alpha) * self.futVolI4[futPos]
        # pass bar data to all strategies
        for strat in self.stratList:
            strat.processBarData(curDate, curTime, row)
        # pass bar data to all book futures
        for futName, fut in self.bookFutDict.items():
            fut.processBarData(row[futName+'_mid'],row[futName+'_hs'],row[futName+'_index'],row[futName+'_bidSize'],row[futName+'_askSize'],self.curTime,0,0)
        # update global margin based on t0
        if not self.runStratsIndependent and (self.curTime[-4:] == '5:00' or self.curTime[-4:] == '0:00'):
            if self.todayT0 / self.targetT0 > self.t0CurveDict[self.curTime]:
                self.globalBookMargin += self.t0MarginInc#min(self.t0MarginInc, self.todayT0 / self.targetT0 - self.t0CurveDict[self.curTime])
            else:
                self.globalBookMargin = max(1.0, self.globalBookMargin - self.t0MarginInc)
        # update market volume for impact calc if we are using impact
        if self.impactLimitOn:
            for fut in self.allFutNames:
                for impSec, impAlpha in zip(self.impactEwms,self.impactEwmAlphas):
                    self.mktVolEwmDict[fut+'_mkt_imp_'+str(impSec)] = ((self.curTime >= self.startTime) * impAlpha * row[fut+'_vol']) + (1 - impAlpha) * self.mktVolEwmDict[fut+'_mkt_imp_'+str(impSec)]
                    self.bookVolEwmDict[fut+'_mkt_imp_'+str(impSec)] = (1 - impAlpha) * self.bookVolEwmDict[fut+'_mkt_imp_'+str(impSec)]
        # update portfolio stats
        self.calcBookStats()
        # compute beta pnl
        self.bookBetaPnl += (self.curBookBeta * self.avgRet)
        #print(self.curTime,self.avgRet,self.curBookBeta,self.bookBetaPnl,' '.join(str(x.lastValid['Mid']) for x in self.bookFutDict.values()))

    def checkAndProcessBookTrades(self):
        if not self.useOrderAggregator:
            mapTimeToSec = lambda x: 3600 * int(x.split(':')[0]) + 60 * int(x.split(':')[1]) + int(x.split(':')[2])
            start_index = mapTimeToSec(self.curTime) % self.numStrats
            for idx in range(self.numStrats):
                self.stratList[(idx + start_index) % self.numStrats].checkAndProcessStrategyTrades()
        #
        else:
            for strat in self.stratList:
                strat.checkAndProcessStrategyTrades()
            self.aggregateAndAllocateTrades()

    def aggregateAndAllocateTrades(self):
        # aggregate trade interests
        aggTrades = np.zeros((self.numFuts, len(self.orderAggFeatures)))
        stratOrders = np.zeros((self.numFuts, self.numStrats))
        stratCount = 0
        for strat in self.stratList:
            for fut, tradeSize, tradePrice in strat.tradeList:
                futRow = self.futRowDict[fut.name]
                if tradeSize > 0:
                    aggTrades[futRow, self.orderAggFeatures.index('BuyQty')] += tradeSize
                    aggTrades[futRow, self.orderAggFeatures.index('BuyPx')] = tradePrice
                    aggTrades[futRow, self.orderAggFeatures.index('AskSize')] = fut.curValues['AskSize']
                else:
                    aggTrades[futRow, self.orderAggFeatures.index('SellQty')] -= tradeSize
                    aggTrades[futRow, self.orderAggFeatures.index('SellPx')] = tradePrice
                    aggTrades[futRow, self.orderAggFeatures.index('BidSize')] = -1 * fut.curValues['BidSize']
                stratOrders[futRow, stratCount] = tradeSize
            #
            stratCount += 1
        # return if no orders
        if np.sum(aggTrades[:,0:2]) == 0:
            return
        # do internal calcs
        aggTrades[:,5] = np.minimum(aggTrades[:,5], [20, 12, 16, 8, 20])  # limit buy size
        aggTrades[:,4] = np.maximum(aggTrades[:,4], [-20, -12, -16, -8, -20]) # limit sell size
        #
        aggTrades[:,2] = aggTrades[:,0] - aggTrades[:,1]  # net quantity
        aggTrades[:,3] = (aggTrades[:,1] < aggTrades[:,0]) * aggTrades[:,1] - (aggTrades[:,1] > aggTrades[:,0]) * aggTrades[:,0] # internalized quantity
        aggTrades[:,6] = ((aggTrades[:,2] > 0) * aggTrades[:,[2,5]].min(axis=1) + (aggTrades[:,2] < 0) * aggTrades[:,[2,4]].max(axis=1)) * (np.random.rand(self.numFuts) <= self.execThresh) # execution qty
        aggTrades[:,9] = aggTrades[:,[7,8]].mean(axis=1) # internal px
        aggTrades[:,10] = (aggTrades[:,2] > 0) * aggTrades[:,7] + (aggTrades[:,2] < 0) * aggTrades[:,8]  # exec px
        # do assignment
        rowNum = 0
        for row in aggTrades:
            totalExec = int(row[3] + row[6])
            if totalExec == 0:
                rowNum += 1
                continue
            totalOrder = row[0] if totalExec > 0 else -1 * row[1]
            alloc = [int(x * totalExec / totalOrder) if x * totalExec > 0 else int(x) for x in stratOrders[rowNum]]
            balance = totalExec - row[3] - sum(alloc)
            balanceDir = int(sign(balance))
            curIdx = np.random.randint(0, self.numStrats)
            while balance != 0: #for i in range(balance):
                if stratOrders[rowNum, curIdx] * balance > 0:
                    alloc[curIdx] += balanceDir
                    balance -= balanceDir
                curIdx = (curIdx + 1) % self.numStrats
            #
            avgPx = (row[3] * row[9] + row[6] * row[10]) / totalExec
            for i in range(self.numStrats):
                if alloc[i] == 0:
                    continue
                strat = self.stratList[i]
                if strat.futList[(rowNum % len(strat.futList))].name == self.allFutNames[rowNum]:
                    strat.processStrategyTrade(strat.futList[rowNum], alloc[i], avgPx if alloc[i] * totalExec > 0 else row[9])
                else:
                    for fut in strat.futList:
                        if fut.name == self.allFutNames[rowNum]:
                            break
                    strat.processStrategyTrade(fut, alloc[i], avgPx if alloc[i] * totalExec > 0 else row[9])
            #
            rowNum += 1

    def checkForBookTrade(self,futName,tradeSize,curDatetime):
        if self.runStratsIndependent:
            return tradeSize, ''
        #
        fut = self.bookFutDict[futName]
        # can't buy if we violate book level limits
        if tradeSize > 0:
            if fut.currentPosition * fut.lastValid['Mid'] * fut.mult * fut.todayCcy >= self.maxBookPos[futName] * self.volPosMult * fut.sizePosMult:
                return 0, 'MaxBookPos'
            elif self.curBookBeta >= self.maxBookBeta * self.maxBetaMult:
                return 0, 'MaxBookBeta'
            elif self.curBookGmv * sign(fut.currentPosition) >= self.maxBookGmv * self.maxGmvMult:
                return 0, 'MaxBookGmv'
            elif self.curBookRf * sign(self.riskFactors[futName]) >= self.maxBookRiskFactor * self.maxRfMult:
                return 0, 'MaxBookRiskFactor'
            for k, v in self.barraFactors.items():
                factorPos = self.barraFactorNameToPos[k]
                if self.curBookBf[factorPos] * sign(v.get(futName,0)) >= self.maxBookBarraFactor[factorPos] * self.maxBarraMult[factorPos]:
                    return 0, 'MaxBookRiskFactor'
        # can't sell if we violate book level limits
        if tradeSize < 0:
            if fut.currentPosition * fut.lastValid['Mid'] * fut.mult * fut.todayCcy <= -1*self.maxBookPos[futName] * self.volPosMult * fut.sizePosMult:
                return 0, 'MaxBookPos'
            elif self.curBookBeta <= -1*self.maxBookBeta * self.maxBetaMult:
                return 0, 'MaxBookBeta'
            elif self.curBookGmv * sign(fut.currentPosition) <= -1*self.maxBookGmv * self.maxGmvMult:
                return 0, 'MaxBookGmv'
            elif self.curBookRf * sign(self.riskFactors[futName]) <= -1*self.maxBookRiskFactor * self.maxRfMult:
                return 0, 'MaxBookRiskFactor'
            for k, v in self.barraFactors.items():
                factorPos = self.barraFactorNameToPos[k]
                if self.curBookBf[factorPos] * sign(v.get(futName,0)) <= -1*self.maxBookBarraFactor[factorPos] * self.maxBarraMult[factorPos]:
                    return 0, 'MaxBookRiskFactor'
        # can't trade if any strategy traded this name recently
        if (curDatetime - fut.lastValid['TradeTime']).total_seconds() < self.timeBetweenTrades:
            return 0, 'TradeCooloff'
        #
        tradeSize, impMsg = self.checkBookImpactLimits(futName,tradeSize)
        return tradeSize, impMsg

    def checkBookImpactLimits(self,futName,tradeSize):
        if self.runStratsIndependent:
            return tradeSize, ''
        #
        for impSec, impLimit in zip(self.impactEwms,self.impactLimits):
            if self.bookVolEwmDict[futName+'_mkt_imp_'+str(impSec)] / self.mktVolEwmDict[futName+'_mkt_imp_'+str(impSec)] > impLimit:
                return 0, 'Impact'
        #
        return tradeSize, ''

    def processBookTrade(self,futName,tradeSize,tradePrice,curDatetime,isExecuted,onlyCountExecs):
        if isExecuted:
            self.bookTradeList.addTradeToday(self.curDate, self.curTime, futName, tradeSize, tradePrice)
            # update book level future position and trade time
            self.bookFutDict[futName].currentPosition += tradeSize
        # update trade cool-off / t0 / impact if we trade or counting all sent orders
        if not onlyCountExecs or isExecuted:
            self.bookFutDict[futName].lastValid['TradeTime'] = curDatetime
            # update book level info after trade
            #self.calcBookStats()
            # update today t0
            self.todayT0 += (abs(tradeSize) * tradePrice * self.bookFutDict[futName].mult * self.bookFutDict[futName].todayCcy)
            # update volume impact
            if self.impactLimitOn:
                for impSec, impAlpha in zip(self.impactEwms,self.impactEwmAlphas):
                    self.bookVolEwmDict[futName+'_mkt_imp_'+str(impSec)] = impAlpha * abs(tradeSize) + self.bookVolEwmDict[futName+'_mkt_imp_'+str(impSec)]

    def progressForwardOneDay(self,curDate):
        # reset volume dictionaries
        self.mktVolEwmDict = {fut+'_mkt_imp_'+str(x):0.01 for x in self.impactEwms for fut in self.allFutNames}
        self.bookVolEwmDict = {fut+'_mkt_imp_'+str(x):0 for x in self.impactEwms for fut in self.allFutNames}
        self.futVolI4 = [0.0] * len(self.allFutNames)
        # update strategy for next day
        for strat in self.stratList:
            strat.progressForwardOneDay(curDate)
        # update data for each future
        for fut in self.bookFutDict.values():
            fut.progressForwardOneDay(curDate,self.closePxUnadj.loc[curDate,fut.name],self.adjFactor.loc[curDate,fut.name])
        # compute trading stats and reset trade list for today
        self.bookTradeList.progressForwardOneDay(curDate,self.bookFutDict.values(),self.bookBetaPnl)
        # update book level info at start of new day
        self.curBookBeta = sum([strat.computeStratBeta() for strat in self.stratList])
        self.curBookGmv = self.computeBookGmv()
        self.curBookRf = self.computeBookRf()
        self.curBookBf = self.computeBookBf()

    def calcBookStats(self):
        tempBeta = sum([strat.curStratBeta for strat in self.stratList])
        if tempBeta == tempBeta:
            self.curBookBeta = tempBeta
        tempGmv = self.computeBookGmv()
        if tempGmv == tempGmv:
            self.curBookGmv = tempGmv
        tempRf = self.computeBookRf()
        if tempRf == tempRf:
            self.curBookRf = tempRf
        tempBf = self.computeBookBf()
        if tempBf == tempBf:
            self.curBookBf = tempBf

    def computeBookGmv(self):
        return sum([abs(fut.currentPosition) * fut.lastValid['Mid'] * fut.mult * fut.todayCcy for fut in self.bookFutDict.values()])

    def computeBookRf(self):
        return sum([fut.currentPosition * fut.lastValid['Mid'] * fut.mult * fut.todayCcy * self.riskFactors[fut.name] for fut in self.bookFutDict.values()])

    def computeBookBf(self):
        return [sum([fut.currentPosition * fut.lastValid['Mid'] * fut.mult * fut.todayCcy * self.barraFactors[factor].get(fut.name,0) for fut in self.bookFutDict.values()]) for factor in self.barraFactorNameList]

    def computeBetaMatrix(self,futNameList,closePx):
        closePx = closePx / closePx.ffill().shift(1) - 1.0
        closePx['Bench'] = closePx[futNameList].mean(axis=1)
        beta = pd.DataFrame()
        for col in futNameList:
            beta[col] = 0.3 + 0.7 * 0.5 * ((closePx[col].rolling(window=252,min_periods=63).cov(closePx['Bench']) / closePx['Bench'].rolling(window=252,min_periods=63).var()) + (closePx[col].rolling(window=63,min_periods=21).cov(closePx['Bench']) / closePx['Bench'].rolling(window=63,min_periods=21).var()))#closePx[col].rolling(window=252,min_periods=61).cov(closePx['Bench']) / closePx['Bench'].rolling(window=252,min_periods=61).var()
            beta[col] = beta[col].apply(lambda x: min(max(x,0.7),1.3))
        #
        beta = beta.shift(1).bfill()
        beta['Month'] = [x[0:7] for x in beta.index]
        beta = beta.drop_duplicates('Month',keep='first')
        beta = beta.reindex(closePx.index,method='ffill').drop('Month',axis=1)
        # assign beta dictionaries to each future
        for futName, fut in self.bookFutDict.items():
            fut.betaDict = beta[futName].to_dict()

    def getOutput(self):
        # output book level stats
        stats = self.bookTradeList.printStatisticsAndOutputFiles('Book_'+str(self.sampleInterval)+'S')
        # output strat level stats
        for strat in self.stratList:
            stats = pd.concat([stats, strat.getOutput()],axis=1)
        #
        stats.columns = ['Book'] + [x.stratName for x in self.stratList]
        print(stats)
        return stats

