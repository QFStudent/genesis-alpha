'''
Created on Oct 7, 2019

module for individual Future class

@author: john.casale
'''
import numpy as np
import pandas as pd
import datetime, pytz

class Future:
    multDict = {'ES':50,'NQ':20,'YM':5,'DM':100,'TFS':50,'SXF':200,\
                'JNM':100, 'JTI':10000, 'YAP':25, 'HCEI':50, 'HSI':50, 'KS':250000, 'STW':100, 'TX':200, 'SSG':100, 'SFC':1,\
                'STXE':10, 'FFI':10, 'FDX':25, 'FSMI':10, 'IFS':5, 'FXXP':50, 'FESB':50, 'FDXM':5, 'AEX':200, 'FCE':10, 'MFXI':10, 'OMX30':100, 'TOP40':10}
    ccyDict = {'ES':'USD','NQ':'USD','YM':'USD','DM':'USD','TFS':'USD','SXF':'CAD',\
               'JNM':'JPY', 'JTI':'JPY', 'YAP':'AUD', 'HCEI':'HKD', 'HSI':'HKD', 'KS':'KRW', 'STW':'USD', 'TX':'TWD', 'SSG':'SGD', 'SFC':'USD',\
               'STXE':'EUR', 'FFI':'GBP', 'FDX':'EUR', 'FSMI':'CHF', 'IFS':'EUR', 'FXXP':'EUR', 'FESB':'EUR', 'FDXM':'EUR', 'AEX':'EUR', 'FCE':'EUR', 'MFXI':'EUR', 'OMX30':'SEK', 'TOP40':'ZAR'}

    def __init__(self, strat, name, basisDict, params, startDate, fx):
        # future specific parameters
        self.name = name
        self.basisDict = basisDict
        self.mult = self.multDict[name]
        #
        self.setFutureParams(strat, params)
        # initialize other parameters
        self.lastValid = {'Mid':0,'TradeTime':datetime.datetime(2016,1,1)}
        self.sizePosMult = 1.0
        self.posDir = 0
        self.posTime = 0
        self.currentPosition = 0
        self.todayEodPosition = 0
        self.yestEodPosition = 0
        self.todayClosePx = 0
        self.yestClosePx = 0
        self.todayAdjFactor = np.NaN
        self.yestAdjFactor = np.NaN
        self.todayCumLimitCount = {'Impact':0,'TradeCooloff':0,'MaxStratPos':0,'MaxBookPos':0,'MaxStratBeta':0,'MaxBookBeta':0,'MaxBookGmv':0,'MaxStratRiskFactor':0,'MaxBookRiskFactor':0}
        #
        self.fx = fx
        self.ccyPair = self.ccyDict[name]+'USD'
        self.yestCcy = self.todayCcy = self.fx.loc[startDate, self.ccyPair]
        #
        self.closeOutPositionPct = 0

    def setFutureParams(self, strat, params):
        if len(params.keys()) > 0:
            self.indexName = params['indexName']
            if strat.useProportionalMaxPos:
                self.maxStratPos = strat.maxBookConsumption * strat.book.maxBookPos[self.name]
            else:
                self.maxStratPos = params['maxStratPos']
            self.margin = params['margin']
            self.maxTradeSize = params['maxTradeSize']
            self.i2Weight = params['i2Weight']
            self.i1Weight = params['i1Weight']
            self.muI1 = params['muI1']
            self.capI1 = params['capI1']
            self.muI3 = params['muI3']
            self.riskFactor = strat.book.riskFactors[self.name]
            self.barraFactor = {x: 0 for x in strat.book.barraFactorNameList}
            self.muI5 = params.get('muI5', 0)
            self.stdI5 = params.get('stdI5', 0.10)
            self.iExtMult = params.get('iExt', 0.0)
        # future agnostic parameters
        if not strat is None:
            self.i1EwmAlpha = 1.0 - np.exp(-2.0 * strat.book.sampleInterval / strat.book.i1EwmSecs)
            self.i4EwmAlpha = strat.book.ewmI4Alpha
            if strat.book.region == 'US':
                self.i2AvgSignalPeriod = strat.rtySignalMin * 60
                self.i2AvgSignalTime = (datetime.datetime(2000,1,1,9,30) + datetime.timedelta(seconds=self.i2AvgSignalPeriod)).strftime('%H:%M:%S')
            else:
                self.i2AvgSignalPeriod = np.NaN
                self.i2AvgSignalTime = '23:59:59'
            self.bookLevel = False
            self.useLiveBasis = strat.useLiveBasis
            if self.useLiveBasis:
                self.basisLiveDict = strat.basisLive[self.name].to_dict()
            #
            self.startTimeLocal = params.get('startTime', strat.startTime)
            self.endTimeLocal = params.get('endTime', strat.endTime)
            self.timezoneLocal = pytz.timezone(params.get('timezone', strat.timezone))
            self.timezoneGlobal = pytz.timezone(strat.book.timezone)
        else:
            self.i1EwmAlpha = np.NaN
            self.i4EwmAlpha = np.NaN
            self.i2AvgSignalPeriod = np.NaN
            self.i2AvgSignalTime = '23:59:59'
            self.bookLevel = True
            self.useLiveBasis = False

    def updateParams(self, strat, params):
        self.setFutureParams(strat, params)

    def initializeTodayValues(self,curDate,barraFactor,riskFactor):
        self.riskFactor = riskFactor
        #self.sizePosMult = sizePosMult #if (self.name in ['DM', 'TFS']) else 1.0
        self.barraFactor = barraFactor
        self.beta = self.betaDict[curDate]
        try:
            self.basis = self.basisDict[curDate]
        except:
            self.basis = np.NaN
        #
        self.firstMidNotFound = True
        self.curValues = {'Mid':np.NaN,'Hs':np.NaN,'Index':np.NaN,'MidEwma':np.NaN,'MidI4':np.NaN,'i4':0, 'i2AvgSum':0, 'i2AvgCount':0, 'iExt':0}
        self.todayCumLimitCount = {'Impact':0,'TradeCooloff':0,'MaxStratPos':0,'MaxBookPos':0,'MaxStratBeta':0,'MaxBookBeta':0,'MaxBookGmv':0,'MaxStratRiskFactor':0,'MaxBookRiskFactor':0}
        #
        try:
            self.startTime = (self.timezoneLocal.localize(pd.to_datetime(curDate+' '+self.startTimeLocal))).astimezone(self.timezoneGlobal).strftime('%H:%M:%S')
            self.endTime = (self.timezoneLocal.localize(pd.to_datetime(curDate+' '+self.endTimeLocal))).astimezone(self.timezoneGlobal).strftime('%H:%M:%S')
            self.posTime = 3600 * int(self.startTime[0:2]) + 60 * int(self.startTime[3:5])
        except:
            self.startTime = '8:00:00'
            self.endTime = '18:00:00'
            self.posTime = 0

    def updateBasisWithLiveValue(self,curDate,curTime):
        try:
            curBasis = self.basisLiveDict[(curDate,curTime)]
            if (curBasis == curBasis):
                self.basis = curBasis
        except:
            return

    def processBarData(self,mid,hs,index,bidSize,askSize,curTime,ratio,iExt):
        # update last valid values if mid is a valid number
        if (mid == mid):
            self.lastValid['Mid'] = mid
        # if book level future, don't need to update anything else
        if self.bookLevel:
            return
        # update current computations
        if self.firstMidNotFound and (mid == mid):
            self.curValues['MidEwma'] = mid
            self.curValues['MidI4'] = mid
            self.firstMidNotFound = False
        elif (mid == mid):
            self.curValues['MidEwma'] = self.i1EwmAlpha * mid + (1.0 - self.i1EwmAlpha) * self.curValues['MidEwma']
            self.curValues['MidI4'] = self.i4EwmAlpha * mid + (1.0 - self.i4EwmAlpha) * self.curValues['MidI4']
        # update current values
        self.curValues['Mid'] = mid
        self.curValues['Hs'] = hs
        self.curValues['Index'] = index
        self.curValues['BidSize'] = bidSize
        self.curValues['AskSize'] = askSize
        if (iExt == iExt):
            self.curValues['iExt'] = iExt
        #
        if (curTime < self.startTime or curTime > self.endTime):
            self.curValues['i2Raw'] = np.NaN
            self.curValues['i1Raw'] = np.NaN
            self.curValues['i3'] = np.NaN
            self.curValues['i4Raw'] = np.NaN
            self.curValues['i5'] = np.NaN
        else:
            self.curValues['i2Raw'] = (mid - index - self.basis) / (index + self.basis)
            self.curValues['i1Raw'] = mid / self.curValues['MidEwma'] - 1.0
            self.curValues['i4Raw'] = mid / self.curValues['MidI4'] - 1.0
            try:
                self.curValues['i3'] = (bidSize - askSize) / (bidSize + askSize)
            except:
                self.curValues['i3'] = 0
            self.curValues['i5'] = self.muI5 * min(max(ratio / self.stdI5 / 4.0, -1.0), 1.0)
            #
            if curTime <= self.i2AvgSignalTime:
                if (self.curValues['i2Raw'] == self.curValues['i2Raw']):
                    self.curValues['i2AvgSum'] += self.curValues['i2Raw']
                    self.curValues['i2AvgCount'] += 1.0

    def progressForwardOneDay(self,curDate,closePx,adjFact):
        if self.closeOutPositionPct > 0:
            self.currentPosition = int(self.currentPosition * (1.0 - self.closeOutPositionPct))
        if self.name == 'YAP' and curDate == '2020-07-27':
            self.currentPosition = 0
        self.yestClosePx = self.todayClosePx
        self.todayClosePx = closePx
        self.yestEodPosition = self.todayEodPosition
        self.todayEodPosition = self.currentPosition
        self.yestAdjFactor = self.todayAdjFactor if not np.isnan(self.yestAdjFactor) else adjFact
        self.todayAdjFactor = adjFact
        self.yestCcy = self.todayCcy
        self.todayCcy = self.fx.loc[curDate, self.ccyPair]

    def futureAuditString(self,retHeaders=False):
        pairs = [
                ('Name',self.name),\
                ('PosContracts',str(self.currentPosition)),\
                ('PosUsd',str(self.currentPosition*self.lastValid['Mid']*self.mult*self.todayCcy)),\
                ('Fx',str(self.todayCcy)),\
                ('Index',str(self.curValues['Index'])),\
                ('IndexUnknown',str(0)),\
                ('Mid',str(self.curValues['Mid'])),\
                ('Hs',str(self.curValues['Hs'])),\
                ('BidSize',str(self.curValues['BidSize'])),\
                ('AskSize',str(self.curValues['AskSize'])),\
                ('Basis',str(self.basis)),\
                ('Beta',str(self.beta)),\
                ('BetaUsd',str(self.currentPosition*self.lastValid['Mid']*self.mult*self.beta)),\
                ('RiskFactor',str(self.riskFactor)),\
                ('RiskFactorUsd',str(self.currentPosition*self.lastValid['Mid']*self.mult*self.riskFactor)),\
                ('MidEwma',str(self.curValues['MidEwma'])),\
                ('i2Raw',str(self.curValues['i2Raw'])),\
                ('i2',str(self.curValues['i2'])),\
                ('i1Raw',str(self.curValues['i1Raw'])),\
                ('i1',str(self.curValues['i1'])),\
                ('i3',str(self.curValues['i3'])),\
                ('tajAlpha',str(self.curValues['tajAlpha'])),\
                ('taj',str(self.curValues['taj'])),\
                ('Margin',str(self.margin)),\
                ('AdjMargin',str(self.curValues['AdjMargin'])),\
                ('Cost',str(self.curValues['Cost'])),\
                ('LastTradeTime',str(self.lastValid['TradeTime'])),\
                ('ImpactLimit',str(self.todayCumLimitCount['Impact'])),\
                ('TradeCooloffLimit',str(self.todayCumLimitCount['TradeCooloff'])),\
                ('StratPosLimit',str(self.todayCumLimitCount['MaxStratPos'])),\
                ('BookPosLimit',str(self.todayCumLimitCount['MaxBookPos'])),\
                ('StratBetaLimit',str(self.todayCumLimitCount['MaxStratBeta'])),\
                ('BookBetaLimit',str(self.todayCumLimitCount['MaxBookBeta'])),\
                ('StratRiskFactorLimit',str(self.todayCumLimitCount['MaxStratRiskFactor'])),\
                ('BookRiskFactorLimit',str(self.todayCumLimitCount['MaxBookRiskFactor'])),\
                ('StratBarraFactorLimit','0'),\
                ('BookBarraFactorLimit','0'),\
                ('BookGmvLimit',str(self.todayCumLimitCount['MaxBookGmv'])),\
                ('CumVolume',str(0))
                ]
        zipPairs = list(zip(*pairs))
        if retHeaders:
            return ','.join([self.name+'_'+x for x in zipPairs[0]])
        return ','.join(zipPairs[1])

