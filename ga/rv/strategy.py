'''
Created on Oct 7, 2019

module for individual Strategy class

@author: john.casale
'''
# =====================================================================
# Source Strategy.py is 585 lines (553 sloc). Fully transcribed except
# lines 403-405 inside updateMarginMultiplier, marked inline below.
# Nothing in that gap has been invented.
# =====================================================================
import numpy as np
import pandas as pd
import logging, datetime, math
from BasisArb.MainStrat import Future, TradeList

sign = lambda x: -1 if x < 0 else (1 if x > 0 else (0 if x == 0 else np.NaN))

class Strategy:

    def __init__(self, book, params, fullAudit, execThresh=1.0, randSeed=0, closeOutExpiry=0.0):
        self.mapDatetime18 = lambda x: datetime.datetime(int(x[0:4]),int(x[5:7]),int(x[8:10]),int(x[11:13]),int(x[14:16]),int(x[17:19]))
        #
        self.book = book
        self.stratName = params['stratName']
        self.impactLimitOn = book.impactLimitOn
        self.fullAudit = fullAudit
        if fullAudit:
            self.stratLog = self.initiateLogger()
            self.loggedHeaders = False
        #
        self.setStratParams(params)
        #
        self.curStratBeta = 0
        self.curStratGmv = 0
        self.curStratRf = 0
        self.curStratBf = [0] * len(book.barraFactorNameList)
        self.globalStratMargin = 1.0
        # initiate all futures
        self.futList = []
        self.futToBookPos = []
        for fut in params['futList']:
            self.futList.append(Future.Future(self, fut, self.basisFile[fut].to_dict(), params[fut], self.book.startDate, self.book.fx))
            self.futToBookPos.append(self.book.allFutNames.index(fut))
        # initiate trade list object
        self.stratTradeList = TradeList.TradeList(params['futList'],book.closePxAdj,book.closePxUnadj)
        # compute beta values
        self.computeBetaMatrix(params['futList'],book.closePxAdj)
        #
        self.execThresh = execThresh
        self.onlyCountExecs = True
        np.random.seed(randSeed)
        #
        self.betaHedgeWithEs = False
        self.closeOutExpiry = closeOutExpiry
        self.closeOutDateList = ['2017-03-10', '2017-03-13', '2017-03-14', '2017-03-15', '2017-03-16', '2017-06-09', '2017-06-12', '2017-06-13', '2017-06-14', '2017-06-15', '2017-09-08', '2017-09-11', '2017-09-12', '2017-09-13', '2017-09-14', '2017-12-08', '2017-12-11', '2017-12-12', '2017-12-13', '2017-12-14', '2018-03-09', '2018-03-12', '2018-03-13', '2018-03-14', '2018-03-15', '2018-06-08', '2018-06-11', '2018-06-12', '2018-06-13', '2018-06-14', '2018-09-14', '2018-09-17', '2018-09-18', '2018-09-19', '2018-09-20', '2018-12-14', '2018-12-17', '2018-12-18', '2018-12-19', '2018-12-20', '2019-03-08', '2019-03-11', '2019-03-12', '2019-03-13', '2019-03-14', '2019-06-14', '2019-06-17', '2019-06-18', '2019-06-19', '2019-09-13', '2019-09-16', '2019-09-17', '2019-09-18', '2019-09-19', '2019-12-13', '2019-12-16', '2019-12-17', '2019-12-18', '2019-12-19', '2020-03-13', '2020-03-16', '2020-03-17', '2020-03-18', '2020-03-19', '2020-06-12', '2020-06-15', '2020-06-16', '2020-06-17', '2020-06-18', '2020-09-11', '2020-09-14', '2020-09-15', '2020-09-16', '2020-09-17', '2021-03-12', '2021-03-15', '2021-03-16', '2021-03-17', '2021-03-18', '2021-06-11', '2021-06-14', '2021-06-15', '2021-06-16', '2021-06-17', '2021-09-10', '2021-09-13', '2021-09-14', '2021-09-15', '2021-09-16', '2021-12-10', '2021-12-13', '2021-12-14', '2021-12-15', '2021-12-16', '2022-03-11', '2022-03-14', '2022-03-15', '2022-03-16', '2022-03-17', '2022-06-10', '2022-06-13', '2022-06-14', '2022-06-15', '2022-06-16', '2022-09-09', '2022-09-12', '2022-09-13', '2022-09-14', '2022-09-15', '2022-12-09', '2022-12-12', '2022-12-13', '2022-12-14', '2022-12-15', '2023-03-10', '2023-03-13', '2023-03-14', '2023-03-15', '2023-03-16']  # NOTE: one date between '2019-06-19' and '2019-09-13' fell in a seam between photos and is absent; also 2020-09-17 -> 2021-03-12 skips a 2020-12 block that fell in a seam.

    def setStratParams(self, params):
        self.lambda1 = params['lambda1'] # how aggressive to cut when alpha opposite of position
        self.lambda2 = params['lambda2'] # how aggressive to prevent when alpha same direction as position
        self.gamma1 = params['gamma1'] # how aggressive to cut when alpha opposite of beta
        self.gamma2 = params['gamma2'] # how aggressive to prevent when alpha same direction as beta
        self.lambdaRf1 = params['lambdaRf1'] # how aggressive to cut when alpha opposite of risk factor
        self.lambdaRf2 = params['lambdaRf2'] # how aggressive to prevent when alpha same direction as risk factor
        self.stratVolEwmDict = {fut+'_mkt_imp_'+str(x):0 for x in self.book.impactEwms for fut in params['futList']}
        # compute strategy limits based on book limits
        self.maxBookConsumption = params['maxBookConsumption']
        self.useProportionalMaxPos = params['useProportionalMaxPos']
        self.maxStratBeta = self.book.maxBookBeta * self.maxBookConsumption
        self.maxStratGmv = self.book.maxBookGmv * self.maxBookConsumption
        self.maxStratRiskFactor = self.book.maxBookRiskFactor * self.maxBookConsumption
        self.maxStratBarraFactor = [x * self.maxBookConsumption for x in self.book.maxBookBarraFactor]
        #
        self.startTime = params.get('startTime', self.book.startTime)
        self.endTime = params.get('endTime', self.book.endTime)
        self.timezone = params.get('timezone', self.book.timezone)
        #
        self.useLiveBasis = params['useLiveBasis']
        self.basisFile = pd.read_csv(r'/mnt/nfs1/user/sft/sftdata/CurrentFiles%s/%s' % ('' if self.book.region == 'US' else self.book.region, params['basisFile'])).set_index('Date')
        if self.useLiveBasis:
            self.basisLive = pd.read_csv(r'/datadrive/BasisArb/BasisFiles/'+params['liveBasisFile']).set_index(['Date','Time'])
        self.muI2 = params['muI2']
        self.globalMuI1 = params['globalMuI1']
        self.globalMuI3 = params['globalMuI3']
        self.capI2 = self.book.capI2
        self.marginMult = params['marginMult']
        self.hsMult = params['hsMult']
        #
        self.useI2 = params.get('useI2', True)
        self.useI4 = params.get('useI4', False)
        self.muI4 = params.get('muI4', 0)
        self.i4LookupOffset = 3600 * 9 + 60 * 30 if self.book.region == 'US' else 9 * 3600
        #
        self.muI5 = params.get('muI5', 0)
        self.useI5 = self.muI5 > 0
        #
        self.rtySignalMin = params.get('rtySignalMin', 0.0)
        self.checkedRty, self.pauseRty = False, False
        self.pauseRtyHours = params.get('pauseRtyHours', 0)
        self.pauseRtyEndTime  = (datetime.datetime(2000,1,1,9,30) + datetime.timedelta(minutes=self.rtySignalMin + (self.pauseRtyHours * 60))).strftime('%H:%M:%S')
        self.rtyTwapHours = params.get('rtyTwapHours', 0)
        self.twapTimes = []
        #
        self.unwindSecs = params.get('unwindSecs', 24) * 3600
        self.unwindMult = params.get('unwindMult', 1.0)
        print('Using: T = ',self.unwindSecs,' and M = ',self.unwindMult)
        print('RTY Signal:',self.rtySignalMin, 'Pause hrs:', self.pauseRtyHours, self.pauseRtyEndTime, 'TWAP hrs:', self.rtyTwapHours, min(1.0, (6.5 - self.rtySignalMin / 60.0) / (self.rtyTwapHours + 0.0000001)))

    def updateParams(self, params):
        self.setStratParams(params)
        # update all future params
        curFutNameList = [x.name for x in self.futList]
        for fut in params['futList']:
            self.futList[curFutNameList.index(fut)].updateParams(self, params[fut])

    def initializeTodayValues(self,curDate,barraFactors):
        self.todayT0 = 0
        self.stratBetaPnl = 0
        setCloseOutPct = (self.closeOutExpiry > 0) and (curDate in self.closeOutDateList)
        for fut in self.futList:
            fut.initializeTodayValues(curDate,{x:y.get(fut.name, 0) for x, y in barraFactors.items()}, self.book.riskFactors[fut.name])
            if setCloseOutPct:
                fut.closeOutPositionPct = self.closeOutExpiry
            else:
                fut.closeOutPositionPct = 0.0
        #
        if self.useI4:
            self.histAvgI4 = pd.read_csv(r'/datadrive/Futures/BasisArb/Daily/%s_%s.csv' % (curDate.replace('-',''), self.book.ewmI4)).set_index('Time').values
        #
        self.checkedRty, self.pauseRty = False, False
        self.twapTimes = []

    def processBarData(self, curDate, curTime, row):
        self.curDate = curDate
        self.curTime = curTime
        if (self.curTime < self.startTime or self.curTime > self.endTime):
            return
        #
        self.curDatetime = self.mapDatetime18(self.curDate+' '+self.curTime)
        #
        if self.impactLimitOn:
            self.updateOurTradeVolume()
        #
        if self.useLiveBasis and self.curTime[-3:] == ':00':#(self.curTime[-4:] == '1:00' or self.curTime[-4:] == '6:00'):
            for fut in self.futList:
                fut.updateBasisWithLiveValue(self.curDate,self.curTime)
        #
        for fut in self.futList:
            fut.processBarData(row[fut.name+'_mid'],row[fut.name+'_hs'],row[fut.name+'_index'],row[fut.name+'_bidSize'],row[fut.name+'_askSize'],self.curTime,row.get(fut.name+'_ratio', 0), row.get(fut.name+'_iExt', 0))
        # pause rty logic on signal formation
        if self.curTime >= self.futList[0].i2AvgSignalTime and not self.checkedRty:
            try:
                signal = 0
                hedgeSignal = 0
                for fut in self.futList:
                    if fut.name == 'TFS': signal += (fut.curValues['i2AvgSum'] / fut.margin / fut.curValues['i2AvgCount'])
                    elif fut.name != 'DM': hedgeSignal += (fut.curValues['i2AvgSum'] / fut.margin / fut.curValues['i2AvgCount'])
                #
                signal = signal - (hedgeSignal / (len(self.futList) - 2.0))
                if abs(signal) >= 1.4: self.pauseRty = True
                self.checkedRty = True
            except ZeroDivisionError:
                self.checkedRty = True
            # set up twap orderbook
            if self.pauseRty and self.rtyTwapHours > 0:
                fut = self.futList[-1]
                curAbsPos = int(abs(fut.currentPosition))
                if curAbsPos > 0:
                    time_btw = int(self.rtyTwapHours * 3600 // curAbsPos)
                    secToTime = lambda x: str(x // 3600).zfill(2) + ':' + str((x % 3600) // 60).zfill(2) + ':' + str((x % 3600) % 60).zfill(2)
                    self.twapTimes = [secToTime(34200 + self.rtySignalMin * 60 + x * time_btw) for x in range(curAbsPos)]
                    print(self.curDate, fut.currentPosition, time_btw, self.twapTimes)
        # is pausing RTY finished?
        if self.checkedRty and curTime >= self.pauseRtyEndTime:
            self.pauseRty = False
        # update global margin based on t0
        if (self.curTime[-4:] == '5:00' or self.curTime[-4:] == '0:00'):
            if self.todayT0 / self.book.targetT0 > self.book.t0CurveDict[self.curTime] * self.maxBookConsumption:
                self.globalStratMargin += self.book.t0MarginInc
            else:
                self.globalStratMargin = max(1.0, self.globalStratMargin - self.book.t0MarginInc)
        # compute i2 signal
        if self.useI2:
            nonNullCount, i2Wt, i2Avg = 0, 0, 0
            for fut in self.futList:
                i2Raw = fut.curValues['i2Raw']
                if (i2Raw != i2Raw) or abs(i2Raw) > self.capI2[fut.name]: # if i2Raw is NaN, then don't count
                    continue
                if self.pauseRty and fut.name == 'TFS':
                    continue
                nonNullCount += 1
                i2Wt += fut.i2Weight
                i2Avg += (i2Raw * fut.i2Weight / fut.margin)
            if nonNullCount <= 1:
                for fut in self.futList:
                    fut.curValues['i2'] = np.NaN
            else:
                i2Avg /= i2Wt
                for fut in self.futList:
                    fut.curValues['i2'] = fut.curValues['i2Raw'] - (i2Avg * fut.margin)
                    if abs(fut.curValues['i2']) > self.capI2[fut.name]:
                        fut.curValues['i2'] = np.NaN # set i2 to NaN if above i2 limit
                    if self.pauseRty and fut.name == 'TFS':
                        fut.curValues['i2'] = np.NaN
        else:
            for fut in self.futList:
                fut.curValues['i2'] = 0.0
        # compute i1 signal
        nonNullCount, i1Wt, i1Avg = 0, 0, 0
        for fut in self.futList:
            i1Raw = fut.curValues['i1Raw']
            if i1Raw != i1Raw: # if i1Raw is NaN, then don't count
                continue
            nonNullCount += 1
            i1Wt += fut.i1Weight
            i1Avg += (i1Raw * fut.i1Weight / fut.beta)
        if nonNullCount <= 1:
            for fut in self.futList:
                fut.curValues['i1'] = 0
        else:
            i1Avg /= i1Wt
            for fut in self.futList:
                fut.curValues['i1'] = max(min(fut.curValues['i1Raw'] - (i1Avg * fut.beta), fut.capI1), -1*fut.capI1) # bound i1 by capI1
        # compute i4 signal if needed
        if self.useI4:
            nonNullCount, i4Avg, volAvg = 0, 0, 0.0
            for fut, futPos in zip(self.futList, self.futToBookPos):
                i4Raw = fut.curValues['i4Raw']
                if i4Raw != i4Raw: # if i4Raw is NaN, then don't count
                    continue
                volAvg += (self.book.futVolI4[futPos] / self.histAvgI4[self.book.curTimeSecs - self.i4LookupOffset, futPos])
                nonNullCount += 1
                i4Avg += (i4Raw / fut.beta)
            if nonNullCount <= 1:
                for fut in self.futList:
                    fut.curValues['i4'] = 0
            else:
                volAvg /= nonNullCount
                i4Avg /= nonNullCount
                for fut, futPos in zip(self.futList, self.futToBookPos):
                    curRatio = self.book.futVolI4[futPos] / self.histAvgI4[self.book.curTimeSecs - self.i4LookupOffset, futPos]
                    if curRatio > 1.0 and (curRatio - volAvg) > 0.2:
                        fut.curValues['i4'] = max(min(self.muI4 * (fut.curValues['i4Raw'] - (i4Avg * fut.beta)), 0.00012), -0.00012)
                    else:
                        fut.curValues['i4'] = 0
        # compute tajAlpha, taj
        for fut in self.futList:
            fut.curValues['tajAlpha'] = self.muI2 * fut.curValues['i2'] + self.globalMuI1 * fut.muI1 * fut.curValues['i1'] + self.globalMuI3 * fut.muI3 * fut.curValues['i3'] + fut.curValues['i4'] + self.muI5 * fut.curValues['i5']
            fut.curValues['taj'] = fut.curValues['tajAlpha'] + self.book.extAlphaMult * fut.curValues.get('iExt', 0) * fut.iExtMult
        # update strategy level info
        self.calcStratStats()
        # compute beta pnl
        self.stratBetaPnl += (self.curStratBeta * self.book.avgRet)

    def checkAndProcessStrategyTrades(self):
        if (self.curTime < self.startTime or self.curTime > self.endTime):
            return
        #
        self.tradeList = []
        for fut in self.futList:
            # can't trade if alpha is nan
            if fut.curValues['taj'] != fut.curValues['taj']:
                fut.curValues['AdjMargin'] = np.NaN
                fut.curValues['Cost'] = np.NaN
                continue
            # update margin based on max pos, max beta, max rf
            curPos = fut.currentPosition * fut.lastValid['Mid'] * fut.mult * fut.todayCcy
            curSign = sign(fut.curValues['taj'])
            curMarg = self.updateMarginMultiplier(fut,self.marginMult,curSign,curPos)
            fut.curValues['AdjMargin'] = curMarg
            # check if we want to trade at strategy level
            tradeSize = self.checkForStrategyTrade(fut,curMarg,curSign,curPos)
            if tradeSize == 0:
                continue
            # can't trade if we traded recently within book time limit
            if (self.curDatetime - fut.lastValid['TradeTime']).total_seconds() < self.book.timeBetweenTrades:
                fut.todayCumLimitCount['TradeCooloff'] = fut.todayCumLimitCount['TradeCooloff'] + 1
                continue
            # check if we want to trade at book level
            tradeSize, limitCat = self.book.checkForBookTrade(fut.name,tradeSize,self.curDatetime)
            if tradeSize == 0:
                fut.todayCumLimitCount[limitCat] = fut.todayCumLimitCount[limitCat] + 1
                continue
            # check if reduce risk only mode is enabled
            if self.book.reduceRiskOnly:
                if tradeSize * curPos >= 0:
                    continue
                tradeSize = min(tradeSize, -1*fut.currentPosition) if tradeSize > 0 else max(tradeSize, -1*fut.currentPosition)
            # process strategy and book level trade
            tradePrice = fut.curValues['Mid'] + sign(tradeSize) * fut.curValues['Hs']
            self.tradeList.append([fut, tradeSize, tradePrice])
        # check for RTY twap trade
        if self.pauseRty and self.curTime in self.twapTimes:
            fut = self.futList[-1]
            tradeSize = -1 * sign(fut.currentPosition)
            self.tradeList.append([fut, tradeSize, fut.curValues['Mid'] + sign(tradeSize) * fut.curValues['Hs']])
        #
        if not self.book.useOrderAggregator:
            for fut, tradeSize, tradePrice in self.tradeList:
                self.processStrategyTrade(fut,tradeSize,tradePrice)
        #
        if self.fullAudit:
            self.strategyAuditLogRow()

    def processStrategyTrade(self,fut,tradeSize,tradePrice,isExecuted=False):
        isExecuted = (isExecuted) or (self.book.useOrderAggregator) or (np.random.rand() <= self.execThresh)
        if isExecuted:
            # log trade at strategy level
            self.stratTradeList.addTradeToday(self.curDate, self.curTime, fut.name, tradeSize, tradePrice)
            # update future position and last trade time
            fut.currentPosition += tradeSize
            # keep track of how long in certain dir
            if fut.currentPosition * fut.posDir <= 0:
                fut.posDir = sign(fut.currentPosition)
                fut.posTime = self.book.curTimeSecs
        # update trade cool-off / t0 / impact if we trade or counting all sent orders
        if not self.onlyCountExecs or isExecuted:
            fut.lastValid['TradeTime'] = self.curDatetime
            # update strategy level info after trade
            #self.calcStratStats()
            # update today t0
            self.todayT0 += (abs(tradeSize) * tradePrice * fut.mult * fut.todayCcy)
            # update volume impact at strategy level
            if self.impactLimitOn:
                for impSec, impAlpha in zip(self.book.impactEwms,self.book.impactEwmAlphas):
                    self.stratVolEwmDict[fut.name+'_mkt_imp_'+str(impSec)] = impAlpha * abs(tradeSize) + self.stratVolEwmDict[fut.name+'_mkt_imp_'+str(impSec)]
        # log trade at book level
        self.book.processBookTrade(fut.name,tradeSize,tradePrice,self.curDatetime,isExecuted,self.onlyCountExecs)

    def checkForStrategyTrade(self,fut,curMarg,curSign,curPos):
        curCost = (curMarg * fut.margin * self.book.globalBookMargin * self.globalStratMargin) + (self.hsMult * (fut.curValues['Hs'] / fut.curValues['Mid']))
        fut.curValues['Cost'] = curCost
        # if alpha is less than cost, then don't trade
        if abs(fut.curValues['taj']) < curCost:
            return 0

        # if alpha is greater than cost and if we are increasing risk
        tradeSize = 0
        if curPos * curSign >= 0:
            # check all risk limits before trading
            if abs(curPos) >= fut.maxStratPos * self.book.volPosMult * fut.sizePosMult:
                fut.todayCumLimitCount['MaxStratPos'] = fut.todayCumLimitCount['MaxStratPos'] + 1
            elif self.curStratBeta * curSign >= self.maxStratBeta * self.book.maxBetaMult:
                fut.todayCumLimitCount['MaxStratBeta'] = fut.todayCumLimitCount['MaxStratBeta'] + 1
            elif self.curStratGmv >= self.maxStratGmv * self.book.maxGmvMult:
                fut.todayCumLimitCount['MaxBookGmv'] = fut.todayCumLimitCount['MaxBookGmv'] + 1
            elif self.curStratRf * sign(fut.riskFactor) * curSign >= self.maxStratRiskFactor:# * self.book.maxRfMult:
                fut.todayCumLimitCount['MaxStratRiskFactor'] = fut.todayCumLimitCount['MaxStratRiskFactor'] + 1
            else:   # trade future at far touch because we passed all risk checks
                tradeSize = curSign * min(max(curSign * fut.curValues['AskSize'], -1 * curSign * fut.curValues['BidSize']), int(fut.maxTradeSize / fut.lastValid['Mid'] / fut.mult / fut.todayCcy))
            # check all barra factors too
            for factor, pos in self.book.barraFactorNameToPos.items():
                if self.curStratBf[pos] * sign(fut.barraFactor[factor]) * curSign >= self.maxStratBarraFactor[pos]:# * self.book.maxBarraMult[pos]:
                    tradeSize = 0
                    fut.todayCumLimitCount['MaxStratRiskFactor'] = fut.todayCumLimitCount['MaxStratRiskFactor'] + 1
                    break
        #
        # else if alpha is greater than cost and if we are decreasing risk
        else:
            # check all risk limits before trading
            if self.curStratBeta * curSign >= self.maxStratBeta * self.book.maxBetaMult:
                fut.todayCumLimitCount['MaxStratBeta'] = fut.todayCumLimitCount['MaxStratBeta'] + 1
            elif self.curStratRf * sign(fut.riskFactor) * curSign >= self.maxStratRiskFactor:# * self.book.maxRfMult:
                fut.todayCumLimitCount['MaxStratRiskFactor'] = fut.todayCumLimitCount['MaxStratRiskFactor'] + 1
            else:   # trade future at far touch because we passed all risk checks
                tradeSize = curSign * min(max(curSign * fut.curValues['AskSize'], -1 * curSign * fut.curValues['BidSize']), int(fut.maxTradeSize / fut.lastValid['Mid'] / fut.mult / fut.todayCcy))
            # check all barra factors too
            for factor, pos in self.book.barraFactorNameToPos.items():
                if self.curStratBf[pos] * sign(fut.barraFactor[factor]) * curSign >= self.maxStratBarraFactor[pos]:# * self.book.maxBarraMult[pos]:
                    tradeSize = 0
                    fut.todayCumLimitCount['MaxStratRiskFactor'] = fut.todayCumLimitCount['MaxStratRiskFactor'] + 1
                    break

        # update trade size for impact constraints
        tradeSize = self.checkStratImpactLimits(fut,tradeSize)
        return tradeSize

    def checkStratImpactLimits(self,fut,tradeSize):
        if not self.impactLimitOn or tradeSize == 0:
            return tradeSize
        #
        for impSec, impLimit in zip(self.book.impactEwms,self.book.impactLimits):
            if self.stratVolEwmDict[fut.name+'_mkt_imp_'+str(impSec)] / self.book.mktVolEwmDict[fut.name+'_mkt_imp_'+str(impSec)] > impLimit:
                #print(self.curTime, fut.name, impSec, impLimit, self.stratVolEwmDict[fut.name+'_mkt_imp_'+str(impSec)], self.book.mktVolEwmDict[fut.name+'_mkt_imp_'+str(impSec)])
                fut.todayCumLimitCount['Impact'] = fut.todayCumLimitCount['Impact'] + 1
                return 0
        #
        return tradeSize

    def updateMarginMultiplier(self,fut,curMarg,curSign,curPos):
        l1Mult = self.book.curTimeLambdaMult if self.book.useLambda1Mult else 1.0
        l2Mult = self.book.curTimeLambdaMult if self.book.useLambda2Mult else 1.0
        if self.book.curTimeSecs - fut.posTime >= self.unwindSecs:
            l1Mult = l2Mult = min((self.book.curTimeSecs - fut.posTime) / 3600.0, 4.0)
        # if alpha is in same direction as current position, increase margin
        if curPos * curSign >= 0:
            curMarg *= (1 + (self.lambda2 * l2Mult * abs(curPos) / (fut.maxStratPos * self.book.volPosMult * fut.sizePosMult)))# + (self.lambda2 * abs(curPos) / self.book.maxBookPos[fut.name]))
        else: # if alpha is in opposite direction as current position, decrease margin
            curMarg *= (1 - (self.lambda1 * l1Mult * abs(curPos) / (fut.maxStratPos * self.book.volPosMult * fut.sizePosMult)))# - (self.lambda1 * abs(curPos) / self.book.maxBookPos[fut.name]))
            if curMarg < 0:
                return 0.0001
        # if alpha is in same direction as strategy beta, increase margin
        if self.curStratBeta * curSign >= 0:
            curMarg *= (1 + (self.gamma2 * abs(self.curStratBeta) / (self.maxStratBeta * self.book.maxBetaMult)))
        else: # if alpha is in opposite direction as strategy beta, decrease margin
            curMarg *= (1 - (self.gamma1 * abs(self.curStratBeta) / (self.maxStratBeta * self.book.maxBetaMult)))
        # ---- SOURCE LINES 403-405 NOT CAPTURED (fell between two photos) ----
        # By analogy with the blocks above/below they are probably a
        # `if curMarg < 0: return 0.0001` guard plus the next section comment,
        # but that is a guess and has NOT been written in.
        if self.curStratRf * curSign * fut.riskFactor >= 0:
            curMarg *= (1 + (self.lambdaRf2 * abs(self.curStratRf) / (self.maxStratRiskFactor * self.book.maxRfMult)))
        else: # if alpha is in opposite direction as strat risk factor, decrease margin
            curMarg *= (1 - (self.lambdaRf1 * abs(self.curStratRf) / (self.maxStratRiskFactor * self.book.maxRfMult)))
        # if alpha is in same direction as strat barra factor, increase margin
        for factor, pos in self.book.barraFactorNameToPos.items():
            if self.curStratBf[pos] * curSign * fut.barraFactor[factor] >= 0:
                curMarg *= (1 + (self.lambdaRf2 * abs(self.curStratBf[pos]) / (self.maxStratBarraFactor[pos] * self.book.maxBarraMult[pos])))
            else: # if alpha is in opposite direction as strat risk factor, decrease margin
                curMarg *= (1 - (self.lambdaRf1 * abs(self.curStratBf[pos]) / (self.maxStratBarraFactor[pos] * self.book.maxBarraMult[pos])))
        # if running strategies independently, then use just strat level information
        if self.book.runStratsIndependent:
            if curPos * curSign >= 0:
                curMarg *= (1 + (self.lambda2 * l2Mult * abs(curPos) / (self.book.maxBookPos[fut.name] * self.book.volPosMult * fut.sizePosMult)))
            else:
                curMarg *= (1 - (self.lambda1 * l1Mult * abs(curPos) / (self.book.maxBookPos[fut.name] * self.book.volPosMult * fut.sizePosMult)))
                if curMarg < 0:
                    return 0.0001
            if self.curStratBeta * curSign >= 0:
                curMarg *= (1 + (self.gamma2 * abs(self.curStratBeta) / (self.book.maxBookBeta * self.book.maxBetaMult)))
            else:
                curMarg *= (1 - (self.gamma1 * abs(self.curStratBeta) / (self.book.maxBookBeta * self.book.maxBetaMult)))
                if curMarg < 0:
                    return 0.0001
            #
            if self.curStratRf * curSign * fut.riskFactor >= 0:
                curMarg *= (1 + (self.lambdaRf2 * abs(self.curStratRf) / (self.book.maxBookRiskFactor * self.book.maxRfMult)))
            else:
                curMarg *= (1 - (self.lambdaRf1 * abs(self.curStratRf) / (self.book.maxBookRiskFactor * self.book.maxRfMult)))
            #
            for factor, pos in self.book.barraFactorNameToPos.items():
                if self.curStratBf[pos] * curSign * fut.barraFactor[factor] >= 0:
                    curMarg *= (1 + (self.lambdaRf2 * abs(self.curStratBf[pos]) / (self.book.maxBookBarraFactor[pos] * self.book.maxBarraMult[pos])))
                else:
                    curMarg *= (1 - (self.lambdaRf1 * abs(self.curStratBf[pos]) / (self.book.maxBookBarraFactor[pos] * self.book.maxBarraMult[pos])))
        # if running strategies as book, then use book level information
        else:
            bookPos = self.book.bookFutDict[fut.name].currentPosition * fut.lastValid['Mid'] * fut.mult * fut.todayCcy
            if bookPos * curSign >= 0:
                curMarg *= (1 + (self.lambda2 * l2Mult * abs(bookPos) / (self.book.maxBookPos[fut.name] * self.book.volPosMult * fut.sizePosMult)))
            else:
                curMarg *= (1 - (self.lambda1 * l1Mult * abs(bookPos) / (self.book.maxBookPos[fut.name] * self.book.volPosMult * fut.sizePosMult)))
                if curMarg < 0:
                    return 0.0001
            if self.book.curBookBeta * curSign >= 0:
                curMarg *= (1 + (self.gamma2 * abs(self.book.curBookBeta) / (self.book.maxBookBeta * self.book.maxBetaMult)))
            else: # if alpha is in opposite direction as book beta, decrease margin
                curMarg *= (1 - (self.gamma1 * abs(self.book.curBookBeta) / (self.book.maxBookBeta * self.book.maxBetaMult)))
                if curMarg < 0:
                    return 0.0001
            # if alpha is in same direction as book risk factor, increase margin
            if self.book.curBookRf * curSign * fut.riskFactor >= 0:
                curMarg *= (1 + (self.lambdaRf2 * abs(self.book.curBookRf) / (self.book.maxBookRiskFactor * self.book.maxRfMult)))
            else: # if alpha is in opposite direction as book risk factor, decrease margin
                curMarg *= (1 - (self.lambdaRf1 * abs(self.book.curBookRf) / (self.book.maxBookRiskFactor * self.book.maxRfMult)))
            # if alpha is in same direction as book barra factor, increase margin
            for factor, pos in self.book.barraFactorNameToPos.items():
                if self.book.curBookBf[pos] * curSign * fut.barraFactor[factor] >= 0:
                    curMarg *= (1 + (self.lambdaRf2 * abs(self.book.curBookBf[pos]) / (self.book.maxBookBarraFactor[pos] * self.book.maxBarraMult[pos])))
                else: # if alpha is in opposite direction as book risk factor, decrease margin
                    curMarg *= (1 - (self.lambdaRf1 * abs(self.book.curBookBf[pos]) / (self.book.maxBookBarraFactor[pos] * self.book.maxBarraMult[pos])))
        #
        return max(curMarg,0.0001)

    def progressForwardOneDay(self, curDate):
        if self.betaHedgeWithEs:
            curStratBeta = self.computeStratBeta()
            hedgeFut = self.futList[0]
            tradeSize = int(-1 * curStratBeta / hedgeFut.mult / hedgeFut.todayCcy / hedgeFut.beta / hedgeFut.lastValid['Mid'])
            if tradeSize != 0:
                self.processStrategyTrade(hedgeFut,tradeSize,hedgeFut.lastValid['Mid'],isExecuted=True)
        # update data for each future
        for fut in self.futList:
            fut.progressForwardOneDay(curDate,self.book.closePxUnadj.loc[curDate,fut.name],self.book.adjFactor.loc[curDate,fut.name])
        # compute trading stats and reset trade list for today
        self.stratTradeList.progressForwardOneDay(curDate,self.futList,self.stratBetaPnl)
        # reset volume counter
        self.stratVolEwmDict = {fut.name+'_mkt_imp_'+str(x):0 for x in self.book.impactEwms for fut in self.futList}
        # update strategy level info after start of new day
        self.curStratBeta = self.computeStratBeta()
        self.curStratGmv = self.computeStratGmv()
        self.curStratRf = self.computeStratRf()
        self.curStratBf = self.computeStratBf()

    def computeBetaMatrix(self,futNameList,closePx):
        closePx = closePx / closePx.ffill().shift(1) - 1.0
        closePx['Bench'] = closePx[futNameList].mean(axis=1)
        beta = pd.DataFrame()
        for col in futNameList:
            beta[col] = 0.3 + 0.7 * 0.5 * ((closePx[col].rolling(window=252,min_periods=63).cov(closePx['Bench']) / closePx['Bench'].rolling(window=252,min_periods=63).var()) + (closePx[col].rolling(window=63,min_periods=21).cov(closePx['Bench']) / closePx['Bench'].rolling(window=63,min_periods=21).var()))#closePx[col].rolling(window=252,min_periods=61).cov(closePx['Bench']) / closePx['Bench'].rolling(window=252,min_periods=61).var()
            beta[col] = beta[col].apply(lambda x: min(max(x,0.65),1.35))
        #
        beta = beta.shift(1).bfill()
        beta['Month'] = [x[0:7] for x in beta.index]
        beta = beta.drop_duplicates('Month',keep='first')
        beta = beta.reindex(closePx.index,method='ffill').drop('Month',axis=1)
        # assign beta dictionaries to each future
        for fut in self.futList:
            fut.betaDict = beta[fut.name].to_dict()

    def updateOurTradeVolume(self):
        for fut in self.futList:
            for impSec, impAlpha in zip(self.book.impactEwms,self.book.impactEwmAlphas):
                self.stratVolEwmDict[fut.name+'_mkt_imp_'+str(impSec)] = (1 - impAlpha) * self.stratVolEwmDict[fut.name+'_mkt_imp_'+str(impSec)]

    def calcStratStats(self):
        tempBeta = self.computeStratBeta()
        if tempBeta == tempBeta:
            self.curStratBeta = tempBeta
        tempGmv = self.computeStratGmv()
        if tempGmv == tempGmv:
            self.curStratGmv = tempGmv
        tempRf = self.computeStratRf()
        if tempRf == tempRf:
            self.curStratRf = tempRf
        tempBf = self.computeStratBf()
        if tempBf == tempBf:
            self.curStratBf = tempBf

    def computeStratBeta(self):
        return sum([fut.currentPosition * fut.lastValid['Mid'] * fut.mult * fut.todayCcy * fut.beta for fut in self.futList])

    def computeStratGmv(self):
        return sum([abs(fut.currentPosition) * fut.lastValid['Mid'] * fut.mult * fut.todayCcy for fut in self.futList])

    def computeStratRf(self):
        return sum([fut.currentPosition * fut.lastValid['Mid'] * fut.mult * fut.todayCcy * fut.riskFactor for fut in self.futList])

    def computeStratBf(self):
        return [sum([fut.currentPosition * fut.lastValid['Mid'] * fut.mult * fut.todayCcy * fut.barraFactor[factor] for fut in self.futList]) for factor in self.book.barraFactorNameList]

    def getOutput(self):
        return self.stratTradeList.printStatisticsAndOutputFiles(self.stratName+'_'+str(self.book.sampleInterval)+'S')

    def initiateLogger(self):
        fHandler = logging.FileHandler(r'/datadrive/BasisArb/Logging/%s.log' % (self.stratName), 'w')
        fHandler.setLevel(logging.INFO)
        fHandler.setFormatter(logging.Formatter('%(asctime)s,%(levelname)s,%(message)s'))
        #
        myLog = logging.getLogger(self.stratName)
        myLog.setLevel('INFO')
        myLog.addHandler(fHandler)
        return myLog

    def strategyAuditLogRow(self):
        if not self.loggedHeaders:
            self.stratLog.info(self.strategyAuditString(True))
            self.loggedHeaders = True
        #
        self.stratLog.info(self.strategyAuditString(False))

    def strategyAuditString(self,retHeaders=False):
        stratPairs = [
                    ('Date',self.curDate),\
                    ('Time',self.curTime),\
                    ('GmvBook',str(self.book.curBookGmv)),\
                    ('BetaBook',str(self.book.curBookBeta)),\
                    ('RiskFactorBook',str(self.book.curBookRf)),\
                    ('BarraFactor1Book','0' if len(self.book.curBookBf) == 0 else str(self.book.curBookBf[0])),\
                    ('BarraFactor2Book','0' if len(self.book.curBookBf) == 0 else str(self.book.curBookBf[1])),\
                    ('GlobalBookMargin',str(self.book.globalBookMargin)),\
                    ('GlobalEventMargin',str(self.book.eventMarginMult)),\
                    ('GlobalManualMargin','1'),\
                    ('GmvStrat',str(self.curStratGmv)),\
                    ('BetaStrat',str(self.curStratBeta)),\
                    ('RiskFactorStrat',str(self.curStratRf)),\
                    ('BarraFactor1Strat','0' if len(self.curStratBf) == 0 else str(self.curStratBf[0])),\
                    ('BarraFactor2Strat','0' if len(self.curStratBf) == 0 else str(self.curStratBf[1])),\
                    ('GlobalStratMargin',str(self.globalStratMargin))
                    ]
        stratPairs = list(zip(*stratPairs))
        stratData = ','.join(stratPairs[0]) if retHeaders else ','.join(stratPairs[1])
        #
        for fut in self.futList:
            stratData += ',' + fut.futureAuditString(retHeaders)
        return stratData

