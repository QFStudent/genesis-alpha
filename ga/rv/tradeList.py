'''
Created on Oct 9, 2019

module to save trade level data for futures

@author: john.casale
'''
import pandas as pd
import numpy as np
import os

class TradeList:
    tCostFixed = 2.5
    tCostFixedDict = {'ES':2.5, 'NQ':2.5, 'YM':2.5, 'DM':2.5, 'TFS':2.5, \
                        'JNM':250.0, 'JTI':200.0, 'YAP':2.0, 'HCEI':16.0, 'HSI':16.0, 'KS':3200.0,\
                        'STXE':0.89, 'FFI':0.61, 'FDX':1.24, 'FSMI':0.81, 'IFS':1.09, 'FXXP':0.64, 'FESB':0.57, 'FDXM':0.53, 'AEX':1.2, 'FCE':0.62, 'MFXI':1.59, 'OMX30':10.41}

    def __init__(self, futList, closePxAdj, closePxUnadj):
        self.futList = futList
        self.rawTradeList = []
        self.closePxAdj = closePxAdj
        self.closePxUnadj = closePxUnadj
        self.adjFactor = self.closePxAdj / self.closePxUnadj
        self.retAdj = (closePxAdj / closePxAdj.shift(1).ffill() - 1.0)
        self.resetTodayTradeDicts()
        #
        self.pnlIntraDict = {}
        self.pnlExtraDict = {}
        self.pnlExtraOldDict = {}
        self.pnlBetaDict = {}
        self.pnlAllDict = {}
        self.gmvDict = {}
        self.deltaDict = {}
        self.t0Dict = {}
        self.varDict = {}
        self.tradeLimitDict = {}
        self.tCostDict = {}

    def addTradeToday(self,curDate,curTime,futName,tradeSize,tradePrice):
        self.rawTradeList.append([curDate,curTime,futName,tradeSize,tradePrice])
        if tradeSize > 0:
            self.todayAvgBuyPxDict[futName] = ((self.todayAvgBuyPxDict[futName]*self.todayBuyVolDict[futName]) + (tradeSize*tradePrice)) / (self.todayBuyVolDict[futName]+tradeSize)
            self.todayBuyVolDict[futName] += tradeSize
        elif tradeSize < 0:
            self.todayAvgSellPxDict[futName] = ((self.todayAvgSellPxDict[futName]*self.todaySellVolDict[futName]) - (tradeSize*tradePrice)) / (self.todaySellVolDict[futName]-tradeSize)
            self.todaySellVolDict[futName] -= tradeSize
        #
        self.todayTCostDict[futName] += abs(tradeSize) * self.tCostFixedDict[futName]

    def progressForwardOneDay(self,curDate,futObjList,betaPnl):
        self.computeAndSaveTodayT0Gmv(curDate,futObjList)
        self.computeAndSaveTodayPnl(curDate,futObjList)
        self.computeAndSaveTodayVar(curDate,futObjList)
        self.saveCumTradeLimits(curDate,futObjList)
        self.pnlBetaDict[curDate] = betaPnl
        #
        self.resetTodayTradeDicts()

    def resetTodayTradeDicts(self):
        self.todayBuyVolDict = {x:0 for x in self.futList}
        self.todaySellVolDict = {x:0 for x in self.futList}
        self.todayAvgBuyPxDict = {x:0 for x in self.futList}
        self.todayAvgSellPxDict = {x:0 for x in self.futList}
        self.todayTCostDict = {x:0 for x in self.futList}

    def computeAndSaveTodayPnl(self,curDate,futObjList):
        # compute & save intraday, extraday pnl
        pnlIntra, pnlExtra, pnlExtraOld, pnlAll = {}, {}, {}, {}
        for fut in futObjList:
            todayClosePx, yestClosePx = fut.todayClosePx, fut.yestClosePx
            todayAdj, yestAdj = fut.todayAdjFactor, fut.yestAdjFactor
            self.todayTCostDict[fut.name] *= fut.todayCcy
            pnlIntra[fut.name] = ((todayClosePx - self.todayAvgBuyPxDict[fut.name]) * self.todayBuyVolDict[fut.name] +
                                (self.todayAvgSellPxDict[fut.name] - todayClosePx) * self.todaySellVolDict[fut.name]) * fut.mult * fut.todayCcy - self.todayTCostDict[fut.name]
            pnlExtra[fut.name] = (todayClosePx - (yestClosePx * (yestAdj / todayAdj))) * fut.todayCcy * fut.yestEodPosition * fut.mult
            pnlExtraOld[fut.name] = (todayClosePx * fut.todayCcy - (yestClosePx * fut.yestCcy * (yestAdj / todayAdj))) * fut.yestEodPosition * fut.mult
            pnlAll[fut.name] = pnlIntra[fut.name] + pnlExtra[fut.name]
        #
        self.pnlIntraDict[curDate] = pnlIntra
        self.pnlExtraDict[curDate] = pnlExtra
        self.pnlExtraOldDict[curDate] = pnlExtraOld
        self.pnlAllDict[curDate] = pnlAll
        self.tCostDict[curDate] = self.todayTCostDict

    def computeAndSaveTodayT0Gmv(self,curDate,futObjList):
        # compute & save turnover, GMV
        t0, gmv, delta = {}, {}, {}
        for fut in futObjList:
            t0[fut.name] = (self.todayBuyVolDict[fut.name] + self.todaySellVolDict[fut.name]) * fut.todayClosePx * fut.mult * fut.todayCcy
            delta[fut.name] = fut.currentPosition * fut.todayClosePx * fut.mult * fut.todayCcy
            gmv[fut.name] = abs(delta[fut.name])
        #
        self.t0Dict[curDate] = t0
        self.gmvDict[curDate] = gmv
        self.deltaDict[curDate] = delta

    def computeAndSaveTodayVar(self,curDate,futObjList):
        posDict = {}
        tempFutList = []
        for fut in futObjList:
            posDict[fut.name] = fut.currentPosition * fut.todayClosePx * fut.mult * fut.todayCcy
            tempFutList.append(fut.name)
        var = self.retAdj[tempFutList].dot(pd.Series(posDict))
        dateInd = var.index.get_loc(curDate)
        try:
            self.varDict[curDate] = var.iloc[max(0,dateInd-252):dateInd].nsmallest(2)[-1]
        except:
            self.varDict[curDate] = 0

    def saveCumTradeLimits(self,curDate,futObjList):
        tempLimitDict = {fut.name: fut.todayCumLimitCount for fut in futObjList}
        self.tradeLimitDict[curDate] = tempLimitDict

    def printStatisticsAndOutputFiles(self,fname):
        self.pnlAll = self.createDataframeFromDict(self.pnlAllDict)
        self.pnlIntra = self.createDataframeFromDict(self.pnlIntraDict)
        self.pnlExtra = self.createDataframeFromDict(self.pnlExtraDict)
        self.pnlExtraOld = self.createDataframeFromDict(self.pnlExtraOldDict)
        self.tCost = self.createDataframeFromDict(self.tCostDict)
        self.t0 = self.createDataframeFromDict(self.t0Dict)
        self.gmv = self.createDataframeFromDict(self.gmvDict)
        self.delta = self.createDataframeFromDict(self.deltaDict)
        self.pnlBeta = pd.Series(self.pnlBetaDict)
        # output all the raw file data
        self.outputAllTradeFiles(fname)
        # compute statics
        return self.printStatistics(fname)

    def printStatistics(self,fname):
        self.pnlAll['DD'] = self.pnlAll['Total'].cumsum() - self.pnlAll['Total'].cumsum().cummax().apply(lambda x: max(x,0))
        self.pnlAll['DDDays'] = (self.pnlAll['DD'] == 0).cumsum()
        self.pnlAll['DDDays'] = self.pnlAll.groupby('DDDays')['DDDays'].cumcount() + 1
        try:
            sortino1 = self.pnlAll['Total'].mean()/np.sqrt(self.pnlAll['Total'].where(self.pnlAll['Total']<0,0).pow(2).sum() / len(self.pnlAll['Total']))*np.sqrt(252)
            sortino2 = self.pnlAll['Total'].mean()/self.pnlAll[self.pnlAll['Total']<0]['Total'].std()*np.sqrt(252)
        except:
            sortino1 = np.NaN
            sortino2 = np.NaN
        try:
            betaSharpe = self.pnlBeta.mean() / self.pnlBeta.std() * np.sqrt(252)
        except:
            betaSharpe = np.NaN
        #
        try:
            statPairs = [
                        ('AvgPnl',self.pnlAll['Total'].mean()),\
                        ('StdPnl',self.pnlAll['Total'].std()),\
                        ('Sharpe',self.pnlAll['Total'].mean()/self.pnlAll['Total'].std()*np.sqrt(252)),\
                        ('BetaSharpe',betaSharpe),\
                        ('SpreadSharpe',(self.pnlAll['Total']-self.pnlBeta).mean() / (self.pnlAll['Total']-self.pnlBeta).std() * np.sqrt(252)),\
                        ('Sortino1',sortino1),\
                        ('Sortino2',sortino2),\
                        ('AvgYearPnl',self.pnlAll['Total'].mean() * 252.0),\
                        ('TotalPnl',self.pnlAll['Total'].sum()),\
                        ('TotalBetaPnl',self.pnlBeta.sum()),\
                        ('TotalSpreadPnl',(self.pnlAll['Total']-self.pnlBeta).sum()),\
                        ('TotalTCost',self.tCost['Total'].sum()),\
                        ('AvgT0',self.t0['Total'].mean()),\
                        ('AvgGmv',self.gmv['Total'].mean()),\
                        ('Pnl/T0',self.pnlAll['Total'].sum()/self.t0['Total'].sum()),\
                        ('Pnl/Gmv',self.pnlAll['Total'].sum()/(self.gmv['Total'].sum()+1)*252.0),\
                        ('AvgVar',np.mean(list(self.varDict.values()))),\
                        ('MaxVar',min(self.varDict.values())),\
                        ('MaxDraw_5%',self.pnlAll['DD'].quantile(0.05, 'linear')),\
                        ('MaxDraw_1%',self.pnlAll['DD'].quantile(0.01, 'linear')),\
                        ('MaxDraw',self.pnlAll['DD'].min()),\
                        ('MaxDraw/AvgGmv',self.pnlAll['DD'].min()/(self.gmv['Total'].mean()+1)),\
                        ('MaxDrawDays',self.pnlAll['DDDays'].max()),\
                        ('HP',2.0*self.gmv['Total'].sum()/self.t0['Total'].sum()),\
                        ('NumDays',len(self.pnlAll))
                        ]
            df = pd.DataFrame(statPairs,columns=['Stat','Value'])
            #print(df.set_index('Stat'))
            df.to_csv(r'/datadrive/BasisArb/Logging/%s_pnl_summary.csv' % (fname),index=False)
            return df.set_index('Stat')
        except:
            return pd.Series(np.NaN, index=['AvgPnl','StdPnl','Sharpe','Sortino1','Sortino2','TotalPnl','TotalIntraPnl',\
                            'TotalExtraPnl','TotalTCost','AvgT0','AvgGmv','Pnl/T0','Pnl/Gmv','AvgVar',\
                            'MaxVar','MaxDraw','MaxDraw/AvgGmv','MaxDrawDays','HP','NumDays'])

    def outputAllTradeFiles(self,fname):
        if not os.path.exists(r'/datadrive/BasisArb/Logging/Detail_%s/' % (fname.split('_')[1])):
            os.makedirs(r'/datadrive/BasisArb/Logging/Detail_%s/' % (fname.split('_')[1]))
        #
        self.pnlAll.to_csv(r'/datadrive/BasisArb/Logging/Detail_%s/%s_pnl.csv' % (fname.split('_')[1], fname),index=True)
        self.pnlIntra.to_csv(r'/datadrive/BasisArb/Logging/Detail_%s/%s_pnlIntra.csv' % (fname.split('_')[1], fname),index=True)
        self.pnlExtra.to_csv(r'/datadrive/BasisArb/Logging/Detail_%s/%s_pnlExtra.csv' % (fname.split('_')[1], fname),index=True)
        self.pnlExtraOld.to_csv(r'/datadrive/BasisArb/Logging/Detail_%s/%s_pnlExtraOld.csv' % (fname.split('_')[1], fname),index=True)
        self.t0.to_csv(r'/datadrive/BasisArb/Logging/Detail_%s/%s_t0.csv' % (fname.split('_')[1], fname),index=True)
        self.gmv.to_csv(r'/datadrive/BasisArb/Logging/Detail_%s/%s_gmv.csv' % (fname.split('_')[1], fname),index=True)
        #
        df = pd.DataFrame(self.rawTradeList,columns=['Date','Time','Name','Size','Price'])
        df.to_csv(r'/datadrive/BasisArb/Logging/Detail_%s/%s_tradeList.csv' % (fname.split('_')[1], fname),index=False)

    def createDataframeFromDict(self,curDict):
        df = pd.DataFrame.from_dict(curDict, orient='index')
        df['Total'] = df.sum(axis=1)
        df.index.name = 'Date'
        return df

