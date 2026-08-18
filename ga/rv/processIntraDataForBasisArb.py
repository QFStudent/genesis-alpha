'''
Created on Oct 10, 2019

module to manipulate historical data into Backtester readable form

@author: john.casale
'''
import sys
sys.path.append('/home/fut_strat/git/research')

import numpy as np
import pandas as pd
import datetime, os, gc
from utils import FutureRoll, BusinessDays
from utils.funcs import *

# update barra factor exposures and BTIC data first
exec(open(r'/home/fut_strat/git/research/BasisArb/ProcessData/ProcessBticBarData.py').read())
#exec(open(r'/home/fut_strat/git/research/Equity/Barra/BarraEtfExposures.py').read())

# start custom variables
inPathFuture = r'/mnt/sftdata/cc_data15/sft_fut_bar/'#r'/datadrive/Bars/sft_fut_bar/'#
inPathIndex = r'/mnt/data/Datasets/cc_data17/index_1s_bar/'#r'/datadrive/Bars/index_bar_1s_v2_NDX/'#
outPath = r'/datadrive/Futures/BasisArb/'
dailyOutPath = r'/datadrive/Futures/BasisArb/Daily/'
fileNameExtra = ''

futList = ['ES','NQ','YM','DM','TFS']
futToIndex = {'ES':'SPX','NQ':'NDX','YM':'INDU','DM':'MID','TFS':'RUT'}
ricToFut = {'.SPX':'ES','.NDX':'NQ','.DJI':'YM','.MID':'DM','.RUT':'TFS'}
indNaNThresh = {'ES':0.002,'NQ':0.002,'YM':0.002,'DM':0.004,'TFS':0.004}
indNaNExcep = {'2017-10-03':{'TFS':0.0046},\
                '2017-06-01':{'DM':0.0050},'2017-09-13':{'DM':0.0047},'2017-09-14':{'DM':0.0047},'2017-09-15':{'DM':0.0047},\
                '2018-05-29':{'DM':0.0050},'2018-05-30':{'DM':0.0050},'2018-06-07':{'DM':0.0045},'2018-06-08':{'DM':0.0045},\
                '2017-07-25':{'ES':0.0028},'2018-06-15':{'ES':0.0035},'2018-06-18':{'ES':0.0035},'2018-06-19':{'ES':0.0035},\
                '2018-11-29':{'ES':0.0041},'2018-11-30':{'ES':0.0041},'2018-12-21':{'ES':0.0028,'NQ':0.0077},\
                '2019-05-03':{'DM':0.0061},'2019-05-06':{'DM':0.0061},'2019-07-30':{'TFS':0.0055},'2020-04-15':{'DM':0.0067},\
                '2020-05-11':{'ES':0.003}}

rollDays = 1

outputExtra = True
outputIntra = True
extraStartDate = datetime.date(2011,1,3)#(2021,5,21)#(2019,9,20)#
intraStartDate = datetime.date(2017,1,1)#(2021,5,21)#(2019,10,1)#
endDate = datetime.date.today()#datetime.date(2019,10,31)#datetime.date.today()#
startTime = '09:30:00'
endTime = '16:00:00'
endTime2 = '16:15:00'
halfDayEndTime = '13:00:00'
halfDayEndTime2 = '13:15:00'
# end custom variables

# initialize objects
bDays = BusinessDays.BusinessDays('NYT',False)
dateList = [d.date() for d in bDays.getDateListBetween(extraStartDate, endDate)]
halfDayList = [d.date() for d in bDays.getHalfdayList()]

futRollDict = {x: FutureRoll.FutureRoll(x,rollDays-1,bDays) for x in futList}
rollDateDict = {k: [bDays.getNBDays(x, -1*rollDays).date() for x in v.generateExpiryDates() if extraStartDate <= x <= endDate] for k,v in futRollDict.items()}
expSym = pd.DataFrame(index=dateList)
for fut in futList:
    expSym[fut] = futRollDict[fut].getFutureSymbolForDateList(dateList)

def computePivotTable(df,valCol,extName):
    mid = df.pivot_table(index=['Date','Time'],columns='Future',values=valCol,dropna=False)#.ffill()
    mid.columns = [x+extName for x in mid.columns]
    for fut in futList:
        if fut+extName not in mid.columns:
            mid[fut+extName] = np.NaN
    #
    isNullDict = mid.isnull().sum().to_dict()
    return mid[[x+extName for x in futList]].ffill(), isNullDict

def getNaNThresh(dateStr):
    if not dateStr in indNaNExcep.keys():
        return indNaNThresh
    #
    tempNaNThresh = indNaNThresh.copy()
    for k,v in indNaNExcep[dateStr].items():
        tempNaNThresh[k] = v
    return tempNaNThresh

def processTodayData(df, curDate, ind, extraOnly=False):
    futToSym = expSym.loc[curDate].to_dict()
    symToFut = {v:k for k,v in futToSym.items()}
    # filter unneeded data
    df = df[df['SYMBOL'].isin(symToFut.keys())]
    df['Date'] = curDate.strftime('%Y-%m-%d')
    df['Time'] = df['INTERVAL_END'].apply(lambda x: x.split(' ')[1])
    df['Future'] = df['SYMBOL'].replace(symToFut)
    df['Mid'] = 0.5 * (df['LAST_BID'] + df['LAST_ASK'])
    df['Hs'] = (0.5 * (df['LAST_ASK'] - df['LAST_BID']))
    #
    if curDate in halfDayList:
        close2 = df[df['INTERVAL_END'] <= curDate.strftime('%Y%m%d ')+halfDayEndTime2]
        df = df[(df['INTERVAL_END'] >= curDate.strftime('%Y%m%d ')+startTime) & (df['INTERVAL_END'] <= curDate.strftime('%Y%m%d ')+halfDayEndTime)]
        ind = ind[(ind['INTERVAL_END'] >= curDate.strftime('%Y%m%d ')+startTime) & (ind['INTERVAL_END'] <= curDate.strftime('%Y%m%d ')+halfDayEndTime)]
    else:
        close2 = df[df['INTERVAL_END'] <= curDate.strftime('%Y%m%d ')+endTime2]
        df = df[(df['INTERVAL_END'] >= curDate.strftime('%Y%m%d ')+startTime) & (df['INTERVAL_END'] <= curDate.strftime('%Y%m%d ')+endTime)]
        ind = ind[(ind['INTERVAL_END'] >= curDate.strftime('%Y%m%d ')+startTime) & (ind['INTERVAL_END'] <= curDate.strftime('%Y%m%d ')+endTime)]
    #
    close2, close2Null = computePivotTable(close2,'Mid','')
    close2 = close2.iloc[-1].to_dict()
    #
    mid, midNull = computePivotTable(df,'Mid','_mid')
    if extraOnly:
        return mid, futToSym, close2
    #
    hs, hsNull = computePivotTable(df,'Hs','_hs')
    bidS, bidSNull = computePivotTable(df,'LAST_BID_SIZE','_bidSize')
    askS, askSNull = computePivotTable(df,'LAST_ASK_SIZE','_askSize')
    vol, volNull = computePivotTable(df,'TRADE_VOLUME','_vol')
    #
    if len(ind) > 0:
        ind['Date'] = curDate.strftime('%Y-%m-%d')
        ind['Time'] = ind['INTERVAL_END'].apply(lambda x: x.split(' ')[1])
        ind['Future'] = ind['INDEX_RIC'].replace(ricToFut)
        tempNaNThresh = getNaNThresh(curDate.strftime('%Y-%m-%d'))
        ind['AboveThresh'] = ind['Future'].apply(lambda x: tempNaNThresh[x])
        ind['AboveThresh'] = (ind['UNKNOWN_RATIO'] > ind['AboveThresh'])
        ind.loc[ind['AboveThresh'],'INDEX_VALUE'] = np.NaN
        ind['AboveThresh'] *= 1.0
        print(ind.groupby('Future')[['AboveThresh']].sum().T)
        indVal, indValNull = computePivotTable(ind,'INDEX_VALUE','_index')
    else:
        indVal = pd.DataFrame(index=mid.index,columns=[x+'_index' for x in futList])
        indValNull = {}
    #
    nullDict = {**midNull, **hsNull, **bidSNull, **askSNull, **volNull, **indValNull}
    return pd.concat([mid,hs,bidS,askS,vol,indVal],axis=1), nullDict, futToSym, close2

def checkAndComputeRollRatio(df, curDate):
    # loop through futures to check if any need to be rolled on curDate
    for k,v in rollDateDict.items():
        if not curDate in v:
            continue
        #
        frontCont = expSym.loc[curDate,k]
        nextCont = expSym.shift(-1*rollDays).loc[curDate,k]
        if frontCont == nextCont:
            print('Error: front and next contract are matching:',curDate,frontCont,nextCont)
            #input('Press enter to continue...')
        #
        temp = df[df['SYMBOL'].isin([frontCont,nextCont])]
        if curDate in halfDayList:
            temp = temp[(temp['INTERVAL_END'] >= curDate.strftime('%Y%m%d ')+startTime) & (temp['INTERVAL_END'] <= curDate.strftime('%Y%m%d ')+halfDayEndTime)]
        else:
            temp = temp[(temp['INTERVAL_END'] >= curDate.strftime('%Y%m%d ')+startTime) & (temp['INTERVAL_END'] <= curDate.strftime('%Y%m%d ')+endTime)]
        temp['Mid'] = 0.5*(temp['LAST_BID']+temp['LAST_ASK'])
        #
        temp = temp.pivot('INTERVAL_END','SYMBOL','Mid').ffill()
        midDict = temp.rolling(window=30*60).mean().iloc[-1]
        adjFactor.at[curDate.strftime('%Y-%m-%d'),k] = midDict[nextCont] / midDict[frontCont]
    #


allData = pd.DataFrame()
openPx = pd.DataFrame()
closePx = pd.DataFrame()
adjFactor = pd.DataFrame()
futSymDict = {}
closePx2 = {}
nullDict = {}
# loop through all days
for root, dir, files in os.walk(inPathFuture):
    for fname in files:
        if fname[-9:] != '1s.csv.gz':
            continue
        #
        curDate = mapDateStr8(fname)
        if curDate < extraStartDate or curDate > endDate or curDate not in dateList:
            continue
        # switch to T-4 roll in September 2020
        if curDate == datetime.date(2020, 9, 1):
            rollDays = 4
            #
            bDays = BusinessDays.BusinessDays('NYT',False)
            dateList = [d.date() for d in bDays.getDateListBetween(extraStartDate, endDate)]
            halfDayList = [d.date() for d in bDays.getHalfdayList()]
            #
            futRollDict = {x: FutureRoll.FutureRoll(x,rollDays-1,bDays) for x in futList}
            rollDateDict = {k: [bDays.getNBDays(x, -1*rollDays).date() for x in v.generateExpiryDates() if extraStartDate <= x <= endDate] for k,v in futRollDict.items()}
            expSym = pd.DataFrame(index=dateList)
            for fut in futList:
                expSym[fut] = futRollDict[fut].getFutureSymbolForDateList(dateList)
        #
        df = pd.read_csv(root+'/'+fname)
        # if we want intra, append intra to allData then handle extra
        if intraStartDate <= curDate <= endDate:
            try:
                ind = pd.read_csv(inPathIndex+curDate.strftime('%Y')+'/'+curDate.strftime('%m')+'/'+fname.replace('.gz',''))
            except:
                ind = pd.DataFrame(columns=['INTERVAL_END'])
                print('No index data found on date:',curDate)
            try:
                dfOut, tempNullDict, futToSym, close2 = processTodayData(df, curDate, ind)
            except Exception as e:
                print('Error on date',curDate,':',str(type(e)),str(e))
                continue
            allData = allData.append(dfOut)
            openPx = openPx.append(dfOut.groupby(level=0)[[x+'_mid' for x in futList]].first())
            closePx = closePx.append(dfOut.groupby(level=0)[[x+'_mid' for x in futList]].last())
            closePx2[curDate.strftime('%Y-%m-%d')] = close2
            nullDict[curDate.strftime('%Y-%m-%d')] = tempNullDict
        #
        else:
            try:
                dfOut, futToSym, close2 = processTodayData(df, curDate, pd.DataFrame(columns=['INTERVAL_END']), extraOnly=True)
            except Exception as e:
                print('Error on date',curDate,':',str(type(e)),str(e))
                continue
            openPx = openPx.append(dfOut.groupby(level=0)[[x+'_mid' for x in futList]].first())
            closePx = closePx.append(dfOut.groupby(level=0)[[x+'_mid' for x in futList]].last())
            closePx2[curDate.strftime('%Y-%m-%d')] = close2
        #
        checkAndComputeRollRatio(df, curDate)
        futSymDict[curDate.strftime('%Y-%m-%d')] = futToSym
        print('Finished processing',fname)


# output close price adj, unadjusted, and adjustment factor
if outputExtra:
    #fileNameExtra = ''
    openPx.columns = [x.split('_')[0] for x in openPx.columns]
    closePx.columns = [x.split('_')[0] for x in closePx.columns]
    outputFileAndCreateDir(openPx,outPath,'US_futures_open_unadj%s.csv' % (fileNameExtra),incIndex=True)
    outputFileAndCreateDir(closePx,outPath,'US_futures_close_unadj%s.csv' % (fileNameExtra),incIndex=True)
    #
    closePx2 = pd.DataFrame(closePx2).T
    closePx2.index.name = 'Date'
    outputFileAndCreateDir(closePx2,outPath,'US_futures_close2_unadj%s.csv' % (fileNameExtra),incIndex=True)
    #
    adjFactor = adjFactor.reindex(closePx.index).fillna(1.0)
    adjFactor = adjFactor.sort_index(ascending=False)
    adjFactor = adjFactor.cumprod()
    adjFactor = adjFactor.sort_index(ascending=True)
    outputFileAndCreateDir(adjFactor,outPath,'US_futures_adj_factor%s.csv' % (fileNameExtra),incIndex=True)
    #
    closePxAdj = closePx.copy()
    for col in closePxAdj.columns:
        closePxAdj[col] *= adjFactor[col]
    outputFileAndCreateDir(closePxAdj,outPath,'US_futures_close_adj%s.csv' % (fileNameExtra),incIndex=True)

if outputIntra:
    allData = allData.sort_index()
    #allData = allData.fillna('NaN')
    # output 1 second intra second
    allData.to_csv(outPath + 'US_futures_1sec_bars%s.csv' % (fileNameExtra), na_rep='NaN')
    # output daily volatility curves
    vol = allData[[x+'_vol' for x in futList]]
    vol = vol.groupby(level=0)[[x+'_vol' for x in futList]].transform(lambda x: x.ewm(span=60,adjust=False).mean())
    vol = vol.groupby(level=1).transform(lambda x: x.rolling(window=252,min_periods=10).mean().shift(1))
    vol.columns = [x.split('_')[0] for x in vol.columns]
    for gr in vol.groupby(level=0):
        if gr[0] >= '2017-01-01' and not os.path.exists(dailyOutPath+gr[0].replace('-','')+'_60.csv'):
            temp = gr[1].reset_index().drop('Date',axis=1)
            outputFileAndCreateDir(temp, dailyOutPath, gr[0].replace('-','')+'_60.csv', incIndex=False)
    #
    vol = allData[[x+'_vol' for x in futList]]
    vol = vol.groupby(level=0)[[x+'_vol' for x in futList]].transform(lambda x: x.ewm(span=1200,adjust=False).mean())
    vol = vol.groupby(level=1).transform(lambda x: x.rolling(window=252,min_periods=10).mean().shift(1))
    vol.columns = [x.split('_')[0] for x in vol.columns]
    for gr in vol.groupby(level=0):
        if gr[0] >= '2017-01-01' and not os.path.exists(dailyOutPath+gr[0].replace('-','')+'_1200.csv'):
            temp = gr[1].reset_index().drop('Date',axis=1)
            outputFileAndCreateDir(temp, dailyOutPath, gr[0].replace('-','')+'_1200.csv', incIndex=False)
    #
    # extrapolate to 5 second snaps and output
    vol = allData[[x+'_vol' for x in futList]]
    vol['Reset'] = [(x[-1:] == '1') | (x[-1:] == '6') for x in vol.index.get_level_values(1)]
    vol['Reset'] = vol['Reset'].cumsum()
    vol = vol.groupby('Reset').cumsum()
    allDataCopy = allData.copy()
    for col in vol.columns:
        allDataCopy[col] = vol[col]
    allDataCopy[[(x[-1:] == '5') | (x[-1:] == '0') for x in allDataCopy.index.get_level_values(1)]].to_csv(outPath + 'US_futures_5sec_bars%s.csv' % (fileNameExtra), na_rep='NaN')
    #
    # extrapolate to 10 second snaps and output
    vol = allData[[x+'_vol' for x in futList]]
    vol['Reset'] = [(x[-1:] == '1') for x in vol.index.get_level_values(1)]
    vol['Reset'] = vol['Reset'].cumsum()
    vol = vol.groupby('Reset').cumsum()
    allDataCopy = allData.copy()
    for col in vol.columns:
        allDataCopy[col] = vol[col]
    allDataCopy[[(x[-1:] == '0') for x in allDataCopy.index.get_level_values(1)]].to_csv(outPath + 'US_futures_10sec_bars%s.csv' % (fileNameExtra), na_rep='NaN')
    #
    del allDataCopy
    gc.collect()
    #
    # output null data tabulation
    nullVal = pd.DataFrame(nullDict).T
    nullVal.index.name = 'Date'
    nullValPct = nullVal / 23400.0
    outputFileAndCreateDir(nullVal,outPath,'US_futures_null_value_count.csv',incIndex=True)

#####################################################
# output barra returns
# barraPath = r'/mnt/data/Common/barra_usmed/daily/'
# yearPath = barraPath + max(os.listdir(barraPath)) + '/'
# monthPath = yearPath + max(os.listdir(yearPath)) + '/'
# datePath = monthPath + max(os.listdir(monthPath)) + '/'
# f = [x for x in os.listdir(datePath) if 'DlyFacRet' in x][0]
# print('Reading Daily Factor Return Data')
# df = pd.read_csv(datePath+f,sep='|',skiprows=2,skipfooter=1)
# df.columns = [x.replace('!','').replace('%','') for x in df.columns]
# df = df[df['DataDate'] >= 20060101]
# df['DataDate'] = df['DataDate'].astype(str).apply(lambda x: x[0:4]+'-'+x[4:6]+'-'+x[6:8])
# df['Factor'] = df['Factor'].apply(lambda x: x.split('_')[1])
# df = df.pivot('DataDate', 'Factor', 'DlyReturn')
# outputFileAndCreateDir(df,r'/datadrive/Equity/Barra/USMED/FactorReturn/','BarraFactorReturn.csv',incIndex=True)


# analyze allData afterwards
# import matplotlib.pyplot as plt
# indToFut = {v:k for k,v in futToIndex.items()}
# comb['Future'] = comb['Index'].replace(indToFut)
# comb['Mid'] = comb[['Bid','Ask']].mean(axis=1)
# btic = comb.pivot('Date','Future','Mid')
# for fut in futList:
#     allData[fut+'_basis'] = btic.lookup(allData.index.get_level_values(0),[fut]*len(allData))
#     allData[fut+'_i2'] = (allData[fut+'_mid'] - (allData[fut+'_index'] + allData[fut+'_basis'])) / allData[fut+'_mid']
#
# allData.groupby(level=0)[[fut+'_i2' for fut in futList]].quantile(0.99) - allData.groupby(level=0)[[fut+'_i2' for fut in futList]].quantile(0.01)
#
# for fut in futList:
#     allData[[fut+'_i2']].plot(rot=45)
#     plt.tight_layout()
# plt.show()
#
#
# dateList = allData.index.get_level_values(0).unique()
# monthList = list(set([x[0:7] for x in dateList]))
# monthList.sort()
#
# m = monthList[31]
# for d in dateList:
#     if d[0:7] == m:
#         allData[allData.index.get_level_values(0) == d].droplevel(level=0)[[fut+'_i2' for fut in futList]].plot(rot=45,title=d,alpha=75)
#         plt.tight_layout()
# plt.show()
