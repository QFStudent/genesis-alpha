'''
Created on Feb 4, 2020

module to compute pnl and statistics from Alex's sim output

@author: john.casale
'''
import pandas as pd
import numpy as np
import datetime, os

###### start custom variables
execPath = r'/mnt/tradingtech/data/archive/logs/sftprod/%s/%s/%s.execs.csv' #r'/home/tradingtech/dev/cts/logs/FinalBt/%s/%s/%s.execs.csv'
eodPath = r'/mnt/tradingtech/data/archive/logs/sftprod/%s/%s/%s.eod.posn.csv' #r'/home/tradingtech/dev/cts/logs/FinalBt/%s/%s/%s.eod.posn.csv'
startDate = datetime.date(2021,5,21)#(2019,10,11)
endDate = datetime.date(2021,7,16)#(2017,1,3)#(2019,10,14)

futList = ['ES','NQ','YM','DM','TFS']
multDict = {'ES': 50, 'NQ': 20, 'YM': 5, 'DM': 100, 'TFS': 50}
multSeries = pd.Series(multDict)
futToRoot = {'ES': 'ES', 'NQ': 'NQ', 'YM': 'DM', 'DM': 'FA', 'TFS': 'RTY'}
rootToFut = {v:k for k,v in futToRoot.items()}

closePriceAdjFile = '/mnt/nfs1/user/sft/sftdata/CurrentFiles/US_futures_close_adj.csv'
closePriceUnadjFile = '/mnt/nfs1/user/sft/sftdata/CurrentFiles/US_futures_close_unadj.csv'

tCostFixed = 2.5
###### end custom variables

def checkIfTodayIsRollDay(cur_date):
    if cur_date.month % 3 != 0:
        return False
    if cur_date.day < 14 or cur_date.day > 20:
        return False
    return cur_date.weekday() == 3

def matchEodPosToFile(eodPos,curDate):
    alexPos = pd.read_csv(eodPath % (curDate.strftime('%Y'), curDate.strftime('%Y%m%d'), curDate.strftime('%Y%m%d')))
    alexPos = alexPos[alexPos['book'] == 'main']
    alexPos['Future'] = alexPos['ticker'].apply(lambda x: rootToFut[x[0:len(x)-2]])
    comp = pd.concat([eodPos,alexPos.set_index('Future')['posn']],axis=1).fillna(0)
    comp.columns = ['John','Alex']
    return (comp['John'] - comp['Alex']).abs().sum() < 1

close_px_adj = pd.read_csv(closePriceAdjFile).set_index('Date')
close_px_unadj = pd.read_csv(closePriceUnadjFile).set_index('Date')
adjFactor = close_px_adj / close_px_unadj
retAdj = (close_px_adj / close_px_adj.shift(1).ffill() - 1.0)

pnlIntra, pnlExtra, tCost, gmv, t0, var = {}, {}, {}, {}, {}, {}
eodPos = pd.Series([79, -38, -66, -18, 81], index=futList).fillna(0)
prevDate = startDate - datetime.timedelta(days=1)
curDate = startDate
while(curDate <= endDate):
    if not os.path.exists(execPath % (curDate.strftime('%Y'), curDate.strftime('%Y%m%d'), curDate.strftime('%Y%m%d'))):
        #print('No exec found on:',curDate)
        curDate += datetime.timedelta(days=1)
        continue
    #
    df = pd.read_csv(execPath % (curDate.strftime('%Y'), curDate.strftime('%Y%m%d'), curDate.strftime('%Y%m%d')))
    df['Future'] = df['symbol'].apply(lambda x: rootToFut[x[:-2]])
    df['PriceQty'] = df['price'] * df['qty']
    #if checkIfTodayIsRollDay(curDate):
    #    df = df[df['time'] < '15:55:00']
    #
    buyQty = df[df['side'] == 'BUY'].groupby('Future')['qty'].sum().reindex(futList).fillna(0)
    sellQty = df[df['side'] == 'SELL'].groupby('Future')['qty'].sum().reindex(futList).fillna(0)
    buyPx = (df[df['side'] == 'BUY'].groupby('Future')['PriceQty'].sum().reindex(futList) / buyQty).fillna(0)
    sellPx = (df[df['side'] == 'SELL'].groupby('Future')['PriceQty'].sum().reindex(futList) / sellQty).fillna(0)
    # t0, gmv, pnl stats
    curDateStr = curDate.strftime('%Y-%m-%d')
    todayClose, yestClose = close_px_unadj.loc[curDateStr], close_px_unadj.loc[prevDate.strftime('%Y-%m-%d')]
    todayAdj, yestAdj = adjFactor.loc[curDateStr], adjFactor.loc[prevDate.strftime('%Y-%m-%d')]
    t0[curDateStr] = (buyQty + sellQty) * todayClose * multSeries #(buyQty * buyPx * multSeries) + (sellQty * sellPx * multSeries)
    gmv[curDateStr] = (eodPos + buyQty - sellQty).abs() * todayClose * multSeries
    pnlIntra[curDateStr] = ((todayClose - buyPx) * buyQty + (sellPx - todayClose) * sellQty) * multSeries
    pnlExtra[curDateStr] = (todayClose - (yestClose * (yestAdj / todayAdj))) * eodPos * multSeries
    tCost[curDateStr] = (buyQty + sellQty) * tCostFixed
    #
    eodPos = (eodPos + buyQty - sellQty)
    # var stats
    temp_var = retAdj[futList].dot((eodPos * todayClose * multSeries))
    dateInd = temp_var.index.get_loc(curDateStr)
    var[curDate] = temp_var.iloc[max(0,dateInd-504):dateInd].nsmallest(5)[-1]
    #
    if not matchEodPosToFile(eodPos,curDate):
        print('Error matching eod file on',curDate)
        break
    prevDate = curDate
    curDate += datetime.timedelta(days=1)

t0 = pd.DataFrame(t0).T
gmv = pd.DataFrame(gmv).T
pnlIntra = pd.DataFrame(pnlIntra).T
pnlExtra = pd.DataFrame(pnlExtra).T
tCost = pd.DataFrame(tCost).T
pnl = pnlIntra + pnlExtra - tCost
var = pd.Series(var)
var.index = [x.strftime('%Y-%m-%d') for x in var.index]

t0['Total'] = t0.sum(axis=1)
gmv['Total'] = gmv.sum(axis=1)
pnlIntra['Total'] = pnlIntra.sum(axis=1)
pnlExtra['Total'] = pnlExtra.sum(axis=1)
tCost['Total'] = tCost.sum(axis=1)
pnl['Total'] = pnl.sum(axis=1)

# compute stats
# for year in range(2017,2021):
#     yearStr = str(year)
#     pnl = pJohn[pJohn.index.str.contains(yearStr)]#[pJohn['vix'] >= 15]#
#     pnlIntra = pIJohn[pIJohn.index.str.contains(yearStr)]#[pIJohn['vix'] >= 15]#
#     pnlExtra = pEJohn[pEJohn.index.str.contains(yearStr)]#[pEJohn['vix'] >= 15]#
#     t0 = tJohn[tJohn.index.str.contains(yearStr)]#[tJohn['vix'] >= 15]#
#     gmv = gJohn[gJohn.index.str.contains(yearStr)]#[gJohn['vix'] >= 15]#
#     tCost = tCJohn[tCJohn.index.str.contains(yearStr)]#[tCJohn['vix'] >= 15]#
#     var = vJohn[vJohn.index.str.contains(yearStr)]#[vJohn['vix'] >= 15]#
pnl['DD'] = pnl['Total'].cumsum() - pnl['Total'].cumsum().cummax().apply(lambda x: max(0,x))
pnl['DDDays'] = (pnl['DD'] == 0).cumsum()
pnl['DDDays'] = pnl.groupby('DDDays')['DDDays'].cumcount() + 1
try:
    sortino1 = pnl['Total'].mean()/np.sqrt(pnl['Total'].where(pnl['Total']<0,0).pow(2).sum() / len(pnl['Total']))*np.sqrt(252)
    sortino2 = pnl['Total'].mean()/pnl[pnl['Total']<0]['Total'].std()*np.sqrt(252)
except:
    sortino1 = np.NaN
    sortino2 = np.NaN

statPairs = [
            ('AvgPnl',pnl['Total'].mean()),\
            ('StdPnl',pnl['Total'].std()),\
            ('Sharpe',pnl['Total'].mean()/pnl['Total'].std()*np.sqrt(252)),\
            ('Sortino1',sortino1),\
            ('Sortino2',sortino2),\
            ('TotalPnl',pnl['Total'].sum()),\
            ('TotalIntraPnl',pnlIntra['Total'].sum()),\
            ('TotalExtraPnl',pnlExtra['Total'].sum()),\
            ('TotalTCost',tCost['Total'].sum()),\
            ('AvgT0',t0['Total'].mean()),\
            ('AvgGmv',gmv['Total'].mean()),\
            ('Pnl/T0',pnl['Total'].sum()/t0['Total'].sum()),\
            ('Pnl/Gmv',pnl['Total'].sum()/gmv['Total'].sum()*252.0),\
            ('AvgVar',np.mean(var)),\
            ('MaxVar',min(var)),\
            ('MaxDraw',pnl['DD'].min()),\
            ('MaxDraw/AvgGmv',pnl['DD'].min()/gmv['Total'].mean()),\
            ('MaxDrawDays',pnl['DDDays'].max()),\
            ('HP',2.0*gmv['Total'].sum()/t0['Total'].sum()),\
            ('NumDays',len(pnl))
            ]
stats = pd.DataFrame(statPairs,columns=['Stat','Value']).set_index('Stat')
print(stats)


pJohn = bt.getBookLevelPnlAll()#pnl.copy()#
pIJohn = bt.getBookLevelPnlIntra()#pnlIntra.copy()#
pEJohn = bt.getBookLevelPnlExtra()#pnlExtra.copy()#
gJohn = bt.getBookLevelGmv()#gmv.copy()#
tJohn = bt.getBookLevelT0()#t0.copy()#
tCJohn = bt.book.bookTradeList.tCost#tCost.copy()#
vJohn = pd.Series(bt.book.bookTradeList.varDict)#var.copy()#

for col in pnl.columns:
    try:
        print(col,pnl[col].corr(pJohn[col]))
    except:
        continue

for col in gmv.columns:
    try:
        print(col,gmv[col].corr(gJohn[col]))
    except:
        continue

for col in t0.columns:
    try:
        print(col,t0[col].corr(tJohn[col]))
    except:
        continue



futToRoot = {'ES': 'ES', 'NQ': 'NQ', 'YM': 'DM', 'DM': 'FA', 'TFS': 'RTY'}
a = pd.read_csv(r'/home/tradingtech/dev/cts/logs/%s/%s/%s.execs.csv' % (startDate.strftime('%Y'),startDate.strftime('%Y%m%d'),startDate.strftime('%Y%m%d')))
a = pd.read_csv(r'/home/tradingtech/dev/cts/logs/test/2019/20191011/20191011.execs.csv')
a['time'] = a['time'].apply(lambda x: x.split('.')[0])
b = pd.DataFrame(bt.book.bookTradeList.rawTradeList, columns=['Date','time','Future','Size','Price'])
b['symbol'] = b['Future'].apply(lambda x: futToRoot[x]+'Z9')
c = pd.merge(a[['time','symbol','side','qty','price']],b[['time','symbol','Size','Price']],on=['time','symbol'],how='outer')
c = c.sort_values('time')

d = pd.read_csv(r'/home/fut_strat/media/disk/_cygdrive_C_JOHNCA1_temp/Ju1_old.log')#john.log')
e = pd.read_csv(r'/home/fut_strat/media/disk/_cygdrive_C_JOHNCA1_temp/Ju1_opt.log')#alex.log')
for col in d.columns:
    try:
        print(col,d[col].corr(e[col]))
    except:
        continue

