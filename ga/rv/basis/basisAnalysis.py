'''
Created on Apr 23, 2020

@author: fut_strat
'''
import pandas as pd
import numpy as np
import datetime, os, smtplib, sys
import matplotlib.pyplot as plt
from email.message import EmailMessage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage

# start custom variables
inPath = r'/cc/dat/tech/trading/cts/sftprod/%s/%s/%s/log/%s.S1.main.log'
#/mnt/tradingtech/data/archive/logs/sftprod/%s/%s/%s.S1.main.log'

futList = ['ES','NQ','YM','DM','TFS']
# end custom variables

def sendEmail(msgBody):
    msg = MIMEMultipart('related')
    msg.attach(MIMEText(msgBody.replace('\n','<br>') + r'<br><img src="cid:image1">', 'html'))
    #
    fp = open(r'/home/fut_strat/US_Basis.png', 'rb')
    msgImage = MIMEImage(fp.read())
    fp.close()
    msgImage.add_header('Content-ID', '<image1>')
    msg.attach(msgImage)
    #
    msg['Subject'] = 'Basis report for trade date: %s' % (curDate.strftime('%Y-%m-%d'))
    # ------------------------------------------------------------------
    # TRANSCRIPTION GAP: original lines 33-38 were not visible in the
    # source images. Not reconstructed -- fill in from the original file.
    # Judging by the surrounding code they set msg['From'] / msg['To'],
    # opened the smtplib.SMTP connection as `s`, and sent the message,
    # since `s` is referenced by s.quit() below and is otherwise undefined.
    # ------------------------------------------------------------------
    s.quit()

def showBasisStats(df, printGraphs):
    outStr = ''
    df = df[(df['Time'] >= '09:30:00.000') & (df['Time'] <= '16:00:00.900')]
    outStr += 'Observed Basis:\n'
    for fut in futList:
        outStr += fut+' = '+str(np.round((df[fut+'_Mid'] - df[fut+'_Index']).median(),2))+'\n'
    #
    outStr += '\nBasis Quote:\n'
    for fut in futList:
        outStr += fut+' = '+str(np.round(df[(df['Time'] >= '13:30:00.000') & (df['Time'] <= '15:30:00.900')][fut+'_Basis'].median(),2))+'\n'
    #
    print(outStr)
    #
    if True:
        fig, axes = plt.subplots(2, 3, figsize=[12.8, 9.6])
        c = 0
        for fut in futList:
            df[fut+'_i2RawAvg'] = df[fut+'_i2Raw'].rolling(window=300, min_periods=200).mean()
            (10000*df[[fut+'_i2Raw', fut+'_i2RawAvg']]).plot(label=fut,rot=45,ax=axes[int(c/3), c%3])
            axes[int(c/3), c%3].axhline(0,color='k')
            c += 1
        fig.tight_layout()
        fig.savefig(r'/home/fut_strat/US_Basis.png')
    #
    sendEmail(outStr)
    #
    if printGraphs:
        for fut in futList:
            df[fut+'_i2RawAvg'] = df[fut+'_i2Raw'].rolling(window=300, min_periods=200).mean()
            (10000*df[[fut+'_i2Raw', fut+'_i2RawAvg']]).plot(label=fut,rot=45)
            plt.axhline(0,color='k')
            plt.tight_layout()
            plt.legend()
        plt.show()


if __name__ == '__main__':
    printGraphs = (len(sys.argv) <= 1) or (sys.argv[1].lower() == 'true')
    curDate = datetime.date.today()#-datetime.timedelta(days=1)
    while curDate >= datetime.date(2016,1,1):
        curDateStr = curDate.strftime('%Y%m%d')
        try:
            df = pd.read_csv(inPath % (curDateStr[0:4], curDateStr[4:6], curDateStr[6:8], curDateStr))
            break
        except:
            print('No File found on:',curDate)
            curDate = curDate - datetime.timedelta(days=1)
    #
    print('Showing basis stats for date:',curDate)
    showBasisStats(df, printGraphs)
