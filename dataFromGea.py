#! /usr/bin/python3.5

import xmlrpc.client as xmlClient
import numpy as np
import sys
from datetime import datetime, timedelta

geaKeys = {
    'ag':'ag',
    'aom':'aom',
    'cr':'crcs',
    'ec':'ecs',
    'gc':'gcal',
    'gis':'gis',
    'lis':'lis',
    'mc':'mcs',
    'pr':'pr',
    'm1':'pcs',
    'm2':'scs',
    'tcs':'tcs',
    'bto':'bto',
    'pwfs1':'wfs',
    'pwfs2':'wfs',
    'gmoi':'wfs',
}

def geaExtractor(record, startDate, endDate):
    geaData = []
    archiver = geaKeys[record.split(':')[0]]
    ts = datetime.strptime(startDate, '%Y-%m-%d %H:%M:%S')
    te = datetime.strptime(endDate, '%Y-%m-%d %H:%M:%S')
    td = timedelta(hours=3)
    tq = timedelta(minutes=15)
    tsc = ts + td
    tec = te + td
    timeTotal = (tec - tsc).total_seconds()
    tw = tsc
    gea = xmlClient.ServerProxy('http://geasouth.cl.gemini.edu/run/ArchiveDataServer.cgi')
    geaDict = {al['name']:al['key'] for al in gea.archiver.archives()}
    sys.stdout.write('\r' + 'Progress: 0.00%')
    while tec > tw:
        tss = (tw - datetime.utcfromtimestamp(0)).total_seconds()
        tes = ((tw + tq) - datetime.utcfromtimestamp(0)).total_seconds()
        # print('mcs archiver:', geaDict['mcs'])
        geaRecord = gea.archiver.values(geaDict[archiver], [record],
                                        int(tss), int(tss - int(tss))*1000000000,
                                        int(tes), int(tes - int(tes))*1000000000,
                                        10000,
                                        0)
        geaDataAux = [[datetime.fromtimestamp(val['secs'] + val['nano']/1000000000),
                    val['value'][0]] for val in geaRecord[0]['values']]
        geaData = geaData + geaDataAux
        tw = tw + tq
        timeSpan = (tw - tsc).total_seconds()
        progress = round((timeSpan / timeTotal) * 100, 2)
        sys.stdout.write('\r' + 'Progress: ' + str(progress) + '% ')
    sys.stdout.write('\r' + 'Progress: 100.00%\n')
    print('size of data array:', len(geaData))
    print('First timestamp:', geaData[0][0])
    print('Last timestamp:', geaData[-1][0])
    return geaData

# geaData = []
# ts = datetime.strptime('2019-10-18 20:15:00', '%Y-%m-%d %H:%M:%S')
# te = datetime.strptime('2019-10-19 05:15:00', '%Y-%m-%d %H:%M:%S')
# td = timedelta(hours=3)
# tq = timedelta(minutes=15)
# tsc = ts + td
# tec = te + td
# timeTotal = (tec - tsc).total_seconds()
# tw = tsc
# gea = xmlClient.ServerProxy('http://geasouth.cl.gemini.edu/run/ArchiveDataServer.cgi')
# geaDict = {al['name']:al['key'] for al in gea.archiver.archives()}
# sys.stdout.write('\r' + 'Progress: 0.00%')
# while tec > tw:
    # tss = (tw - datetime.utcfromtimestamp(0)).total_seconds()
    # tes = ((tw + tq) - datetime.utcfromtimestamp(0)).total_seconds()
    # # print('mcs archiver:', geaDict['mcs'])
    # geaRecord = gea.archiver.values(geaDict['wfs'], ['pwfs1:dc:ttf.VALA'],
                                    # int(tss), int(tss - int(tss))*1000000000,
                                    # int(tes), int(tes - int(tes))*1000000000,
                                    # 10000,
                                    # 0)
    # geaDataAux = [[datetime.fromtimestamp(val['secs'] + val['nano']/1000000000),
                   # val['value'][0]] for val in geaRecord[0]['values']]
    # geaData = geaData + geaDataAux
    # tw = tw + tq
    # timeSpan = (tw - tsc).total_seconds()
    # progress = round((timeSpan / timeTotal) * 100, 2)
    # sys.stdout.write('\r' + 'Progress: ' + str(progress) + '% ')
# sys.stdout.write('\r' + 'Progress: 100.00%\n')
# print('size of data array:', len(geaData))
# print('First timestamp:', geaData[0][0])
# print('Last timestamp:', geaData[-1][0])

if __name__ == '__main__':
    geaDataArray = geaExtractor('mc:elCurrentPos','2019-10-18 20:15:00',
                                '2019-10-19 05:15:00')

