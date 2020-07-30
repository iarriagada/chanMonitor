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
    '''
    This method is used to extract data from GEA.
    It uses the name of the record to determine the archiver to be used
    '''
    geaData = []
    # Get IOC top from the name of the record, then look in the GEA keyword
    # dictionary
    archiver = geaKeys[record.split(':')[0]]
    # Format the start and end date as datetime objects, with the specified
    # format
    ts = datetime.strptime(startDate, '%Y-%m-%d %H:%M:%S')
    te = datetime.strptime(endDate, '%Y-%m-%d %H:%M:%S')
    # Define the timedelta to switch date to UTC and define a delta of 15
    # minutes to get data from GEA
    td = timedelta(hours=3)
    tq = timedelta(minutes=15)
    tsc = ts + td
    tec = te + td
    # Determine the total amount of time to be extracted
    timeTotal = (tec - tsc).total_seconds()
    # Define the start of the time window
    tw = tsc
    # Connect to GEA xml server
    gea = xmlClient.ServerProxy(
        'http://geasouth.cl.gemini.edu/run/ArchiveDataServer.cgi'
    )
    # Get the list of archivers available
    geaDict = {al['name']:al['key'] for al in gea.archiver.archives()}
    sys.stdout.write('\r' + 'Progress: 0.00%')
    while tec > tw:
        # Get a window of 15 minutes of data from GEA
        tss = (tw - datetime.utcfromtimestamp(0)).total_seconds()
        tes = ((tw + tq) - datetime.utcfromtimestamp(0)).total_seconds()
        geaRecord = gea.archiver.values(geaDict[archiver], [record],
                                        int(tss), int(tss - int(tss))*1000000000,
                                        int(tes), int(tes - int(tes))*1000000000,
                                        10000,
                                        0)
        # print(geaRecord)
        # Auxiliary array for the 15 min of data
        # geaDataAux = [[datetime.fromtimestamp(val['secs'] + val['nano']/1000000000),
                    # val['value'][0]] for val in geaRecord[0]['values']]
        if not(geaRecord[0]['type']):
            geaDataAux = [[val['secs'] + val['nano']/1000000000,
                        np.array(val['value'][0], dtype='S32')] for val in geaRecord[0]['values']]
        else:
            geaDataAux = [[val['secs'] + val['nano']/1000000000,
                           val['value'][0]] for val in geaRecord[0]['values']]
        # Final array to be returned
        geaData = geaData + geaDataAux
        # Advance time window
        tw = tw + tq
        # Calculate and display progress
        timeSpan = (tw - tsc).total_seconds()
        progress = round((timeSpan / timeTotal) * 100, 2)
        sys.stdout.write('\r' + 'Progress: ' + str(progress) + '% ')
    sys.stdout.write('\r' + 'Progress: 100.00%\n')
    print('size of data array:', len(geaData))
    print('First timestamp:', datetime.fromtimestamp(geaData[0][0]))
    print('Last timestamp:', datetime.fromtimestamp(geaData[-1][0]))
    return geaData

if __name__ == '__main__':
    geaDataArray = geaExtractor('mc:azCurrentPos','2019-10-07 11:00:00',
                                '2019-10-07 13:00:00')

