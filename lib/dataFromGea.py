#! /usr/bin/python3

import xmlrpc.client as xmlClient
import numpy as np
import sys
from datetime import datetime, timedelta, timezone, tzinfo

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
    'm1p':'pcsact',
    'm2':'scs',
    'tcs':'tcs',
    'bto':'bto',
    'pwfs1':'wfs',
    'pwfs2':'wfs',
    'gmoi':'wfs',
    'oiwfs':'wfs',
    'gm':'gmos',
    'ta':'ta',
    'ws':'gws'
}

gea_xml_server = {
    'gs':'http://geasouth.cl.gemini.edu/run/ArchiveDataServer.cgi',
    'gn':'http://geanorth.hi.gemini.edu/run/ArchiveDataServer.cgi'
}

time_zone = {
    'gs':timezone(timedelta(hours=-4), 'CLT'),
    'gs_st':timezone(timedelta(hours=-3), 'CLST'),
    'gn':timezone(timedelta(hours=-10), 'HST'),
    'utc':timezone.utc,
    'no':None
}

def utc2site_time(timestamp, site):
    time_dt_site = datetime.fromtimestamp(timestamp, tz=time_zone[site])
    return time_dt_site

def geaExtractor(record, tsc, tec, site='gs'):
    '''
    This method is used to extract data from GEA.
    It uses the name of the record to determine the archiver to be used
    '''
    geaData = {'timestamp':[], 'value':[]}
    # Get IOC top from the name of the record, then look in the GEA keyword
    # dictionary
    top = record.split(':')
    archiver = geaKeys[top[0]]
    # TODO: This is just for the especial case of a weird PCS archiver, it needs
    # to be fixed
    if top[0] == 'm1p':
        record = 'm1:' + ':'.join(top[1:])
    # Define the number of data points to be harvested in each iteration
    data_points = 10000
    # Determine the total amount of time to be extracted
    timeTotal = (tec - tsc).total_seconds()
    # Define the start and end of the time window, in seconds and nano seconds
    # to be handled by the archiver
    ts = tsc.timestamp()
    te = tec.timestamp()
    tss = int(ts)
    tsn = int((ts - tss) * 1000000000)
    tes = int(te)
    ten = int((te - tes) * 1000000000)
    # Connect to GEA xml server
    gea = xmlClient.ServerProxy(gea_xml_server[site.split('_')[0]])
    # Get the list of archivers available
    geaDict = {al['name']:al['key'] for al in gea.archiver.archives()}
    # Check to see if the record is being archived, if not, return False
    try:
        if not(gea.archiver.names(geaDict[archiver], record)):
            print('Channels {} is not being archived'.format(record))
            return False
    except Exception as error:
        print("\nConnection with server {} failed".format(gea_xml_server[site]))
        print("with error: {}".format(error))
        return False

    print('Extracting', record)
    sys.stdout.write('\r' + 'Progress: 0.00%')
    try:
        # Start a while loop to extract data. The loop will run until the number
        # of data points extracted is less than 10000
        while data_points >= 10000:
            # Get the data from the archiver. The archiver will extract data
            # points between the timestamps ts and te, or until it extracts
            # 10000 points, whatever happens first.
            geaRecord = gea.archiver.values(geaDict[archiver], [record],
                                            tss,
                                            tsn,
                                            tes,
                                            ten,
                                            100000,
                                            0)
            # Check the amount of data points harvested
            data_points = len(geaRecord[0]['values'])
            # Define a new start to the time window to harvest
            tss = geaRecord[0]['values'][-1]['secs']
            tsn = geaRecord[0]['values'][-1]['nano']
            ts = tss + tsn/1000000000
            timestamps = [val['secs'] + val['nano']/1000000000\
                            for val in geaRecord[0]['values']]
            geaData['timestamp'] += timestamps[:-1]
            # Check to see if the record value is a string, and create a numpy
            # array with the proper type
            if not(geaRecord[0]['type']):
                values = np.array([val['value'][0]
                                   for val in geaRecord[0]['values']],
                                  dtype='S32')
                geaData['value'] = np.concatenate((geaData['value'], values[:-1]))
            else:
                values = [val['value']
                          if (len(val['value'])>1)
                          else val['value'][0]
                          for val in geaRecord[0]['values']]
                geaData['value'] += values[:-1]
            # Calculate and display progress
            timeSpan = 1 - (te - ts)/timeTotal
            progress = round(timeSpan * 100, 2)
            sys.stdout.write('\r' + 'Progress: ' + str(progress) + '% ')
        sys.stdout.write('\r' + 'Progress: 100.00%\n')
        print('size of data array:', len(geaData['timestamp']))
        print('First timestamp:',
            utc2site_time(geaData['timestamp'][0], site))
        print('Last timestamp:',
            utc2site_time(geaData['timestamp'][-1], site))
        return geaData
    except Exception as error:
        print("\nChannel {} extraction failed".format(record))
        print("with error: {}".format(error))
        return False

if __name__ == '__main__':
    ts = '210107T0100'
    te = '210107T0115'
    tsd = datetime.strptime(ts, '%y%m%dT%H%M')
    ted = datetime.strptime(te, '%y%m%dT%H%M')

    geaDataArray = geaExtractor('mc:azCurrentPos',tsd,
                                ted, 'gn')

