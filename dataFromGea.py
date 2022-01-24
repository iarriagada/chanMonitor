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
    'm2':'scs',
    'tcs':'tcs',
    'bto':'bto',
    'pwfs1':'wfs',
    'pwfs2':'wfs',
    'gmoi':'wfs',
    'oiwfs':'wfs',
    'gm':'gmos',
    'ta':'temporary',
    'ws':'gws'
}

gea_xml_server = {
    'gs':'http://geasouth.cl.gemini.edu/run/ArchiveDataServer.cgi',
    'gn':'http://geanorth.hi.gemini.edu/run/ArchiveDataServer.cgi'
}

time_zone = {
    'gs':timezone(timedelta(hours=-4), 'CLT'),
    'gn':timezone(timedelta(hours=-10), 'HST'),
    'utc':timezone.utc
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
    archiver = geaKeys[record.split(':')[0]]
    # Define the timedelta to switch date to UTC and define a delta of 5
    # minutes to get data from GEA
    tq = timedelta(minutes=5)
    # Determine the total amount of time to be extracted
    timeTotal = (tec - tsc).total_seconds()
    # Define the start of the time window
    tw = tsc
    # Connect to GEA xml server
    gea = xmlClient.ServerProxy(gea_xml_server[site])
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
        while tec > tw:
            # Get a window of 5 minutes of data from GEA. For some reason
            # umknown to me, you can only extract 15 min worth of data at a
            # time
            tss = (tw).timestamp()
            tes = (tw + tq).timestamp()
            geaRecord = gea.archiver.values(geaDict[archiver], [record],
                                            int(tss),
                                            int(tss - int(tss))*1000000000,
                                            int(tes),
                                            int(tes - int(tes))*1000000000,
                                            10000,
                                            0)
            # Auxiliary array for the 15 min of data
            # Check to see if the record value is a string, and create a numpy
            # array with the proper type
            timestamps = [val['secs'] + val['nano']/1000000000\
                            for val in geaRecord[0]['values']]

            geaData['timestamp'] += timestamps
            if not(geaRecord[0]['type']):
                values = np.array([val['value'][0]
                                   for val in geaRecord[0]['values']],
                                  dtype='S32')
                geaData['value'] = np.concatenate((geaData['value'], values))
            else:
                values = [val['value']
                          if (len(val['value'])>1)
                          else val['value'][0]
                          for val in geaRecord[0]['values']]
                geaData['value'] += values
            # Advance time window
            tw = tw + tq
            # Calculate and display progress
            timeSpan = (tw - tsc).total_seconds()
            progress = round((timeSpan / timeTotal) * 100, 2)
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

