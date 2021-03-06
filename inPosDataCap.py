#!/usr/bin/env python3.5

import os
import argparse
import epics
import time
import h5py
import re
import numpy as np

from datetime import datetime, timedelta
from dataFromGea import geaKeys, geaExtractor

statusStr = 'not connected'
# set the the IP Address for the systems from which you need to read
# os.environ["EPICS_CA_ADDR_LIST"] = "172.17.2.36 172.17.2.32"
os.environ["EPICS_CA_ADDR_LIST"] = "172.17.2.255 172.16.71.11"

def parse_args():
    '''
    This routines parses the arguments used when running this script
    '''
    parser = argparse.ArgumentParser(
        description='Use this script to capture data from different\
        telescope subsystems')

    subparser = parser.add_subparsers(
        help='Capture options: ca for realtime, gea for GEA extraction'
    )
    # Define gea data extraction option
    parser_gea = subparser.add_parser('gea',
                                      help='Use this option to extract data\
                                      from gea')

    parser_gea.set_defaults(func=geaExtraction)

    parser_gea.add_argument('recFileG',
                            metavar='REC-FILE',
                            help='path to text file with records to be monitored')

    parser_gea.add_argument('tw',
                            nargs=2,
                            help='Time window for data extraction\
                            e.g.: "2019-11-21 02:00:00"\
                            "2019-11-22 03:00:00"')

    parser_gea.add_argument('-rn',
                            '--recname',
                            dest='rn',
                            action='store_true',
                            help='Use this option to get data from single Epics\
                            record')

    # Define real time data capture using Channel access
    parser_ca = subparser.add_parser('ca',
                                     help='Use this option to extract data\
                                     from channel access')

    parser_ca.set_defaults(func=caRealTimeCap)

    parser_ca.add_argument('recFile',
                           metavar='REC-FILE',
                           help='path to text file with records to be monitored')

    parser_ca.add_argument('-d',
                           '--days',
                           dest='rtDays',
                           default=0,
                           help='Number of days for data capture\
                           e.g.: -d 3')

    parser_ca.add_argument('-hr',
                           '--hrs',
                           dest='rtHrs',
                           default=0,
                           help='Number of hours for data capture\
                           e.g.: -hr 3')

    parser_ca.add_argument('-m',
                           '--min',
                           dest='rtMin',
                           default=1,
                           help='Number of minutes for data capture\
                           e.g.: -m 60')

    parser_ca.add_argument('-st',
                           '--starttime',
                           dest='sttime',
                           default='',
                           help='User defined start time of data capture')

    parser_ca.add_argument('-tc',
                           '--timecapture',
                           dest='tcap',
                           default='time',
                           help='Where to get timestamp, i.e "time", "native"')

    args = parser.parse_args()
    args.func(args)
    return args

def monChan(chanNames, frm):
    print(os.environ['EPICS_CA_ADDR_LIST'])
    # Initialize empty arrays
    chanList = []
    cnameList = []
    chanString = []
    # Go through EPICS channel names in array and initializes PV object
    for cn in chanNames:
        chan = epics.PV(cn, form=frm)
        time.sleep(0.25) # needed to give the PV object time to connect
        chanSt = repr(chan)
        # Check if PV object connected successfully, if not, do not add to
        # final array and print message
        if re.search(statusStr,chanSt) == None:
            cnameList.append(cn)
            chanList.append(chan)
        else:
            print("{0} not connected".format(cn))
    for ch in chanList:
        if re.search('string', ch.info):
            chanString.append(1)
        else:
            chanString.append(0)
    # Return array with PV object and channel names
    return [chanList, cnameList, chanString]

def createH5F(fname, rInfo, rDic):
    # Create and initialize .h5 file
    recHF = h5py.File(fname, 'w')
    groupsHF = [[recHF.create_group(name),name] for name in rInfo]
    for g in groupsHF:
        g[0].create_dataset('timestamp', data=rDic[g[1]][0])
        g[0].create_dataset('value', data=rDic[g[1]][1])


def caRealTimeCap(args):
    # args = parse_args() # capture the input arguments
    startTime = datetime.now() # starting time of the capture
    if not(args.sttime == ''):
        startTime = datetime.strptime(args.sttime, '%Y-%m-%d %H:%M:%S')
    while datetime.now() < startTime:
        time.sleep(1)
    startDateStr = datetime.strftime(startTime, '%Y%m%dT%H%M%S')
    startDateP = datetime.strftime(startTime, '%Y-%m-%d %H:%M:%S')
    print('Starting data capture at', startDateP)
    fileName = 'recMonCA-'+startDateStr+'.h5' # define file name
    currTime = datetime.now()
    dataCapDur = timedelta(days=int(args.rtDays), hours=int(args.rtHrs),
                        minutes=int(args.rtMin))
    with open(args.recFile, 'r') as f:
        recList = f.read().splitlines()
    recInfo = monChan(recList, args.tcap)
    # Create Dictionary for each EPICS Record data
    recDic = {name:[[],[],chan, isStr] for chan,name,isStr in np.array(recInfo).T}
    firstPass = True
    loopcnt = 0
    while ((currTime - startTime) < dataCapDur):
        startWhile = datetime.now()
        for name in recInfo[1]:
            if firstPass or not(float(recDic[name][2].timestamp)
                                == float(recDic[name][0][-1])):
                if (recDic[name][3]):
                    recDic[name][1] = np.append(recDic[name][1],
                                                np.array(recDic[name][2].value,
                                                        dtype='S32'))
                else:
                    recDic[name][1].append(recDic[name][2].value)
                recDic[name][0].append(recDic[name][2].timestamp)
        loopcnt += 1
        currTime = datetime.now()
        loopTime = (currTime - startWhile).total_seconds()
        waitTime = 0.1 - loopTime
        firstPass = False
        if loopTime > 0.1:
            print("Loop {1} took too long: {0} [s]".format(loopTime, loopcnt))
            continue
        time.sleep(waitTime)

    createH5F(fileName, recInfo[1], recDic)

def geaExtraction(args):
    # args = parse_args() # capture the input arguments
    startTime = datetime.now() # starting time of the capture
    startDateStr = datetime.strftime(startTime, '%Y%m%dT%H%M%S')
    fileName = 'recDataGea-'+startDateStr+'.h5' # define file name
    print(args.tw)
    if not(args.rn):
        with open(args.recFileG, 'r') as f:
            recList = f.read().splitlines()
    else:
        recList = [args.recFileG]
    recDic = {name:[[],[]] for name in recList}
    for recname in recList:
        print('Extracting', recname)
        recDataT = geaExtractor(recname, args.tw[0], args.tw[1])
        recData = np.array(recDataT).T
        recDic[recname][0] = recData[0]
        recDic[recname][1] = recData[1]

    createH5F(fileName, recList, recDic)


    # pass


if __name__ == '__main__':
    args = parse_args() # capture the input arguments


