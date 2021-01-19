#!/usr/bin/env python3

import os
import sys
import argparse
import epics
import time
import h5py
import re
import numpy as np
import _thread

from datetime import datetime, timedelta
from dataFromGea import geaKeys, geaExtractor

STATUS_STR = 'not connected'

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
    print('EPICS_CA_ADDR_LIST = {}'.format(os.environ['EPICS_CA_ADDR_LIST']))
    # Initialize empty arrays
    chanList = []
    cnameList = []
    chanString = []
    # Go through EPICS channel names in array and initializes PV object
    for cn in chanNames:
        conn_ctr = 1 # Connection counter
        print("Connecting to {}...".format(cn))
        # Go into connection loop
        while True:
            chan = epics.PV(cn, form=frm)
            time.sleep(0.25) # needed to give the PV object time to connect
            chanSt = repr(chan)
            # Check if PV object connected successfully and exit loop
            if re.search(STATUS_STR,chanSt) == None:
                cnameList.append(cn)
                chanList.append(chan)
                print("Connected!")
                break
            conn_ctr += 1
            # If channel doesn't connect after 8 retries (2 sec), break from
            # loop
            if conn_ctr == 9:
                print("{0} didn't connect after {1} tries".format(cn, conn_ctr-1))
                break

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

def on_press_thread(run_flag):
    key_press = input()
    if key_press == 'x':
        print('Stopping capture')
        run_flag[0] = False

def caRealTimeCap(args):
    run_flag = [True]
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
    # Read file with record names, ignore comments
    with open(args.recFile, 'r') as f:
        recList = [l.split('#')[0].strip()
                   for l in f.read().splitlines() if l.split('#')[0]]
    # Connect to the EPICS channels
    recInfo = monChan(recList, args.tcap)
    # If no channel connected, abort the program
    if not(recInfo[0]):
        sys.exit('No channels connected, aborting')
    # Create Dictionary for each EPICS Record data
    recDic = {name:[[],[],chan, isStr] for chan,name,isStr in np.array(recInfo).T}
    firstPass = True # First data value capture flag
    loopcnt = 0
    # Start "abort data capture" thread
    _thread.start_new_thread(on_press_thread, (run_flag,))
    print('To stop capture now: x + [Enter]')
    while (((currTime - startTime) < dataCapDur) and run_flag[0]):
        startWhile = datetime.now()
        for name in recInfo[1]:
            if firstPass or not(float(recDic[name][2].timestamp)
                                == float(recDic[name][0][-1])):
                # Check if value is a string and store it in array with proper
                # data type
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
        # Set loop time to 100 ms minus the time it takes the loop to process
        waitTime = 0.1 - loopTime
        firstPass = False
        if loopTime > 0.1:
            print("Loop {1} took too long: {0} [s]".format(loopTime, loopcnt))
            continue
        time.sleep(waitTime)

    createH5F(fileName, recInfo[1], recDic)
    print('Data capture complete with file: {}'.format(fileName))

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

if __name__ == '__main__':
    args = parse_args() # capture the input arguments


