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
from dataFromGea import geaExtractor

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

    parser_gea.add_argument('site',
                            metavar='SITE',
                            help='Site from which to extract the data. {gs, gn}')

    parser_gea.add_argument('recFileG',
                            metavar='REC-FILE',
                            help='path to text file with records to be monitored')

    parser_gea.add_argument('tw',
                            nargs=2,
                            help='Time window for data extraction\
                            e.g.: "191121T0200"\
                            "191122T0300"')

    parser_gea.add_argument('-rn',
                            '--recname',
                            dest='rn',
                            action='store_true',
                            help='Replace file path with the name of a single\
                            EPICS record')

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
                           help='User defined start time of data capture.\
                           Format yymmddThhmm')

    args = parser.parse_args()
    args.func(args)
    return args

def monChan(chanNames):
    print('EPICS_CA_ADDR_LIST = {}'.format(os.environ['EPICS_CA_ADDR_LIST']))
    # Initialize empty dictionary
    chan_dict = {}
    # Go through EPICS channel names in array and initializes PV object
    for cn in chanNames:
        if cn in chan_dict.keys():
            print('{} repeated in the list'.format(cn))
            continue
        conn_ctr = 1 # Connection counter
        print("Connecting to {}...".format(cn))
        # Go into connection loop
        while True:
            chan = epics.PV(cn)
            time.sleep(0.25) # needed to give the PV object time to connect
            chanSt = repr(chan)
            # Check if PV object connected successfully and exit loop
            if re.search(STATUS_STR,chanSt) == None:
                chan_dict[cn] = {}
                chan_dict[cn]['chan'] = chan
                chan_dict[cn]['value'] = []
                chan_dict[cn]['timestamp'] = []
                chan_dict[cn]['string'] = False
                print("Connected!")
                break
            conn_ctr += 1
            # If channel doesn't connect after 8 retries (2 sec), break from
            # loop
            if conn_ctr == 9:
                print("{0} didn't connect after {1} tries".format(cn, conn_ctr-1))
                break

    for ch in chan_dict:
        if re.search('string', chan_dict[ch]['chan'].info):
            chan_dict[ch]['string'] = True
    # Return dictionary with connected channels
    return chan_dict

def createH5F(fname, rDic):
    # Create and initialize .h5 file
    recHF = h5py.File(fname, 'w')
    groupsHF = [[recHF.create_group(name),name] for name in rDic]
    for g in groupsHF:
        g[0].create_dataset('timestamp', data=rDic[g[1]]['timestamp'])
        g[0].create_dataset('value', data=rDic[g[1]]['value'])
    print('Data capture complete with file: {}'.format(fname))

def on_change(pvname=None, value=None, timestamp=None, **kw):
    kw['rdict'][pvname]['timestamp'].append(timestamp)
    if kw['rdict'][pvname]['string']:
        kw['rdict'][pvname]['value'] = np.append(kw['rdict'][pvname]['value'],
                                              np.array(value, dtype='S32'))
    else:
        kw['rdict'][pvname]['value'].append(value)

def on_press_thread(run_flag):
    key_press = input()
    if key_press == 'x':
        print('Stopping capture')
        run_flag[0] = False

def caRealTimeCap(args):
    run_flag = [True]
    startTime = datetime.now() # starting time of the capture
    if not(args.sttime == ''):
        startTime = datetime.strptime(args.sttime, '%y%m%dT%H%M%S')
    while datetime.now() < startTime:
        time.sleep(1)
    startDateStr = datetime.strftime(startTime, '%Y%m%dT%H%M%S')
    startDateP = datetime.strftime(startTime, '%Y-%m-%d %H:%M:%S')
    print('Starting data capture at', startDateP)
    fileName = 'recMonCA-'+startDateStr+'.h5' # define file name
    dataCapDur = timedelta(days=int(args.rtDays), hours=int(args.rtHrs),
                        minutes=int(args.rtMin))
    # Read file with record names, ignore comments
    with open(args.recFile, 'r') as f:
        recList = [l.split('#')[0].strip()
                   for l in f.read().splitlines() if l.split('#')[0]]
    # Connect to the EPICS channels
    rec_dict = monChan(recList)
    # If no channel connected, abort the program
    if not(rec_dict):
        sys.exit('No channels connected, aborting')
    # Start the monitors
    for cname in rec_dict:
        print("Starting monitor for {}".format(cname))
        # Add the on_change routine and pass the record dictionary as argument
        rec_dict[cname]['chan'].add_callback(on_change, rdict=rec_dict)
    # Start "abort data capture" thread
    _thread.start_new_thread(on_press_thread, (run_flag,))
    print('To stop capture now: x + [Enter]')
    wait_time = 0.005 # Are 5 ms enough for while loop not to hog the CPU?
    init_time = datetime.now()
    curr_time = init_time
    # Enter loop and wait for the callbacks to fill up the record dict
    while (((curr_time - init_time) < dataCapDur) and run_flag[0]):
        time.sleep(wait_time)
        curr_time = datetime.now()
    # Once data capture ends, stop the monitors. I don't want to continue
    # capturing data at this point
    print("Stopping monitors")
    for cname in rec_dict:
        rec_dict[cname]['chan'].clear_callbacks()
    print("Done!")

    createH5F(fileName, rec_dict)

def geaExtraction(args):
    startTime = datetime.now() # starting time of the capture
    startDateStr = datetime.strftime(startTime, '%Y%m%dT%H%M%S')
    fileName = 'recDataGea-'+startDateStr+'.h5' # define file name
    print(args.tw)
    if not(args.rn):
        with open(args.recFileG, 'r') as f:
            recList = [l.split('#')[0].strip()
                       for l in f.read().splitlines() if l.split('#')[0]]
    else:
        recList = [args.recFileG]
        recDic = {name:{'value':[], 'timestamp':[]} for name in recList}
    chan_count = 0
    # Format the start and end date as datetime objects, with the specified
    # format
    try:
        startDate = datetime.strptime(args.tw[0], '%y%m%dT%H%M')
        endDate = datetime.strptime(args.tw[1], '%y%m%dT%H%M')
    except ValueError as err:
        sys.exit("ValueError timewindow: {}".format(err))
    for recname in recList:
        recDataT = geaExtractor(recname, startDate, endDate, args.site)
        if not(recDataT):
            continue
        recData = np.array(recDataT).T
        recDic[recname]['timestamp'] = recData[0]
        recDic[recname]['value'] = recData[1]
        chan_count += 1

    if not(chan_count):
        sys.exit('No channels were extracted')
    createH5F(fileName, recDic)

if __name__ == '__main__':
    args = parse_args() # capture the input arguments


