#!/usr/bin/env python3.5

import os
import argparse
from datetime import datetime, timedelta
import epics
import time
import h5py
import re
import numpy as np

statusStr = 'not connected'
# set the the IP Address for the systems from which you need to read
# os.environ["EPICS_CA_ADDR_LIST"] = "172.17.2.36 172.17.2.32"
os.environ["EPICS_CA_ADDR_LIST"] = "172.17.2.255"

def parse_args():
    '''
    This routines parses the arguments used when running this script
    '''
    parser = argparse.ArgumentParser(
        description='Use this script to capture Time data from different\
        telescope subsystems')

    parser.add_argument('recFile',
                        metavar='REC-FILE',
                        help='path to text file with records to be monitored')

    parser.add_argument('-d',
                        '--days',
                        dest='rtDays',
                        default=0,
                        help='Number of days for data capture\
                        e.g.: -d 3')

    parser.add_argument('-hr',
                        '--hrs',
                        dest='rtHrs',
                        default=0,
                        help='Number of hours for data capture\
                        e.g.: -h 3')

    parser.add_argument('-m',
                        '--min',
                        dest='rtMin',
                        default=1,
                        help='Number of minutes for data capture\
                        e.g.: -m 60')

    parser.add_argument('-st',
                        '--status',
                        dest='rtStatus',
                        action='store_true',
                        help='Use this option to get the status of channels')

    args = parser.parse_args()
    return args

def monChan(chanNames):
    print(os.environ['EPICS_CA_ADDR_LIST'])
    # Initialize empty arrays
    chanList = []
    cnameList = []
    chanString = []
    # Go through EPICS channel names in array and initializes PV object
    for cn in chanNames:
        chan = epics.PV(cn)
        time.sleep(0.2) # needed to give the PV object time to connect
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

if __name__ == '__main__':
    args = parse_args() # capture the input arguments
    prevTS = 0
    startTime = datetime.now() # starting time of the capture
    startDateStr = datetime.strftime(startTime, '%Y%m%dT%H%M%S')
    fileName = 'recMon-'+startDateStr+'.h5' # define file name
    currTime = datetime.now()
    dataCapDur = timedelta(days=int(args.rtDays), hours=int(args.rtHrs),
                        minutes=int(args.rtMin))
    with open(args.recFile, 'r') as f:
        recList = f.read().splitlines()
    recInfo = monChan(recList)
    # Create Dictionary for each EPICS Record data
    recDic = {name:[[],[],chan, isStr] for chan,name,isStr in np.array(recInfo).T}
    firstPass = True
    loopcnt = 0
    while ((currTime - startTime) < dataCapDur):
        startWhile = datetime.now()
        for name in recInfo[1]:
            # if firstPass:
            if (recDic[name][3]):
                recDic[name][1] = np.append(recDic[name][1],
                                            np.array(recDic[name][2].value,
                                                     dtype='S10'))
            else:
                recDic[name][1].append(recDic[name][2].value)
                # continue
            # if not(float(recDic[name][2].timestamp) == float(recDic[name][0][-1])):
                # recDic[name][1].append(recDic[name][2].value)
                # recDic[name][0].append(recDic[name][2].timestamp)
            # else:
                # print('Identical Timestamps for {2}: {0} = {1}'.format(recDic[name][2].timestamp, recDic[name][0][-1], name))
        loopcnt += 1
        currTime = datetime.now()
        loopTime = (currTime - startWhile).total_seconds()
        waitTime = 0.1 - loopTime
        firstPass = False
        if loopTime > 0.1:
            print("Loop {1} took too long: {0} [s]".format(loopTime, loopcnt))
            continue
        time.sleep(waitTime)
    # Create and initialize .h5 file
    recHF = h5py.File(fileName, 'w')
    groupsHF = [[recHF.create_group(name),name] for name in recInfo[1]]
    for g in groupsHF:
        g[0].create_dataset('timestamp', data=recDic[g[1]][0])
        g[0].create_dataset('value', data=recDic[g[1]][1])




