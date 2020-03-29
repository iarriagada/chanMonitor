#!/usr/bin/env python3.5

from datetime import datetime, timedelta
from plotEpicsData import DataAx, DataAxesPlotter

import plotEpicsData as ped
import collections
import argparse
import numpy as np
import h5py
import re


def parse_args():
    '''
    This routines parses the arguments used when running this script
    '''
    parser = argparse.ArgumentParser(
        description='Use this script to plot captured Time data from different\
        telescope subsystems from a pickled file')

    parser.add_argument('hdf5File',
                        metavar='H5-FILE',
                        nargs='*',
                        help='path to hdf5 file with data to be plotted')

    # parser.add_argument('-d',
                        # '--days',
                        # dest='rtDays',
                        # default=0,
                        # help='Number of days for data capture\
                        # e.g.: -d 3')

    # parser.add_argument('-hr',
                        # '--hrs',
                        # dest='rtHrs',
                        # default=0,
                        # help='Number of hours for data capture\
                        # e.g.: -h 3')

    # parser.add_argument('-m',
                        # '--min',
                        # dest='rtMin',
                        # default=1,
                        # help='Number of minutes for data capture\
                        # e.g.: -m 60')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    # Create empty arrays for TCS in position time stamp and value
    recGroups = []
    print(args.hdf5File)
    # Read h5 file
    posFile = [h5py.File(datafile, 'r') for datafile in args.hdf5File]
    # Extract record names and create an array with group object and name of
    # group
    for seqData in posFile:
        recGroups += [[seqData.get(recname),recname] for recname in seqData]
    print(np.array(recGroups).T[1])
    # Create dictionary with time and process value data
    print('Generating data dictionary')
    recData = {name:[np.array(g.get('timestamp')),np.array(g.get('value'))]\
               for g,name in recGroups}
    # Reformat timestamp data
    print('Reformating timestamps')
    for n in recData:
        print(n)
        if (re.search('ips:', n)):
            print('Fixing 37sec offset for ', n)
            recTime = [datetime.fromtimestamp(ts+37) for ts in recData[n][0]]
        else:
            recTime = [datetime.fromtimestamp(ts) for ts in recData[n][0]]
        timeStampS = [datetime.strftime(rt, '%m/%d/%Y %H:%M:%S.%f')\
                      for rt in recTime]
        timeStamp = [datetime.strptime(ts, '%m/%d/%Y %H:%M:%S.%f')\
                     for ts in timeStampS]
        recData[n][0] = timeStamp

    print('Performing calculations on data')
    azdP = ped.diffData(recData['mc:azCurrentPos'])
    azdTD = ped.diffData(recData['tcs:demandAz'])
    azdD = ped.diffData(recData['mc:azDemandPos'])
    # azErr = diffData(recData['mc:azPosError'])
    eldP = ped.diffData(recData['mc:elCurrentPos'])
    eldTD = ped.diffData(recData['tcs:demandEl'])
    eldD = ped.diffData(recData['mc:elDemandPos'])
    # elErr = diffData(recData['mc:elPosError'])
    # crdP = diffData(recData['cr:crCurrentPos'])
    # crdV = derivData(recData['cr:crCurrentVel'])
    crdErr = ped.derivData(recData['cr:crPosError'])
    avgCRdErr = ped.avgData(crdErr,3)
    # avgCRdV = crdV

    # inFalsePos = inPosDur(tcsipData)
    inFalsePos = ped.inPosDur(recData['tcs:inPosCombine'])
    # print(inFalsePos)
    ifpMCS = ped.inPosDur(recData['mc:inPosition'])
    ifpCRCS = ped.inPosDur(recData['cr:inPosition'])
    ifpP1 = ped.inPosDur(recData['ag:p1:inPosition.VAL'])
    ifpP2 = ped.inPosDur(recData['ag:p2:inPosition.VAL'])
    ifpAG = ped.inPosDur(recData['ag:inPosition'])
    ipcCRCS = ped.inPosCorr(recData['cr:crPosError'])
    # ifpAGOI = inPosDur(recData['ag:oi:inPosition'])
    ipZones = ped.ipZones(recData['tcs:inPosCombine'], 'r')
    mcZones = ped.ipZones(recData['mc:inPosition'], 'g')
    crZones = ped.ipZones(recData['cr:inPosition'], 'b')
    crCZones = ped.ipZones(ipcCRCS, 'y')
    # offTotalZones = offPZones

    print('Size of InPos False Flag duration array: ', len(inFalsePos))

    yLims = (-0.1, 0.1)
    # ylimsIP = (0.45, 0.55)
    ylimsIP = (0.0, 5.55)
    ylimsOFF = (-0.5, 1.5)
    ax_lst = list()

    plts = collections.OrderedDict()
    pC1 = collections.OrderedDict()
    pC2 = collections.OrderedDict()
    pC1['tcsIPdur'] = DataAx(inFalsePos,
                             'r*', 'In Pos\nTime\n[sec]',
                             label="TCS",
                             ylims=ylimsIP)
    # pC1['pqAOff'] = DataAx(recData['tcs:offsetPoA1.VALA'],
                           # 'b', "PQ AB\noffset\n[s]",
                           # label = 'PQ Off A',
                           # drawstyle='steps-post')
    # pC1['pqBOff'] = DataAx(recData['tcs:offsetPoA1.VALB'],
                           # 'r',
                           # shax='pqAOff',
                           # label = 'PQ Off B',
                           # drawstyle='steps-post')
    pC1['tcsIP'] = DataAx(recData['tcs:inPosCombine'],
                          'k', "In Pos\n[bool]",
                          label='TCS',
                          ylims=ylimsOFF,
                          zone=[[crCZones,'CRCS !InPos']],
                          drawstyle='steps-post')
    # pC1['mcsIP'] = DataAx(recData['mc:inPosition'],
                          # 'k', "In Pos\n[bool]",
                          # label='MCS',
                          # ylims=ylimsOFF,
                          # drawstyle='steps-post')
    pC1['crIP'] = DataAx(recData['cr:inPosition'],
                         'k', "In Pos\n[bool]",
                         label='CRCS',
                         ylims=ylimsOFF,
                         drawstyle='steps-post')
    pC1['crIPc'] = DataAx(ipcCRCS,
                          'g', "In Pos\n[bool]",
                          label='CRCS Corr',
                          shax='crIP',
                          linewidth=2.50,
                          ylims=ylimsOFF,
                          drawstyle='steps-post')
    if 'ips:inPosCalc.VALC' in recData.keys():
        pC1['crIPs'] = DataAx(recData['ips:inPosCalc.VALC'],
                            'r', "In Pos\n[bool]",
                            label='CRCS Sim',
                            linewidth=2.50,
                            ylims=ylimsOFF,
                            drawstyle='steps-post')
    if 'ips:inPosCalc.VALD' in recData.keys():
        pC1['crIPOs'] = DataAx(recData['ips:inPosCalc.VALD'],
                            'g', "In Pos\n[bool]",
                            label='CRCS Sim Op',
                            shax='crIPs',
                            linewidth=1.50,
                            ylims=ylimsOFF,
                            drawstyle='steps-post')
    # pC1['agIP'] = DataAx(recData['ag:inPosition'],
                         # 'k', "AG\nIn Pos\n[bool]",
                         # ylims=ylimsOFF,
                         # drawstyle='steps-post')
    # pC1['oiIP'] = DataAx(recData['ag:oi:inPosition'],
                         # 'k', "OI\nIn Pos\n[bool]",
                         # ylims=ylimsOFF,
                         # drawstyle='steps-post')
    # pC1['p1IP'] = DataAx(recData['ag:p1:inPosition.VAL'], 'k', "P1S In Position\n[bool]",
                              # ylims=ylimsOFF, drawstyle='steps-post')
    # pC1['p2IP'] = DataAx(recData['ag:p2:inPosition.VAL'], 'k', "P2S In Position\n[bool]",
                              # ylims=ylimsOFF, drawstyle='steps-post')

    # pC2['mcAzPos'] = DataAx(recData['mc:azCurrentPos'],
                            # 'b-', "MCS Azimuth\nPosition\n[deg]",
                            # height=2)
    # pC2['mcElPos'] = DataAx(recData['mc:elCurrentPos'],
                            # 'g-', "MCS Elevation\nPosition\n[deg]",
                            # height=2)
    # pC2['mcAzdD'] = DataAx(azdD,
                           # 'b-',
                            # label='Az Dmd',
                           # ylims=[-0.05,0.05],
                           # errorline=[-0.004165,0.004165],
                           # zone=mcZones,
                           # zlabel='MCS !InPos',
                           # height=2)
    # pC2['mcAzdP'] = DataAx(azdP,
                           # 'g-',
                            # label='Az Pos',
                           # shax='mcAzdD')
    # pC2['mcAzdTD'] = DataAx(azdTD,
                           # 'c-',
                            # label='Az TCS Dmd',
                           # shax='mcAzdD')
    # pC2['mcEldD'] = DataAx(eldD,
                           # 'b-',
                           # label='El Dmd',
                           # ylims=[-0.05,0.05],
                           # errorline=[-0.004165,0.004165],
                           # zone=mcZones,
                           # zlabel='MCS !InPos',
                           # height=2)
    # pC2['mcEldP'] = DataAx(eldP,
                           # 'g-',
                            # label='El Pos',
                           # shax='mcEldD')
    # pC2['mcEldTD'] = DataAx(eldTD,
                           # 'c-',
                            # label='El TCS Dmd',
                           # shax='mcEldD')
    # pC2['mcAzPErr'] = DataAx(recData['mc:azPosError'],
                             # 'r-', "MCS Az\nPosition Error\n[deg]",
                            # label='Az Err',
                           # shax='mcAzdD')
                             # # ylims=yLims)
    # pC2['mcElPErr'] = DataAx(recData['mc:elPosError'],
                               # 'r-', "MCS El\nPosition Error\n[deg]",
                            # label='El Err',
                           # shax='mcEldD')
                               # # ylims=yLims)
    pC2['crPErr'] = DataAx(recData['cr:crPosError'],
                           'g-', "CRCS\nPos Error\n[deg]",
                           errorline=[-0.0068, 0.0068],
                           label='CRCS Pos Error',
                          zone=[[crZones,'CRCS IP'],
                                [crCZones,'CRCS IP corr']],
                           ylims=yLims)
    pC2['crVel'] = DataAx(recData['cr:crCurrentVel'],
                          'r-', "CRCS\nSpeed\nand Error Speed\n[arcsec/sec]",
                          label='CRCS Velocity',
                          zone=[[crZones,'CRCS IP'],
                                [crCZones,'CRCS IP corr']],
                          ylims=yLims)
    pC2['crdErr'] = DataAx(avgCRdErr,
                           'b-',
                           errorline=[-0.004, 0.004],
                           label='CRCS Err Speed',
                           shax='crVel',
                           ylims=yLims)
    if 'ips:inPosCalc.VALA' in recData.keys():
        pC2['crdErrSim'] = DataAx(recData['ips:inPosCalc.VALA'],
                            'g-', "CRCS\nError Speed Sim\n[arcsec/sec]",
                           errorline=[-0.004, 0.004],
                            label='CRCS Er Speed Sim',
                            ylims=yLims)
    if 'ips:inPosCalc.VALB' in recData.keys():
        pC2['crErrSim'] = DataAx(recData['ips:inPosCalc.VALB'],
                            'b-',
                            errorline=[-0.0068, 0.0068],
                            label='CRCS Err Sim',
                            shax='crdErrSim',
                            ylims=yLims)
    # pC2['crIPo'] = DataAx(recData['cr:inPosition'],
                          # 'k', "CRCS\nIn Pos\n[bool]",
                          # ylims=ylimsOFF,
                          # drawstyle='steps-post')
    # plt.suptitle('In Position Characterization of Gemini South Telescope Systems')
    # plt.show()
    plts['c1'] = pC1
    plts['c2'] = pC2
    pltAx = DataAxesPlotter(plts)
    pltAx.positionPlot()
    pltAx.plotConfig()
