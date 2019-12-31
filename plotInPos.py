#!/usr/bin/env python3.5

from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from matplotlib import dates
import matplotlib

import pickle
import argparse
import numpy as np
import h5py
import re

PLOT_AREA = (2,2)

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

    args = parser.parse_args()
    return args

def rmsChan(dataSet):
    time, data = zip(*dataSet)
    ts = time[0]
    i = 0
    s = i
    drms = 0
    dataRMS = list()
    for t in time:
        if (((t - ts) < timedelta(minutes=60)) and not(t == time[len(time)-1])):
            drms += data[i]**2
            i += 1
        else:
            drms = np.sqrt(drms) / (i - s)
            auxRMS = [x * drms for x in np.ones((i-s))]
            dataRMS.extend(auxRMS)
            ts = time[i]
            drms = 0
            s = i - 1
    return list(zip(time,dataRMS))

def inPosCorr(dataA):
    dataCR = derivData(dataA)
    dataCRA = avgData(dataCR, 3)
    if (abs(dataCRA[1][0])>0.003) and (abs(dataA[1][3])>0.0068):
        currIP = 0
    else:
        currIP = 1
    corrInPos = [currIP]
    for y,dy in np.array([dataA[1][4:],dataCRA[1][1:]]).T:
        if ((abs(y)<0.0068) and (abs(dy)<0.004) and not(currIP)):
            currIP = 1
        elif ((abs(y)>0.0068) and (abs(dy)>0.004) and currIP):
            currIP = 0
        corrInPos.append(currIP)
    # corrInPos = [1 if ((abs(y)<0.0068) and (abs(dy)<0.004)) else 0 for y,dy in np.array([dataA[1][3:],dataCRA[1]]).T]
    corrInPosData = [dataA[0][3:],corrInPos]
    return corrInPosData

def diffData(dataA):
    '''
    This function is used to calculate the differential of the values of an
    array between consecutive samples
    '''
    diffPV = [y - x \
               for x,y in np.array([dataA[1][:-1],
                                    dataA[1][1:]]).T]
    diffA = [dataA[0][1:], diffPV]
    return diffA

def avgData(dataA, avg):
    '''
    This function is used to calculate the moving average of a data array
    '''
    # avg = 3
    avgPV = [np.average(a)\
             for a in np.array([dataA[1][(avg-1)-x:len(dataA[1])-x] for x in range(avg)]).T]
    avgA = [dataA[0][(avg-1):], avgPV]
    return avgA

def derivData(dataA):
    '''
    This function is used to calculate the first derivative of data array
    between consecutive values
    '''
    diffPV = [(x - xa)/(y - ya).total_seconds() \
              for ya,y,xa,x in np.array([dataA[0][:-1], dataA[0][1:],
                                   dataA[1][:-1], dataA[1][1:]]).T]
    diffA = [dataA[0][1:], diffPV]
    return diffA

def inPosDur(ipArray):
    inPos = False
    ipd = list()
    previpt = ipArray[0][0]
    for t,ip in np.array(ipArray).T:
        if not(ip) and inPos:
            inPos = False
            ipd.append([t, (t - previpt).total_seconds()])
        if ip and not(inPos):
            inPos = True
            previpt = t
    return ipd

def offsetZones(offsetArray, zcolor):
    zeroFlag = False
    offsetPA = []
    offsetNA = []
    for t,ov in np.array(offsetArray).T:
        if not(ov) and zeroFlag:
            zeroFlag = False
            offsetPA.append(t)
        if ov and not(zeroFlag):
            zeroFlag = True
            offsetNA.append(t)
    # offsZones = [[row[i] for row in [offsetPA,offsetNA]]
                 # for i in range(len(offsetPA))]
    if len(offsetNA) > len(offsetPA):
        offsetNA = offsetNA[:-1]
    zonecolor = [zcolor for i in range(len(offsetPA))]
    print('size array PA: ', len(offsetPA), 'size array NA: ', len(offsetNA))
    offsetZonesArray = np.array([offsetNA,offsetPA,zonecolor]).T
    return offsetZonesArray

def plotConfig(ax_lst):
    '''
    This function is used to:
        - Disable tick labels for plots other than the bottom positions
        - Set the date labels on the bottom plots
        - Update Plot Event Zones
    '''
    for ax in ax_lst:
        ax[0].ax.grid(True)
        ax[0].ax.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
        # ax[0].ax.xaxis_date()
        if not(ax[1]):
            plt.setp(ax[0].ax.get_xticklabels(), fontsize=9, visible=False)
        else:
            ax[0].ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
            ax[0].ax.xaxis.set_major_formatter(
                matplotlib.dates.DateFormatter("%d/%m %H:%M:%S.%f"))
            ax[0].ax.xaxis.set_minor_locator(ticker.MaxNLocator(200))
            plt.setp(ax[0].ax.get_xticklabels(), fontsize=9,
                     rotation=30, ha='right')
            # # ax[0].ax.set_xlim(lowX,highX)
            # ax[0].ax.set_xlim(time[5],time[len(time)-5])
            # h, l = ax[0].ax.get_legend_handles_labels()
            # ax[0].ax.legend(h, l, loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')

class AxPlt:
    '''
    This class defines a plot object
    Attributes
        - plotArea: Figure grid area for the plots
        - gs: Grid size, related to plotArea
        - sharedAx, shareFlag: Flags used to determine axis behavior
        - linestyle: Color, line type and marker for each plot
        - ylabel: Y-axis label
        - autoyax: used when Y-axis limits are not pre-defined
        - ylims: Y-axis limits
        - ax: the pyplot.ax object to be plotted
        - dataT: Time-axis data
        - dataY: Y-axis data

    Methods
        - plot_ax: Plots axis object. Checks and adds the shared-x feature.
          Adds horizontal lines (to visualize limits.
        - plot_twin: Plot a second set of data on an additional axis to
          the right of the plot box.
        - add_plot: Adds additional plots to the existing axis. This needs to
          be called after plot_ax as been defined on
          the corresponding object.

    '''
    # Attributes common to every instance
    plotArea = None
    gs = None
    sharedAx = []
    shareFlag = []

    def __init__(self, linestyle, ylabel,
                 ylims=(0,0), autoyax=True):
        '''
        Creates the plot obejct
        '''
        self.linestyle = linestyle
        self.drawstyle = None
        self.label = None
        self.ylabel = ylabel
        self.autoyax = autoyax
        self.ylims = ylims
        self.ax = None
        self.dataT = None
        self.dataY = None
        self.alpha = 1.0

    def plot_ax(self, ax_lst, position, height, width, data, offZone = [],
                bottomPlot=False, errorLine=False, zonesPlot=True,
                ds=None, plotLabel=None, alphaT=1.0, stemp = False):
        '''
        Plot data in corresponding axes
        '''
        self.dataT, self.dataY = zip(*data)
        self.label = plotLabel
        self.drawstyle = ds
        self.alpha = alphaT
        # The position and spans are collected this way because it was adjusted
        # to a previous version of this method :P
        yp = position[0]
        xp = position[1]
        yl = yp + height
        xl = xp + width

        # Checks if the first plot has been defined, if not, plot the current
        # data set and set the shareFlag for the rest of the plots.
        if not(self.shareFlag):
            # Defines the plot box
            self.ax = plt.subplot(self.gs[yp:yl,xp:xl])
            self.sharedAx.append(self.ax)
            self.shareFlag.append(True)
        else:
            # Defines the plot box and main axis with which Y is shared
            self.ax = plt.subplot(self.gs[yp:yl,xp:xl],
                                  sharex=self.sharedAx[0])

        # Plot the data
        if not(self.drawstyle):
            self.ax.plot(self.dataT, self.dataY, self.linestyle,
                         label=self.label, alpha=self.alpha)
        else:
            self.ax.plot(self.dataT, self.dataY, c=self.linestyle,
                         drawstyle=self.drawstyle, label=self.label, alpha=self.alpha)

        if errorLine:
            self.ax.axhline(y=0.0068, linestyle='-.', linewidth=1.25, color='crimson')
            # self.ax.axhline(y=0.006, linestyle='-.', linewidth=1.25, color='darkorange')
            self.ax.axhline(y=-0.0068, linestyle='-.', linewidth=1.25, color='crimson')
            # self.ax.axhline(y=-0.006, linestyle='-.', linewidth=1.25, color='darkorange')
            self.ax.fill_between(self.dataT, self.dataY, 0.0068,
                                 where=[True if y>0.0068 else False for y in self.dataY], color='y')
            self.ax.fill_between(self.dataT, self.dataY, -0.0068,
                                 where=[True if y<-0.0068 else False for y in self.dataY], color='y')

        if len(offZone):
            for oz in offZone:
                self.ax.axvspan(oz[0], oz[1], facecolor=oz[2], alpha=0.4)

        self.ax.grid(True)
        self.ax.tick_params("y", colors="b")
        self.ax.set_ylabel(self.ylabel, color="b")
        if not(self.autoyax):
            self.ax.set_ylim(self.ylims[0], self.ylims[1])
        ax_lst.append([self, bottomPlot, zonesPlot])

    def plot_twin(self, ax_lst, ax_prime, data,
                  bottomPlot=False, errorLine=False, zonesPlot=True):
        '''
        Plot data sharing plot box with different Y axis
        '''
        # Set the twin axis
        self.ax = ax_prime.twinx()
        self.dataT, self.dataY = zip(*data)
        self.ax.plot(self.dataT, self.dataY, self.linestyle)
        self.ax.tick_params("y", colors="b")
        # This line is important if you want to zoom both axis to the same
        # values at the same time.
        self.ax.get_shared_y_axes().join(self.ax, ax_prime)
        if not(self.autoyax):
            self.ax.set_ylim(self.ylims[0], self.ylims[1])
        ax_lst.append([self, bottomPlot, zonesPlot])

    def add_plot(self, ax_lst, data, ls=None, ds=None, plotLabel=None,
                 bottomPlot=False, errorLine=False, zonesPlot=False, alphaT=1.0):
        '''
        Add plots to the same box
        '''
        if not(ls):
            ls = self.linestyle
        self.dataT, self.dataY = zip(*data)
        self.ax.plot(self.dataT, self.dataY, ls, drawstyle=ds,
                     label=plotLabel, alpha=alphaT)
        self.ax.tick_params("y", colors="b")
        # self.ax.get_shared_y_axes().join(self.ax, ax_prime)
        if not(self.autoyax):
            self.ax.set_ylim(self.ylims[0], self.ylims[1])
        ax_lst.append([self, bottomPlot, zonesPlot])

# axPlt object end

if __name__ == '__main__':
    args = parse_args()
    # Read TCS in position from camonitor txt file
    # with open('testinpos.txt', 'r') as f:
        # tcsIPD = f.read().splitlines()

    # Create empty arrays for TCS in position time stamp and value
    dataT = list()
    dataV = list()
    # for l in tcsIPD:
        # # Search for tcs:inPosition record
        # if re.search('inPosition', l):
            # # Append time stamp of inPosition instance
            # dataT.append(re.search('2019-10-07 [012][0-9]:[0-9]{2}:[0-9]{2}.[0-9]+', l).group(0))
            # # Append 1 if value is TRUE string and 0 if FALSE
            # if re.search('TRUE', l):
                # dataV.append(1)
            # elif re.search('FALSE', l):
                # dataV.append(0)
    # # Format time stamp to datetime object
    # dtT = [datetime.strptime(t, '%Y-%m-%d %H:%M:%S.%f') for t in dataT]
    # Create final array
    # tcsipData = [dtT,dataV]

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
    recData = {name:[np.array(g.get('timestamp')),np.array(g.get('value'))]\
               for g,name in recGroups}
    # Reformat timestamp data
    for n in recData:
        recTime = [datetime.fromtimestamp(ts) for ts in recData[n][0]]
        timeStampS = [datetime.strftime(rt, '%m/%d/%Y %H:%M:%S.%f')\
                      for rt in recTime]
        timeStamp = [datetime.strptime(ts, '%m/%d/%Y %H:%M:%S.%f')\
                     for ts in timeStampS]
        recData[n][0] = timeStamp

    # azdP = diffData(recData['mc:azCurrentPos'])
    # azdD = diffData(recData['tcs:demandAz'])
    # azErr = diffData(recData['mc:azPosError'])
    # eldP = diffData(recData['mc:elCurrentPos'])
    # eldD = diffData(recData['tcs:demandEl'])
    # elErr = diffData(recData['mc:elPosError'])
    # crdP = diffData(recData['cr:crCurrentPos'])
    # crdV = derivData(recData['cr:crCurrentVel'])
    crdV = derivData(recData['cr:crPosError'])
    avgCRdV = avgData(crdV,3)
    # avgCRdV = crdV

    # inFalsePos = inPosDur(tcsipData)
    inFalsePos = inPosDur(recData['tcs:inPosCombine'])
    ifpMCS = inPosDur(recData['mc:inPosition'])
    ifpCRCS = inPosDur(recData['cr:inPosition'])
    ifpP1 = inPosDur(recData['ag:p1:inPosition.VAL'])
    ifpP2 = inPosDur(recData['ag:p2:inPosition.VAL'])
    ifpAG = inPosDur(recData['ag:inPosition'])
    ipcCRCS = inPosCorr(recData['cr:crPosError'])
    # ifpAGOI = inPosDur(recData['ag:oi:inPosition'])
    # offPZones = offsetZones(recData['tcs:offsetPoA1.VALA'], 'r')
    # offTotalZones = offPZones

    print('Size of InPos False Flag duration array: ', len(inFalsePos))

    ylims = (-0.1, 0.1)
    # ylimsIP = (0.45, 0.55)
    ylimsIP = (0.0, 5.55)
    ylimsOFF = (-1, 2)
    ax_lst = list()
    AxPlt.plotArea = (10,2)
    AxPlt.gs = gridspec.GridSpec(AxPlt.plotArea[0], AxPlt.plotArea[1])
    # print('size az dmd array: ', len(recData['tcs:demandAz'][1]),
          # ' size az pos array ', len(recData['mc:azCurrentPos'][1]))

    plot10 = AxPlt('r*', "TCS In Position Duration\n[sec]", ylimsIP, False)
    plot11 = AxPlt('b', "PQ offset\n[s]", False)
    plot12 = AxPlt('k', "TCS\nIn Position\n[bool]", ylimsOFF, False)
    plot13 = AxPlt('k', "CRCS\nIn Position\n[bool]", ylimsOFF, False)
    plot14 = AxPlt('k', "MCS\nIn Position\n[bool]", ylimsOFF, False)
    plot15 = AxPlt('k', "AG\nIn position\n[bool]", ylimsOFF, False)
    plot16 = AxPlt('k', "OI\nIn position\n[bool]", ylimsOFF, False)
    plot17 = AxPlt('k', "P1\nIn Position\n[bool]", ylimsOFF, False)
    plot18 = AxPlt('k', "P2\nIn Position\n[bool]", ylimsOFF, False)
    # plot18 = AxPlt('r*', "In Position Flag Duration\n[s]", ylimsIP, False)
    # plot19 = AxPlt('r*', "In Position Flag Duration\n[s]", ylimsIP, False)

    plot20 = AxPlt('b-', "MCS Azimuth\nPosition\n[deg]", False)
    plot21 = AxPlt('g-', "MCS Elevation\nPosition\n[deg]", False)
    plot22 = AxPlt('b-', "MCS Azimuth\nPosition Error\n[deg]", ylims, False)
    plot23 = AxPlt('g-', "MCS Elevation\nPosition Error\n[deg]", ylims, False)
    plot24 = AxPlt('g-*', "CRCS\nPosition Error\n[deg]", ylims, False)
    plot25 = AxPlt('r-', "CRCS\nVelocity\n[bool]", ylims, False)
    plot26 = AxPlt('b-', "CRCS\ndiff Velocity\n[bool]", ylims, False)
    plot27 = AxPlt('k', "CRCS\nIn Position\n[bool]", ylimsOFF, False)
    # plot24 = AxPlt('c', "PQ offset\n[s]", False)
    # plot25 = AxPlt('c', "PQ offset\n[s]", False)
    # plot26 = AxPlt('k', "MCS in position\n[s]", ylimsOFF, False)
    # plot27 = AxPlt('k', "MCS in position\n[s]", ylimsOFF, False)
    # plot28 = AxPlt('r*', "In Position Flag Duration\n[s]", ylimsIP, False)
    # plot29 = AxPlt('r*', "In Position Flag Duration\n[s]", ylimsIP, False)

    plot10Array = inFalsePos
    plot11Array = np.array(recData['tcs:offsetPoA1.VALA']).T
    plot11bArray = np.array(recData['tcs:offsetPoA1.VALB']).T
    # plot12Array = np.array(tcsipData).T
    plot12Array = np.array(recData['tcs:inPosCombine']).T
    plot13Array = np.array(recData['cr:inPosition']).T
    plot14Array = np.array(ipcCRCS).T
    # plot14Array = np.array(recData['mc:inPosition']).T
    plot15Array = np.array(recData['ag:inPosition']).T
    plot16Array = np.array(recData['ag:oi:inPosition']).T
    plot17Array = np.array(recData['ag:p1:inPosition.VAL']).T
    plot18Array = np.array(recData['ag:p2:inPosition.VAL']).T

    # plot20Array = np.array(azdP).T
    # plot20bArray = np.array(azdD).T
    # plot21Array = np.array(eldP).T
    # plot21bArray = np.array(eldD).T
    # plot22Array = np.array(recData['cr:crPmacPosError']).T
    plot20Array = np.array(recData['mc:azCurrentPos']).T
    plot21Array = np.array(recData['mc:elCurrentPos']).T
    plot22Array = np.array(recData['mc:azPosError']).T
    plot23Array = np.array(recData['mc:elPosError']).T
    plot24Array = np.array(recData['cr:crPosError']).T
    plot25Array = np.array(recData['cr:crCurrentVel']).T
    # plot24Array = np.array(crdP).T
    # plot24bArray = np.array(crdD).T
    plot26Array = np.array(avgCRdV).T
    plot27Array = plot13Array
    # plot26Array = np.array(tcsipData).T
    # plot13Array = np.array(recData['mc:inPosition']).T
    # plot14Array = np.array(recData['pwfs1:inPosition']).T
    # plot15Array = np.array(recData['pwfs2:inPosition']).T

    plot10.plot_ax(ax_lst, (0,0), 1, 1,
                     plot10Array,
                     plotLabel='InPos Duration')
                     # plotLabel='InPos Duration',
                     # offZone = offTotalZones)
    plot11.plot_ax(ax_lst, (1,0), 1, 1,
                      plot11Array,
                      plotLabel='TCS InPos', ds='steps-post')
                      # plotLabel='TCS InPos', ds='steps-post',
                      # offZone = offTotalZones)
    plot11.add_plot(ax_lst, plot11bArray, ls='r',
                      plotLabel='TCS InPos', ds='steps-post')
    plot12.plot_ax(ax_lst, (2,0), 1, 1,
                      plot12Array,
                      plotLabel='TCS InPos', ds='steps-post')
                      # plotLabel='TCS InPos', ds='steps-post',
                      # offZone = offTotalZones)
    plot13.plot_ax(ax_lst, (3,0), 1, 1,
                      plot13Array,
                      plotLabel='CRCS InPos', ds='steps-post')
                      # plotLabel='CRCS InPos', ds='steps-post',
                      # offZone = offTotalZones)
    plot14.plot_ax(ax_lst, (4,0), 1, 1,
                      plot14Array,
                      plotLabel='MCS InPos', ds='steps-post')
                      # plotLabel='MCS InPos', ds='steps-post',
                      # offZone = offTotalZones)
    plot15.plot_ax(ax_lst, (5,0), 1, 1,
                      plot15Array,
                      plotLabel='AG InPos', ds='steps-post')
                      # plotLabel='AG InPos', ds='steps-post',
                      # offZone = offTotalZones)
    plot16.plot_ax(ax_lst, (6,0), 1, 1,
                      plot16Array,
                      plotLabel='AG-OI InPos', ds='steps-post')
                      # plotLabel='AG-OI InPos', ds='steps-post',
                      # offZone = offTotalZones)
    plot17.plot_ax(ax_lst, (7,0), 1, 1,
                      plot17Array,
                      plotLabel='P1 InPos', ds='steps-post')
                      # plotLabel='P1 InPos', ds='steps-post',
                      # offZone = offTotalZones)
    plot18.plot_ax(ax_lst, (8,0), 1, 1,
                      plot18Array,
                      bottomPlot=True,
                      plotLabel='P2 InPos', ds='steps-post')
                      # plotLabel='P2 InPos', ds='steps-post',
                      # offZone = offTotalZones)

    plot20.plot_ax(ax_lst, (0,1), 2, 1, plot20Array,
                   plotLabel='Az Position')
                   # plotLabel='Az Position',
                   # offZone = offTotalZones)
    # plot20.add_plot(ax_lst, plot20bArray, ls='r-',
                    # plotLabel='Az Diff Dmd', alphaT=0.40)
    plot21.plot_ax(ax_lst, (2,1), 2, 1, plot21Array,
                   plotLabel='El Position')
                   # plotLabel='El Position',
                   # offZone = offTotalZones)
    # plot21.add_plot(ax_lst, plot21bArray, ls='r-',
                    # plotLabel='El Diff Dmd', alphaT=0.40)
    plot22.plot_ax(ax_lst, (4,1), 1, 1,
                   plot22Array,
                   plotLabel='Az Pos Error')
                   # plotLabel='Az Pos Error',
                   # offZone = offTotalZones)
    plot23.plot_ax(ax_lst, (5,1), 1, 1,
                   plot23Array,
                   plotLabel='El Pos Error')
                   # plotLabel='El Pos Error',
                   # offZone = offTotalZones)
    plot24.plot_ax(ax_lst, (6,1), 1, 1,
                   plot24Array,
                   errorLine=True,
                   plotLabel='CRCS Error')
                   # plotLabel='CRCS Error',
                   # offZone = offTotalZones)
    # plot24.add_plot(ax_lst, plot24bArray, ls='r-',
                    # plotLabel='El Diff Dmd', alphaT=0.40)
    plot25.plot_ax(ax_lst, (7,1), 1, 1,
                   plot25Array,
                   plotLabel='CRCS Vel')
                   # plotLabel='CRCS InPos', ds='steps-post',
                   # offZone = offTotalZones)
    plot26.plot_ax(ax_lst, (8,1), 1, 1,
                   plot26Array,
                   plotLabel='CRCS diff Vel')
                   # plotLabel='TCS InPos', ds='steps-post',
                   # offZone = offTotalZones)
    plot27.plot_ax(ax_lst, (9,1), 1, 1,
                   plot27Array,
                   bottomPlot=True,
                   plotLabel='CRCS InPos', ds='steps-post')
                   # plotLabel='TCS InPos', ds='steps-post',
                   # offZone = offTotalZones)

    plotConfig(ax_lst)

    plt.suptitle('In Position Characterization of Gemini South Telescope Systems')
    plt.show()
