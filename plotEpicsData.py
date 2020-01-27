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
             for a in np.array([dataA[1][(avg-1)-x:len(dataA[1])-x]
                                for x in range(avg)]).T]
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

def gcd(a,b):
    c = a % b
    if not(c):
        return b
    b = gcd(b,c)
    return b

def lcm(a,b):
    c = gcd(a,b)
    d = (a/c) * b
    return d

def lcmArray(a):
    if len(a) == 1:
        return int(a[0])
    r = lcm(a[0],a[1])
    if len(a) == 2:
        return int(r)
    b = [r] + a[2:]
    r = lcmArray(b)
    return int(r)

class DataAx:
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

    def __init__(self, data, linestyle, ylabel, ylims=[], label=None,
                 drawstyle=None, pos=[1,1], height=1,
                 masterax=False, errorline=[], offzone=[], alpha=1.0):
        '''
        Creates the plot object
        '''
        self.linestyle = linestyle
        self.drawstyle = drawstyle
        self.label = label
        self.ylabel = ylabel
        self.ylims = ylims
        self.ax = None
        self.alpha = alpha
        self.data = data
        self.errorline = errorline
        self.offzone = offzone
        self.pos = pos
        self.height = height
        self.mstax = masterax
        self.ys = None
        self.ye = None
        self.xs = pos[1] - 1
        self.xe = pos[1]
        self.bottomax = False


    def plot_ax(self, gs, masterax=None):
        '''
        Plot data in corresponding axes
        '''
        # dataT, dataY = zip(*self.data)
        dataT, dataY = self.data

        # Checks if the first plot has been defined, if not, plot the current
        # data set and set the shareFlag for the rest of the plots.
        if masterax and not(self.mstax) and not(masterax.pos == self.pos):
            # Defines the plot box and main axis with which Y is shared
            self.ax = plt.subplot(gs[self.ys:self.ye,self.xs:self.xe],
                                  sharex=masterax.ax)
        else:
            # Defines the plot box
            self.ax = plt.subplot(gs[self.ys:self.ye,self.xs:self.xe])

        # Plot the data
        if not(self.drawstyle):
            self.ax.plot(dataT, dataY, self.linestyle,
                         label=self.label, alpha=self.alpha)
        else:
            self.ax.plot(dataT, dataY, c=self.linestyle,
                         drawstyle=self.drawstyle, label=self.label,
                         alpha=self.alpha)

        if self.errorline:
            self.ax.axhline(y=self.errorline[1], linestyle='-.',
                            linewidth=1.25, color='crimson')
            self.ax.axhline(y=self.errorline[0], linestyle='-.',
                            linewidth=1.25, color='crimson')
            self.ax.fill_between(dataT, dataY, self.errorline[1],
                                 where=[True if y>self.errorline[1]
                                        else False for y in dataY], color='y')
            self.ax.fill_between(dataT, dataY, self.errorline[0],
                                 where=[True if y<self.errorline[0]
                                        else False for y in dataY], color='y')

        if self.offzone:
            for oz in self.offzone:
                self.ax.axvspan(oz[0], oz[1], facecolor=oz[2], alpha=0.4)

        self.ax.grid(True)
        self.ax.tick_params("y", colors="b")
        self.ax.set_ylabel(self.ylabel, color="b")
        if self.ylims:
            self.ax.set_ylim(self.ylims[0], self.ylims[1])

# axPlt object end

class DataAxesPlotter:
    def __init__(self, Axes):
        self.dictAxes = Axes
        self.masterax = None
        # self.bottomax = []
        self.gs = None
        self.fig = plt.figure()

    def addDataAxes(self):
        # This method starts by constructing a dictionary for the position of
        # each plot within a matrix. The idea is to be able to calculate the
        # proper coordinates and dimentions for each plot
        # The dictionary contains a subdict for each column as well as a key
        # with the total number of columns for the plot
        Dax = self.dictAxes
        Da = {'nCols':0} # Initialize dictionary
        for n in Dax:
            # Initialize each subdict
            if not(str(Dax[n].pos[1]) in Da):
                Da[str(Dax[n].pos[1])] = {'nRows':0, 'nCells':0}
            # Add the height and name of each plot
            Da[str(Dax[n].pos[1])][str(Dax[n].pos[0])] = [[0,Dax[n].height], n]
            # Update the number of cols of the plot area
            if Dax[n].pos[1] > Da['nCols']:
                Da['nCols'] = Dax[n].pos[1]
            # Update the number of rows in each column
            if Dax[n].pos[0] > Da[str(Dax[n].pos[1])]['nRows']:
                Da[str(Dax[n].pos[1])]['nRows'] = Dax[n].pos[0]
        # Complete each column with unit spaces for non defined plots. These
        # units spaces are place holders for plots with height 1
        for i in range(1,Da['nCols']+1):
            for j in range(1,Da[str(i)]['nRows']+1):
                if not(str(j) in Da[str(i)].keys()):
                    Da[str(i)][str(j)] = [[0, 1], '']
                # Calculate the number of rows in each column
                Da[str(i)]['nCells'] += Da[str(i)][str(j)][0][1]
        # print(Da)
        # Array with the number of rows in each column. This array is used to
        # calculate the least common multiple for these numbers
        rows = [Da[i]['nCells'] for i in Da.keys() - ['nCols','nCells']]
        # print(rows)
        totCells = lcmArray(rows)
        self.gs = gridspec.GridSpec(totCells,Da['nCols'])
        # Go through the indexing dictionary to calculate the final coordinates
        # for each plot
        for i in range(1,Da['nCols']+1):
            for j in range(1,Da[str(i)]['nRows']+1):
                plotname = Da[str(i)][str(j)][1]
                if plotname in Dax.keys():
                    Dax[plotname].ys = int(Da[str(i)][str(j)][0][0]
                                           * (totCells/Da[str(i)]['nCells']))
                    Dax[plotname].ye = int(Da[str(i)][str(j)][0][1]
                                           * (totCells/Da[str(i)]['nCells']))
                    # Define the bottom plot for each column
                    if j == Da[str(i)]['nRows']:
                        Dax[plotname].bottomax = True
                    if Dax[plotname].mstax:
                        self.masterax = plotname
                        Dax[plotname].plot_ax(self.gs)
                if (j < Da[str(i)]['nRows']):
                    # calculate the starting and ending X-coords
                    ys = Da[str(i)][str(j)][0][1]
                    ye = Da[str(i)][str(j+1)][0][1] + ys
                    Da[str(i)][str(j+1)][0] =  [ys, ye]
        # print(Da)
        for n in Dax.keys()-[self.masterax]:
            Dax[n].plot_ax(self.gs, Dax[self.masterax])

    def plotConfig(self):
        '''
        This function is used to:
            - Disable tick labels for plots other than the bottom positions
            - Set the date labels on the bottom plots
            - Update Plot Event Zones
        '''
        Dax = self.dictAxes
        for n in Dax:
            Dax[n].ax.grid(True)
            Dax[n].ax.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
            # ax[0].ax.xaxis_date()
            if not(Dax[n].bottomax):
                plt.setp(Dax[n].ax.get_xticklabels(), fontsize=9, visible=False)
            else:
                # Dax[n].ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
                # Dax[n].ax.xaxis.set_major_formatter(
                    # matplotlib.dates.DateFormatter("%d/%m %H:%M:%S.%f"))
                # Dax[n].ax.xaxis.set_minor_locator(ticker.MaxNLocator(200))
                plt.setp(Dax[n].ax.get_xticklabels(), fontsize=9,
                        rotation=30, ha='right')
                # # ax[0].ax.set_xlim(lowX,highX)
                # ax[0].ax.set_xlim(time[5],time[len(time)-5])
                # h, l = ax[0].ax.get_legend_handles_labels()
                # ax[0].ax.legend(h, l, loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
        plt.show()



if __name__ == '__main__':
    # plotConfig(ax_lst)
    x = np.arange(21)
    y1 = x**2
    y2 = 50 * (np.sin((3.1415/2) * x))
    y3 = 100 / x
    plts = {}
    # plts['g1'] = DataAx('', '', '', pos=[2,1], height=1)
    # plts['g2'] = DataAx('', '', '', pos=[1,2], height=2)
    # plts['g3'] = DataAx('', '', '', pos=[3,1], height=4)
    plts['g4'] = DataAx([x,y1], 'r-', 't1', pos=[1,2], height=2, masterax=True)
    # plts['g5'] = DataAx([x,y2], 'b-', 't2', pos=[1,1], height=5, masterax = True)
    plts['g2'] = DataAx([x,((x-10)**2)+40], 'g-', 't4', pos=[2,1], height=2)
    plts['g6'] = DataAx([x,y3], 'k-', 't3', pos=[2,2], height=1)
    pltAx = DataAxesPlotter(plts)
    pltAx.addDataAxes()
    for p in plts:
        print('Plot {2} Span in col: {0} - {1}'.format(plts[p].ys, plts[p].ye, p))
    pltAx.plotConfig()

    # plt.suptitle('In Position Characterization of Gemini South Telescope Systems')
    # plt.show()
