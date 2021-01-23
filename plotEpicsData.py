#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec

from matplotlib import dates
from datetime import datetime, timedelta

import collections
import numpy as np
import h5py
import sys
import re

def extract_h5df(hdf5File, stime, etime, listonly=False):
    # Create empty arrays for TCS in position time stamp and value
    recGroups = []
    recData = {}
    print(hdf5File)
    # Read h5 file
    data_files = [h5py.File(datafile, 'r') for datafile in hdf5File]
    # Extract record names and create an array with group object and name of
    # group
    for data_group in data_files:
        for rn in data_group:
            timestamps = data_group.get(rn).get('timestamp')[0:]
            ts_index = np.where(np.logical_and(stime<timestamps,
                                               etime>timestamps))
            recData[rn] =[timestamps[ts_index],
                          data_group.get(rn).get('value')[ts_index]]
    print(recData.keys())
    if listonly:
        sys.exit()
    return recData

def rmsChan(dataSet):
    '''
    Calculate the RMS value of a dataset

    Parameters
    ----------
    dataSet : array
        dataSet is a list with two arrays, Process Value data and timestamps of
        an EPICS channel.

    Returns
    -------
    (time, dataRMS) : list with two arrays array with timestamps and RMS value
    of the EPICS channel process value data
    '''
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
    if (abs(dataCRA[1][0])>0.004) and (abs(dataA[1][3])>0.0068):
        currIP = 0
    else:
        currIP = 1
    corrInPos = [currIP]
    for y,dy in np.array([dataA[1][4:],dataCRA[1][1:]]).T:
        # if ((abs(y)<0.0068) and (abs(dy)<0.004) and not(currIP)):
        if ((abs(y)<0.0068) and (abs(dy)<0.004)):
            currIP = 1
        # elif ((abs(y)>0.0068) and (abs(dy)>0.004) and currIP):
        elif ((abs(y)>0.0068) and (abs(dy)>0.004)):
            currIP = 0
        corrInPos.append(currIP)
    corrInPosData = [dataA[0][3:],corrInPos]
    return corrInPosData

def inPosCorrH(dataA):
    dataCR = derivData(dataA)
    dataCRA = avgData(dataCR, 3)
    if (abs(dataCRA[1][0])>0.004) and (abs(dataA[1][3])>0.0068):
        currIP = 0
    else:
        currIP = 1
    corrInPos = [currIP]
    for y,dy in np.array([dataA[1][4:],dataCRA[1][1:]]).T:
        # if ((abs(y)<0.0068) and (abs(dy)<0.004) and not(currIP)):
        if ((abs(y)<0.006) and (abs(dy)<0.0004)):
            currIP = 1
        # elif ((abs(y)>0.0068) and (abs(dy)>0.004) and currIP):
        elif ((abs(y)>0.0068) and (abs(dy)>0.004)):
            currIP = 0
        corrInPos.append(currIP)
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
    return np.array(ipd).T

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
    if len(offsetNA) > len(offsetPA):
        offsetNA = offsetNA[:-1]
    zonecolor = [zcolor for i in range(len(offsetPA))]
    print('size array PA: ', len(offsetPA), 'size array NA: ', len(offsetNA))
    offsetZonesArray = np.array([offsetNA,offsetPA,zonecolor]).T
    return offsetZonesArray

def ipZones(offsetArray, zcolor):
    zeroFlag = False
    offsetPA = []
    offsetNA = []
    for t,ov in np.array(offsetArray).T:
        if ov and zeroFlag:
            zeroFlag = False
            offsetPA.append(t)
        if not(ov) and not(zeroFlag):
            zeroFlag = True
            offsetNA.append(t)
    if len(offsetNA) > len(offsetPA):
        offsetNA = offsetNA[:-1]
    zonecolor = [zcolor for i in range(len(offsetPA))]
    print('size array PA: ', len(offsetPA), 'size array NA: ', len(offsetNA))
    offsetZonesArray = np.array([offsetNA,offsetPA,zonecolor]).T
    return offsetZonesArray

def gcd(a,b):
    '''
    gcd() calculates the greatest common denominator
    '''
    c = a % b
    if not(c):
        return b
    b = gcd(b,c)
    return b

def lcm(a,b):
    '''
    lcm() calculates the least common multiple between two numbers
    '''
    c = gcd(a,b)
    d = (a/c) * b
    return d

def lcmArray(a):
    '''
    lcm() calculates the least common multiple between two or more numbers
    '''
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
        - ax: the pyplot.ax handling each pair of axe
        - shaxname: Name of DataAx object handling the common axe
        - shax: the pyplot.ax handler of the shaxname object
        - ylabel: Y-axis label
        - ylims: Y-axis limits
        - errorline: Y coord for an error line
        - linestyle: Color, line type and marker for each plot
        - autoyax: used when Y-axis limits are not pre-defined
        - dataT: Time-axis data
        - color: Line color
        - lwidth: Line width
        - label: Line label, displayed on legend box
        - ylabel: Y axis label
        - ylims: Y axis limits
        - alpha: Line transparency
        - zone: Define the limits of a colored zone
        - pos: Plot position within the figure
        - height: Plot height within the figure
        - mstax: Define one of the plots as the master axis
        - ys: Plot vertical starting point within the figure
        - ye: Plot vertical end point within the figure
        - xs: Plot horizontal starting point within the figure
        - xe: Plot horizontal end point within the figure
        - bottomax: Bottom axis flag for each plot
        - sharedax: Shared axis flag
        - histbins: Number of bins for a histogram
        - limsbins: Limits for the span of the bins
        - bindata: Output of the plot.hist function
        - bins: Output of the plot.hist function
        - patches: Output of the plot.hist function

    Methods
        - plot_ax: Plots axis object. Checks and adds the shared-x feature.
          Adds horizontal lines (to visualize limits).

    '''

    def __init__(self, data, color, linestyle='-', marker=None,
                 drawstyle='default', label=None, ylabel='',
                 xlabel='', ylims=[], shax=None, height=1,
                 errorline=[], zone=[], alpha=1.0,
                 linewidth=1.25, histbins=False, limsbins=None):
        '''
        Creates the plot object
        '''
        # Axe attributes
        self.ax = None
        self.shaxname = shax
        self.shax = None
        self.ylabel = ylabel
        self.xlabel = xlabel
        self.ylims = ylims
        self.errorline = errorline
        self.zone = zone
        self.height = height
        self.ys = None
        self.ye = None
        self.xs = None
        self.xe = None
        self.bottomax = False
        # Plot attributes
        self.data = data
        self.color = color
        self.linestyle = linestyle
        self.marker = marker
        self.drawstyle = drawstyle
        self.lwidth = linewidth
        self.alpha = alpha
        self.label = label
        self.histbins = histbins
        self.limsbins = limsbins
        self.bindata = False
        self.bins = False
        self.patches = False


    def plot_ax(self, gs, tcells, rows, masterax):
        ''' plot_ax

        This method defines each axe to be plotted. It will prepare each ax
        depending if is a regular plot or a histogram
        It places each ax within the gridspace, and defines all its
        characteristics
        '''
        if not(self.histbins):
            shared_x_axis = masterax.ax
        else:
            shared_x_axis = None

        if self.shax:
            self.ax = self.shax.ax
            self.ylabel = self.shax.ylabel
            self.bottomax = self.shax.bottomax
        else:
            ys = int(self.ys * (tcells/rows))
            ye = int(self.ye * (tcells/rows))
            self.ax = plt.subplot(gs[ys:ye,self.xs:self.xe],
                                  sharex=shared_x_axis)
        # Check if plot is a histogram. If it's not, prepares the x & y axis for
        # plotting. If it is prepares an array with the x values.
        if not(self.histbins):
            dataT, dataY = self.data
            dataT = [datetime.fromtimestamp(ts) for ts in dataT]
        else:
            dataX = self.data

        # Checks if the first plot has been defined, if not, plot the current
        # data set and set the shareFlag for the rest of the plots.

        # Configures each axe for plotting
        if not(self.histbins):
            self.ax.plot(dataT, dataY, linestyle=self.linestyle,
                         linewidth=self.lwidth, color=self.color,
                         marker=self.marker, drawstyle=self.drawstyle,
                         label=self.label, alpha=self.alpha)
            # Define a two horizontal line in the ax, if errorline has been set
            if self.errorline:
                self.ax.axhline(y=self.errorline[1], linestyle='-.',
                                linewidth=1.55, color='crimson',
                                label='|error|={0}'.format(self.errorline[1]))
                self.ax.axhline(y=self.errorline[0], linestyle='-.',
                                linewidth=1.55, color='crimson')
                self.ax.fill_between(dataT, dataY, self.errorline[1],
                                    where=[True if y>self.errorline[1]
                                            else False for y in dataY],
                                     color='y')
                self.ax.fill_between(dataT, dataY, self.errorline[0],
                                    where=[True if y<self.errorline[0]
                                            else False for y in dataY],
                                     color='y')

            # Define an area of the plot to be shaded
            if len(self.zone):
                for zs in self.zone:
                    self.ax.axvspan(zs[0][0][0], zs[0][0][1],
                                    facecolor=zs[0][0][2], alpha=0.15,
                                    label=zs[1])
                    for oz in zs[0][1:]:
                        self.ax.axvspan(oz[0], oz[1], facecolor=oz[2],
                                        alpha=0.15)
                self.ax.legend(loc='upper right',
                               bbox_to_anchor=(1, 1),
                               fontsize = 'small')

            if not(self.bottomax):
                plt.setp(self.ax.get_xticklabels(), fontsize=9, visible=False)
            else:
                self.ax.grid(True)
                self.ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
                self.ax.xaxis.set_major_formatter(
                    matplotlib.dates.DateFormatter("%d/%m %H:%M:%S.%f"))
                self.ax.xaxis.set_minor_locator(ticker.MaxNLocator(200))
                plt.setp(self.ax.get_xticklabels(), fontsize=8,
                         rotation=25, ha='right')
        # Configure an histogram ax
        else:
            self.bindata, self.bins, self.patches = \
                self.ax.hist(dataX, bins=self.histbins, label=self.label,
                             range=self.limsbins, edgecolor='black',
                             alpha=self.alpha, color=self.color)

            if self.bottomax:
                self.ax.set_xlabel(self.xlabel, fontsize=9)
            plt.xticks(self.bins)
            plt.setp(self.ax.get_xticklabels(), fontsize=8,
                    rotation=35, ha='right')

        # Define the plot area style for the ax
        self.ax.grid(True)
        self.ax.tick_params("y", colors="k")
        self.ax.set_ylabel(self.ylabel, color="k", fontsize=9)
        if self.ylims:
            self.ax.set_ylim(self.ylims[0], self.ylims[1])
        self.ax.legend(loc='upper right',
                       bbox_to_anchor=(1, 1),
                       fontsize='small')
        print('Setup done')

    def set_ax(self, gs, tcells, rows, masterax):
        if self.shax:
            self.ax = self.shax.ax
            self.ylabel = self.shax.ylabel
            self.bottomax = self.shax.bottomax
        else:
            ys = int(self.ys * (tcells/rows))
            ye = int(self.ye * (tcells/rows))
            self.ax = plt.subplot(gs[ys:ye,self.xs:self.xe],
                                  sharex=masterax.ax)

# axPlt object end

class DataAxePlotter:
    def __init__(self, ncols=1):
        self.Axe = collections.OrderedDict()
        for n in range(ncols):
            self.Axe['c'+str(n+1)] = collections.OrderedDict()
        self.masterax = None
        self.gs = None

    def positionPlot(self):
        # This method starts by constructing a dictionary for the position of
        # each plot within a matrix. The idea is to be able to calculate the
        # proper coordinates and dimentions for each plot
        # The dictionary contains a subdict for each column as well as a key
        # with the total number of columns for the plot
        Dax = self.Axe
        ric = {}
        i = 0
        for c in Dax:
            j = 0
            for n in Dax[c]:
                shaxname = Dax[c][n].shaxname
                if not(shaxname):
                    Dax[c][n].xs = i
                    Dax[c][n].xe = i + 1
                    Dax[c][n].ys = j
                    Dax[c][n].ye = Dax[c][n].height + j
                    j = Dax[c][n].ye
                    if not(self.masterax):
                        self.masterax = Dax[c][n]
                    bax = n
                else:
                    Dax[c][n].shax = Dax[c][shaxname]
            ric[c] = j
            Dax[c][bax].bottomax = True
            i += 1
        totCells = lcmArray([ric[c] for c in ric])
        # print('totCells', totCells)
        # print('totCols', i)
        self.gs = gridspec.GridSpec(totCells,i)
        for c in Dax:
            for n in Dax[c]:
                print('Setting up plot {}'.format(n))
                Dax[c][n].plot_ax(self.gs, totCells, ric[c], self.masterax)

    def plotConfig(self, title=None, xlabels=None):
        '''
        This function is used to:
            - Disable tick labels for plots other than the bottom positions
            - Set the date labels on the bottom plots
            - Update Plot Event Zones
        '''
        Dax = self.Axe
        if xlabels:
            for c in Dax:
                for n in Dax[c]:
                    if Dax[c][n].bottomax:
                        Dax[c][n].ax.set_xlabel(xlabels[c])
        plt.suptitle(title)
        plt.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.95,
                            hspace=0.20, wspace=0.2)
        plt.show()

if __name__ == '__main__':
    x = np.arange(21)
    y1 = x**2
    y2 = 50 * (np.sin((3.1415/6) * x))
    y3 = 500 / (x + 1)
    plts = DataAxePlotter(2)
    plts.Axe['c1']['g4'] = DataAx([x,y1], 'r', ylabel='t1', height=5, label='g4')
    plts.Axe['c1']['g6'] = DataAx([x,y3], 'k', ylabel='t3', label='g6', height=1, shax='g4')
    plts.Axe['c1']['g5'] = DataAx([x,y2], 'b', linestyle='--', label='g5', shax='g4')
    plts.Axe['c2']['g2'] = DataAx([x,((x-10)**2)+40], 'g', ylabel='t4', label='g2', height=2)
    plts.positionPlot()
    plts.plotConfig()

