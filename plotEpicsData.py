#!/usr/bin/env python3.5

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec

from matplotlib import dates
from datetime import datetime, timedelta

import collections
import numpy as np
import h5py

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
        - color: Line color
        - lwidth: Line width
        - label: Line label, displayed on legend box
        - ylabel: Y axis label
        - ylims: Y axis limits
        - alpha: Line transparency
        - errorline: Y coord for errorline
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
          Adds horizontal lines (to visualize limits.

    '''

    def __init__(self, data, linestyle='False', ylabel='', ylims=[], label=None,
                 drawstyle=None, shax=None, height=1,
                 masterax=False, errorline=[], zone=[], alpha=1.0,
                 linewidth=1.25, histbins=False, limsbins=None, color=None):
        '''
        Creates the plot object
        '''
        # Axes attributes
        self.ylabel = ylabel
        self.ylims = ylims
        self.ax = None
        self.errorline = errorline
        self.zone = zone
        self.height = height
        self.mstax = masterax
        self.ys = None
        self.ye = None
        self.xs = None
        self.xe = None
        self.bottomax = False
        self.shaxname = shax
        self.shax = None
        # Plot attributes
        self.data = data
        self.label = label
        self.linestyle = linestyle
        self.label = label
        self.drawstyle = drawstyle
        self.color = color
        self.lwidth = linewidth
        self.alpha = alpha
        self.histbins = histbins
        self.limsbins = limsbins
        self.bindata = False
        self.bins = False
        self.patches = False


    def plot_ax(self, gs, masterax=None):
        ''' plot_ax

        This method defines each axes to be plotted. It will prepare each ax
        depending if is a regular plot or a histogram
        It places each ax within the gridspace, and defines all its
        characteristics
        '''

        # Check if plot is a histogram. If it's not, prepares the x & y axis for
        # plotting. If it is prepares an array with the x values.
        if not(self.histbins):
            dataT, dataY = self.data
        else:
            dataX = self.data

        # Checks if the first plot has been defined, if not, plot the current
        # data set and set the shareFlag for the rest of the plots.
        if not(self.mstax):
            # Defines the plot box and main axis with which X-axis is shared
            print('creating', self.label)
            self.ax = plt.subplot(gs[self.ys:self.ye,self.xs:self.xe],
                                  sharex=masterax.ax)
            print('ax created')
        else:
            # Defines the plot box
            print('creating', self.label)
            self.ax = plt.subplot(gs[self.ys:self.ye,self.xs:self.xe])
            print('ax created')

        # Configures each axes for plotting
        if not(self.histbins):
            # Check to see if drawstyle has been defined. This is used when
            # plotting something like a boolean plot
            if not(self.drawstyle):
                self.ax.plot(dataT, dataY, self.linestyle, linewidth=self.lwidth,
                            label=self.label, alpha=self.alpha)
            else:
                self.ax.plot(dataT, dataY, c=self.linestyle, linewidth=self.lwidth,
                            drawstyle=self.drawstyle, label=self.label,
                            alpha=self.alpha)

            self.ax.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
            # Define a two horizontal line in the ax, if errorline has been set
            if self.errorline:
                self.ax.axhline(y=self.errorline[1], linestyle='-.',
                                linewidth=1.55, color='crimson',
                                label='|error|={0}'.format(self.errorline[1]))
                self.ax.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
                self.ax.axhline(y=self.errorline[0], linestyle='-.',
                                linewidth=1.55, color='crimson')
                self.ax.fill_between(dataT, dataY, self.errorline[1],
                                    where=[True if y>self.errorline[1]
                                            else False for y in dataY], color='y')
                self.ax.fill_between(dataT, dataY, self.errorline[0],
                                    where=[True if y<self.errorline[0]
                                            else False for y in dataY], color='y')

            # Define an area of the plot to be shaded
            if len(self.zone):
                for zs in self.zone:
                    # print(zs)
                    self.ax.axvspan(zs[0][0][0], zs[0][0][1],
                                    facecolor=zs[0][0][2], alpha=0.15,
                                    label=zs[1])
                    for oz in zs[0][1:]:
                        self.ax.axvspan(oz[0], oz[1], facecolor=oz[2], alpha=0.15)
                self.ax.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
        # Configure an histogram ax
        else:
            self.bindata, self.bins, self.patches = \
                self.ax.hist(dataX, bins=self.histbins, label=self.label,
                             range=self.limsbins, edgecolor='black',
                             alpha=self.alpha, color=self.color)
            self.ax.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
            # print(self.bins)

        # Define the plot area style for the ax
        self.ax.grid(True)
        self.ax.tick_params("y", colors="b")
        self.ax.set_ylabel(self.ylabel, color="b")
        if self.ylims:
            self.ax.set_ylim(self.ylims[0], self.ylims[1])

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

class DataAxesPlotter:
    def __init__(self, ncols=1):
        # self.dictAxes = Axes
        self.Axes = collections.OrderedDict()
        for n in ncols:
            self.Axes['c'+str(n+1)] = collections.OrderedDict()
        self.masterax = None
        self.bottomax = []
        self.gs = None
        self.fig = plt.figure()
        self.maFlag = False

    def positionPlot(self):
        # This method starts by constructing a dictionary for the position of
        # each plot within a matrix. The idea is to be able to calculate the
        # proper coordinates and dimentions for each plot
        # The dictionary contains a subdict for each column as well as a key
        # with the total number of columns for the plot
        Dax = self.Axes
        Dshax = {}
        ric = {}
        i = 0
        for c in Dax:
            Dshax[c] = {}
            j = 0
            for n in Dax[c]:
                shaxname = Dax[c][n].shaxname
                if not(shaxname):
                    Dax[c][n].xs = i
                    Dax[c][n].xe = i + 1
                    Dax[c][n].ys = j
                    Dax[c][n].ye = Dax[c][n].height + j
                    j = Dax[c][n].ye
                    # if not(self.maFlag):
                    if not(self.masterax):
                        self.masterax = Dax[c][n]
                        Dax[c][n].mstax = True
                        self.maFlag = True
                    bax = n
                else:
                    Dax[c][n].shax = Dax[c][shaxname]
                    # if not(Dax[c][n].sharedax in Dshax[c].keys()):
                        # Dshax[c][Dax[c][n].sharedax] = []
                    # Dshax[c][Dax[c][n].sharedax] += [n]
            ric[c] = j
            Dax[c][bax].bottomax = True
            i += 1
        totCells = lcmArray([ric[c] for c in ric])
        print('totCells', totCells)
        print('totCols', i)
        self.gs = gridspec.GridSpec(totCells,i)
        for c in Dshax:
            for n in Dshax[c]:
                for a in Dshax[c][n]:
                    Dax[c][a].xs = Dax[c][n].xs
                    Dax[c][a].xe = Dax[c][n].xe
                    Dax[c][a].ys = Dax[c][n].ys
                    Dax[c][a].ye = Dax[c][n].ye
                    Dax[c][a].ylabel = Dax[c][n].ylabel
                    Dax[c][a].bottomax = Dax[c][n].bottomax
                    Dax[c][a].mstax = Dax[c][n].mstax

        for c in Dax:
            for n in Dax[c]:
                Dax[c][n].ys = int(Dax[c][n].ys * (totCells/ric[c]))
                Dax[c][n].ye = int(Dax[c][n].ye * (totCells/ric[c]))

        Dax[self.masterax[0]][self.masterax[1]].plot_ax(self.gs)
        for c in Dax:
            for n in Dax[c].keys()-[self.masterax[1]]:
                # if not(Dax[c][n].histbins):
                Dax[c][n].plot_ax(self.gs,
                                Dax[self.masterax[0]][self.masterax[1]])

    def coordPlot(self):
        # This method starts by constructing a dictionary for the position of
        # each plot within a matrix. The idea is to be able to calculate the
        # proper coordinates and dimentions for each plot
        # The dictionary contains a subdict for each column as well as a key
        # with the total number of columns for the plot
        Dax = self.dictAxes
        Da = {'nCols':0} # Initialize dictionary
        Dshax = {}
        for n in Dax:
            if not(Dax[n].sharedax):
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
            else:
                if not(Dax[n].sharedax in Dshax.keys()):
                    Dshax[Dax[n].sharedax] = []
                Dshax[Dax[n].sharedax] += [n]
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
                    xs = i - 1
                    xe = i
                    ys = int(Da[str(i)][str(j)][0][0]
                             * (totCells/Da[str(i)]['nCells']))
                    ye = int(Da[str(i)][str(j)][0][1]
                             * (totCells/Da[str(i)]['nCells']))
                    Dax[plotname].xs = xs
                    Dax[plotname].xe = xe
                    Dax[plotname].ys = ys
                    Dax[plotname].ye = ye
                    # Define the bottom plot for each column
                    if j == Da[str(i)]['nRows']:
                        Dax[plotname].bottomax = True
                    if not(self.maFlag):
                        self.masterax = plotname
                        Dax[plotname].mstax = True
                        self.maFlag = True
                if (j < Da[str(i)]['nRows']):
                    # calculate the starting and ending X-coords
                    ysp = Da[str(i)][str(j)][0][1]
                    yep = Da[str(i)][str(j+1)][0][1] + ysp
                    Da[str(i)][str(j+1)][0] =  [ysp, yep]
        for n in Dshax:
            for a in Dshax[n]:
                Dax[a].xs = Dax[n].xs
                Dax[a].xe = Dax[n].xe
                Dax[a].ys = Dax[n].ys
                Dax[a].ye = Dax[n].ye
                Dax[a].bottomax = Dax[n].bottomax
                Dax[a].mstax = Dax[n].mstax
        Dax[self.masterax[0]][self.masterax[1]].plot_ax(self.gs)
        for c in Dax:
            for n in Dax[c].keys()-[self.masterax[1]]:
                Dax[c][n].plot_ax(self.gs,
                                  Dax[self.masterax[0]][self.masterax[1]])
        self.dictAxes = {'s':Dax}

    def plotConfig(self, title=None, xlabels=None):
        '''
        This function is used to:
            - Disable tick labels for plots other than the bottom positions
            - Set the date labels on the bottom plots
            - Update Plot Event Zones
        '''
        Dax = self.dictAxes
        for c in Dax:
            for n in Dax[c]:
                # Dax[c][n].ax.grid(True)
                # Dax[c][n].ax.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
                # ax[0].ax.xaxis_date()
                if not(Dax[c][n].bottomax):
                    plt.setp(Dax[c][n].ax.get_xticklabels(), fontsize=9, visible=False)
                else:
                    if not(Dax[c][n].histbins):
                        Dax[c][n].ax.grid(True)
                        Dax[c][n].ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
                        Dax[c][n].ax.xaxis.set_major_formatter(
                            matplotlib.dates.DateFormatter("%d/%m %H:%M:%S.%f"))
                        Dax[c][n].ax.xaxis.set_minor_locator(ticker.MaxNLocator(200))
                        plt.setp(Dax[c][n].ax.get_xticklabels(), fontsize=9,
                                rotation=30, ha='right')
                        # # ax[0].ax.set_xlim(lowX,highX)
                        # ax[0].ax.set_xlim(time[5],time[len(time)-5])
                        # h, l = ax[0].ax.get_legend_handles_labels()
                        # ax[0].ax.legend(h, l, loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
                    else:
                        plt.xticks(Dax[c][n].bins)
                        # Dax[c][n].ax.xaxis.set_major_locator(ticker.MaxNLocator(Dax[c][n].histbins))
                        # Dax[c][n].ax.xaxis.set_major_locator(
                            # ticker.MaxNLocator(nbins=11, symmetric=True))
                        # Dax[c][n].ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=2))
                        plt.setp(Dax[c][n].ax.get_xticklabels(), fontsize=9,
                                rotation=30, ha='right')
                        if not(xlabels):
                            continue
                        Dax[c][n].ax.set_xlabel(xlabels[c])
        plt.suptitle(title)
        plt.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.95,
                            hspace=0.20, wspace=0.2)
        plt.show()

    @staticmethod
    def plot_attributes(data, linestyle='False', ylabel='', ylims=[], label=None,
                 drawstyle=None, shax=None, height=1,
                 masterax=False, errorline=[], zone=[], alpha=1.0,
                 linewidth=1.25, histbins=False, limsbins=None, color=None):
        '''
        Creates the plot object
        '''
        return {
        # Axes attributes
        'ylabel':ylabel,
        'ylims':ylims,
        'ax':None,
        'errorline':errorline,
        'zone':zone,
        'height':height,
        'mstax':masterax,
        'ys':None,
        'ye':None,
        'xs':None,
        'xe':None,
        'bottomax':False,
        'sharedax':shax,
        # Plot attributes
        'plots':[],
        'data':data,
        'label':label,
        'linestyle':linestyle,
        'label':label,
        'drawstyle':drawstyle,
        'color':color,
        'lwidth':linewidth,
        'alpha':alpha,
        'histbins':histbins,
        'limsbins':limsbins,
        'bindata':False,
        'bins':False,
        'patches':False
        }


if __name__ == '__main__':
    # plotConfig(ax_lst)
    x = np.arange(21)
    y1 = x**2
    y2 = 50 * (np.sin((3.1415/6) * x))
    y3 = 500 / (x + 1)
    # plts = collections.OrderedDict()
    plts = collections.OrderedDict()
    pC1 = collections.OrderedDict()
    pC2 = collections.OrderedDict()
    # plts['g1'] = DataAx('', '', '', pos=[2,1], height=1)
    # plts['g2'] = DataAx('', '', '', pos=[1,2], height=2)
    # plts['g3'] = DataAx('', '', '', pos=[3,1], height=4)
    pC1['g4'] = DataAx([x,y1], 'r-', 't1', height=5, label='g4')
    pC1['g5'] = DataAx([x,y2], 'b--', label='g5', shax='g4')
    pC2['g2'] = DataAx([x,((x-10)**2)+40], 'g-', 't4', label='g2', height=2)
    pC2['g6'] = DataAx([x,y3], 'k-', 't3', label='g6', height=1)
    # plts['c1']['g4'] = DataAx([x,y1], 'r-', 't1', shax='g2', pos=[1,2], height=2, label='g4')
    # plts['c1']['g5'] = DataAx([x,y2], 'b--', label='g5', shax='g2')
    # plts['c2']['g2'] = DataAx([x,((x-10)**2)+40], 'g-', 't4', pos=[1,1], label='g2', height=2)
    # plts['c2']['g6'] = DataAx([x,y3], 'k-', 't3', shax='g2', pos=[2,2], label='g6', height=1)
    plts['c1'] = pC1
    plts['c2'] = pC2
    pltAx = DataAxesPlotter(plts)
    pltAx.positionPlot()
    # for plts in [pltC1, pltC2]:
        # for p in plts:
            # print('Plot {2} Span in col: {0} - {1}'.format(plts[p].ys,
                                                           # plts[p].ye, p))
    pltAx.plotConfig()

    # plt.suptitle('In Position Characterization of Gemini South Telescope Systems')
    # plt.show()
