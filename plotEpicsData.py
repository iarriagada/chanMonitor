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
    print(offsetZonesArray[1:10])
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

    def __init__(self, data, linestyle, ylabel='', ylims=[], label=None,
                 drawstyle=None, shax=None, pos=[1,1], height=1,
                 masterax=False, errorline=[], zone=[], alpha=1.0,
                 linewidth=1.25):
        '''
        Creates the plot object
        '''
        self.linestyle = linestyle
        self.drawstyle = drawstyle
        self.lwidth = linewidth
        self.label = label
        self.ylabel = ylabel
        # self.zlabel = zlabel
        self.ylims = ylims
        self.ax = None
        self.alpha = alpha
        self.data = data
        self.errorline = errorline
        self.zone = zone
        self.pos = pos
        self.height = height
        self.mstax = masterax
        self.ys = None
        self.ye = None
        self.xs = None
        self.xe = None
        self.bottomax = False
        self.sharedax = shax


    def plot_ax(self, gs, masterax=None):
        '''
        Plot data in corresponding axes
        '''
        # dataT, dataY = zip(*self.data)
        dataT, dataY = self.data

        # Checks if the first plot has been defined, if not, plot the current
        # data set and set the shareFlag for the rest of the plots.
        # if masterax and not(self.mstax) and not(masterax.pos == self.pos):
        if not(self.mstax):
            # Defines the plot box and main axis with which Y is shared
            self.ax = plt.subplot(gs[self.ys:self.ye,self.xs:self.xe],
                                  sharex=masterax.ax)
        else:
            # Defines the plot box
            self.ax = plt.subplot(gs[self.ys:self.ye,self.xs:self.xe])

        # Plot the data
        if not(self.drawstyle):
            self.ax.plot(dataT, dataY, self.linestyle, linewidth=self.lwidth,
                         label=self.label, alpha=self.alpha)
        else:
            self.ax.plot(dataT, dataY, c=self.linestyle, linewidth=self.lwidth,
                         drawstyle=self.drawstyle, label=self.label,
                         alpha=self.alpha)

        if self.errorline:
            self.ax.axhline(y=self.errorline[1], linestyle='-.',
                            linewidth=1.55, color='crimson',
                            label='|error|={0}'.format(self.errorline[1]))
            self.ax.axhline(y=self.errorline[0], linestyle='-.',
                            linewidth=1.55, color='crimson')
            self.ax.fill_between(dataT, dataY, self.errorline[1],
                                 where=[True if y>self.errorline[1]
                                        else False for y in dataY], color='y')
            self.ax.fill_between(dataT, dataY, self.errorline[0],
                                 where=[True if y<self.errorline[0]
                                        else False for y in dataY], color='y')

        if len(self.zone):
            for zs in self.zone:
                # print(zs)
                self.ax.axvspan(zs[0][0][0], zs[0][0][1],
                                facecolor=zs[0][0][2], alpha=0.15,
                                label=zs[1])
                for oz in zs[0][1:]:
                    self.ax.axvspan(oz[0], oz[1], facecolor=oz[2], alpha=0.15)

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
        Dax = self.dictAxes
        Dshax = {}
        ric = {}
        i = 0
        for c in Dax:
            Dshax[c] = {}
            j = 0
            for n in Dax[c]:
                if not(Dax[c][n].sharedax):
                    Dax[c][n].xs = i
                    Dax[c][n].xe = i + 1
                    Dax[c][n].ys = j
                    Dax[c][n].ye = Dax[c][n].height + j
                    j = Dax[c][n].ye
                    if not(self.maFlag):
                        self.masterax = [c,n]
                        Dax[c][n].mstax = True
                        self.maFlag = True
                    bax = n
                else:
                    if not(Dax[c][n].sharedax in Dshax[c].keys()):
                        Dshax[c][Dax[c][n].sharedax] = []
                    Dshax[c][Dax[c][n].sharedax] += [n]
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
                    # if Dax[plotname].mstax:
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

    def plotConfig(self):
        '''
        This function is used to:
            - Disable tick labels for plots other than the bottom positions
            - Set the date labels on the bottom plots
            - Update Plot Event Zones
        '''
        Dax = self.dictAxes
        for c in Dax:
            for n in Dax[c]:
                Dax[c][n].ax.grid(True)
                Dax[c][n].ax.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize = 'small')
                # ax[0].ax.xaxis_date()
                if not(Dax[c][n].bottomax):
                    plt.setp(Dax[c][n].ax.get_xticklabels(), fontsize=9, visible=False)
                else:
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
        plt.subplots_adjust(top=0.90, bottom=0.08, left=0.06, right=0.95,
                            hspace=0.20, wspace=0.15)
        plt.show()



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
