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

def diffData(dataA):
    diffPV = [y - x \
               for x,y in np.array([dataA[1][:-1],
                                    dataA[1][1:]]).T]
    diffA = [dataA[0][1:], diffPV]
    return diffA

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

    def plot_ax(self, ax_lst, position, height, width, data,
                bottomPlot=False, errorLine=False, zonesPlot=True,
                ds=None, plotLabel=None, alphaT=1.0):
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
            self.ax.axhline(y=1.0, linestyle='-.', linewidth=1.25, color='crimson')
            self.ax.axhline(y=0.5, linestyle='-.', linewidth=1.25, color='darkorange')
            self.ax.axhline(y=-1.0, linestyle='-.', linewidth=1.25, color='crimson')
            self.ax.axhline(y=-0.5, linestyle='-.', linewidth=1.25, color='darkorange')
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

    posFile = h5py.File(args.hdf5File, 'r')
    recGroups = [[posFile.get(recname),recname] for recname in posFile]
    recData = {name:[np.array(g.get('timestamp')),np.array(g.get('value'))]\
               for g,name in recGroups}
    for n in recData:
        recTime = [datetime.fromtimestamp(ts) for ts in recData[n][0]]
        timeStampS = [datetime.strftime(rt, '%m/%d/%Y %H:%M:%S.%f')\
                      for rt in recTime]
        timeStamp = [datetime.strptime(ts, '%m/%d/%Y %H:%M:%S.%f')\
                     for ts in timeStampS]
        recData[n][0] = timeStamp

    azdP = diffData(recData['mc:azCurrentPos'])
    azdD = diffData(recData['tcs:demandAz'])
    eldP = diffData(recData['mc:elCurrentPos'])
    eldD = diffData(recData['tcs:demandEl'])

    ylims = (-1, 7)
    ylimsOFF = (-1, 2)
    ax_lst = list()
    AxPlt.plotArea = (6,2)
    AxPlt.gs = gridspec.GridSpec(AxPlt.plotArea[0], AxPlt.plotArea[1])
    print('size az dmd array: ', len(recData['tcs:demandAz'][1]),
          ' size az pos array ', len(recData['mc:azCurrentPos'][1]))

    azPlot = AxPlt('b-', "MCS Azimuth Position & Demand\n[deg]", False)
    azDPlot = AxPlt('b-', "MCS Azimuth Position Differential\n[deg]", False)
    # azDmd = AxPlt('b.', "MCS Azimuth Demand\n[ms]", ylims, False)
    elPlot = AxPlt('g-', "MCS Elevation Position & Demand\n[deg]", False)
    elDPlot = AxPlt('g-', "MCS Elevation Position Differential\n[deg]", False)
    # elDmd = AxPlt('b', "MCS Elevation Demand\n[ms]", ylims, False)
    offsPlot = AxPlt('c', "PQ offset\n[s]", False)
    offsPlot2 = AxPlt('c', "PQ offset\n[s]", False)
    inposPlot = AxPlt('k', "MCS in position\n[s]", ylimsOFF, False)
    inposPlot2 = AxPlt('k', "MCS in position\n[s]", ylimsOFF, False)

    azPlot.plot_ax(ax_lst, (0,0), 2, 1, np.array(recData['mc:azCurrentPos']).T,
                   plotLabel='Az Pos')
    azPlot.add_plot(ax_lst, np.array(recData['tcs:demandAz']).T, ls='r-',
                    plotLabel='Az Dmd', alphaT=0.40)
    azDPlot.plot_ax(ax_lst, (2,0), 2, 1, np.array(azdP).T,
                   plotLabel='Az Pos')
    azDPlot.add_plot(ax_lst, np.array(azdD).T, ls='r-',
                    plotLabel='Az Dmd', alphaT=0.40)
    elPlot.plot_ax(ax_lst, (0,1), 2, 1, np.array(recData['mc:elCurrentPos']).T,
                   plotLabel='El Pos')
    elPlot.add_plot(ax_lst, np.array(recData['tcs:demandEl']).T, ls='r-',
                    plotLabel='El Dmd', alphaT=0.40)
    elDPlot.plot_ax(ax_lst, (2,1), 2, 1, np.array(eldP).T,
                   plotLabel='El Pos')
    elDPlot.add_plot(ax_lst, np.array(eldD).T, ls='r-',
                    plotLabel='El Dmd', alphaT=0.40)
    # tdLA1.add_plot(ax_lst, diffMCS, ls='g.', plotLabel='MCS-LA1', alphaT=0.30)
    # tdLA1.add_plot(ax_lst, diffSCS, ls='y.', plotLabel='SCS-LA1', alphaT=0.20)
    # tdLA1.add_plot(ax_lst, diffCRCS, ls='c.', plotLabel='CRCS-LA1', alphaT=0.10)

    offsPlot.plot_ax(ax_lst, (4,0), 1, 1,
                     np.array(recData['tcs:offsetPoA1.VALA']).T,
                     plotLabel='P Offset', ds='steps-post')
    offsPlot.add_plot(ax_lst,
                      np.array(recData['tcs:offsetPoA1.VALB']).T,
                      ls='r', plotLabel='Q Offset',alphaT=0.40, ds='steps-post')
    offsPlot2.plot_ax(ax_lst, (4,1), 1, 1,
                     np.array(recData['tcs:offsetPoA1.VALA']).T,
                     plotLabel='P Offset', ds='steps-post')
    offsPlot2.add_plot(ax_lst,
                      np.array(recData['tcs:offsetPoA1.VALB']).T,
                      ls='r', plotLabel='Q Offset',alphaT=0.40, ds='steps-post')

    inposPlot.plot_ax(ax_lst, (5,0), 1, 1,
                      np.array(recData['mc:inPosition']).T,
                      bottomPlot=True,
                      plotLabel='InPos', ds='steps-post')
    inposPlot2.plot_ax(ax_lst, (5,1), 1, 1,
                      np.array(recData['mc:inPosition']).T,
                      bottomPlot=True,
                      plotLabel='InPos', ds='steps-post')
    # elPlot.plot_ax(ax_lst, (2,0), 2, 2, np.array(recData['mc:elCurrentPos']).T, plotLabel='El Pos')
    # elPlot.add_plot(ax_lst, np.array(recData['tcs:demandEl']).T, ls='r-', plotLabel='El Dmd',alphaT=0.40)
    # tdTC1.plot_ax(ax_lst, (2,0), 1, 2, diffTC1, plotLabel='TC1-LA1')
    # deltaTLA1.plot_ax(ax_lst, (3,0), 1, 2, dtLA1, plotLabel='delta LA1')

    # rmsTDLA1.plot_ax(ax_lst, (4,0), 1, 2, rmsTC1, bottomPlot=True, ds='steps-post')
    # rmsTDLA1.add_plot(ax_lst, rmsTCS, ls='r', bottomPlot=True, ds='steps-post', alphaT=0.90)
    # rmsTDLA1.add_plot(ax_lst, rmsMCS, ls='g', bottomPlot=True, ds='steps-post', alphaT=0.70)
    # rmsTDLA1.add_plot(ax_lst, rmsSCS, ls='y', bottomPlot=True, ds='steps-post', alphaT=0.50)
    # rmsTDLA1.add_plot(ax_lst, rmsCRCS, ls='c', bottomPlot=True, ds='steps-post', alphaT=0.30)

    plotConfig(ax_lst)

    plt.suptitle('Time synchronization performance of Gemini Subsystems')

    # tsDTC1, dTC1 = zip(*diffTC1)
    # tsDTCS, dTCS = zip(*diffTCS)
    # tsDMCS, dMCS = zip(*diffMCS)
    # tsDSCS, dSCS = zip(*diffSCS)
    # tsDCRCS, dCRCS = zip(*diffCRCS)

    # fig, ax1 = plt.subplots()
    # plt.title("pvload+pvsave time on simple soft-ioc running on sbfrtdev-lv1 and r/w to cportprd-lv1 1Hz")
    # ax1.plot(tsDTC1, dTC1, "b.", label=dataLabels[1])
    # ax1.plot(tsDTCS, dTCS, "r.", label=dataLabels[2])
    # ax1.plot(tsDMCS, dMCS, "g.", label=dataLabels[3])
    # ax1.plot(tsDSCS, dSCS, "y.", label=dataLabels[4])
    # ax1.plot(tsDCRCS, dCRCS, "k.", label=dataLabels[5])
    # ax1.grid(True)
    # ax1.set_ylabel("Milliseconds")
    # ax1.set_ylim(-1, 7)
    # plt.gcf().autofmt_xdate()

    # plt.legend()
    plt.show()
