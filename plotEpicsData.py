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
                 drawstyle='default', label=None, ylabel='', rawx=False,
                 xlabel='', ylims=[], shax=None, height=1, marksize=5,
                 errorline=[], zone={}, alpha=1.0, ticklabels=False,
                 standalone=False, linewidth=1.25, histbins=False,
                 limsbins=None, timezone=None):
        '''
        Creates the plot object
        '''
        # Axe attributes
        self.ax = None
        self.shaxname = shax
        self.shax = None
        self.mstax = None
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
        self.rawx = rawx
        self.standalone = standalone
        # Plot attributes
        self.data = data
        self.color = color
        self.linestyle = linestyle
        self.marker = marker
        self.marksize = marksize
        self.drawstyle = drawstyle
        self.lwidth = linewidth
        self.alpha = alpha
        self.label = label
        self.histbins = histbins
        self.limsbins = limsbins
        self.ticklabels = ticklabels
        self.bindata = False
        self.bins = False
        self.patches = False
        self.timezone = timezone


    def plot_ax(self, gs, tcells, rows, masterax):
        ''' plot_ax

        This method defines each pair of axe to be plotted. It will prepare each
        axe set depending if is a regular plot or a histogram.
        It places each axe pair within the gridspace, and set all its
        characteristics
        '''
        shared_x_axis = None
        if not(self.histbins) and not(self.standalone):
            shared_x_axis = masterax.ax
        if self.mstax:
            shared_x_axis = self.mstax.ax

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
            if not(self.rawx):
                dataT = [datetime.fromtimestamp(ts, tz=self.timezone) for ts in dataT]
        if self.histbins:
            dataX = self.data

        # Checks if the first plot has been defined, if not, plot the current
        # data set and set the shareFlag for the rest of the plots.

        # Configures each axe for plotting
        if not(self.histbins):
            l = self.ax.plot(dataT, dataY, linestyle=self.linestyle,
                         linewidth=self.lwidth, color=self.color,
                         marker=self.marker, drawstyle=self.drawstyle,
                         label=self.label, alpha=self.alpha,
                         markersize=self.marksize)
            if not(self.linestyle):
                l[0].set_markerfacecolor(
                    matplotlib.colors.to_rgba(self.color, self.alpha))
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
            if self.zone:
                for k,zs in self.zone['zones'].items():

                    zn_start = datetime.strptime(zs['timestamps'][0][0],
                                                 '%y%m%dT%H%M')
                    zn_end = datetime.strptime(zs['timestamps'][0][1],
                                                 '%y%m%dT%H%M')
                    self.ax.axvspan(zn_start, zn_end,
                                    facecolor=zs['color'], alpha=0.15,
                                    label=zs['label'])
                    for oz in zs['timestamps'][1:]:
                        oz_start = datetime.strptime(oz[0], '%y%m%dT%H%M')
                        oz_end = datetime.strptime(oz[1], '%y%m%dT%H%M')
                        self.ax.axvspan(oz_start, oz_end, facecolor=zs['color'],
                                        alpha=0.15)
                self.ax.legend(loc='upper right',
                               bbox_to_anchor=(1, 1),
                               fontsize = 'small')

            # if len(self.zone):
                # for zs in self.zone:
                    # zn_start = datetime.fromtimestamp(zs[0][0][0], tz=self.timezone)
                    # zn_end = datetime.fromtimestamp(zs[0][0][1], tz=self.timezone)
                    # self.ax.axvspan(zn_start, zn_end,
                                    # facecolor=zs[0][0][2], alpha=0.15,
                                    # label=zs[1])
                    # for oz in zs[0][1:]:
                        # oz_start = datetime.fromtimestamp(oz[0], tz=self.timezone)
                        # oz_end = datetime.fromtimestamp(oz[1], tz=self.timezone)
                        # self.ax.axvspan(oz_start, oz_end, facecolor=oz[2],
                                        # alpha=0.15)
                # self.ax.legend(loc='upper right',
                               # bbox_to_anchor=(1, 1),
                               # fontsize = 'small')

            if self.bottomax and not(self.rawx):
                self.ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
                self.ax.xaxis.set_major_formatter(
                    matplotlib.dates.DateFormatter("%d/%m %H:%M:%S.%f"))
                self.ax.xaxis.set_minor_locator(ticker.MaxNLocator(200))
            if self.bottomax:
                plt.setp(self.ax.get_xticklabels(), fontsize=8,
                         rotation=25, ha='right')
            if not(self.bottomax):
                plt.setp(self.ax.get_xticklabels(), fontsize=9, visible=False)
        # Configure an histogram ax
        else:
            self.bindata, self.bins, self.patches = \
                self.ax.hist(dataX, bins=self.histbins, label=self.label,
                             range=self.limsbins, edgecolor='black',
                             alpha=self.alpha, color=self.color)

            if not(self.bottomax) and not(self.ticklabels):
                plt.setp(self.ax.get_xticklabels(), fontsize=8, visible=False)
            if self.bottomax or self.ticklabels:
                print('Bottom ax!!!')
                plt.setp(self.ax.get_xticklabels(), fontsize=8,
                         rotation=45, ha='right', visible=True)

            plt.xticks(self.bins)
        # Set x label if defined by user
        if self.bottomax:
            self.ax.set_xlabel(self.xlabel, color='k', fontsize=9, visible=True)
            print('setting x label!')

        # Define the plot area style for the ax
        self.ax.grid(True)
        self.ax.tick_params("y", colors="k")
        self.ax.set_ylabel(self.ylabel, color="k", fontsize=9)
        if self.ylims:
            self.ax.set_ylim(**self.ylims)
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

    @staticmethod
    def update_axe(data_axe, **kwargs):
        for attr in kwargs:
            setattr(data_axe, attr, kwargs[attr])
        return data_axe

# DataAx class end

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
            - Set the date labels on the bottom plots
        '''
        Dax = self.Axe
        if xlabels:
            for c in Dax:
                for n in Dax[c]:
                    if Dax[c][n].bottomax:
                        Dax[c][n].ax.set_xlabel(xlabels[c], fontsize=9)
        plt.suptitle(title)
        plt.subplots_adjust(top=0.95, bottom=0.1, left=0.085, right=0.95,
                            hspace=0.225, wspace=0.225)
        plt.show()

# DataAxePlotter class end

def extract_hdf5(hdf5File, start_time=None,
                 end_time=None, channel_mask=None):
    # Handle the start and end times
    if start_time:
        try:
            stime_dt = datetime.strptime(start_time, '%y%m%dT%H%M')
            stime = datetime.timestamp(stime_dt)

        except ValueError as err:
            print("ValueError --starttime: {}".format(err))
            sys.exit()
    else:
        stime = 0.0
    if end_time:
        try:
            etime_dt = datetime.strptime(end_time, '%y%m%dT%H%M')
            etime = datetime.timestamp(etime_dt)

        except ValueError as err:
            print("ValueError --endtime: {}".format(err))
            sys.exit()
    else:
        etime = datetime.timestamp(datetime.now())

    if etime < stime:
        sys.exit("End time can't be earlier than start time")
    # Create empty arrays for TCS in position time stamp and value
    recGroups = []
    recData = {}
    channel_list = []
    print(hdf5File)
    # Read h5 file
    data_files = [h5py.File(datafile, 'r') for datafile in hdf5File]
    # If channel mask not defined, make the mask be every channel in the data
    # set

    # Generate list with all channels in dataset
    for data_group in data_files:
        for cn in data_group:
            channel_list.append(cn)
    # If there's a channel mask defined, check if channels are in data set.
    if channel_mask:
        mask = [chan in channel_list for chan in channel_mask]
        if not(all(mask)):
            missing_chan = [not(m) for m in mask]
            # Print all missing channels and abort execution
            print('Channels missing from dataset!!!')
            for ch in np.array(channel_mask)[missing_chan]:
                print(ch)
            sys.exit('Aborting')
    else:
        # If there's no mask defined, make it the data channel list.
        channel_mask = channel_list

    # Extract record names and create an array with group object and name of
    # group
    aux_mask = []
    for data_group in data_files:
        for rn in channel_mask:
            # Check to see if data name is in channel mask, if not, skip it
            if not(rn in data_group.keys()):
                # Create array with leftover channels
                aux_mask.append(rn)
                continue
            timestamps = data_group.get(rn).get('timestamp')[0:]
            values = data_group.get(rn).get('value')[0:]
            # Define a mask for timestamps inside data range
            filt_mask = (stime<timestamps) & (timestamps<etime)
            try:
                recData[rn] =[timestamps[filt_mask],
                              values[filt_mask]]
            except IndexError as error:
                print('Correcting {0} indexing error: {1}'.format(rn, error))
                timestamps = timestamps[1:]
                recData[rn] =[timestamps[filt_mask[1:]],
                              values[filt_mask[1:]]]
        channel_mask = aux_mask
    return recData

def list_hdf5(hdf5File):
    # Read h5 file
    data_files = [h5py.File(datafile, 'r') for datafile in hdf5File]
    # Extract record names and create an array with group object and name of
    # group
    print('Channels in Dataset:')
    for data_group in data_files:
        for rn in data_group:
            print(rn)

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
    diffA = [dataA[0][1:], np.array(diffPV)]
    return diffA

def avgData(dataA, avg):
    '''
    This function is used to calculate the moving average of a data array
    '''
    # avg = 3
    avgPV = [np.average(a)\
             for a in np.array([dataA[1][(avg-1)-x:len(dataA[1])-x]
                                for x in range(avg)]).T]
    avgA = [dataA[0][(avg-1):], np.array(avgPV)]
    return avgA

def derivData(dataA):
    '''
    This function is used to calculate the first derivative of data array
    between consecutive values
    '''
    diffPV = [(x - xa)/(y - ya)\
              for ya,y,xa,x in np.array([dataA[0][:-1], dataA[0][1:],
                                   dataA[1][:-1], dataA[1][1:]]).T]
    diffA = [dataA[0][1:], np.array(diffPV)]
    return diffA

def inPosDur(ipArray):
    inPos = False
    ipd = list()
    previpt = ipArray[0][0]
    for t,ip in np.array(ipArray).T:
        if not(ip) and inPos:
            inPos = False
            ipd.append([t, (t - previpt)])
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

def ipZonesT(offsetArray, zcolor='k'):
    '''
    ipZonesT() Generates a data matrix with timestamp start/end pairs for in
    position timeframes. It adds a color string to each zone for plotting
    purposes
    '''
    # Return an empty array if offset array is empty
    if (len(offsetArray[0]) < 2):
        return []
    # If first sample of in position flag is 1, set it as start time, if it's 0
    # assume is not in position.
    if offsetArray[1][0]:
        zeroFlag = True
    else:
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
    if offsetArray[1][-1] and not(zeroFlag):
        offsetNA.append(offsetArray[0][-1])
    if len(offsetPA) > len(offsetNA):
        offsetPA = offsetPA[:-1]
    zonecolor = [zcolor for i in range(len(offsetPA))]
    print('size array PA: ', len(offsetPA), 'size array NA: ', len(offsetNA))
    offsetZonesArray = np.array([offsetPA,offsetNA,zonecolor]).T
    if len(offsetZonesArray) == 1 and (offsetZonesArray[0][0] == offsetZonesArray[0][1]):
        return []
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

def filter_outliers(data, low_lim, high_lim):
    '''
    filter_outliers() filters a numpy array according to low and high limits.
    It generates two arrays, one with data within the limits, the other with the
    outliers.
    '''
    # Create mask based on low and high values. Filter data array accordingly
    filt_mask = (low_lim<data[1]) & (data[1]<high_lim)
    filtered_data = [data[0][filt_mask], data[1][filt_mask]]
    outlier_data = [data[0][~filt_mask], data[1][~filt_mask]]
    return filtered_data, outlier_data

def tracking_filter(data, ip_data):
    '''
    tracking_filter() filters a data array based on an in-position array,
    returning a the data within in position timestamps
    '''
    # Generate the in position start - end timestamps
    in_pos_data = ipZonesT(ip_data)
    # If there's no in position timestamps, return the data as-is. I'm
    # shamelessly assuming the data is in position. Sue me
    if not(len(in_pos_data)):
        return data
    # Lots of support arrays. Not sure how "best practices" or pythonic this is
    ts_aux = data[0]
    val_aux = data[1]
    ts_array = []
    val_array = []
    # Go through all the zone timestamps pairs
    for zone in in_pos_data:
        # Create a filter mask
        zone_mask = (float(zone[0])<ts_aux) & (ts_aux<float(zone[1]))
        # Apply filter mask. I'm also assuming deleting trailing and leading 5
        # points of data will filter out the huge error seen the first moments
        # after a slew. Fudge factors are getting on my nerves
        ts_array = np.append(ts_array[:-5], ts_aux[zone_mask])
        val_array = np.append(val_array[:-5], val_aux[zone_mask])
        # Trying to optimize operations, eliminate the filtered data from the
        # source array
        val_aux = val_aux[ts_aux>float(zone[1])]
        ts_aux = ts_aux[ts_aux>float(zone[1])]
    return [ts_array, val_array]

def fft_generator(data):
    '''
    fft_generator() Calculates the Fast Fourier Transform of a data set. It
    formats the resulting array so it makes sense when it's plotted (ignores
    the theoretical FFT analysis mumbo jumbo and goes straight to something that
    we sw folk can make sense of.
    '''
    # Generate an array of timedeltas for every sample and then calculate the
    # average to determine the sampling interval
    dt_array = data[0][1:] - data[0][:-1]
    dt = np.average(dt_array)
    # Generate the FFT and the array with the frequencies
    raw_fft = np.fft.fft(data[1])
    raw_fft_freq = np.fft.fftfreq(len(data[1]), dt)
    # Filter the weird math stuff, leave the real world info
    freq_mask = raw_fft_freq >= 0
    # Scale the Y-axis according to the number of samples
    data_fft_freq = raw_fft_freq[freq_mask]
    # data_fft = abs(raw_fft[freq_mask])*(2/len(data[1]))
    data_fft = abs(raw_fft[freq_mask])*(1/len(data[1]))
    data_fft_norm = np.concatenate(([data_fft[0]],data_fft[1:]*2))
    return [data_fft_freq, data_fft_norm]
    # return [data_fft_freq, data_fft]

def lost_dmd(tx_data, rx_data):
    '''
    lost_dmd() creates indicator array for lost demand pkgs between 2 subsystems
    '''
    # Extract timestamp from tracking arrays
    tx_stream = [x[0] for x in tx_data[1]]
    rx_stream = [x[0] for x in rx_data[1]]
    # Calculate the mask for data in Tx stream that is on in Rx stream
    mask = ~np.isin(tx_stream, rx_stream)
    # Filter Tx stream
    tx_indicators = np.ones(len(tx_stream))[mask]
    tx_ts = np.array(tx_data[0])[mask]
    return [tx_ts, tx_indicators]

def lost_dmd_diff(lost_pkg, diff_sec=0, diff_min=0):
    '''
    lost_dmd_diff() generates differentials of width "diff_width" over the time
    span of an array
    '''
    lwts = lost_pkg[0][0] # Assign first timestamp
    # Calculate differential window
    diff_width = diff_sec + diff_min*60
    uwts = lwts + diff_width # Differential window upper limit
    # Initialize auxiliary arrays and return array
    aux_ts, aux_count= lost_pkg
    diff_accum = []
    diff_ts = []

    # While the lower diff window timestamp is less than the last timestamp in
    # array
    while lwts <= lost_pkg[0][-1]:
        # Define a mask for all data pairs with timestamp less than upper lim
        mask = aux_ts < uwts
        # Sum all lost pkg counts in the mask, append result and ts in arrays
        diff_sum = np.sum(aux_count[mask])
        diff_accum.append(diff_sum)
        diff_ts.append(lwts)
        # Update aux arrays and diff window
        aux_ts = aux_ts[~mask]
        aux_count = aux_count[~mask]
        lwts = uwts
        uwts = lwts + diff_width

    return [diff_ts, diff_accum]

if __name__ == '__main__':
    x = np.arange(6000,step=0.1)
    x2 = np.arange(6000,step=0.01)
    y1 = (x/6000)**2
    y2 = (0.26*(np.sin((2*3.1415*0.040)*x)) + 0.27*(np.sin((2*3.1415*0.080) * x)) +
          0.22*(np.sin((2*3.1415*0.12)*x)) + 0.24*(np.sin((2*3.1415*0.395)*x)) +
          0.12*(np.sin((2*3.1415*0.786)*x)) + 0.1*(np.sin((2*3.1415*1.173)*x)) +
          0.07*(np.sin((2*3.1415*1.539)*x)) + 0.08*(np.sin((2*3.1415*2.100)*x)) +
          0.07*(np.sin((2*3.1415*2.800)*x)) + 0.07*(np.sin((2*3.1415*2.900)*x)) +
          0.07*(np.sin((2*3.1415*3.100)*x)) + 5)
          # 0.06*(np.sin((2*3.1415*1.539)*x)) + 2*x)
    print(np.average(y2))
    y5 = 50*(np.sin((2*3.1415*10)*x2))
    y3 = 500 / (x + 1)
    y4 = fft_generator([x,y2])
    plts = DataAxePlotter(2)
    plts.Axe['c1']['g4'] = DataAx([x,y1], 'r', ylabel='t1', height=5, label='g4', rawx=True)
    plts.Axe['c1']['g6'] = DataAx([x,y3], 'k', ylabel='t3', label='g6', height=1, shax='g4', rawx=True)
    plts.Axe['c1']['g5'] = DataAx([x,y2], 'b', linestyle='-', label='g5', shax='g4', rawx=True)
    plts.Axe['c1']['g7'] = DataAx([x2,y5], 'k', linestyle='--', label='g7', rawx=True)
    plts.Axe['c2']['g1'] = DataAx(y4, 'b', ylabel='t5', label='g2', height=2, rawx=True, standalone=True)
    # plts.Axe['c2']['g2'] = DataAx([x,((x-10)**2)+40], 'g', ylabel='t4', label='g2', height=2, rawx=True)
    plts.positionPlot()
    plts.plotConfig()

