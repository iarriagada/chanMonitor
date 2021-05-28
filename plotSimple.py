#!/usr/bin/env python3

from plotEpicsData import DataAx, DataAxePlotter
from datetime import datetime
import numpy as np
import plotEpicsData as ped
import argparse
import sys

CHAN_LIST = ['pwfs2:dc:fgDiag1P2.VALQ']

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

    parser.add_argument('-st',
                           '--starttime',
                           dest='stime',
                           default=None,
                           help='Plot range start time, eg. 210111T2130')

    parser.add_argument('-et',
                           '--endtime',
                           dest='etime',
                           default=None,
                           help='Plot range end time, eg. 210111T2315')

    parser.add_argument('-lc',
                        '--list-channels',
                        dest='list_chan',
                        action='store_true',
                        help='List channel names contained in data file')

    parser.add_argument('-hst',
                        '--histogram',
                        dest='hst',
                        action='store_true',
                        help='Use this option to plot In Pos Duration Histogram')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if args.list_chan:
        ped.list_hdf5(args.hdf5File)
        sys.exit()
    # Read h5 file
    recData = ped.extract_hdf5(args.hdf5File,
                               args.stime,
                               args.etime,
                               CHAN_LIST)

    p2_counts = DataAx(recData['pwfs2:dc:fgDiag1P2.VALQ'],
                      'xkcd:grass green',
                      label='P2 Counts',
                      ylabel='Counts [un]',
                      linewidth=1.5)

    plts = DataAxePlotter(ncols=1)

    plts.Axe['c1']['pwfs2_counts'] = p2_counts

    plts.positionPlot()
    plts.plotConfig('FR-41226')
