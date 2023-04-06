#!/usr/bin/env python3

from plotEpicsData import DataAx, DataAxePlotter
from datetime import datetime
import numpy as np
import plotEpicsData as ped
import argparse
import sys

# CHAN_LIST = ['tcs:drives:driveMCS.VALI',
             # 'mc:azCurrentPos',
             # 'mc:azDemandPos',
             # 'mc:azPmacDemandPos',
             # 'mc:elCurrentPos',
             # 'mc:elDemandPos',
             # 'mc:elPmacDemandPos',
             # 'mc:Tracking.VALC',
             # 'mc:azMPHandshake',
             # 'mc:elMPHandshake']

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

    parser.add_argument('-cf',
                        '--channels-file',
                        dest='chan_file',
                        default=None,
                        help='Specify channel file with list of names')

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

    if args.chan_file:
        with open(args.chan_file, 'r') as f:
            CHAN_LIST = [c.strip('\n') for c in f.readlines()]

    # Read h5 file
    recData = ped.extract_hdf5(args.hdf5File,
                               args.stime,
                               args.etime)

    lim_dict = lambda x,y:{'bottom':x,'top':y}

    ylim_ta = lim_dict(-0.0001, 0.008)
    ylim_ta_out = lim_dict(-1.1, 0.1)
    ylim_dr_mcs = lim_dict(0.07, 0.13)
    ylim_dr_crcs = lim_dict(0.07, 0.21)

    # mcs_aztrack_err = ped.tracking_filter(recData['mc:azPosError'],
                                        # recData['mc:inPositionAz'])

    # mc_azterror_fft = ped.fft_generator(mcs_aztrack_err)
    mc_azterror_fft = ped.fft_generator(recData['mc:azPosError'])

    # mcs_eltrack_err = ped.tracking_filter(recData['mc:elPosError'],
                                        # recData['mc:inPositionEl'])

    # mc_elterror_fft = ped.fft_generator(mcs_eltrack_err)
    mc_elterror_fft = ped.fft_generator(recData['mc:elPosError'])

    mc_azDmd = DataAx(recData['mc:azDemandPos'],
                      'xkcd:grass green',
                      label='mc:azDemandPos',
                      ylabel='Position [deg]',
                      linewidth=1.5)

    mc_azPos = DataAx(recData['mc:azCurrentPos'],
                      'xkcd:bright blue',
                      label='mc:azCurrentPos',
                      ylabel='Position [deg]',
                      linewidth=1.25)

    mc_azPmacDmd = DataAx(recData['mc:azPmacDemandPos'],
                      'xkcd:diarrhea',
                      label='mc:azPmacDemandPos',
                      ylabel='Position [deg]',
                      linewidth=1.25)

    mc_elDmd = DataAx(recData['mc:elDemandPos'],
                      'xkcd:grass green',
                      label='mc:elDemandPos',
                      ylabel='Position [deg]',
                      linewidth=1.5)

    mc_elPos = DataAx(recData['mc:elCurrentPos'],
                      'xkcd:bright blue',
                      label='mc:elCurrentPos',
                      ylabel='Position [deg]',
                      linewidth=1.25)

    mc_elPmacDmd = DataAx(recData['mc:elPmacDemandPos'],
                      'xkcd:diarrhea',
                      label='mc:elPmacDemandPos',
                      ylabel='Position [deg]',
                      linewidth=1.25)

    mc_azErr = DataAx(recData['mc:azPosError'],
                      'xkcd:brick red',
                      label='mc:azPosError',
                      ylabel='Position [deg]',
                      linewidth=1.25)

    mc_azPmacErr = DataAx(recData['mc:azPmacPosError'],
                      'xkcd:plum',
                      label='mc:azPmacPosError',
                      ylabel='Position [deg]',
                      linewidth=1.25)

    mc_elErr = DataAx(recData['mc:elPosError'],
                      'xkcd:red orange',
                      label='mc:elPosError',
                      ylabel='Position [deg]',
                      linewidth=1.25)

    mc_elPmacErr = DataAx(recData['mc:elPmacPosError'],
                      'xkcd:brownish orange',
                      label='mc:elPmacPosError',
                      ylabel='Position [deg]',
                      linewidth=1.25)

    mc_trk_azerr_fft = DataAx(mc_azterror_fft,
                      'xkcd:drab green',
                      label='mc:azPosError FFT',
                      ylabel='Intensity [dB]',
                            rawx=True,
                      linewidth=1.25)
                            # standalone=True,

    mc_trk_elerr_fft = DataAx(mc_elterror_fft,
                      'xkcd:dark periwinkle',
                      label='mc:elPosError FFT',
                      ylabel='Intensity [dB]',
                            rawx=True,
                      linewidth=1.25)
                            # standalone=True,

    plts = DataAxePlotter(ncols=1)

    # plts.Axe['c1']['mc_azDmd'] = mc_azDmd
    # plts.Axe['c1']['mc_azPos'] = DataAx.update_axe(mc_azPos,
                                                   # shaxname='mc_azDmd')
    # plts.Axe['c1']['mc_azPmacDmd'] = mc_azPmacDmd
    # plts.Axe['c1']['mc_azErr'] = mc_azErr
    # plts.Axe['c1']['mc_azPmacErr'] = mc_azPmacErr

    # plts.Axe['c1']['mc_elDmd'] = mc_elDmd
    # plts.Axe['c1']['mc_elPos'] = DataAx.update_axe(mc_elPos,
                                                   # shaxname='mc_elDmd')
    # plts.Axe['c1']['mc_elPmacDmd'] = mc_elPmacDmd
    # plts.Axe['c1']['mc_elErr'] = mc_elErr
    # plts.Axe['c1']['mc_elPmacErr'] = mc_elPmacErr

    plts.Axe['c1']['mc_trk_azerr_fft'] = mc_trk_azerr_fft
    plts.Axe['c1']['mc_trk_elerr_fft'] = mc_trk_elerr_fft
    plts.positionPlot()
    plts.plotConfig('Fast Track Analysis')
