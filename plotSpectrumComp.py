#!/usr/bin/env python3

from plotEpicsData import DataAx, DataAxePlotter
from datetime import datetime
from dataFromGea import utc2site_time
import numpy as np
import plotEpicsData as ped
import argparse
import sys
import re

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

    parser.add_argument('args_file',
                        metavar='ARGS-FILE',
                        nargs='*',
                        help='path to file with file names and time limits')

    # parser.add_argument('-st',
                           # '--starttime',
                           # dest='stime',
                           # default=None,
                           # help='Plot range start time, eg. 210111T2130')

    # parser.add_argument('-et',
                           # '--endtime',
                           # dest='etime',
                           # default=None,
                           # help='Plot range end time, eg. 210111T2315')

    # parser.add_argument('-lc',
                        # '--list-channels',
                        # dest='list_chan',
                        # action='store_true',
                        # help='List channel names contained in data file')

    # parser.add_argument('-cf',
                        # '--channels-file',
                        # dest='chan_file',
                        # default=None,
                        # help='Specify channel file with list of names')

    # parser.add_argument('-hst',
                        # '--histogram',
                        # dest='hst',
                        # action='store_true',
                        # help='Use this option to plot In Pos Duration Histogram')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    args_dict = {}
    data_dict = {}
    ax_dict = {}
    i = 1
    # if args.list_chan:
        # ped.list_hdf5(args.hdf5File)
        # sys.exit()

    # if args.chan_file:
        # with open(args.chan_file, 'r') as f:
            # CHAN_LIST = [c.strip('\n') for c in f.readlines()]

    # Read h5 file
    with open(args.args_file[0], 'r') as f:
        arguments = f.readlines()

    for l in arguments:
        if re.search(r'^([ ]*#|\n)',l):
            continue
        args_dict[f'Set_{i}'] = {}
        args_dict[f'Set_{i}']['files'] = l.split("-t")[0].strip().split(' ')
        args_dict[f'Set_{i}']['times'] = l.split("-t")[1].strip().split(' ')
        i += 1

    for s in args_dict:
        data_dict[s] = ped.extract_hdf5(args_dict[s]['files'],
                                        args_dict[s]['times'][0],
                                        args_dict[s]['times'][1])

    print(data_dict.keys())

    lim_dict = lambda x,y:{'bottom':x,'top':y}

    ylim_ta = lim_dict(-0.0001, 0.008)
    ylim_ta_out = lim_dict(-1.1, 0.1)
    ylim_dr_mcs = lim_dict(0.07, 0.13)
    ylim_dr_crcs = lim_dict(0.07, 0.21)
    # ylim_fft_el = lim_dict(-0.25e-5, 1.0e-4)
    # ylim_fft_az = lim_dict(-0.25e-5, 15e-5)
    ylim_fft_az = lim_dict(-0.01, 1.00)
    ylim_fft_el = lim_dict(-0.01, 0.25)

    # mcs_aztrack_err = ped.tracking_filter(recData['mc:azPosError'],
                                        # recData['mc:inPositionAz'])

    # az_err = 'mc:azPosError'
    # el_err = 'mc:elPosError'
    az_err = 'mc:azPmacPosError'
    el_err = 'mc:elPmacPosError'

    # mc_azterror_ffst = ped.fft_generator(mcs_aztrack_err)
    for s in data_dict:
        fft_az = ped.fft_generator(data_dict[s][az_err])
        fft_el = ped.fft_generator(data_dict[s][el_err])
        start_ts = utc2site_time(data_dict[s][az_err][0][0],
                                 'gs_st')
        date_str = datetime.strftime(start_ts, "%y/%m/%d %H:%M")
        ax_dict[s] = {'az_error':[fft_az[0], fft_az[1]*3600],
                      'el_error':[fft_el[0], fft_el[1]*3600],
                      'line_lgnd': date_str}
        ave_az = np.average(data_dict[s][az_err][1]*3600)
        no_samples_az = len(data_dict[s][az_err][1])
        ave_el = np.average(data_dict[s][el_err][1]*3600)
        no_samples_el = len(data_dict[s][el_err][1])
        print(f'Azimuth error average: {ave_az}')
        print(f'Number of samples: {no_samples_az}')
        print(f'Elevation error average: {ave_el}')
        print(f'Number of samples: {no_samples_el}')


    for s in ax_dict:
        ax_dict[s]['ax_azerr'] = DataAx(ax_dict[s]['az_error'],
                                        'xkcd:drab green',
                                        label='Az @ ' + ax_dict[s]['line_lgnd'],
                                        ylabel='Pos Error [arcsec]',
                                        ylims = ylim_fft_az,
                                        rawx=True,
                                        linewidth=1.25)
        # standalone=True,
        ax_dict[s]['ax_elerr'] = DataAx(ax_dict[s]['el_error'],
                                        'xkcd:dark periwinkle',
                                        label='El @ ' + ax_dict[s]['line_lgnd'],
                                        ylabel='Pos Error [arcsec]',
                                        xlabel='Frequency [Hz]',
                                        ylims = ylim_fft_el,
                                        rawx=True,
                                        linewidth=1.25)
                                        # label=f'{s} mc:elPosError FFT',


    plts = DataAxePlotter(ncols=len(ax_dict.keys()))
    xlabels = {}
    j = 1

    for c in plts.Axe:
        plts.Axe[c]['mc_azerr_fft'] = ax_dict[f'Set_{j}']['ax_azerr']
        plts.Axe[c]['mc_elerr_fft'] = ax_dict[f'Set_{j}']['ax_elerr']
        xlabels[c] = 'Frequency [Hz]'
        j += 1
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

    # plts.Axe['c1']['mc_trk_azerr_fft'] = mc_trk_azerr_fft
    # plts.Axe['c1']['mc_trk_elerr_fft'] = mc_trk_elerr_fft
    plts.positionPlot()
    # plts.plotConfig('Fast Track Analysis', xlabels)
    plts.plotConfig('Fast Track Analysis')
