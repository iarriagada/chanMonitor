#!/usr/bin/env python3

from chanmonitor.lib.plotEpicsData import DataAx, DataAxePlotter
import chanmonitor.lib.plotEpicsData as ped
from datetime import datetime
import numpy as np
import argparse
import json
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

    parser.add_argument('-zn',
                        '--zones',
                        dest='zones_file',
                        default=None,
                        help='Specify file with list of zones to be drawn')

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

    zones_dict = {}
    if args.zones_file:
        with open(args.zones_file, 'r') as f:
            zones_dict = json.load(f)

    # Read h5 file
    recData = ped.extract_hdf5(args.hdf5File,
                               args.stime,
                               args.etime)

    lim_dict = lambda x,y:{'bottom':x,'top':y}

    ylim_ta = lim_dict(-0.0001, 0.008)
    ylim_ta_out = lim_dict(-1.1, 0.1)
    ylim_dr_mcs = lim_dict(-0.01, 0.01)
    ylim_dr_crcs = lim_dict(0.07, 0.21)
    # zone_upgrd = [[[1680649260.0, 1680655500.0, 'xkcd:vermillion'],
                 # [1680665400.0, 1680674400.0, 'xkcd:vermillion']],'Upgrd Sys']
    # zone_ops = [[[1680655500.0, 1680665400.0, 'xkcd:cerulean blue'],
                 # [1680674400.0, 1680678000.0, 'xkcd:cerulean blue']], 'Ops Sys']

    # mcs_aztrack_err = ped.tracking_filter(recData['mc:azPosError'],
                                        # recData['mc:inPositionAz'])

    # mc_azterror_fft = ped.fft_generator(mcs_aztrack_err)

    # mcs_eltrack_err = ped.tracking_filter(recData['mc:elPosError'],
                                        # recData['mc:inPositionEl'])

    # mc_elterror_fft = ped.fft_generator(mcs_eltrack_err)
    tcs_lost = ped.lost_dmd(recData['tcs:drives:driveMCS.VALI'],
                            recData['mc:followA.J'])
    tcs_total = len(recData['tcs:drives:driveMCS.VALI'][0])
    tcs_total_lost = len(tcs_lost[0])
    lost_prcnt = (tcs_total_lost/tcs_total) * 100
    tcs_total_accum = [tcs_lost[0], list(range(1,tcs_total_lost+1))]
    diff_window = 1
    lost_pkg_diff = ped.lost_dmd_diff(tcs_lost, diff_min=diff_window)
    tcs_dmd_val = [d[3] for d in recData['tcs:drives:driveMCS.VALA'][1]]
    tcs_dmd_array = [recData['tcs:drives:driveMCS.VALA'][0], tcs_dmd_val]
    mcs_fllw_val = [d[3] for d in recData['mc:followA.J'][1]]
    mcs_fllw_array = [recData['mc:followA.J'][0], mcs_fllw_val]

    mc_azDmd = DataAx(recData['mc:azDemandPos'],
                      'xkcd:grass green',
                      label='mc:azDemandPos',
                      ylabel='Position [deg]',
                          marker='o',
                      marksize=5,
                      zone=zones_dict,
                      linewidth=1.5)
                          # zone=[zone_upgrd, zone_ops],

    mc_azPos = DataAx(recData['mc:azCurrentPos'],
                      'xkcd:bright blue',
                      label='mc:azCurrentPos',
                      ylabel='Position [deg]',
                          marker='o',
                      marksize=5,
                      linewidth=1.25)

    mc_azPos_tcs = DataAx(tcs_dmd_array,
                      'xkcd:hot pink',
                          label='tcs:drives:driveMCS.VALI',
                      ylabel='Position [deg]',
                          marker='o',
                      marksize=5,
                      linewidth=1.25)

    mc_azPos_fllw = DataAx(mcs_fllw_array,
                      'xkcd:neon blue',
                          label='mc:followA.J',
                      ylabel='Position [deg]',
                          marker='o',
                      marksize=5,
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

    mc_azErr = DataAx([recData['mc:azPosError'][0],
                       recData['mc:azPosError'][1]],
                      'xkcd:brick red',
                      label='mc:azPosError',
                      ylabel='Position Error [deg]',
                      ylims = ylim_dr_mcs,
                      height = 2,
                      linewidth=1.25)

    mc_azPmacErr = DataAx(recData['mc:azPmacPosError'],
                      'xkcd:plum',
                      label='mc:azPmacPosError',
                      ylabel='Position Error [deg]',
                          ylims = ylim_dr_mcs,
                      height = 2,
                      linewidth=1.25)

    mc_elErr = DataAx([recData['mc:elPosError'][0],
                       recData['mc:elPosError'][1]],
                      'xkcd:red orange',
                      label='mc:elPosError',
                      ylabel='Position Error [deg]',
                          ylims = ylim_dr_mcs,
                      linewidth=1.25)

    mc_elPmacErr = DataAx(recData['mc:elPmacPosError'],
                      'xkcd:brownish orange',
                      label='mc:elPmacPosError',
                      ylabel='Position Error [deg]',
                          ylims = ylim_dr_mcs,
                      linewidth=1.25)

    tcs_lost_dmd = DataAx(tcs_lost,
                          'xkcd:cherry',
                          linestyle='',
                          marker='o',
                          label=f'Lost Dmd: {lost_prcnt:.2f}% {tcs_total_lost}/{tcs_total}',
                          ylabel='Lost Pkg [bool]')
                          # zone=[zone_upgrd, zone_ops],

    # tcs_lost_dmd = DataAx(tcs_total_accum,
                          # 'xkcd:cherry',
                          # linestyle='',
                          # marker='o',
                          # zone=[zone_upgrd, zone_ops],
                          # label=f'Lost TCS Dmd: {tcs_total_lost}/{tcs_total}',
                          # ylabel='Lost Pkg [bool]')

    tcs_lost_diff = DataAx(lost_pkg_diff,
                          'xkcd:cherry',
                          marker='o',
                          drawstyle='steps-post',
                      height = 2,
                           label=f'Acc Lost Dmd, Diff={diff_window} [min] ',
                          ylabel=f'Acc Lost Pkg [count]')
                          # zone=[zone_upgrd, zone_ops],

    plts = DataAxePlotter(ncols=2)

    plts.Axe['c2']['mc_azDmd'] = mc_azDmd
    plts.Axe['c2']['mc_azPos'] = DataAx.update_axe(mc_azPos,
                                                   shaxname='mc_azDmd')
    plts.Axe['c2']['mc_azPos_tcs'] = DataAx.update_axe(mc_azPos_tcs,
                                                   shaxname='mc_azDmd')
    plts.Axe['c2']['mc_azPos_fllw'] = DataAx.update_axe(mc_azPos_fllw,
                                                   shaxname='mc_azDmd')
    # plts.Axe['c1']['mc_azPmacDmd'] = mc_azPmacDmd
    plts.Axe['c1']['mc_azErr'] = mc_azErr
    plts.Axe['c1']['mc_azPmacErr'] = mc_azPmacErr
    plts.Axe['c1']['tcs_lost_dmd'] = tcs_lost_dmd
    plts.Axe['c1']['tcs_lost_diff'] = tcs_lost_diff

    # plts.Axe['c2']['mc_elDmd'] = mc_elDmd
    # plts.Axe['c2']['mc_elPos'] = DataAx.update_axe(mc_elPos,
                                                   # shaxname='mc_elDmd')
    # plts.Axe['c2']['mc_elPmacDmd'] = mc_elPmacDmd
    # plts.Axe['c2']['mc_elErr'] = mc_elErr
    # plts.Axe['c2']['mc_elPmacErr'] = mc_elPmacErr

    # plts.Axe['c3']['mc_trk_azerr_fft'] = mc_trk_azerr_fft
    # plts.Axe['c3']['mc_trk_elerr_fft'] = mc_trk_elerr_fft
    plts.positionPlot()
    plts.plotConfig('FR-42809 Analysis')
