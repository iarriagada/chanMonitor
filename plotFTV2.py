#!/usr/bin/env python3

from plotEpicsData import DataAx, DataAxePlotter
from datetime import datetime
import plotEpicsData as ped
import argparse
import sys


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
                           default=False,
                           help='Plot range start time, eg. 210111T2130')

    parser.add_argument('-et',
                           '--endtime',
                           dest='etime',
                           default=False,
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
    if args.stime:
        try:
            stime_dt = datetime.strptime(args.stime, '%y%m%dT%H%M')
            stime = datetime.timestamp(stime_dt)

        except ValueError as err:
            print("ValueError --starttime: {}".format(err))
            sys.exit()
    else:
        stime = 0
    if args.etime:
        try:
            etime_dt = datetime.strptime(args.etime, '%y%m%dT%H%M')
            etime = datetime.timestamp(etime_dt)

        except ValueError as err:
            print("ValueError --endtime: {}".format(err))
            sys.exit()
    else:
        etime = datetime.timestamp(datetime.now())

    if etime < stime:
        sys.exit("End time can't be earlier than start time")
    # Read h5 file
    recData = ped.extract_h5df(args.hdf5File, stime, etime,  args.list_chan)

    tc1_VALI_4 = [recData['tcs:drives:driveMCS.VALI'][0],
                 [e[4] for e in recData['tcs:drives:driveMCS.VALI'][1]]]
    tc1_VALA_0 = [recData['tcs:drives:driveCRS.VALA'][0],
                 [e[0] for e in recData['tcs:drives:driveCRS.VALA'][1]]]
    mcs_j_VALA_0 = [recData['mc:followA.J'][0],
                 [e[0] for e in recData['mc:followA.J'][1]]]
    mcs_dmd_diff = ped.diffData(mcs_j_VALA_0)
    cr_dmd_val0_diff = ped.diffData(tc1_VALA_0)


    ylim_ta = (-0.0001, 0.001)
    ylim_dr_mcs = (0.07, 0.13)
    ylim_dr_crcs = (0.07, 0.21)

    plts = DataAxePlotter(ncols=2)

    plts.Axe['c1']['mcsA[0]diff'] = DataAx(tc1_VALI_4,
                                        'g',
                                            linestyle='',
                                            marker='o',
                                           ylims=ylim_dr_mcs,
                                        label='tc1 -> mcs loop time diff',
                                        ylabel="Difference\n[sec]",
                                        linewidth=2.50)

    plts.Axe['c1']['crcs[0]diff'] = DataAx(cr_dmd_val0_diff,
                                        'b',
                                            linestyle='',
                                            marker='o',
                                           ylims=ylim_dr_crcs,
                                        label='tc1 -> cr loop time diff',
                                        ylabel="Difference\n[sec]",
                                        linewidth=2.50)

    plts.Axe['c1']['mc_dmd_diff'] = DataAx(mcs_dmd_diff,
                                        'r',
                                            linestyle='',
                                            marker='o',
                                           ylims=ylim_dr_mcs,
                                        label='mcs reception time diff',
                                        ylabel="Difference\n[sec]",
                                        linewidth=2.50)

    plts.Axe['c2']['mcs_diff_hist'] = DataAx(tc1_VALI_4[1],
                                        'g',
                                        label='tc1 -> mcs histogram time diff',
                                        ylabel="No Samples\n[un]",
                                        xlabel="Time diff [sec]",
                                              histbins=30,
                                        linewidth=2.50)

    plts.Axe['c2']['crcs_diff_hist'] = DataAx(cr_dmd_val0_diff[1],
                                        'b',
                                        label='tc1 -> cr histogram time diff',
                                        ylabel="No Samples\n[un]",
                                        xlabel="Time diff [sec]",
                                              histbins=30,
                                        linewidth=2.50)

    plts.Axe['c2']['mcs_dmd_diff_hist'] = DataAx(mcs_dmd_diff[1],
                                        'r',
                                        label='mcs reception histogram time diff',
                                        ylabel="No Samples\n[un]",
                                        xlabel="Time diff [sec]",
                                              histbins=30,
                                        linewidth=2.50)

    # plts.Axe['c1']['ta_tc1'] = DataAx(recData['ta:tcs:diff'],
                                        # 'r',
                                      # label='ta:tcs:diff',
                                      # ylims=ylim_ta,
                                        # ylabel="Difference\n[sec]",
                                        # linewidth=1.50)

    # plts.Axe['c1']['ta_cr1'] = DataAx(recData['ta:cr:diff'],
                                        # 'g',
                                      # label='ta:cr:diff',
                                      # shax='ta_tc1',
                                        # ylabel="Difference\n[sec]",
                                        # linewidth=1.50)

    # # pC1['topsDPos'] = DataAx(tsCorrDmd,
                        # # 'r', "Position\n[deg]",
                        # # label='Top Shut Demand',
                        # # drawstyle='steps-post',
                        # # shax='botsDPos')

    # # pC1['topsPos'] = DataAx(tsCorrPos,
                        # # 'm-', "Position\n[deg]",
                        # # label='Top Shut Position',
                        # # linewidth=2.50,
                        # # shax='botsDPos')

    # # # pC1['botShtrHLS'] = DataAx(recData['ec:botShtrmotor.HLS'],
                        # # # 'r', "Accum State\n[un]",
                        # # # label='BS HL Switch',
                        # # # linewidth=2.50,
                        # # # ylims=(-0.5,3),
                        # # # drawstyle='steps-post')

    # plts.Axe['c1']['botsPos'] = DataAx(recData['ec:fault:word16:bits.B6'],
                                        # 'k',
                                        # drawstyle='steps-post',
                                        # label='Slow Speed, Bit 6',
                                        # ylabel="Switch State\n[bool]",
                                        # linewidth=2.50,
                                        # ylims=ylimsOFF)

    # # pC1['LSbit7'] = DataAx(recData['ec:fault:word16:bits.B7'],
                        # # 'g', "Switch State\n[bool]",
                        # # label='End-of-travel, Bit 7',
                        # # linewidth=2.50,
                        # # ylims=ylimsOFF,
                        # # drawstyle='steps-post')


    plts.positionPlot()
    plts.plotConfig('Fast Track Analysis')
