#!/usr/bin/env python3

from plotEpicsData import DataAx, DataAxePlotter
from datetime import datetime
import numpy as np
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
    if args.list_chan:
        ped.list_hdf5(args.hdf5File)
        sys.exit()
    # Read h5 file
    recData = ped.extract_hdf5(args.hdf5File,
                               args.stime,
                               args.etime,
                               args.list_chan)

    tc1_VALI_4 = [recData['tcs:drives:driveMCS.VALI'][0],
                 np.array([e[4] for e in recData['tcs:drives:driveMCS.VALI'][1]])]
    tc1_VALA_0 = [recData['tcs:drives:driveCRS.VALA'][0],
                 [e[0] for e in recData['tcs:drives:driveCRS.VALA'][1]]]
    mcs_j_VALA_0 = [recData['mc:followA.J'][0],
                 [e[0] for e in recData['mc:followA.J'][1]]]

    mcs_dmd_diff = ped.diffData(mcs_j_VALA_0)
    cr_dmd_val0_diff = ped.diffData(tc1_VALA_0)
    tc_diff_norm, tc_diff_outl = ped.filter_outliers(tc1_VALI_4,
                                                     0.00,
                                                     0.22)
    mc_diff_norm, mc_diff_outl = ped.filter_outliers(mcs_dmd_diff,
                                                     0.00,
                                                     0.22)
    ta_tcs_norm, ta_tcs_outl = ped.filter_outliers(recData['ta:tcs:diff'],
                                                     0.00,
                                                     0.006)
    ta_mc_norm, ta_mc_outl = ped.filter_outliers(recData['ta:mc:diff'],
                                                     0.00,
                                                     0.006)
    ta_cr_norm, ta_cr_outl = ped.filter_outliers(recData['ta:cr:diff'],
                                                     0.00,
                                                     0.006)
    ta_m2_norm, ta_m2_outl = ped.filter_outliers(recData['ta:m2:diff'],
                                                     0.00,
                                                     0.006)
    ta_max_array = [np.amax(a) for a in [ta_tcs_outl[1],
                                         ta_mc_outl[1],
                                         ta_cr_outl[1],
                                         ta_m2_outl[1]] if len(a)]
    ta_min_array = [np.amin(a) for a in [ta_tcs_outl[1],
                                         ta_mc_outl[1],
                                         ta_cr_outl[1],
                                         ta_m2_outl[1]] if len(a)]


    ylim_ta_hist = (0, len(recData['ta:tcs:diff'][1]))
    ylim_ta_out_tcs = (0, len(ta_tcs_outl))
    ylim_ta_out_mcs = (0, len(ta_mc_outl))
    ylim_ta_out_crcs = (0, len(ta_cr_outl))
    ylim_ta_out_scs = (0, len(ta_m2_outl))
    ylim_ta = (-0.0001, 0.008)
    ylim_ta_out = (-1.1, 0.1)
    ylim_dr_mcs = (0.07, 0.13)
    ylim_dr_crcs = (0.07, 0.21)
    lim_bin_ta = (0, 6)
    lim_bin_ta_out = (np.amin(ta_min_array), np.amax(ta_max_array))

    # if 'mc:azDemandPos' in recData.keys():
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

    tcs_drvMCSvI4 = DataAx(tc1_VALI_4,
                            'g',
                            height=2,
                            linestyle='',
                            marker='o',
                            ylims=ylim_dr_mcs,
                            label='tc1 -> mcs loop time diff',
                            ylabel="Difference\n[sec]",
                            linewidth=2.50)

    tcs_drvCRSvA0 = DataAx(cr_dmd_val0_diff,
                            'b',
                            height=2,
                            linestyle='',
                            marker='o',
                            ylims=ylim_dr_crcs,
                            label='tc1 -> cr loop time diff',
                            ylabel="Difference\n[sec]",
                            linewidth=2.50)

    tcs_drvMCSvI4_norm = DataAx(tc_diff_norm,
                               'g',
                               linestyle='',
                               marker='o',
                               label='tc/mc tx time diff',
                               ylabel="Difference\n[sec]",
                               linewidth=2.50)
                                # ylims=ylim_dr_mcs,

    tcs_drvMCSvI4_outl = DataAx(tc_diff_outl,
                                'r',
                                linestyle='',
                                marker='o',
                                label='tc/mc tx time diff outl',
                                ylabel="Difference\n[sec]",
                                linewidth=2.50)
                                # ylims=ylim_dr_mcs,

    mcs_followAJ0_norm = DataAx(mc_diff_norm,
                                'b',
                                linestyle='',
                                marker='o',
                                label='mcs reception time diff',
                                ylabel="Difference\n[sec]",
                                linewidth=2.50)
    # ylims=ylim_dr_mcs,

    mcs_followAJ0_outl = DataAx(mc_diff_outl,
                                'r',
                                linestyle='',
                                marker='o',
                                label='mcs rx time diff outl',
                                ylabel="Difference\n[sec]",
                                linewidth=2.50)
    # ylims=ylim_dr_mcs,

    tcs_drvMCSvI4_hist = DataAx(tc1_VALI_4[1],
                                'g',
                                label='tc1 -> mcs histogram time diff',
                                ylabel="No Samples\n[un]",
                                xlabel="Time diff [sec]",
                                histbins=30,
                                linewidth=2.50)

    tcs_drvCRSvA0_hist = DataAx(cr_dmd_val0_diff[1],
                                'b',
                                label='tc1 -> cr histogram time diff',
                                ylabel="No Samples\n[un]",
                                xlabel="Time diff [sec]",
                                histbins=30,
                                linewidth=2.50)

    tcs_drvMCSvI4_nrm_hist = DataAx(tc_diff_norm[1],
                                    'g',
                                    label='tc/mc tx time diff hist',
                                    ylabel="No Samples\n[un]",
                                    xlabel="Time diff [sec]",
                                    histbins=30,
                                    linewidth=2.50)

    tcs_drvMCSvI4_out_hist = DataAx(tc_diff_outl[1],
                                    'r',
                                    label='tc/mc tx time diff outl hist',
                                    ylabel="No Samples\n[un]",
                                    xlabel="Time diff [sec]",
                                    histbins=30,
                                    linewidth=2.50)

    mcs_followAJ0_norm_hist = DataAx(mc_diff_norm[1],
                                     'b',
                                     label='mcs rx time diff hist',
                                     ylabel="No Samples\n[un]",
                                     xlabel="Time diff [sec]",
                                     histbins=30,
                                     linewidth=2.50)

    mcs_followAJ0_out_dist = DataAx(mc_diff_outl[1],
                                     'r',
                                     label='mcs rx time diff outl hist',
                                     ylabel="No Samples\n[un]",
                                     xlabel="Time diff [sec]",
                                     histbins=30,
                                     linewidth=2.50)

    timeaudit_tcs_norm = DataAx(ta_tcs_norm,
                          'r',
                           linestyle='',
                           marker='o',
                          label='ta:tcs:diff',
                          ylims=ylim_ta,
                          ylabel="Difference\n[sec]",
                           marksize=4,
                                height=2,
                          linewidth=1.0)

    timeaudit_tcs_outl = DataAx(ta_tcs_outl,
                          'r',
                           linestyle='',
                           marker='o',
                          label='ta:tcs:diff outlier',
                          ylabel="Difference\n[sec]",
                           marksize=3,
                                height=2,
                          linewidth=1.0)
                          # ylims=ylim_ta_out,

    timeaudit_crcs_norm = DataAx(ta_cr_norm,
                            'g',
                           linestyle='',
                           marker='o',
                            label='ta:cr:diff',
                            ylabel="Difference\n[sec]",
                           marksize=3,
                            linewidth=1.0)

    timeaudit_crcs_outl = DataAx(ta_cr_outl,
                            'g',
                           linestyle='',
                           marker='o',
                            label='ta:cr:diff outlier',
                            ylabel="Difference\n[sec]",
                           marksize=3,
                            linewidth=1.0)

    timeaudit_mcs_norm = DataAx(ta_mc_norm,
                            'b',
                           linestyle='',
                           marker='o',
                            label='ta:mc:diff',
                            ylabel="Difference\n[sec]",
                           marksize=2,
                            linewidth=1.0)

    timeaudit_mcs_outl = DataAx(ta_mc_outl,
                            'b',
                           linestyle='',
                           marker='o',
                            label='ta:mc:diff outlier',
                            ylabel="Difference\n[sec]",
                           marksize=3,
                            linewidth=1.0)

    timeaudit_scs_norm = DataAx(ta_m2_norm,
                           'xkcd:tangerine',
                           linestyle='',
                           marker='o',
                            label='ta:m2:diff',
                            ylabel="Difference\n[sec]",
                           marksize=1,
                           alpha=0.25,
                            linewidth=0.5)

    timeaudit_scs_outl = DataAx(ta_m2_outl,
                           'xkcd:tangerine',
                           linestyle='',
                           marker='o',
                            label='ta:m2:diff outlier',
                            ylabel="Difference\n[sec]",
                           marksize=3,
                            linewidth=0.5)

    ta_tcs_norm_hist = DataAx(ta_tcs_norm[1]*1000,
                           'r',
                           label='ta:tcs:diff',
                           ylabel="No samples [un]",
                           xlabel="Time diff [ms]",
                           histbins=30,
                                 ylims=ylim_ta_hist,
                                limsbins=lim_bin_ta,
                           linewidth=1.50)

    ta_tcs_outl_hist = DataAx(ta_tcs_outl[1],
                           'r',
                           label='ta:tcs:diff outlier',
                           ylabel="No samples [un]",
                           xlabel="Time diff [s]",
                           histbins=30,
                                limsbins=lim_bin_ta_out,
                           linewidth=1.50)
                                 # ylims=ylim_ta_out_tcs,

    ta_crcs_norm_hist = DataAx(ta_cr_norm[1]*1000,
                            'g',
                            label='ta:cr:diff',
                            ylabel="No samples [un]",
                            xlabel="Time diff [ms]",
                            histbins=30,
                                 ylims=ylim_ta_hist,
                                limsbins=lim_bin_ta,
                            linewidth=1.50)

    ta_crcs_outl_hist = DataAx(ta_cr_outl[1],
                            'g',
                            label='ta:cr:diff outlier',
                            ylabel="No samples [un]",
                            xlabel="Time diff [s]",
                            histbins=30,
                                limsbins=lim_bin_ta_out,
                            linewidth=1.50)
                                 # ylims=ylim_ta_out_crcs,

    ta_mcs_norm_hist = DataAx(ta_mc_norm[1]*1000,
                            'b',
                            label='ta:mc:diff',
                            ylabel="No samples [un]",
                            xlabel="Time diff [ms]",
                            histbins=30,
                                 ylims=ylim_ta_hist,
                                limsbins=lim_bin_ta,
                            linewidth=1.50)

    ta_mcs_outl_hist = DataAx(ta_mc_outl[1],
                            'b',
                            label='ta:mc:diff outlier',
                            ylabel="No samples [un]",
                            xlabel="Time diff [s]",
                            histbins=30,
                                limsbins=lim_bin_ta_out,
                            linewidth=1.50)
                                 # ylims=ylim_ta_out_mcs,

    ta_scs_norm_hist = DataAx(ta_m2_norm[1]*1000,
                           'xkcd:tangerine',
                            label='ta:m2:diff',
                            ylabel="No samples [un]",
                            xlabel="Time diff [ms]",
                            histbins=30,
                                 ylims=ylim_ta_hist,
                                limsbins=lim_bin_ta,
                            linewidth=1.50)

    ta_scs_outl_hist = DataAx(ta_m2_outl[1],
                           'xkcd:tangerine',
                            label='ta:m2:diff outlier',
                            ylabel="No samples [un]",
                            xlabel="Time diff [s]",
                            histbins=30,
                                limsbins=lim_bin_ta_out,
                            linewidth=1.50)
                                 # ylims=ylim_ta_out_scs,

    plts = DataAxePlotter(ncols=3)

    plts.Axe['c1']['mc_azDmd'] = mc_azDmd
    plts.Axe['c1']['mc_azPos'] = DataAx.update_axe(mc_azPos,
                                                   shaxname='mc_azDmd')
    plts.Axe['c1']['ta_tcs_norm'] = timeaudit_tcs_norm
    plts.Axe['c1']['ta_crcs_norm'] = DataAx.update_axe(timeaudit_crcs_norm,
                                                              shaxname='ta_tcs_norm')
    plts.Axe['c1']['ta_mcs_norm'] = DataAx.update_axe(timeaudit_mcs_norm,
                                                             shaxname='ta_tcs_norm')
    plts.Axe['c1']['ta_scs_norm'] = DataAx.update_axe(timeaudit_scs_norm,
                                                             shaxname='ta_tcs_norm')
    plts.Axe['c1']['ta_tcs_outl'] = timeaudit_tcs_outl
    plts.Axe['c1']['ta_crcs_outl'] = DataAx.update_axe(timeaudit_crcs_outl,
                                                              shaxname='ta_tcs_outl')
    plts.Axe['c1']['ta_mcs_outl'] = DataAx.update_axe(timeaudit_mcs_outl,
                                                             shaxname='ta_tcs_outl')
    plts.Axe['c1']['ta_scs_outl'] = DataAx.update_axe(timeaudit_scs_outl,
                                                             shaxname='ta_tcs_outl')
    plts.Axe['c2']['ta_tcs_norm_hist'] = ta_tcs_norm_hist
    plts.Axe['c2']['ta_crcs_norm_hist'] = ta_crcs_norm_hist
    plts.Axe['c2']['ta_mcs_norm_hist'] = ta_mcs_norm_hist
    plts.Axe['c2']['ta_scs_norm_hist'] = ta_scs_norm_hist

    plts.Axe['c3']['ta_tcs_outl_hist'] = ta_tcs_outl_hist
    plts.Axe['c3']['ta_crcs_outl_hist'] = ta_crcs_outl_hist
    plts.Axe['c3']['ta_mcs_outl_hist'] = ta_mcs_outl_hist
    plts.Axe['c3']['ta_scs_outl_hist'] = ta_scs_outl_hist

    plts.positionPlot()
    plts.plotConfig('Fast Track Analysis')
