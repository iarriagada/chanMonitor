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
                               args.etime)

    tc1_VALI_4 = [recData['tcs:drives:driveMCS.VALI'][0],
                 np.array([e[4] for e in recData['tcs:drives:driveMCS.VALI'][1]])]
    tc1_VALA_0 = [recData['tcs:drives:driveCRS.VALA'][0],
                 np.array([e[0] for e in recData['tcs:drives:driveCRS.VALA'][1]])]
    mcs_j_VALA_0 = [recData['mc:followA.J'][0],
                 [e[0] for e in recData['mc:followA.J'][1]]]

    mcs_dmd_diff = ped.diffData(mcs_j_VALA_0)
    cr_dmd_val0_diff = ped.diffData(tc1_VALA_0)
    tc_diff_norm, tc_diff_outl = ped.filter_outliers(tc1_VALI_4,
                                                     0.00,
                                                     0.15)
    mc_diff_norm, mc_diff_outl = ped.filter_outliers(mcs_dmd_diff,
                                                     0.00,
                                                     0.15)
    cr_diff_norm, cr_diff_outl = ped.filter_outliers(cr_dmd_val0_diff,
                                                     0.00,
                                                     0.15)

    ylim_dr_mcs = (0.07, 0.13)
    ylim_dr_crcs = (0.07, 0.21)
    lim_bin_mcs = (0, 0.15)
    # lim_bin_mcs_out = (0.15, np.amax(mc_diff_outl[1]))
    lim_bin_mcs_out = (0.15, 1)

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
                                label='mcs rx time diff',
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

    crs_followAJ0_norm = DataAx(cr_diff_norm,
                                'xkcd:orange',
                                linestyle='',
                                marker='o',
                                label='crcs tx time diff',
                                ylabel="Difference\n[sec]",
                                linewidth=2.50)
    # ylims=ylim_dr_mcs,

    crs_followAJ0_outl = DataAx(cr_diff_outl,
                                'r',
                                linestyle='',
                                marker='o',
                                label='crcs tx time diff outl',
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

    mcs_followAJ0_out_hist = DataAx(mc_diff_outl[1],
                                     'r',
                                     label='mcs rx time diff outl hist',
                                     ylabel="No Samples\n[un]",
                                     xlabel="Time diff [sec]",
                                     histbins=30,
                                     linewidth=2.50)

    tcs_drvCRSvA0_nrm_hist = DataAx(cr_diff_norm[1],
                                    'xkcd:orange',
                                    label='tc/cr tx time diff hist',
                                    ylabel="No Samples\n[un]",
                                    xlabel="Time diff [sec]",
                                    histbins=30,
                                    linewidth=2.50)

    tcs_drvCRSvA0_out_hist = DataAx(cr_diff_outl[1],
                                    'r',
                                    label='tc/cr tx time diff outl hist',
                                    ylabel="No Samples\n[un]",
                                    xlabel="Time diff [sec]",
                                    histbins=30,
                                    linewidth=2.50)


    plts = DataAxePlotter(ncols=3)

    plts.Axe['c1']['az_dmd'] = mc_azDmd
    plts.Axe['c1']['az_pos'] = DataAx.update_axe(mc_azPos, shaxname='az_dmd')
    plts.Axe['c1']['mcs_dmd_tx'] = DataAx.update_axe(tcs_drvMCSvI4_norm,
                                                         marksize=2)
    plts.Axe['c1']['mcs_dmd_tx_out'] = DataAx.update_axe(tcs_drvMCSvI4_outl,
                                                         marksize=2)
    plts.Axe['c1']['mcs_dmd_rx'] = DataAx.update_axe(mcs_followAJ0_norm,
                                                         marksize=2)
    plts.Axe['c1']['mcs_dmd_rx_out'] = DataAx.update_axe(mcs_followAJ0_outl,
                                                         marksize=2)
    plts.Axe['c1']['cr_dmd_rx'] = DataAx.update_axe(crs_followAJ0_norm,
                                                         marksize=2)
    plts.Axe['c1']['cr_dmd_rx_out'] = DataAx.update_axe(crs_followAJ0_outl,
                                                         marksize=2)
    plts.Axe['c2']['mcs_dtx_hist'] = DataAx.update_axe(tcs_drvMCSvI4_nrm_hist,
                                                       limsbins=lim_bin_mcs)
    norm_hist_mstax = plts.Axe['c2']['mcs_dtx_hist']
    plts.Axe['c2']['mcs_drx_hist'] = DataAx.update_axe(mcs_followAJ0_norm_hist,
                                                       mstax=norm_hist_mstax,
                                                       limsbins=lim_bin_mcs)
    plts.Axe['c2']['crs_dtx_hist'] = DataAx.update_axe(tcs_drvCRSvA0_nrm_hist,
                                                       mstax=norm_hist_mstax,
                                                       limsbins=lim_bin_mcs)
    plts.Axe['c3']['mcs_dtx_out_hist'] = DataAx.update_axe(tcs_drvMCSvI4_out_hist,
                                                       limsbins=lim_bin_mcs_out)
    outl_hist_mstax = plts.Axe['c3']['mcs_dtx_out_hist']
    plts.Axe['c3']['mcs_drx_out_hist'] = DataAx.update_axe(mcs_followAJ0_out_hist,
                                                       mstax=outl_hist_mstax,
                                                       limsbins=lim_bin_mcs_out)
    plts.Axe['c3']['crs_dtx_out_hist'] = DataAx.update_axe(tcs_drvCRSvA0_out_hist,
                                                       mstax=outl_hist_mstax,
                                                       limsbins=lim_bin_mcs_out)

    plts.positionPlot()
    plts.plotConfig('Fast Track Analysis')
