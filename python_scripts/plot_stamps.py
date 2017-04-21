"""
Filename:     plot_stamps.py
Author:       Stephan Rasp, s.rasp@lmu.de
Description:  Plot stamps of precipitation

"""

import argparse
from helpers import get_config, save_fig_and_log, \
    make_datelist, get_and_crop_radar_fobj
from datetime import timedelta
from cosmo_utils.plot import ax_contourf
from cosmo_utils.pyncdf import getfobj_ncdf, getfobj_ncdf_ens
from cosmo_utils.helpers import ddhhmmss, yyyymmddhh_strtotime,\
    make_timelist
import matplotlib.pyplot as plt
import numpy as np


################################################################################
# PLOTTING FUNCTIONS
################################################################################
def plot_prec_stamps(inargs):
    """
    Plots precipitation stamps of obs, det and ensemble every hour for each
    date and time specified

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    """

    # TODO: Update colors
    cmPrec = ((1, 1, 1),
              (0, 0.627, 1),
              (0.137, 0.235, 0.98),
              (0.392, 0, 0.627),
              (0.784, 0, 0.627))
    # (0.1  , 0.1   , 0.784),
    levelsPrec = [0, 1, 3, 10, 30, 100.]

    # Loop over dates
    for idate, date in enumerate(make_datelist(inargs)):
        print('Date ' + date)
        # Loop over times
        for t in make_timelist(timedelta(hours=inargs.time_start),
                               timedelta(hours=inargs.time_end),
                               timedelta(hours=1)):
            print('Time ' + str(t))
            # Load all CU objects
            fobjlist = []
            titlelist = []

            # 1st: Radar
            fobjlist.append(get_and_crop_radar_fobj(inargs, date, t))
            titlelist.append('Radar')

            # 2nd: det
            ncdffn_pref = (get_config(inargs, 'paths', 'raw_data') +
                           date + '/deout_ceu_pspens/' + 'det' +
                           '/OUTPUT/lfff')
            fobjlist.append(getfobj_ncdf(ncdffn_pref + ddhhmmss(t) +
                                         '.nc_30m_surf',
                                         'PREC_ACCUM'))
            titlelist.append('Det')

            # 3rd: ens
            ncdffn = 'lfff' + ddhhmmss(t) + '.nc_30m_surf'
            date_dir = (get_config(inargs, 'paths', 'raw_data') +
                        date + '/deout_ceu_pspens/')
            fobjlist.extend(getfobj_ncdf_ens(date_dir, 'sub', inargs.nens,
                                             ncdffn, dir_suffix='/OUTPUT/',
                                             fieldn='PREC_ACCUM', nfill=1))
            titlelist.extend(['Mem ' + str(i + 1) for i in range(inargs.nens)])

            # Now plot
            n_panels = len(fobjlist)
            n_cols = 4
            n_rows = int(np.ceil(float(n_panels) / n_cols))
            fig, axmat = plt.subplots(n_rows, n_cols,
                                      figsize=(10, 3.5 * n_rows))
            axflat = np.ravel(axmat)

            for i in range(len(fobjlist)):
                plt.sca(axflat[i])
                cf, tmp = ax_contourf(axflat[i], fobjlist[i], colors=cmPrec,
                                      pllevels=levelsPrec,
                                      ji0=(50 + inargs.zoom_lat1,
                                           50 + inargs.zoom_lon1),
                                      ji1=(357 - 51 + inargs.zoom_lat2,
                                           357 - 51 + inargs.zoom_lon2),
                                      sp_title=titlelist[i],
                                      Basemap_drawrivers=False,
                                      npars=0, nmers=0)
            cb = fig.colorbar(cf, cax=fig.add_axes([0.4, 0.1, 0.2, 0.02]),
                              orientation='horizontal')
            cb.set_label('Accumulation [mm/h]')
            titlestr = (yyyymmddhh_strtotime(date).strftime(
                get_config(inargs, 'plotting', 'date_fmt')) +
                        ' ' + str(t.seconds / 3600).zfill(2) + 'UTC')
            fig.suptitle(titlestr)
            plt.tight_layout(rect=[0, 0.1, 1, 0.93])

            # Save figure and log
            save_fig_and_log(fig, None, inargs, 'prec_stamps',
                             date=date, time=str(t.seconds / 3600).zfill(2))


def main(inargs):
    """
    Runs the main program

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    """

    plot_prec_stamps(inargs)


if __name__ == '__main__':
    description = __doc__

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--date_start',
                        type=str,
                        help='Start date of analysis in yyyymmddhh')
    parser.add_argument('--date_end',
                        type=str,
                        help='End date of analysis in yyyymmddhh')
    parser.add_argument('--time_start',
                        type=int,
                        default=1,
                        help='Analysis start time in hrs [including]. \
                              Default = 1')
    parser.add_argument('--time_end',
                        type=int,
                        default=24,
                        help='Analysis end time in hrs [including]. \
                              Default = 24')
    parser.add_argument('--time_inc',
                        type=float,
                        default=1,
                        help='Analysis increment in hrs. Default = 1')
    parser.add_argument('--nens',
                        type=int,
                        help='Number of ensemble members')
    parser.add_argument('--config_file',
                        type=str,
                        default='config.yml',
                        help='Config file in relative directory ../config. \
                              Default = config.yml')
    parser.add_argument('--zoom_lon1',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lon1')
    parser.add_argument('--zoom_lon2',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lon2')
    parser.add_argument('--zoom_lat1',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lat1')
    parser.add_argument('--zoom_lat2',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lat2')
    parser.add_argument('--sub_dir',
                        type=str,
                        default='prec_stamps',
                        help='Sub-directory for figures and pp_files')
    parser.add_argument('--plot_name',
                        type=str,
                        default='',
                        help='Custom plot name.')

    args = parser.parse_args()

    main(args)