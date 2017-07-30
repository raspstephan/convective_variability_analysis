"""
Filename:     plot_stamps.py
Author:       Stephan Rasp, s.rasp@lmu.de
Description:  Plot stamps of precipitation

"""

import argparse
from helpers import get_config, save_fig_and_log, \
    make_datelist, identify_clouds, get_domain_limits, plot_stamp
from datetime import timedelta
from cosmo_utils.plot import ax_contourf
from cosmo_utils.pyncdf import getfobj_ncdf, getfobj_ncdf_ens, getfield_ncdf
from cosmo_utils.helpers import ddhhmmss, yyyymmddhh_strtotime,\
    make_timelist, yymmddhhmm
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap

# TODO: Update colors
cmPrec = ((1, 1, 1),
          (0, 0.627, 1),
          (0.137, 0.235, 0.98),
          (0.392, 0, 0.627),
          (0.784, 0, 0.627))
# (0.1  , 0.1   , 0.784),
levelsPrec = [0, 1, 3, 10, 30, 100.]


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
            radarpref =(get_config(inargs, 'paths', 'radar_data') +
                        get_config(inargs, 'paths', 'radar_prefx'))
            radarsufx = get_config(inargs, 'paths', 'radar_sufix')
            dtradar = timedelta(minutes=10)

            radartime = yymmddhhmm(yyyymmddhh_strtotime(date) + t - dtradar)
            radarfn = radarpref + radartime + radarsufx
            radar_fobj = getfobj_ncdf(radarfn, fieldn='pr', dwdradar=True)
            # Crop data
            l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad = \
                get_domain_limits(inargs)
            l11_diff = l11_rad - l11
            l12_diff = l12_rad - l12
            l21_diff = l21_rad - l21
            l22_diff = l22_rad - l22
            radar_fobj.data = radar_fobj.data[l11_diff:l12_diff,
                                              l21_diff:l22_diff]
            radar_fobj.lats = radar_fobj.lats[l11_diff:l12_diff,
                              l21_diff:l22_diff]
            radar_fobj.lons = radar_fobj.lons[l11_diff:l12_diff,
                              l21_diff:l22_diff]
            fobjlist.append(radar_fobj)
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
                cf = plot_stamp(inargs, fobjlist[i], cmPrec, levelsPrec,
                                axflat[i], 'PREC_ACCUM')
                axflat[i].set_title(titlelist[i])
            cb = fig.colorbar(cf, cax=fig.add_axes([0.4, 0.15, 0.2, 0.02]),
                              orientation='horizontal')
            cb.set_label('Accumulation [mm/h]')
            titlestr = ((yyyymmddhh_strtotime(date) + t).
                        strftime('%d %b - %H UTC'))
            fig.suptitle(titlestr)
            plt.tight_layout(rect=[0, 0.1, 1, 0.95])

            # Save figure and log
            save_fig_and_log(fig, None, inargs, 'prec_stamps',
                             datestr=((yyyymmddhh_strtotime(date) + t).
                                      strftime('%Y%m%d_%H')),
                             tight=True)


def plot_individual(inargs):
    """
    
    Parameters
    ----------
    inargs

    Returns
    -------

    """

    # Get data
    date_dir = (get_config(inargs, 'paths', 'raw_data') + inargs.date_start +
                '/deout_ceu_pspens/')
    t = timedelta(hours=inargs.time_start)

    if inargs.ind_var == 'PREC_ACCUM':
        ncdffn = 'lfff' + ddhhmmss(t) + '.nc_30m_surf'
        fobj = getfobj_ncdf(date_dir + str(inargs.ind_ens) + '/OUTPUT/' +
                            ncdffn, fieldn='PREC_ACCUM')
        cmap = None
        colors = cmPrec
        levels = levelsPrec
    elif inargs.ind_var == 'radar':
        radarpref = (get_config(inargs, 'paths', 'radar_data') +
                     get_config(inargs, 'paths', 'radar_prefx'))
        radarsufx = get_config(inargs, 'paths', 'radar_sufix')
        dtradar = timedelta(minutes=10)
        radartime = yymmddhhmm(yyyymmddhh_strtotime(inargs.date_start) + t -
                               dtradar)
        radarfn = radarpref + radartime + radarsufx
        radar_fobj = getfobj_ncdf(radarfn, fieldn='pr', dwdradar=True)
        # Crop data
        l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad = \
            get_domain_limits(inargs)
        l11_diff = l11_rad - l11
        l12_diff = l12_rad - l12
        l21_diff = l21_rad - l21
        l22_diff = l22_rad - l22
        radar_fobj.data = radar_fobj.data[l11_diff:l12_diff,
                          l21_diff:l22_diff]
        radar_fobj.lats = radar_fobj.lats[l11_diff:l12_diff,
                          l21_diff:l22_diff]
        radar_fobj.lons = radar_fobj.lons[l11_diff:l12_diff,
                          l21_diff:l22_diff]
        fobj = radar_fobj
    elif inargs.ind_var in ['obj_m', 'obj_prec']:
        if inargs.ind_var == 'obj_m':
            lvl = 30
            ncdffn = 'lfff' + ddhhmmss(t) + '.nc_30m'
            fobj = getfobj_ncdf(date_dir + str(inargs.ind_ens) + '/OUTPUT/' +
                                ncdffn, fieldn='W', levs=lvl)
            opt_field = (
                getfield_ncdf(date_dir + str(inargs.ind_ens) + '/OUTPUT/' +
                              ncdffn, fieldn='QC', levs=lvl) +
                getfield_ncdf(date_dir + str(inargs.ind_ens) + '/OUTPUT/' +
                              ncdffn, fieldn='QS', levs=lvl) +
                getfield_ncdf(date_dir + str(inargs.ind_ens) + '/OUTPUT/' +
                              ncdffn, fieldn='QI', levs=lvl)
            )
        else:
            ncdffn = 'lfff' + ddhhmmss(t) + '.nc_30m_surf'
            fobj = getfobj_ncdf(date_dir + str(inargs.ind_ens) + '/OUTPUT/' +
                                ncdffn, fieldn='PREC_ACCUM')
            opt_field = None

        labels, size_list, sum_list = identify_clouds(fobj.data, 1.,
                                                      opt_field=opt_field,
                                                      water=True,
                                                      neighborhood=3,
                                                      opt_thresh=0)
        fobj.data = labels
        cmap = plt.cm.prism
        colors = None
        levels = None
    else:
        raise Exception('Wrong variable!')

    # Set up figure
    pw = get_config(inargs, 'plotting', 'page_width')
    width_fraction = 2./9.
    ratio = 1.
    fig, ax = plt.subplots(1, 1, figsize=(pw * width_fraction,
                                          pw * width_fraction * ratio))

    cf = plot_stamp(inargs, fobj, cmPrec, levelsPrec, ax, 'PREC_ACCUM')

    if inargs.ind_colorbar:
        cb = fig.colorbar(cf, orientation='horizontal', shrink=0.6)
        cb.set_label('Accumulation [mm/h]')

    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

    save_fig_and_log(fig, None, inargs, 'prec_individual', tight=False)


################################################################################
# MAIN FUNCTION
################################################################################
def main(inargs):
    """
    Runs the main program

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    """

    if inargs.plot_type == 'stamps':
        plot_prec_stamps(inargs)
    if inargs.plot_type == 'individual':
        plot_individual(inargs)


if __name__ == '__main__':
    description = __doc__

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--date_start',
                        type=str,
                        help='Start date of analysis in yyyymmddhh')
    parser.add_argument('--date_end',
                        type=str,
                        default='',
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
                        default=3,
                        help='Number of ensemble members')
    parser.add_argument('--config_file',
                        type=str,
                        default='config.yml',
                        help='Config file in relative directory ../config. \
                              Default = config.yml')

    # Plotting arguments
    parser.add_argument('--plot_type',
                        type=str,
                        default='',
                        help='Plot type [stamps, individual]')
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
    parser.add_argument('--plot_format',
                        type=str,
                        default='pdf',
                        help='Which format for figure file.')


    # Individual plot arguments
    parser.add_argument('--ind_var',
                        type=str,
                        default='PREC_ACCUM',
                        help='Which field to plot [PREC_ACCUM, obj]')
    parser.add_argument('--ind_ens',
                        type=str,
                        default='',
                        help='which ensemble member')
    parser.add_argument('--ind_scale_pos',
                        type=float,
                        nargs='+',
                        default=[7., 48.4],
                        help='Where to put the scale [lon, lat]')
    parser.add_argument('--ind_scale',
                        dest='ind_scale',
                        action='store_true',
                        help='If given, plots box.')
    parser.set_defaults(ind_scale=False)
    parser.add_argument('--ind_scale_len',
                        type=float,
                        default=200,
                        help='Length of scale in km.')
    parser.add_argument('--ind_box',
                        dest='ind_box',
                        action='store_true',
                        help='If given, plots box.')
    parser.set_defaults(ind_box=False)
    parser.add_argument('--ind_box_lon1',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lon1')
    parser.add_argument('--ind_box_lon2',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lon2')
    parser.add_argument('--ind_box_lat1',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lat1')
    parser.add_argument('--ind_box_lat2',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lat2')
    parser.add_argument('--ind_colorbar',
                        dest='ind_colorbar',
                        action='store_true',
                        help='If given, plots colorbar')
    parser.set_defaults(ind_colorbar=False)

    # Other arguments
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