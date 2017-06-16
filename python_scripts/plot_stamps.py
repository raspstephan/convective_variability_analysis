"""
Filename:     plot_stamps.py
Author:       Stephan Rasp, s.rasp@lmu.de
Description:  Plot stamps of precipitation

"""

import argparse
from helpers import get_config, save_fig_and_log, \
    make_datelist#, get_and_crop_radar_fobj
from datetime import timedelta
from cosmo_utils.plot import ax_contourf
from cosmo_utils.pyncdf import getfobj_ncdf, getfobj_ncdf_ens
from cosmo_utils.helpers import ddhhmmss, yyyymmddhh_strtotime,\
    make_timelist
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
    ncdffn = 'lfff' + ddhhmmss(t) + '.nc_30m_surf'

    fobj = getfobj_ncdf(date_dir + str(inargs.individual_ens) + '/OUTPUT/' +
                        ncdffn, fieldn='PREC_ACCUM')

    # Set up figure
    pw = get_config(inargs, 'plotting', 'page_width')
    fig, ax = plt.subplots(1, 1, figsize=(pw / 2., pw / 2.))

    jpl0, jpl1, ipl0, ipl1 = (50 + inargs.zoom_lat1,
                              357 - 51 + inargs.zoom_lat2,
                              50 + inargs.zoom_lon1,
                              357 - 51 + inargs.zoom_lon2)
    print jpl0, jpl1, ipl0, ipl1
    data = fobj.data
    lats = fobj.lats[jpl0:jpl1, ipl0:ipl1]
    lons = fobj.lons[jpl0:jpl1, ipl0:ipl1]
    polelat = fobj.polelat
    polelon = fobj.polelon
    lllat = lats[0, 0]
    lllon = lons[0, 0]
    urlat = lats[-1, -1]
    urlon = lons[-1, -1]
    Basemap_kwargs = { \
        "projection": "lcc",
    # the projection "lcc" lets the domain appear as rectangular on the 2D plot
        "lon_0": polelon,
        "lat_0": 90. + polelat,
        "llcrnrlat": lllat,
        "urcrnrlat": urlat,
        "llcrnrlon": lllon,
        "urcrnrlon": urlon,
        "resolution": 'h',
    }
    m = Basemap(**Basemap_kwargs)
    x, y = m(lons, lats)

    cfplot = m.contourf(x, y, data[jpl0:jpl1, ipl0:ipl1], levels=levelsPrec,
                        colors=cmPrec, ax=ax)

    m.drawcoastlines()  # linewidth=0.1, antialiased=0)
    m.drawcountries(linewidth=0.1, antialiased=0)
    if inargs.ind_scale:
        m.drawmapscale(inargs.ind_scale_pos[0], inargs.ind_scale_pos[1],
                       inargs.ind_scale_pos[0], inargs.ind_scale_pos[1],
                       inargs.ind_scale_len, barstyle='fancy', fontsize=7)

    if inargs.ind_box:
        jpl0_box, jpl1_box, ipl0_box, ipl1_box = (
            50 + inargs.ind_box_lat1,
            357 - 51 + inargs.ind_box_lat2,
            50 + inargs.ind_box_lon1,
            357 - 51 + inargs.ind_box_lon2)

        lon_bl = fobj.lons[jpl0_box, ipl0_box]
        lon_br = fobj.lons[jpl0_box, ipl1_box]
        lon_tl = fobj.lons[jpl1_box, ipl0_box]
        lon_tr = fobj.lons[jpl1_box, ipl1_box]

        lat_bl = fobj.lats[jpl0_box, ipl0_box]
        lat_br = fobj.lats[jpl0_box, ipl1_box]
        lat_tl = fobj.lats[jpl1_box, ipl0_box]
        lat_tr = fobj.lats[jpl1_box, ipl1_box]

        m.plot((lon_bl, lon_br), (lat_bl, lat_br), latlon=True, color='k')
        m.plot((lon_br, lon_tr), (lat_br, lat_tr), latlon=True, color='k')
        m.plot((lon_tr, lon_tl), (lat_tr, lat_tl), latlon=True, color='k')
        m.plot((lon_tl, lon_bl), (lat_tl, lat_bl), latlon=True, color='k')


# cf, tmp = ax_contourf(ax, fobj, colors=cmPrec,
#                           pllevels=levelsPrec,
#                           ji0=(50 + inargs.zoom_lat1,
#                                50 + inargs.zoom_lon1),
#                           ji1=(357 - 51 + inargs.zoom_lat2,
#                                357 - 51 + inargs.zoom_lon2),
#                           sp_title='',
#                           Basemap_drawrivers=False,
#                           npars=0, nmers=0)
    cb = fig.colorbar(cfplot, orientation='horizontal', shrink=0.6)
    cb.set_label('Accumulation [mm/h]')
    plt.tight_layout()

    save_fig_and_log(fig, None, inargs, 'prec_individual',
                     date=inargs.date_start, time=str(t.seconds / 3600).zfill(2),
                     tight=True)


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

    # Individual plot arguments
    parser.add_argument('--individual_ens',
                        type=str,
                        default='',
                        help='Plot type [stamps, individual]')
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