"""
Filename:     plotting.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains functions to read pre-processed NetCDF files, do some
final analysis and plot the results.

"""

# Import modules
from helpers import read_netcdf_dataset, get_config, save_fig_and_log, \
    make_datelist
from datetime import datetime, timedelta
from cosmo_utils.plot import ax_contourf
from cosmo_utils.pyncdf import getfobj_ncdf, getfobj_ncdf_ens
from cosmo_utils.helpers import yymmddhhmm, ddhhmmss, yyyymmddhh_strtotime,\
    make_timelist
import matplotlib.pyplot as plt
import numpy as np


# Define functions
def plot_domain_mean_timeseries_individual(inargs, plot_type):
    """
    Function to plot time series of domain mean precipitation for each day

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    plot_type : str
      Type of plot. Must be 'precipitation' or 'cape_tauc'

    """

    assert plot_type in ['precipitation', 'cape_tauc'], \
        'Type must be precipitation or cape_tauc'

    # Read pre-processed data
    rootgroup = read_netcdf_dataset(inargs)
    n_days = rootgroup.dimensions['date'].size

    # Set up figure
    n_cols = 4
    n_rows = int(np.ceil(float(n_days) / n_cols))

    fig, axmat = plt.subplots(n_rows, n_cols, sharex=True, sharey=True,
                              figsize=(10, 3 * n_rows))
    axflat = np.ravel(axmat)

    # Loop over axes / days
    for iday in range(n_days):
        dateobj = (timedelta(seconds=int(rootgroup.variables['date'][iday])) +
                   datetime(1, 1, 1))
        datestr = dateobj.strftime(get_config(inargs, 'plotting', 'date_fmt'))
        axflat[iday].set_title(datestr)
        if iday >= ((n_cols * n_rows) - n_cols):  # Only bottom row
            axflat[iday].set_xlabel('Time [UTC]')

        if plot_type == 'precipitaiton':
            plot_precipitation_panel(inargs, axflat, iday, rootgroup)
        if plot_type == 'cape_tauc':
            plot_cape_tauc_panel(inargs, axflat, iday, rootgroup)

    # Finish figure
    axflat[0].legend(loc=0)

    plt.tight_layout()

    # Save figure
    save_fig_and_log(fig, rootgroup, inargs, plot_type + '_ts_individual')


def plot_precipitation_panel(inargs, axflat, iday, rootgroup):
    """
    Plots precipitation timeseries in each panel.
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    axflat : list
      List of axis objects
    iday : int
      Index for list object
    rootgroup : netCDF object
      NetCDF object with data to plot

    """

    dateobj = (timedelta(seconds=int(rootgroup.variables['date'][iday])) +
               datetime(1, 1, 1))
    datestr = dateobj.strftime(get_config(inargs, 'plotting', 'date_fmt'))
    axflat[iday].set_title(datestr)

    for group in rootgroup.groups:
        # Get data do be plotted
        # prec: det, obs or ens mean
        # prec_lower/upper: ensemble minimum, maximum
        if group == 'ens':
            prec_array = rootgroup.groups[group].variables['PREC_ACCUM'] \
                [iday, :, :]
            prec = np.mean(prec_array, axis=1)
            prec_lower = np.amin(prec_array, axis=1)
            prec_upper = np.amax(prec_array, axis=1)
        else:
            prec = rootgroup.groups[group].variables['PREC_ACCUM'] \
                [iday, :, 0]

        # Plot data
        axflat[iday].plot(rootgroup.variables['time'][:], prec, label=group,
                          c=get_config(inargs, 'colors', group))
        if group == 'ens':
            axflat[iday].fill_between(rootgroup.variables['time'][:],
                                      prec_lower, prec_upper,
                                      where=prec_upper >= prec_lower,
                                      facecolor=get_config(inargs,
                                                           'colors',
                                                           'ens_range'))


def plot_cape_tauc_panel(inargs, axflat, iday, rootgroup):
    """
    Plots cape and tau_c timeseries in each panel.

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    axflat : list
      List of axis objects
    iday : int
      Index for list object
    rootgroup : netCDF object
      NetCDF object with data to plot

    """

    ax1 = axflat[iday]
    ax2 = ax1.twinx()
    ax2.set_ylim(0,20)

    dateobj = (timedelta(seconds=int(rootgroup.variables['date'][iday])) +
               datetime(1, 1, 1))
    datestr = dateobj.strftime(get_config(inargs, 'plotting', 'date_fmt'))
    ax1.set_title(datestr)
    if iday % 4 == 0:  # Only left column
        ax1.set_ylabel('CAPE [J/kg]')
    if iday % 4 == 3:
        ax2.set_ylabel('tau_c [h]')
    else:
        ax2.get_yaxis().set_ticks([])

    for group in ['det', 'ens']:
        for var, ax, ls in zip(['CAPE_ML', 'TAU_C'],
                               [ax1, ax2],
                               ['-', '--']):
            # Get data do be plotted
            # prec: det, obs or ens mean
            # prec_lower/upper: ensemble minimum, maximum
            if group == 'ens':
                array = rootgroup.groups[group].variables[var][iday, :, :]
                mean = np.mean(array, axis=1)
                lower = np.amin(array, axis=1)
                upper = np.amax(array, axis=1)
            else:
                mean = rootgroup.groups[group].variables[var][iday, :, 0]

            # Plot data
            ax.plot(rootgroup.variables['time'][:], mean, label=group,
                    c=get_config(inargs, 'colors', group), ls=ls)
            if group == 'ens':
                ax.fill_between(rootgroup.variables['time'][:],
                                lower, upper,
                                where=upper >= lower,
                                facecolor=get_config(inargs, 'colors',
                                                     'ens_range'))


def plot_domain_mean_timeseries_composite(inargs, plot_type):
    """
    Function to plot time series of domain mean as a composite over 
    all days

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    plot_type : str
      Type of plot. Must be 'precipitation' or 'cape_tauc'

    """

    assert plot_type in ['precipitation', 'cape_tauc'], \
        'Type must be precipitation or cape_tauc'

    # Read pre-processed data
    rootgroup = read_netcdf_dataset(inargs)
    x = rootgroup.variables['time'][:]

    fig, ax1 = plt.subplots(1, 1, figsize=(3, 3))

    if plot_type == 'precipitation':
        ax1.set_ylabel('Accumulation [mm/h]')
        for group in rootgroup.groups:
            array = rootgroup.groups[group].variables['PREC_ACCUM'][:]
            mean = np.mean(array, axis=(0,2))
            ax1.plot(x, mean, label=group,
                    c=get_config(inargs, 'colors', group))

    if plot_type == 'cape_tauc':
        ax1.set_ylabel('CAPE [J/kg]')
        ax2 = ax1.twinx()
        ax2.set_ylabel('tau_c [h]')
        for group in ['det', 'ens']:
            for var, ax, ls in zip(['CAPE_ML', 'TAU_C'],
                                   [ax1, ax2],
                                   ['-', '--']):
                array = rootgroup.groups[group].variables[var][:]
                mean = np.mean(array, axis=(0, 2))
                ax.plot(x, mean, label=group,
                        c=get_config(inargs, 'colors', group), ls=ls)

    ax1.set_xlabel('Time [UTC]')
    dateobj_start = (timedelta(seconds=int(rootgroup.variables['date'][0])) +
                     datetime(1,1,1))
    datestr_start = dateobj_start.strftime(get_config(inargs, 'plotting',
                                                      'date_fmt'))
    dateobj_end = (timedelta(seconds=int(rootgroup.variables['date'][-1])) +
                   datetime(1, 1, 1))
    datestr_end = dateobj_end.strftime(get_config(inargs, 'plotting',
                                                  'date_fmt'))
    comp_str = 'Composite ' + datestr_start + ' - ' + datestr_end

    ax1.set_title(comp_str)
    ax1.legend(loc=0)

    plt.tight_layout()

    # Save figure and log
    save_fig_and_log(fig, rootgroup, inargs, plot_type + '_ts_composite')


def plot_prec_stamps(inargs):

    cmPrec = ((1, 1, 1),
              (0, 0.627, 1),
              (0.137, 0.235, 0.98),
              (0.392, 0, 0.627),
              (0.784, 0, 0.627))
    # (0.1  , 0.1   , 0.784),
    levelsPrec = [0, 1, 3, 10, 30, 100.]

    # Loop over dates
    for idate, date in enumerate(make_datelist(inargs)):

        # Loop over times
        for t in make_timelist(timedelta(hours=inargs.time_start),
                               timedelta(hours=inargs.time_end),
                               timedelta(hours=1)):

            # Load all CU objects
            fobjlist = []
            titlelist = []

            # 1st: Radar
            radarpref = (get_config(inargs, 'paths', 'radar_data') +
                         get_config(inargs, 'paths', 'radar_prefx'))
            radarsufx = get_config(inargs, 'paths', 'radar_sufix')
            dtradar = timedelta(minutes=10)
            t_rad = yyyymmddhh_strtotime(date) + t - dtradar
            RADARobj = getfobj_ncdf(radarpref + yymmddhhmm(t_rad) +
                                    radarsufx, 'pr', dwdradar=True)
            # Crop radar field
            RADARobj.data = RADARobj.data[62:62 + 357, 22:22 + 357]
            RADARobj.lats = RADARobj.lats[62:62 + 357, 22:22 + 357]
            RADARobj.lons = RADARobj.lons[62:62 + 357, 22:22 + 357]
            RADARobj.rlats = RADARobj.rlats[62:62 + 357, 22:22 + 357]
            RADARobj.rlons = RADARobj.rlons[62:62 + 357, 22:22 + 357]
            RADARobj.nx = 357
            RADARobj.ny = 357

            fobjlist.append(RADARobj)
            titlelist.append('Radar')

            # 2nd: det
            ncdffn_pref = (get_config(inargs, 'paths', 'raw_data') +
                           date + '/deout_ceu_pspens/' + 'det' +
                           '/OUTPUT/lfff')
            detfobj = getfobj_ncdf(ncdffn_pref + ddhhmmss(t) + '.nc_30m_surf',
                                   'PREC_ACCUM')

            fobjlist.append(detfobj)
            titlelist.append('Det')

            # 3rd: ens
            ncdffn = 'lfff' + ddhhmmss(t) + '.nc_30m_surf'
            date_dir = (get_config(inargs, 'paths', 'raw_data') +
                           date + '/deout_ceu_pspens/')
            enslist = getfobj_ncdf_ens(date_dir, 'sub', inargs.nens, ncdffn,
                                       dir_suffix='/OUTPUT/',
                                       fieldn = 'PREC_ACCUM', nfill=1)
            fobjlist.extend(enslist)
            titlelist.extend(['Mem ' + str(i+1) for i in range(inargs.nens)])

            # Now plot
            n_panels = len(fobjlist)
            n_cols = 4
            n_rows = int(np.ceil(float(n_panels) / n_cols))
            fig, axmat = plt.subplots(n_rows, n_cols, figsize=(10, 3.5 * n_rows))
            axflat = np.ravel(axmat)

            for i in range(len(fobjlist)):
                plt.sca(axflat[i])
                cf, tmp = ax_contourf(axflat[i], fobjlist[i], colors=cmPrec,
                                      pllevels=levelsPrec,
                                      ji0=(50 + inargs.zoom_lat1,
                                           50 + inargs.zoom_lon1),
                                      ji1=(357-51 + inargs.zoom_lat2,
                                           357-51 + inargs.zoom_lon2),
                                      sp_title=titlelist[i],
                                      Basemap_drawrivers=False,
                                      npars=0, nmers=0)
            cb = fig.colorbar(cf, cax=fig.add_axes([0.4, 0.1, 0.2, 0.02]),
                              orientation='horizontal')
            cb.set_label('Accumulation [mm/h]')
        titlestr = (yyyymmddhh_strtotime(date).strftime(get_config(inargs,
                                                                   'plotting',
                                                                   'date_fmt'))+
                    ' ' + str(t.seconds / 3600).zfill(2) + 'UTC')
        fig.suptitle(titlestr)
        plt.tight_layout(rect=[0, 0.1, 1, 0.93])

        # Save figure and log
        fig.savefig('/home/s/S.Rasp/tmp/test_stamps.pdf')



def plotting(inargs):
    """
    Top-level function called by main.py
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------

    """

    # Call appropriate plotting function
    if inargs.plot == 'weather_ts':
        plot_domain_mean_timeseries_individual(inargs,
                                               plot_type='precipitation')
        plot_domain_mean_timeseries_composite(inargs,
                                              plot_type='precipitation')
        plot_domain_mean_timeseries_individual(inargs,
                                               plot_type='cape_tauc')
        plot_domain_mean_timeseries_composite(inargs,
                                              plot_type='cape_tauc')
    if inargs.plot == 'prec_stamps':
        plot_prec_stamps(inargs)
