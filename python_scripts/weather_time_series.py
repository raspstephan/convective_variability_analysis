"""
Filename:     weather_time_series.py
Author:       Stephan Rasp, s.rasp@lmu.de
Description:  Analyze model output and radar observations and compute 
              domain mean time series of precipitation, CAPE and tau_c

"""

# Import modules
import argparse
from netCDF4 import Dataset
from datetime import datetime, timedelta
from helpers import make_datelist, get_radar_mask, get_pp_fn, \
    get_datalist_radar, create_log_str, get_datalist_model, \
    read_netcdf_dataset, get_config, save_fig_and_log, pp_exists, \
    get_composite_str, load_raw_data
import numpy as np
import matplotlib.pyplot as plt

################################################################################
# PREPROCESSING FUNCTIONS
################################################################################


def create_netcdf(inargs, groups, dimensions, variables,
                  ensemble_dim=False):
    """
    Creates a NetCDF object to store data.

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    groups : list
      List of groups
    dimensions : dict
      Dictionary with dimension name as key and dimension value as value
    variables : dict
      Dictionary with variable name as key and variable dimension list as value
      (must be numpy array). Attention: Type i8 is hardcoded for all dimensions.
    ensemble_dim : bool
      If true, group variables have additional dimension ens_no with size 
      inargs.nens for group 'ens' and 1 for all other groups

    Returns
    -------
    rootgroup : NetCDF object

    """

    pp_fn = get_pp_fn(inargs)

    # Create NetCDF file
    rootgroup = Dataset(pp_fn, 'w', format='NETCDF4')
    rootgroup.log = create_log_str(inargs, 'Preprocessing')

    # Create root dimensions and variables
    for dim_name, dim_val in dimensions.items():
        rootgroup.createDimension(dim_name, dim_val.shape[0])
        tmp_var = rootgroup.createVariable(dim_name, 'i8', dim_name)
        tmp_var[:] = dim_val

    # Create group dimensions and variables
    if ensemble_dim:
        [b.append('ens_no') for a, b in variables.items()]
        dimensions['ens_no'] = 1

    for g in groups:
        rootgroup.createGroup(g)
        if g == 'ens' and ensemble_dim:
            dimensions['ens_no'] = inargs.nens

        # Create dimensions
        for dim_name, dim_len in dimensions.items():
            if type(dim_len) is not int:
                dim_len = dim_len.shape[0]
            rootgroup.groups[g].createDimension(dim_name, dim_len)

        # Create variables
        for var_name, var_dims in variables.items():
            rootgroup.groups[g].createVariable(var_name, 'f8', var_dims)
    return rootgroup


def compute_ts_mean(inargs, idate, date, group, ie, var, rootgroup,
                    radar_mask):
    """
    Compute mean time series and appends it to rootgroup object.

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    idate : int
      index of date
    date : str
      Date to analyze
    group : str
      Group name
    ie  : int 
      Ensemble Member
    var : str 
      COSMO variable to analyze
    rootgroup : NetCDF Dataset object
      NetCDF rootgroup
    radar_mask : np.array
      Total radar mask

    """
    if group in ['det', 'ens']:
        if group == 'det':
            ens_no = 'det'
        else:
            ens_no = ie + 1
        datalist = get_datalist_model(inargs, date, ens_no, var)
    elif group == 'obs':
        if not var == 'PREC_ACCUM':
            return
        datalist = get_datalist_radar(inargs, date)
    else:
        raise Exception('Wrong group.')



    # Compute domain mean and save in NetCDF file
    # Note: Need to loop, because array operation ignores mask
    mean_ts = []
    for data in datalist:
        mean_ts.append(np.mean(data[~radar_mask]))
    rootgroup.groups[group].variables[var][idate, :, ie] = np.array(mean_ts)


def domain_mean_weather_ts(inargs):
    """
    Calculate hourly time-series for domain mean variables:

    - hourly precipitation
    - CAPE
    - convective adjustment timescale
    - boundary layer height
    Precipitation is analyzed for ensemble, deterministic and observations.
    All other values are calculated for the ensemble mean and deterministic.

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    """

    # Define NetCDF parameters and create rootgroup
    groups = ['obs', 'det', 'ens']
    datearray = np.array(make_datelist(inargs, out_format='netcdf'))
    timearray = np.arange(inargs.time_start, inargs.time_end + inargs.time_inc,
                          inargs.time_inc)
    dimensions = {
        'time': timearray,
        'date': datearray,
    }
    if inargs.compute_tauc_hpbl:
        variables = {
            'PREC_ACCUM': ['date', 'time'],
            'CAPE_ML': ['date', 'time'],
            'TAU_C': ['date', 'time'],
            'HPBL': ['date', 'time'],
        }
    else:
        variables = {
            'PREC_ACCUM': ['date', 'time'],
            'CAPE_ML': ['date', 'time'],
        }
    rootgroup = create_netcdf(inargs, groups, dimensions, variables,
                              ensemble_dim=True)

    radar_mask = get_radar_mask(inargs)
    print('Number of masked grid points: ' + str(np.sum(radar_mask)) +
          ' from total grid points: ' + str(radar_mask.size))

    # Load analysis data and store in NetCDF
    for idate, date in enumerate(make_datelist(inargs)):
        print('Computing time series for: ' + date)

        # Determine radar mask
        if inargs.radar_mask == 'total':
            tmp_mask = np.any(radar_mask, axis=(0, 1))
        elif inargs.radar_mask == 'day':
            tmp_mask = np.any(radar_mask[idate], axis=0)
        elif inargs.radar_mask == 'hour':
            tmp_mask = radar_mask[idate]

        for group in rootgroup.groups:
            for ie in range(rootgroup.groups[group].dimensions['ens_no'].size):
                for var in rootgroup.groups[group].variables:
                    compute_ts_mean(inargs, idate, date, group, ie, var,
                                    rootgroup, tmp_mask)

    # Close NetCDF file
    rootgroup.close()


################################################################################
# PLOTTING FUNCTIONS
################################################################################
def plot_domain_mean_timeseries_individual(inargs, plot_var):
    """
    Function to plot time series of domain mean precipitation for each day

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    plot_var : str
      Type of plot. Must be 'precipitation' or 'cape_tauc'

    """

    assert plot_var in ['precipitation', 'cape_tauc'], \
        'Type must be precipitation or cape_tauc'

    # Read pre-processed data
    rootgroup = read_netcdf_dataset(inargs)
    n_days = rootgroup.dimensions['date'].size

    # Set up figure
    n_cols = 4
    n_rows = int(np.ceil(float(n_days) / n_cols))

    fig, axmat = plt.subplots(n_rows, n_cols, sharex=True, sharey=True,
                              figsize=(get_config(inargs, 'plotting',
                                                  'page_width'),
                                       2.1 * n_rows))
    axflat = np.ravel(axmat)

    # Loop over axes / days
    for iday in range(n_days):
        ax = axflat[iday]
        dateobj = (timedelta(seconds=int(rootgroup.variables['date'][iday])) +
                   datetime(1, 1, 1))
        datestr = dateobj.strftime(get_config(inargs, 'plotting', 'date_fmt'))
        if not iday == 0:
            ax.text(0.3, 0.9, datestr, fontsize=10, transform = ax.transAxes)
        else:
            ax.text(0.65, 0.9, datestr, fontsize=10, transform = ax.transAxes)
        if iday >= ((n_cols * n_rows) - n_cols):  # Only bottom row
            ax.set_xlabel('Time [UTC]')
            ax.set_xticks([0, 6, 12, 18, 24])

        if plot_var == 'precipitation':
            plot_precipitation_panel(inargs, ax, iday, rootgroup)
        if plot_var == 'cape_tauc':
            plot_cape_tauc_panel(inargs, ax, iday, rootgroup)

        ax.set_ylim(0, inargs.ymax)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # ax.spines['left'].set_position(('outward', 3))
        ax.spines['bottom'].set_position(('outward', 3))

    # Finish figure
    axflat[0].legend(loc=2)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.15, hspace=0.25)

    # Save figure
    save_fig_and_log(fig, rootgroup, inargs, plot_var + '_ts_individual')


def plot_precipitation_panel(inargs, ax, iday, rootgroup):
    """
    Plots precipitation timeseries in each panel.

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    ax : axis object
      Axis objects
    iday : int
      Index for list object
    rootgroup : netCDF object
      NetCDF object with data to plot

    """

    dateobj = (timedelta(seconds=int(rootgroup.variables['date'][iday])) +
               datetime(1, 1, 1))

    if iday % 4 == 0:  # Only left column
        ax.set_ylabel(r'Precip [mm h$^{-1}$]')

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
        ax.plot(rootgroup.variables['time'][:], prec, label=group,
                c=get_config(inargs, 'colors', group), linewidth=2)
        if group == 'ens':
            ax.fill_between(rootgroup.variables['time'][:],
                            prec_lower, prec_upper,
                            where=prec_upper >= prec_lower,
                            facecolor=get_config(inargs,
                                                 'colors',
                                                 'ens_range'))


def plot_cape_tauc_panel(inargs, ax, iday, rootgroup):
    """
    Plots cape and tau_c timeseries in each panel.

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    ax : axis object
      Axis object
    iday : int
      Index for list object
    rootgroup : netCDF object
      NetCDF object with data to plot

    """

    ax1 = ax
    ax2 = ax1.twinx()
    ax2.set_ylim(0, 20)

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


def plot_domain_mean_timeseries_composite(inargs, plot_var):
    """
    Function to plot time series of domain mean as a composite over 
    all days

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    plot_var : str
      Type of plot. Must be 'precipitation' or 'cape_tauc'

    """

    assert plot_var in ['precipitation', 'cape_tauc', 'prec_cape'], \
        'Wrong type!'

    # Read pre-processed data
    rootgroup = read_netcdf_dataset(inargs)
    x = rootgroup.variables['time'][:]

    pw = get_config(inargs, 'plotting', 'page_width')
    ratio = 0.7
    fig, ax1 = plt.subplots(1, 1, figsize=(pw / 2., pw / 2. * ratio))

    if plot_var in ['precipitation', 'prec_cape']:
        ax1.set_ylabel(r'Precip [mm h$^{-1}$]')
        for group in rootgroup.groups:
            array = rootgroup.groups[group].variables['PREC_ACCUM'][:]
            mean = np.mean(array, axis=(0, 2))
            std = np.mean(np.std(array, axis=2, ddof=1), axis=0)
            ax1.plot(x, mean, label=group,
                     c=get_config(inargs, 'colors', group),
                     linewidth=2)
            if group == 'ens':
                lower = mean - std
                upper = mean + std
                ax1.fill_between(x, lower, upper, where=upper >= lower,
                                 facecolor=get_config(inargs, 'colors',
                                                      'ens_range'))
        if plot_var == 'prec_cape':
            ax2 = ax1.twinx()
            ax2.set_ylabel(r'CAPE [J kg$^{-1}$]')

            array = rootgroup.groups['ens'].variables['CAPE_ML'][:]
            mean = np.mean(array, axis=(0, 2))
            ax2.plot(x, mean, c=get_config(inargs, 'colors', 'third'),
                     linestyle='--', linewidth=2)

    if plot_var == 'cape_tauc':
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
    ax1.set_xticks([0, 6, 12, 18, 24])
    ax1.axvline(6, c='lightgray', linewidth=0.5, zorder = 0.001)
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_position(('outward', 3))
        ax.spines['left'].set_position(('outward', 3))
        ax.spines['right'].set_position(('outward', 3))
        ax.set_xlim((0, 24))

    # comp_str = 'Composite ' + get_composite_str(inargs, rootgroup)
    ax1.set_title('Composite precipitation')
    ax1.legend(loc=0, fontsize=8, title='Precip')

    plt.subplots_adjust(left=0.18, right=0.82, bottom=0.2, top=0.9)

    # Save figure and log
    save_fig_and_log(fig, rootgroup, inargs, plot_var + '_ts_composite')


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

    # Check if pre-processed file exists
    if (pp_exists(inargs) is False) or (inargs.recompute is True):
        print('Compute preprocessed file: ' + get_pp_fn(inargs))
        # Call preprocessing routine with arguments
        domain_mean_weather_ts(inargs)
    else:
        print('Found pre-processed file: ' + get_pp_fn(inargs))

    # Call analyzing and plotting routine
    if 'prec_ind' in inargs.plot_type:
        plot_domain_mean_timeseries_individual(inargs,
                                               plot_var='precipitation')
    if 'prec_comp' in inargs.plot_type:
        plot_domain_mean_timeseries_composite(inargs,
                                              plot_var='precipitation')
    if 'cape_tauc_ind' in inargs.plot_type:
        plot_domain_mean_timeseries_individual(inargs,
                                               plot_var='cape_tauc')
    if 'cape_tauc_comp' in inargs.plot_type:
        plot_domain_mean_timeseries_composite(inargs,
                                              plot_var='cape_tauc')
    if 'prec_cape_comp' in inargs.plot_type:
        plot_domain_mean_timeseries_composite(inargs,
                                              plot_var='prec_cape')


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
    parser.add_argument('--radar_mask',
                        type=str,
                        default='hour',
                        help='Radar mask for [hour, day, total]?')
    parser.add_argument('--compute_tauc_hpbl',
                        dest='compute_tauc_hpbl',
                        action='store_true',
                        help='If given, compute mean tau_c and hpbl')
    parser.set_defaults(compute_tauc_hpbl=False)



    # General plotting options
    parser.add_argument('--plot_type',
                        type=str,
                        default='',
                        help='Which plot to plot. [prec_ind, prec_comp, '
                             'cape_tauc_ind, cape_tauc_comp, prec_cape_comp]')
    parser.add_argument('--plot_format',
                        type=str,
                        default='pdf',
                        help='Which format for figure file.')
    parser.add_argument('--ymax',
                        type=float,
                        default=0.9,
                        help='Max for y axis')

    # General settings
    parser.add_argument('--config_file',
                        type=str,
                        default='config.yml',
                        help='Config file in relative directory ../config. \
                                  Default = config.yml')
    parser.add_argument('--sub_dir',
                        type=str,
                        default='weather_time_series',
                        help='Sub-directory for figures and pp_files')
    parser.add_argument('--plot_name',
                        type=str,
                        default='',
                        help='Custom plot name.')
    parser.add_argument('--pp_name',
                        type=str,
                        default='',
                        help='Custom name for preprocessed file.')
    parser.add_argument('--recompute',
                        dest='recompute',
                        action='store_true',
                        help='If True, recompute pre-processed file.')
    parser.set_defaults(recompute=False)

    args = parser.parse_args()

    main(args)