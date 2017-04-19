"""
Filename:     prec_hist.py
Author:       Stephan Rasp, s.rasp@lmu.de
Description:  Compute and plot precipitation histograms of deterministic and 
              ensemble runs and observations

"""

# Import modules
import argparse
from netCDF4 import Dataset
from helpers import make_datelist, get_radar_mask, get_pp_fn, \
    get_datalist_radar, create_log_str, get_datalist_model, \
    read_netcdf_dataset, get_config, save_fig_and_log, pp_exists, \
    get_composite_str
import numpy as np
import matplotlib.pyplot as plt


################################################################################
# PREPROCESSING FUNCTIONS
################################################################################
def create_netcdf(inargs, groups, dimensions, variables):
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
        tmp_var = rootgroup.createVariable(dim_name, 'f8', dim_name)
        tmp_var[:] = dim_val

    # Create group dimensions and variables
    [b.append('ens_no') for a, b in variables.items()]
    dimensions['ens_no'] = 1

    for g in groups:
        rootgroup.createGroup(g)
        if g == 'ens':
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


def prec_stats(inargs):
    """
    Compute and save precipitation amount and cloud size and cloud 
    precipitation histograms.
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments


    """

    # Define bins TODO: Read from config!
    histbinedges = [0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 1000]

    # Make netCDF file
    datearray = np.array(make_datelist(inargs, out_format='netcdf'))
    timearray = np.arange(inargs.time_start, inargs.time_end + inargs.time_inc,
                          inargs.time_inc)
    groups = ['obs', 'det', 'ens']
    dimensions = {
        'time': timearray,
        'date': datearray,
        'bins': np.array(histbinedges[1:]),
    }
    variables = {
        'prec_hist': ['date', 'time', 'bins'],
    }
    rootgroup = create_netcdf(inargs, groups, dimensions, variables)

    # TODO: This is somewhat the same as domain_mean_weather_ts
    radar_mask = get_radar_mask(inargs)
    print('Number of masked grid points: ' + str(np.sum(radar_mask)) +
          ' from total grid points: ' + str(radar_mask.size))

    # Load analysis data and store in NetCDF
    for idate, date in enumerate(make_datelist(inargs)):
        print('Computing prec_hist for: ' + date)
        for group in rootgroup.groups:
            for ie in range(rootgroup.groups[group].dimensions['ens_no'].size):
                if group in ['det', 'ens']:
                    if group == 'det':
                        ens_no = 'det'
                    else:
                        ens_no = ie + 1
                    datalist = get_datalist_model(inargs, date, ens_no,
                                                  'PREC_ACCUM', radar_mask)
                elif group == 'obs':
                    datalist = get_datalist_radar(inargs, date, radar_mask)
                else:
                    raise Exception('Wrong group.')

                # Now do the actually new calculation
                for it, data in enumerate(datalist):

                    rootgroup.groups[group].variables['prec_hist']\
                        [idate, it, :, ie] = np.histogram(data, histbinedges)[0]

    # Close NetCDF file
    rootgroup.close()


################################################################################
# PLOTTING FUNCTIONS
################################################################################
def plot_hist(inargs):
    """
    Plot precipitation histogram comparing the three groups
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    """

    # Read pre-processed data
    rootgroup = read_netcdf_dataset(inargs)

    # Set up figure
    fig, ax = plt.subplots(1, 1, figsize=(4, 3.5))
    x = np.arange(rootgroup.variables['bins'][:].shape[0])

    # Loop over groups
    for ig, group in enumerate(rootgroup.groups):
        # Compute mean in all directions but bins
        mean_hist = np.mean(rootgroup.groups[group].variables['prec_hist'][:],
                            axis=(0, 1, 3))
        ax.bar(x[1:] + ig * 0.2, mean_hist[1:], width=0.2,
               color=get_config(inargs, 'colors', group), label=group)

    # Make figure look nice
    ax.legend(loc=0, prop={'size': 10})
    plt.xticks(x[1:], rootgroup.variables['bins'][:-1])
    ax.set_xlabel('Hourly accumulation [mm/h]')
    ax.set_ylabel('Number of grid points')
    date_str = get_composite_str(inargs, rootgroup)
    ax.set_title(date_str, fontsize=12)

    plt.tight_layout()

    # Save figure and log
    save_fig_and_log(fig, rootgroup, inargs, 'prec_hist')


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
        prec_hist(inargs)
    else:
        print('Found pre-processed file:' + get_pp_fn(inargs))

    # Plotting
    plot_hist(inargs)


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
    parser.add_argument('--sub_dir',
                        type=str,
                        default='prec_hist',
                        help='Sub-directory for figures and pp_files')
    parser.add_argument('--recompute',
                        dest='recompute',
                        action='store_true',
                        help='If True, recompute pre-processed file.')
    parser.set_defaults(recompute=False)

    args = parser.parse_args()

    main(args)