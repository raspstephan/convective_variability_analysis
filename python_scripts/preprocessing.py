"""
Filename:     preprocessing.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains functions to pre-process data and save an intermediate
netCDF File.

"""

# Import modules
from netCDF4 import Dataset
from helpers import make_datelist, get_radar_mask, get_pp_fn, \
    get_datalist_radar, create_log_str, get_datalist_model
import numpy as np


# Define functions
# Weather time series
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
        datalist = get_datalist_model(inargs, date, ens_no, var, radar_mask)
    elif group == 'obs':
        if not var == 'PREC_ACCUM':
            return
        datalist = get_datalist_radar(inargs, date, radar_mask)
    else:
        raise Exception('Wrong group.')

    # Compute domain mean and save in NetCDF file
    # Note: Need to loop, because array operation ignores mask
    mean_ts = []
    for data in datalist:
        mean_ts.append(np.mean(data))
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
    log_str : str
      Log text for NetCDF file

    Returns
    -------

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
    variables = {
        'PREC_ACCUM': ['date', 'time'],
        'CAPE_ML': ['date', 'time'],
        'TAU_C': ['date', 'time'],
        'HPBL': ['date', 'time'],
    }
    rootgroup = create_netcdf(inargs, groups, dimensions, variables,
                                         ensemble_dim=True)

    radar_mask = get_radar_mask(inargs)
    print('Number of masked grid points: ' + str(np.sum(radar_mask)) +
          ' from total grid points: ' + str(radar_mask.size))

    # Load analysis data and store in NetCDF
    for idate, date in enumerate(make_datelist(inargs)):
        print('Computing time series for: ' + date)
        for group in rootgroup.groups:
            for ie in range(rootgroup.groups[group].dimensions['ens_no'].size):
                for var in rootgroup.groups[group].variables:

                    compute_ts_mean(inargs, idate, date, group, ie, var,
                                    rootgroup, radar_mask)

    # Close NetCDF file
    rootgroup.close()


# precipitation histogram
def prec_hist(inargs):

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
    rootgroup = create_netcdf(inargs, groups, dimensions, variables,
                              ensemble_dim=True)

    # TODO: This is somewhat the same as domain_mean_weather_ts
    radar_mask = get_radar_mask(inargs)
    print('Number of masked grid points: ' + str(np.sum(radar_mask)) +
          ' from total grid points: ' + str(radar_mask.size))

    # Load analysis data and store in NetCDF
    for idate, date in enumerate(make_datelist(inargs)):
        print('Computing time series for: ' + date)
        for group in rootgroup.groups:
            # TODO This is copied from compute_ts_mean
            if group in ['det', 'ens']:
                if group == 'det':
                    ens_no = 'det'
                else:
                    ens_no = ie + 1
                datalist = get_datalist_model(inargs, date, ens_no, 'PREC_ACCUM',
                                              radar_mask)
            elif group == 'obs':
                datalist = get_datalist_radar(inargs, date, radar_mask)
            else:
                raise Exception('Wrong group.')

            # Now do the actually new calculation
            for it, data in enumerate(datalist):
                rootgroup.groups[group].variables['prec_hist'][idate, it, :] =\
                asdf


def preprocess(inargs):
    """
    Top-level function called by main.py

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    log_str : str
      Log text for NetCDF file

    Returns
    -------

    """

    # Call analysis function
    if inargs.plot == 'weather_ts':
        domain_mean_weather_ts(inargs)
    elif inargs.plot == 'prec_hist':
        prec_hist(inargs)
    else:
        print('No preprocessing necessary.')
