"""
Filename:     preprocessing.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains functions to pre-process data and save an intermediate
netCDF File.

"""

# Import modules
from netCDF4 import Dataset
from cosmo_utils.pyncdf import getfobj_ncdf_timeseries
from helpers import get_config, make_datelist_yyyymmddhh, get_domain_limits, \
    get_radar_mask, get_pp_fn, get_datalist_radar, yyyymmddhh_strtotime
from datetime import timedelta
import numpy as np
from numpy.ma import masked_array


# Define functions
def create_netcdf_weather_ts(inargs, log_str):
    """
    Creates a NetCDF object to store weather time series data.
    
    3 groups : obs, det, ens
    3 dimensions : date, time, ens_no (1 for det and obs)
    4 variables : mean_prec, mean_cape, mean_tauc, mean_hpbl
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    log_str : str
      Log text for NetCDF file

    Returns
    -------
    rootgroup : NetCDF object

    """

    pp_fn = get_pp_fn(inargs)

    # Create NetCDF file
    rootgroup = Dataset(pp_fn, 'w', format='NETCDF4')
    rootgroup.log = log_str

    groups = ['obs', 'det', 'ens']
    n_days = (yyyymmddhh_strtotime(inargs.date_end) -
              yyyymmddhh_strtotime(inargs.date_start)).days + 1
    dimensions = {
        'time': 24,
        'date': n_days,
    }
    variables = {
        'PREC_ACCUM': ['date', 'time'],
        #'CAPE_ML': ['date', 'time'],
        #'TAU_C': ['date', 'time'],
        #'HPBL': ['date', 'time'],
    }

    # Create root dimensions and variables
    for dim_name, dim_len in dimensions.items():
        rootgroup.createDimension(dim_name, dim_len)

    rootgroup.createVariable('time', 'i4', ('time'))
    rootgroup.createVariable('date', 'i4', ('date'))

    # Create group dimensions and variables
    [b.append('ens_no') for a, b in variables.items()]
    dimensions['ens_no'] = 1

    for g in groups:
        rootgroup.createGroup(g)
        if g == 'ens':
            dimensions['ens_no'] = inargs.nens

        # Create dimensions
        for dim_name, dim_len in dimensions.items():
            rootgroup.groups[g].createDimension(dim_name, dim_len)

        # Create variables
        for var_name, var_dims in variables.items():
            rootgroup.groups[g].createVariable(var_name, 'f8', var_dims)
    return rootgroup


def get_datalist_model(inargs, date, ens_no, var, radar_mask):
    """
    Get data time series for model output.
    Parameters
    ----------
    inargs : : argparse object
      Argparse object with all input arguments
    date : str
      Date in format yyyymmddhh
    ens_no : int or str 
      Ensemble number or str in case of det
    var : str 
      Variable
    radar_mask : 2D numpy array
      Radar mask to create masked arrays

    Returns
    -------
    datalist : list
      List of 2D masked arrays
    """
    # Get file name
    ncdffn_pref = (get_config(inargs, 'paths', 'raw_data') +
                   date + '/deout_ceu_pspens/' + str(ens_no) +
                   '/OUTPUT/lfff')
    datalist = getfobj_ncdf_timeseries(ncdffn_pref,
                                       timedelta(hours=inargs.time_start),
                                       timedelta(hours=inargs.time_end),
                                       timedelta(hours=inargs.time_inc),
                                       ncdffn_sufx='.nc_30m_surf',
                                       return_arrays=True,
                                       fieldn=var)
    # Crop data
    l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad = \
        get_domain_limits(inargs)
    for i, data in enumerate(datalist):
        datalist[i] = masked_array(data[l11:l12, l21:l22],
                                   mask=radar_mask)
    return datalist


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


def domain_mean_weather_ts(inargs, log_str):
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

    rootgroup = create_netcdf_weather_ts(inargs, log_str)

    radar_mask = get_radar_mask(inargs)
    print('Number of masked grid points: ' + str(np.sum(radar_mask)) +
          ' from total grid points: ' + str(radar_mask.size))

    # Load analysis data and store in NetCDF
    for idate, date in enumerate(make_datelist_yyyymmddhh(inargs)):
        print('Computing time series for :' + date)
        for group in rootgroup.groups:
            for ie in range(rootgroup.groups[group].dimensions['ens_no'].size):
                for var in rootgroup.groups[group].variables:

                    compute_ts_mean(inargs, idate, date, group, ie, var,
                                    rootgroup, radar_mask)

    # Close NetCDF file
    rootgroup.close()


def preprocess(inargs, log_str):
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
    domain_mean_weather_ts(inargs, log_str)
