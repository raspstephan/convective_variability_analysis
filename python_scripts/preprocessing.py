"""
Filename:     preprocessing.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains functions to pre-process data and save an intermediate
netCDF File.

"""

# Import modules
from netCDF4 import Dataset
from cosmo_utils.pyncdf import getfobj_ncdf_timeseries
from helpers import get_config, make_datelist_yyyymmddhh
from datetime import timedelta
import numpy as np

# Define functions
def create_netcdf_weather_ts(pp_fn, inargs):
    """
    Creates a NetCDF object to store weather time series data.
    3 groups : obs, det, ens
    3 dimensions : date, time, ens_no (1 for det and obs)
    4 variables : mean_prec, mean_cape, mean_tauc, mean_hpbl
    
    Parameters
    ----------
    pp_fn : str
      Filename with path of pre-processed NetCDF file
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------
    rootgroup : NetCDF object

    """
    # Create NetCDF file
    rootgroup = Dataset(pp_fn, 'w', format='NETCDF4')

    groups = ['obs', 'det', 'ens']
    dimensions = {
        'time': 24,
        'date': int(inargs.date_end) - int(inargs.date_end) + 1,
    }
    variables = {
        'mean_prec': ['date', 'time'],
        'mean_cape': ['date', 'time'],
        'mean_tauc': ['date', 'time'],
        'mean_hpbl': ['date', 'time'],
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


def domain_mean_weather_ts(inargs, pp_fn):
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
    pp_fn : str
      Filename with path of pre-processed NetCDF file

    Returns
    -------

    """

    rootgroup = create_netcdf_weather_ts(pp_fn, inargs)

    # Load analysis data and store in NetCDF
    for id, date in enumerate(make_datelist_yyyymmddhh(inargs)):
        for group in rootgroup.groups:
            for ie in range(rootgroup.groups[group].dimensions['ens_no'].size):
                if group in ['det', 'ens']:
                    ens_str = str(ie + 1)
                    if group == 'det':
                        ens_str = 'det'

                    ncdffn_pref = (get_config(inargs, 'paths', 'raw_data') +
                                   date + '/deout_ceu_pspens/' + ens_str +
                                   '/OUTPUT/lfff')

                    preclist = getfobj_ncdf_timeseries(ncdffn_pref,
                                                       timedelta(hours=inargs.time_start),
                                                       timedelta(hours=inargs.time_end),
                                                       timedelta(hours=inargs.time_inc),
                                                       ncdffn_sufx='.nc_30m_surf',
                                                       return_arrays=True,
                                                       fieldn='TOT_PREC')

                    # Compute domain mean and save in NetCDF file
                    mean_ts = np.mean(preclist, axis=(1, 2))
                    rootgroup.groups[group].variables['mean_prec'][id, :, ie] =\
                        mean_ts

    # Close NetCDF file
    rootgroup.close()



def preprocess(inargs, pp_fn):
    """
    Top-level function called by main.py

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    pp_fn : str
      Filename with path of pre-processed NetCDF file

    Returns
    -------

    """

    # Call analysis function
    domain_mean_weather_ts(inargs, pp_fn)