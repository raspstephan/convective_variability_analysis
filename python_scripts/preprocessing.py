"""
Filename:     preprocessing.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains functions to pre-process data and save an intermediate
netCDF File.

"""

# Import modules
from netCDF4 import Dataset

# Define functions

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
    for g in groups:
        rootgroup.createGroup(g)
        if g == 'ens':
            dimensions['ens_no'] = inargs.nens
            [b.append('ens_no') for a, b in variables.items()]

        # Create dimensions
        for dim_name, dim_len in dimensions.items():
            rootgroup.groups[g].createDimension(dim_name, dim_len)

        # Create variables
        for var_name, var_dims in variables.items():
            rootgroup.groups[g].createVariable(var_name, 'f8', var_dims)
    print rootgroup
    # rootgroup = create_netcdf(pp_fn, dimensions, variables)
    # Load analysis data

    # Compute means

    # Save in NetCDF file
    rootgroup.groups['ens'].variables['mean_prec'][0,0,0] = -99   # TEST

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