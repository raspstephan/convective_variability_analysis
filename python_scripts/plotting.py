"""
Filename:     plotting.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains functions to read pre-processed NetCDF files, do some
final analysis and plot the results.

"""

# Import modules
from helpers import get_pp_fn
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np



def read_netcdf_dataset(inargs):
    """
    Open NetCDF file and return rootgroup object.
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------
    rootgroup : NetCDF object
      NetCDF object
    """
    pp_fn = get_pp_fn(inargs)
    rootgroup = Dataset(pp_fn)
    return rootgroup


def plot_domain_mean_weather_ts(inargs):
    """
    Function to plot time series of domain mean weather.
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------

    """

    # Read pre-processed data
    rootgroup = read_netcdf_dataset(inargs)
    n_days = rootgroup.dimensions['date'].size
    x = range(1, rootgroup.dimensions['time'].size + 1)

    # Set up figure
    n_cols = 4
    n_rows = int(np.ceil(float(n_days) / n_cols))

    fig, axmat = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)   # TODO: figsize
    axflat = np.ravel(axmat)

    # Loop over axes / days
    for iday in range(n_days):
        for group in rootgroup.groups:
            # Get data do be plotted
            # prec: det, obs or ens mean
            # prec_lower/upper: ensemble minimum, maximum
            if group == 'ens':
                prec_array = rootgroup.groups[group].variables['PREC_ACCUM']\
                    [iday, :,:]
                prec = np.mean(prec_array, axis=1)
                prec_lower = np.amin(prec_array, axis=1)
                prec_upper = np.amax(prec_array, axis=1)
            else:
                prec = rootgroup.groups[group].variables['PREC_ACCUM']\
                    [iday, :, 0]

            # Plot data
            axflat[iday].plot(x, prec, label=group)
            if group == 'ens':
                axflat[iday].fill_between(x, prec_lower, prec_upper,
                                          where=prec_upper >= prec_lower)

    # Save figure
    fig.savefig('/home/s/S.Rasp/repositories/convective_variability_analysis/figures/test_ts.pdf')












def plotting(inargs):
    """
    Top-level function called by main.py
    
    Parameters
    ----------
    inargs

    Returns
    -------

    """

    # Call appropriate plotting function
    plot_domain_mean_weather_ts(inargs)