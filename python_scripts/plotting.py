"""
Filename:     plotting.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains functions to read pre-processed NetCDF files, do some
final analysis and plot the results.

"""

# Import modules
from helpers import read_netcdf_dataset, get_config, get_pp_fn
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np


# Define functions
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
    x = rootgroup.variables['time'][:]

    # Set up figure
    n_cols = 4
    n_rows = int(np.ceil(float(n_days) / n_cols))

    fig, axmat = plt.subplots(n_rows, n_cols, sharex=True, sharey=True,
                              figsize=(10,3*n_rows))   # TODO: figsize
    axflat = np.ravel(axmat)

    # Loop over axes / days
    for iday in range(n_days):
        dateobj = (timedelta(seconds=int(rootgroup.variables['date'][iday])) +
                   datetime(1,1,1))
        datestr = dateobj.strftime(get_config(inargs, 'colors', 'date_fmt'))
        axflat[iday].set_title(datestr)
        if iday % 4 == 0:   # Only left column
            axflat[iday].set_ylabel('Accumulation [mm/h]')
        if iday >= ((n_cols * n_rows) - n_cols):   # Only bottom row
            axflat[iday].set_xlabel('Time [UTC]')

        for group in rootgroup.groups:
            # Get data do be plotted
            # prec: det, obs or ens mean
            # prec_lower/upper: ensemble minimum, maximum
            if group == 'ens':
                prec_array = rootgroup.groups[group].variables['PREC_ACCUM']\
                    [iday, :, :]
                prec = np.mean(prec_array, axis=1)
                prec_lower = np.amin(prec_array, axis=1)
                prec_upper = np.amax(prec_array, axis=1)
            else:
                prec = rootgroup.groups[group].variables['PREC_ACCUM'] \
                    [iday, :, 0]

            # Plot data
            axflat[iday].plot(x, prec, label=group,
                              c=get_config(inargs, 'colors', group))
            if group == 'ens':
                axflat[iday].fill_between(x, prec_lower, prec_upper,
                                          where=prec_upper >= prec_lower,
                                          facecolor=get_config(inargs,
                                                               'colors',
                                                               'ens_range'))

    # Finish figure
    axflat[0].legend(loc=0)

    plt.tight_layout()
    # Save figure
    plotfn = (get_config(inargs, 'paths', 'figures') +
               get_pp_fn(inargs, sufx='.pdf', pure_fn=True))
    fig.savefig(plotfn)

    # Save log file
    logfn = (get_config(inargs, 'paths', 'figures') +
               get_pp_fn(inargs, sufx='.log', pure_fn=True))
    logf = open(logfn, 'w+')
    logf.write(rootgroup.log)
    logf.close()


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
    plot_domain_mean_weather_ts(inargs)
