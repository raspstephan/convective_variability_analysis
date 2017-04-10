"""
Filename:     helpers.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains helper functions.

"""


# import modules
import os
import yaml
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, date2num
from datetime import datetime, timedelta
from cosmo_utils.helpers import yyyymmddhh_strtotime
from numpy.ma import masked_array
from cosmo_utils.pyncdf import getfobj_ncdf_timeseries


# Define functions
def get_config(inargs, top_key, bottom_key):
    """
    Reads the config JSON file

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    top_key : str
      Name of top category 
    bottom_key : key
      Name of key in top category

    Returns
    -------
    value : value for specified key pair
    """
    config = yaml.safe_load(open('../config/' + inargs.config_file))
    return config[top_key][bottom_key]


def make_datelist(inargs, out_format='yyyymmddhh'):
    """
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    out_format : str
      Output format. 'yyyymmddhh'[default] or 'netcdf'
 
    Returns
    -------
    datelist_yyyymmddhh : list
      List containing dates as strings. If date_start == date_end, returns list
      with one element.
    """
    dateob_start = datetime(int(inargs.date_start[0:4]),
                            int(inargs.date_start[4:6]),
                            int(inargs.date_start[6:8]),
                            int(inargs.date_start[8:10]))
    dateob_end = datetime(int(inargs.date_end[0:4]),
                          int(inargs.date_end[4:6]),
                          int(inargs.date_end[6:8]),
                          int(inargs.date_start[8:10]))
    date_inc = timedelta(days=1)

    datelist = []
    date = dateob_start
    while date <= dateob_end:
        if out_format == 'yyymmddhh':
            date_str = (str(date.year) + str(date.month).zfill(2) +
                        str(date.day).zfill(2) + str(date.hour).zfill(2))
            datelist.append(date_str)
        elif out_format == 'netcdf':
            datelist.append((date - datetime(1,1,1)).total_seconds())
        date += date_inc
    return datelist


def get_domain_limits(inargs):
    """
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------
    l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad : int
      Start and end indices for model and radar domain
    """
    ie = get_config(inargs, 'domain', 'ie')
    je = get_config(inargs, 'domain', 'je')
    ana_irange = get_config(inargs, 'domain', 'ana_irange')
    ana_jrange = get_config(inargs, 'domain', 'ana_jrange')

    # Calculate analysis limits
    assert (ie % 2 == 1) and (je % 2 == 1), 'Function assumes odd ie and je.'
    l11 = (ie - ana_irange - 1) / 2
    l12 = -(l11 + 1)
    l21 = (je - ana_jrange - 1) / 2
    l22 = -(l21 + 1)

    # Get radar limits. Hard coded as of now
    l11_rad = get_config(inargs, 'domain', 'radar_istart')
    l12_rad = get_config(inargs, 'domain', 'radar_istop')
    l21_rad = get_config(inargs, 'domain', 'radar_jstart')
    l22_rad = get_config(inargs, 'domain', 'radar_jstop')

    return l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad


def get_radar_mask(inargs):
    """
    Returns a radar_mask for all analysis days
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    Returns
    -------
    radar_mask :  2D numpy array 
      Mask where True values are masked (invalid)
    """
    # Check if radar_mask exists for days
    radar_mask_fn = ('../aux_files/radar_tot_mask_' + inargs.date_start +
                     '_' + inargs.date_end + '.npy')
    if (os.path.isfile(radar_mask_fn)) and (inargs.recompute is False):
        print('Found radar mask: ' + radar_mask_fn)
        return np.load(radar_mask_fn)
    else:
        print('Compute radar mask: ' + radar_mask_fn)
        datalist = get_datalist_radar(inargs, 'all')
        mask = np.max(datalist, axis=0) > 100   # TODO This has to go in the paper!

        # Crop mask
        l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad = \
            get_domain_limits(inargs)
        mask = mask[l11_rad:l12_rad, l21_rad:l22_rad]

        # Save figure
        plt.imshow(mask)
        plt.colorbar()
        fig_fn = (get_config(inargs, 'paths', 'figures') +
                    'radar_tot_mask_' + inargs.date_start +
                    '_' + inargs.date_end + '.pdf')
        print('Save radar mask figure as ' + fig_fn)
        plt.savefig(fig_fn)

        np.save(radar_mask_fn, mask)
        return mask


def get_pp_fn(inargs):
    """
    Creates a filename for the pre-processed NetCDF file
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------
    pp_fn : str
      Filename with path of pre-processed NetCDF file

    """
    pp_fn = get_config(inargs, 'paths', 'preproc_data')
    for key, value in vars(inargs).items():
        if not key is 'recompute':
            pp_fn += key + '-' + str(value) + '_'
    pp_fn = pp_fn[:-1] + '.nc'  # remove last '_'
    return pp_fn


def pp_exists(inargs):
    """
    Check whether  preprocessed file exists.
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------
    pp_exists : bool
      True if pre-processed file exitst
    """
    return os.path.isfile(get_pp_fn(inargs))


def get_datalist_radar(inargs, date, radar_mask=False):
    """
    Get data time series for radar observation.
    Parameters
    ----------
    inargs : : argparse object
      Argparse object with all input arguments
    date : str
      Date in format yyyymmddhh. If 'all', radar data for all days are returned
    radar_mask : 2D numpy array
      Radar mask to create masked arrays.

    Returns
    -------
    datalist : list
      List of 2D masked arrays
    """
    # Get file name
    radarpref = (get_config(inargs, 'paths', 'radar_data') +
                 get_config(inargs, 'paths', 'radar_prefx'))
    radarsufx = get_config(inargs, 'paths', 'radar_sufix')
    dtradar = timedelta(minutes=10)
    if date is 'all':
        date_start = (yyyymmddhh_strtotime(inargs.date_start) +
                      timedelta(hours=inargs.time_start) - dtradar)
        date_end = (yyyymmddhh_strtotime(inargs.date_end) +
                    timedelta(hours=inargs.time_end)-dtradar)
    else:
        dateobj = yyyymmddhh_strtotime(date)
        date_start = dateobj + timedelta(hours=inargs.time_start) - dtradar
        date_end = dateobj + timedelta(hours=inargs.time_end) - dtradar
    datalist = getfobj_ncdf_timeseries(radarpref,
                                       date_start,
                                       date_end,
                                       timedelta(hours=inargs.time_inc),
                                       reftime=date_start,
                                       ncdffn_sufx=radarsufx,
                                       fieldn='pr',
                                       abs_datestr='yymmddhhmm',
                                       dwdradar=True,
                                       return_arrays=True)
    # Crop data
    l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad = \
        get_domain_limits(inargs)
    if not radar_mask is False:
        # Apply total mask to data
        for i, data in enumerate(datalist):
            datalist[i] = masked_array(data[l11_rad:l12_rad,
                                            l21_rad:l22_rad],
                                       mask=radar_mask)
    return datalist


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