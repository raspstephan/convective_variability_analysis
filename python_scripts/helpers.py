"""
Filename:     helpers.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains helper functions.

"""


# import modules
import yaml
import numpy as np
from datetime import datetime, timedelta


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


def make_datelist_yyyymmddhh(inargs):
    """
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

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

    datelist_yyyymmddhh = []
    date = dateob_start
    while date <= dateob_end:
        date_str = (str(date.year) + str(date.month).zfill(2) +
                    str(date.day).zfill(2) + str(date.hour).zfill(2))
        datelist_yyyymmddhh.append(date_str)
        date += date_inc
    return datelist_yyyymmddhh


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


def get_radar_mask():
    """
    Returns
    -------
    radar_mask :  2D numpy array 
      Mask where True values are masked (invalid)
    """

    return np.load('../aux_files/radar_tot_mask.npy')





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
        pp_fn += key + '-' + str(value) + '_'
    pp_fn = pp_fn[:-1] + '.nc'  # remove last '_'
    print('Pre-processed file: ' + pp_fn)
    return pp_fn