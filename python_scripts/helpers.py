"""
Filename:     helpers.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains helper functions.

"""


# import modules
import yaml
from datetime import datetime, timedelta


# Define functions
def get_config(inargs, top_key, bottom_key):
    """
    Reads the config JSON file

    Parameters
    ----------
    config_file : str
      Path to config file

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

