"""
Filename:     helpers.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains helper functions.

"""


# import modules
import yaml


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