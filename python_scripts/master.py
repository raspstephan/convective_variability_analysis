"""
Filename:     master.py
Author:       Stephan Rasp, s.rasp@lmu.de

This file is the top-level file for analyzing COSMO output data for my 
convective variability project.

Much of the structure and many of the ideas are directly taken from 
https://github.com/DamienIrving

"""

# Import modules

import sys, os
import yaml
import argparse
import datetime
from git import Repo
from preprocessing import preprocess


# Define functions
def read_config_file(config_file):
    """
    Reads the config JSON file
    
    Parameters
    ----------
    config_file : str 
      Path to config file
      
    Returns
    -------
    config : dict
      Dictionary with all config settings
    """
    config = yaml.safe_load(open('../config/' + config_file))
    return config


def get_pp_fn(config, inargs):
    """
    Creates a filename for the pre-processed NetCDF file
    Parameters
    ----------
    config : dict
      Configuration dictionary
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------
    pp_fn : str
      Filename with path of pre-processed NetCDF file

    """
    pp_fn = config['paths']['preproc_data']
    for key, value in vars(inargs).items():
        pp_fn += key + '-' + str(value) + '_'
    pp_fn = pp_fn[:-1]   # remove last '_'
    print('Pre-processed file: ' + pp_fn)
    return pp_fn



def create_log_str(inargs):
    """
    Function to create a log file tracking all steps from initial call to
    figure.
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    """
    time_stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    exe = sys.executable
    py_vers = sys.version
    args_str = '  '.join(vars(inargs))
    # TODO: Get base dir automatically
    git_hash =  Repo('~/repositories/convective_variability_analysis').heads[0].commit
    
    # TODO: How do I dynamically get the filename?
    log_fn = '/home/s/S.Rasp/repositories/convective_variability_analysis/log_files/test_log.log'
    # TODO: Assert log file with same name does not already exist
    log_file = open(log_fn, 'w+')
    log_file.write("""%s: %s %s %s (Git hash: %s)""" 
                   %(time_stamp, exe, args, py_vers, str(git_hash)[0:7]))
    log_file.close()
    

def main(inargs):
    """
    Runs the main program
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    """

    # Some preliminary setup
    config = read_config_file(inargs.config_file)
    create_log_str(inargs)
    pp_fn = get_pp_fn(config, inargs)


    # Step 2: Call preprocessing routine with arguments
    preprocess(inargs, pp_fn)


if __name__ == '__main__':
    
    extra_info = 'Top-level script to analyze data'
    description = 'Bla bla'
    
    parser = argparse.ArgumentParser(description = description,
                                     epilog = extra_info)
    
    parser.add_argument('--date_start',
                        type = str,
                        help = 'Start date of analysis in yyyymmdd')
    parser.add_argument('--date_end',
                        type = str,
                        help = 'End date of analysis in yyyymmdd')
    parser.add_argument('--nens',
                        type = int,
                        help = 'Number of ensemble members')
    parser.add_argument('--config_file',
                        type = str,
                        default = 'config.yml',
                        help = 'Config file')

    args = parser.parse_args()

    main(args)






