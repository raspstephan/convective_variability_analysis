"""
Filename:     master.py
Author:       Stephan Rasp, s.rasp@lmu.de

This file is the top-level file for analyzing COSMO output data for my 
convective variability project.

Much of the structure and many of the ideas are directly taken from 
https://github.com/DamienIrving

"""

# Import modules
import os
import argparse
from preprocessing import preprocess
from subprocess import check_output
from git import Repo
from datetime import datetime
from helpers import pp_exists, get_pp_fn


# Define functions
def create_log_str(inargs):
    """
    Function to create a log file tracking all steps from initial call to
    figure.
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    """
    time_stamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    conda_info = check_output(['conda', 'info'])
    conda_list = check_output(['conda', 'list'])
    # TODO: Get base dir automatically
    git_dir = '~/repositories/convective_variability_analysis'
    git_hash = Repo(git_dir).heads[0].commit
    pwd = check_output(['pwd'])
    script_name = os.path.basename(__file__)
    args_str = ''
    for arg in vars(inargs):
        args_str += ('--' + arg + ' ' + str(getattr(inargs, arg)) + ' ')

    log_str = ("""
    Preprocessing log\n
    -----------------\n
    %s\n
    %s\n
    %s\n
    Git hash: %s\n
    In directory: %s\n
    %s %s\n
    """ % (time_stamp, conda_info, conda_list, str(git_hash)[0:7], pwd,
           script_name, args_str))
    return log_str


def main(inargs):
    """
    Runs the main program
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    """
    log_str = create_log_str(inargs)

    # Check if pre-processed file exists
    if (pp_exists(inargs) is False) or (inargs.recompute is True):
        # Call preprocessing routine with arguments
        preprocess(inargs, log_str)
    else:
        print('Found pre-processed file:' + get_pp_fn(inargs))


if __name__ == '__main__':
    
    extra_info = 'Top-level script to analyze data'
    description = 'Bla bla'
    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info)
    
    parser.add_argument('--date_start',
                        type=str,
                        help='Start date of analysis in yyyymmddhh')
    parser.add_argument('--date_end',
                        type=str,
                        help='End date of analysis in yyyymmddhh')
    parser.add_argument('--time_start',
                        type=int,
                        default=1,
                        help='Analysis start time in hrs [including]. \
                              Default = 1')
    parser.add_argument('--time_end',
                        type=int,
                        default=24,
                        help='Analysis end time in hrs [including]. \
                              Default = 24')
    parser.add_argument('--time_inc',
                        type=float,
                        default=1,
                        help='Analysis increment in hrs. Default = 1')
    parser.add_argument('--nens',
                        type=int,
                        help='Number of ensemble members')
    parser.add_argument('--config_file',
                        type=str,
                        default='config.yml',
                        help='Config file in relative directory ../config. \
                              Default = config.yml')
    parser.add_argument('--recompute',
                        dest='recompute',
                        action='store_true',
                        help='If True, recompute pre-processed file.')
    parser.set_defaults(recompute=False)

    args = parser.parse_args()

    main(args)
