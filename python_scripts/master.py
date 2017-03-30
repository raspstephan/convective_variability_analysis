"""
Filename:     master.py
Author:       Stephan Rasp, s.rasp@lmu.de

This file is the top-level file for analyzing COSMO output data for my 
convective variability project.

"""

# Import modules

import yaml
import argparse


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
    config = yaml.safe_load(open(config_file))
    return config
    
def main(inargs):
    """
    Runs the main program
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    """
    config = read_config_file(inargs.config_file)


if __name__ == '__main__':
    
    extra_info = 'Top-level script to analyze data'
    description = 'Bla bla'
    
    parser = argparse.ArgumentParser(description = description,
                                     epilog = extra_info)
    
    parser.add_argument('--config_file', 
                        type = str, 
                        default = '../config/config.yml',
                        help = 'Config file')
    
    args = parser.parse_args()
    
    main(args)






