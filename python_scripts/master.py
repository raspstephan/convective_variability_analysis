"""
Filename:     master.py
Author:       Stephan Rasp, s.rasp@lmu.de

This file is the top-level file for analyzing COSMO output data for my 
convective variability project.

Much of the structure and many of the ideas are directly taken from 
https://github.com/DamienIrving

"""

# Import modules
import argparse
from preprocessing import preprocess
from plotting import plotting
from helpers import pp_exists, get_pp_fn, create_log_str


# Define functions
def main(inargs):
    """
    Runs the main program
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    """

    # Check arguments
    assert inargs.plot in ['weather_ts', 'prec_stamps', 'prec_hist'], \
        'Plot not supported.'

    # Check if pre-processed file exists
    if (pp_exists(inargs) is False) or (inargs.recompute is True):
        print('Compute preprocessed file: ' + get_pp_fn(inargs))
        # Call preprocessing routine with arguments
        preprocess(inargs)
    else:
        print('Found pre-processed file:' + get_pp_fn(inargs))

    # Call analyzing and plotting routine
    plotting(inargs)


if __name__ == '__main__':
    
    extra_info = 'Top-level script to analyze data'
    description = __doc__
    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info)

    parser.add_argument('--plot',
                        type=str,
                        help='Which plot? [weather_ts, prec_stamps, prec_hist]')
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
    parser.add_argument('--zoom_lon1',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lon1')
    parser.add_argument('--zoom_lon2',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lon2')
    parser.add_argument('--zoom_lat1',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lat1')
    parser.add_argument('--zoom_lat2',
                        type=int,
                        default=0,
                        help='Zoom index for prec_stamps. Lat2')
    parser.add_argument('--recompute',
                        dest='recompute',
                        action='store_true',
                        help='If True, recompute pre-processed file.')
    parser.set_defaults(recompute=False)

    args = parser.parse_args()

    main(args)
