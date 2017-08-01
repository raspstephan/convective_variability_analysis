"""
Filename:     spectra.py
Author:       Stephan Rasp, s.rasp@lmu.de
Description:  Plot precipitation spectra

"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import timedelta
from helpers import save_fig_and_log, get_config, make_datelist

################################################################################
# PLOTTING FUNCTIONS
################################################################################
def plot_spectra(inargs):
    """
    For now this loads previously computed files.

    Returns
    -------

    """
    savesuf = '_ana-spectra_wat-True_height-3000_nens-50_tstart-3_tend-24_tinc-180_minmem-5_dr-2.nc'
    # Load data
    dke_spec_list = []
    bgke_spec_list = []
    dprec_spec_list = []
    bgprec_spec_list = []
    for d in make_datelist(inargs):
        dataset = Dataset(inargs.preproc_dir + d + savesuf, 'r')
        dke_spec_list.append(dataset.variables['dkespec'][:])
        bgke_spec_list.append(dataset.variables['bgkespec'][:])
        dprec_spec_list.append(dataset.variables['dprecspec'][:])
        bgprec_spec_list.append(dataset.variables['bgprecspec'][:])
    dke_spec = np.nanmean(dke_spec_list, axis = 0)
    bgke_spec = np.nanmean(bgke_spec_list, axis = 0)
    dprec_spec = np.nanmean(dprec_spec_list, axis = 0)
    bgprec_spec = np.nanmean(bgprec_spec_list, axis = 0)
    timelist = [timedelta(seconds=ts) for ts in dataset.variables['time']]
    timelist_plot = [(dt.total_seconds() / 3600) for dt in timelist]

    # Define colors
    cyc = [plt.cm.jet(i) for i in np.linspace(0, 1, len(timelist))]
    cyc = ("#E7A7FF", "#FF84DB", "#EF8974", "#AF9300", "#529324", "#008768",
           "#006C88", "#2D3184")
    speclam = dataset.variables['speclam'][:]

    # Set up figures
    for diff, bg, name in zip([dke_spec, dprec_spec],
                              [bgke_spec, bgprec_spec],
                              ['Kinetic energy', 'Precipitation']):
        pw = get_config(inargs, 'plotting', 'page_width')
        width_fraction = 3. / 9.
        ratio = 1.
        fig, ax = plt.subplots(1, 1, figsize=(pw * width_fraction,
                                              pw * width_fraction * ratio))

        ############# Time loop ##############
        for it, t in enumerate(timelist):
            print 'time: ', t
            # Get ratio
            ratio = diff[it] / bg[it] / 2.
            ax.plot(speclam / 1000., ratio, c=cyc[it],
                       label=str(int(timelist_plot[it])).zfill(2),
                       linewidth=1.5)

        ax.legend(loc=3, ncol=2, fontsize=8, title='Time [UTC]')
        ax.plot([5, 1000.], [1, 1], c='gray', alpha=0.5)
        ax.set_xlabel('Wavelength [km]')
        # ax.set_ylabel('Saturation ratio')
        # ax.set_title("Saturation of KE spectrum")
        ax.set_ylim(1e-2, 1.1)
        ax.set_xlim(5, 1000.)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_position(('outward', 3))
        ax.spines['left'].set_position(('outward', 3))
        plt.yticks(rotation=90)

        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.3, top=0.85)

        save_fig_and_log(fig, None, inargs, name[:4], tight=False)





################################################################################
# MAIN FUNCTION
################################################################################
def main(inargs):
    """
    Runs the main program

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    """

    # Plot spectra
    plot_spectra(inargs)

if __name__ == '__main__':

    description = __doc__

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--date_start',
                        type=str,
                        default='2016052800',
                        help='Start date of analysis in yyyymmddhh')
    parser.add_argument('--date_end',
                        type=str,
                        default='2016060800',
                        help='End date of analysis in yyyymmddhh')
    parser.add_argument('--preproc_dir',
                        type=str,
                        default='/project/meteo/scratch/S.Rasp/convective_variability_data/old_preproc_data/',
                        help='Directory where old preprocessed files are stored.')
    parser.add_argument('--config_file',
                        type=str,
                        default='config.yml',
                        help='Config file in relative directory ../config. \
                              Default = config.yml')

    parser.add_argument('--sub_dir',
                        type=str,
                        default='spectra',
                        help='Sub-directory for figures and pp_files')
    parser.add_argument('--plot_name',
                        type=str,
                        default='',
                        help='Custom plot name.')
    parser.add_argument('--plot_format',
                        type=str,
                        default='pdf',
                        help='Which format for figure file.')

    args = parser.parse_args()

    main(args)