# python_scripts - Data processing and plotting

This directory contains all the Python scripts to process the raw data and produce figures, with the exception of the ECMWF synoptic plots with are the directory `synop_plots`.

The scripts are all designed (with the exception of `helpers.py` which contains additional functions for the other scripts) to be called from the command line with arguments. Type `-h` to show the arguments.

Note that I am using an Anaconda environent. The details are specified in the log file.

The table below lists the files along with the plotting options. Every script processes the raw data in a certain way and creates different plots. Only one plot is produced at one time, depending on the `plot_type` parameter.


| File | `plot_type` | Description |
|------|-------|----|
| `cloud_stats.py` | `freq_hist `|Precipitation frequency histogram |
|                 | `size_hist` | Cloud size and cloud precipitation sum distributions |
|                 | `rdf_individual` | Radial distribution function plotted for each day individually |
|                 | `rdf_composite` | Radial distribution function as a composite over all days |
|                 | `m_evolution` | Timeseries composite of mean m |
| `variability.py` | `r_v` or variant | Plots the requested metric as a function of scale and time of day |
|                 | `std_vs_mean` | Standard deviation against mean on log-log plot |
|                 | `correlation` | Plot correlation between m and N |
| `weather_time_series.py` | `prec_ind` | Plot domain mean timeseries of precipitation for each day individually |
|                 | `prec_comp` | Plot domain mean timeseries of precipitation as a composite over all days |
|                 | `cape_tauc_ind` | Plot domain mean CAPE and tau_c for each day individually |
|                 | `cape_tauc_comp` | Plot domain mean CAPE and tau_c as a composite over all days |
|                 | `prec_cape_comp` | Plot domain mean precipitation  and CAPE as a composite over all days |
| `plot_stamps.py` | `stamps` | Stamps comparing precipitation in observations, deterministic and select ensemble members |
|                 | `individual` | Single plot of one ensemble member. Either precipitation or separated cloud objects. |
| `spectra.py` |  | Plot spectra of precipitation and kinetic energy |


