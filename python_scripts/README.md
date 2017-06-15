# python_scripts - Data processing and plotting

This directory contains all the Python scripts to process the raw data and produce figures.

The scripts are all designed (with the exception of `helpers.py` which contains additional functions for the other scripts) to be called from the command line with arguments. Type `-h` to show the arguments.

Note that I am using an Anaconda environent. The details are specified in the log file.

The table below lists the files plots and options. Every script processes the raw data in a certain way and creates different plots. Only one plot is produced at one time, depending on the `plot_type` parameter.

| File | Plots |
|------|-------|
| `cloud_stats.py` | Precipitation frequency histogram |
|                 | Cloud size and cloud precipitation sum distributions |
|                 | Radial distribution function |
| `variability.py` | Std_vs_mean |
|                 | Cloud size and cloud precipitation sum distributions |
|                 | Precipitation radial distribution function |


