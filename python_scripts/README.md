# README for python_scripts

The main script,`master.py`,calls `preprocessing.py` which reads the model output data, does some data manipulation and saves the manipulated data as a netCDF file. Then `plotting.py` are called to do some final data manipulation and plot the figures.

All scripts call function from `helpers.py`.

`master.py` has the following arguments:

Argument | Description
--- | ---
`--date_start` | Start date of analysis in yyyymmddhh
`--date_end`    | End date of analysis in yyyymmddhh
`--time_start` | Analysis start time in hrs [including]. Default = 1
`--time_end`  |  Analysis end time in hrs [including]. Default = 24
`--time_inc`  | Analysis increment in hrs. Default = 1
`--nens`     |       Number of ensemble members
`--config_file` | Config file in relative directory ../config. Default = config.yml
`--recompute`      |     If True, recompute pre-processed file.

The `config_file` gives additional arguements which do not quantitatively change the results, such as paths, domain size specifics and plotting features.

## Log files
For each plot computed a log_file is saved which gives the exact command executed and all version information.
