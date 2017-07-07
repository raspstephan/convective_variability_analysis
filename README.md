## Release Notes: poster_delft
This is the version of the repository used to make the plots for the poster at the "Future of Cumulus Parameterization" workshop in Delft, Netherlands in July 2017. The figures and log-files for creating the figures on the poster can be found at https://figshare.com/collections/convective_variability_poster_delft_jul_2017/3806629

The Jupyter notebook illustrating cloud separation and the radial distribution function is in the `jupyter_notebooks` directory. If you want to play around with it, go there; if you simply want to look at it you can do that here: http://nbviewer.jupyter.org/github/raspstephan/convective_variability_analysis/blob/master/jupyter_notebooks/cloud_identification_and_rdf.ipynb

# Convective Variability in COSMO ensembles

This is my working repository for my convective organization and variability research which I am doing for my PhD. One main goal of this work is to make the data analysis **reproducible**. Here I am loosely following the footsteps of [Damien Irving](https://github.com/DamienIrving).

The repository is structured as follows:

- `python scripts` contains the main data analysis and plotting scripts. All scripts are designed to be run from the command line. Typing `python <script_name>.py -h` gives a summary of what the script does and which arguments it takes. Each python script produces figures (figure directory is specified in the configuration file) along with a log file which describes the computing setup and evironment and the command with which the script was called. More details on this can be found in the README inside the directory.

- `config` contains the configuration file which is read by the python scripts. Here some basic settings are specified along with all the directories.

- `jupyter_notebooks` contains, as the name suggests, Jupyter notebooks. This can be little test snippets or larger notebooks.

- `cosmo_runscripts` contains the COSMO runscripts with which the raw data were produced. Again, the README inside the directory gives more information.

- `aux_files` contains some auxiliary files, such as the radar masks.

- `synop_plots` contains scripts to plot a synoptic overview from ECMWF data


