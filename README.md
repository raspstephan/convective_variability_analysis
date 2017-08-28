## Release Notes: submission
This is the version of the code repository at the time the paper was submitted. The Jupyter notebooks mentioned in the paper can be downloaded or viewed here:
- http://nbviewer.jupyter.org/github/raspstephan/convective_variability_analysis/blob/master/jupyter_notebooks/cloud_identification_and_rdf.ipynb
- http://nbviewer.jupyter.org/github/raspstephan/convective_variability_analysis/blob/master/jupyter_notebooks/beta_sample_size_dependency.ipynb


# Convective Variability in COSMO ensembles

This is the working repository for my convective variability research. One main goal of this work is to make the data analysis reproducible. I am loosely following the footsteps of [Damien Irving](https://github.com/DamienIrving).

The repository is structured as follows:

- `python scripts` contains the main data analysis and plotting scripts. All scripts are designed to be run from the command line. Typing `python <script_name>.py -h` gives a summary of what the script does and which arguments it takes. Each python script produces figures (figure directory is specified in the configuration file) along with a log file which describes the computing setup and evironment and the command with which the script was called. More details on this can be found in the README inside the directory.

- `config` contains the configuration file which is read by the python scripts. Here some basic settings are specified along with all the directories.

- `jupyter_notebooks` contains, as the name suggests, Jupyter notebooks. This can be little test snippets or larger notebooks.

- `cosmo_runscripts` contains the COSMO runscripts with which the raw data were produced. Again, the README inside the directory gives more information.

- `synop_plots` contains scripts to plot a synoptic overview from ECMWF data.

- `aux_files` contains some auxiliary files, such as the radar masks. *Not used anymore*

- `fortran_scripts` *Not used anymore*




