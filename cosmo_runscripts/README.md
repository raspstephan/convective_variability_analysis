# Description of INT2LM and COSMO runs
**Disclaimer:** Because Deutscher Wetterdienst (DWD) does not provide their models and data as open source, we are not able to provide any actual data or model code. With this limitation in mind, we nonetheless will try to provide a complete description of how our final data was produced. The hope is that this makes it easy for ther researchers with access to DWD data to set up similar experiments. Furthermore, it should provide internal (LMU) users with detailed instructions on the data structure and model setup. For this reason, all the specific paths are provided.


## Initial and Boundary data
The initial and boundary data for these experiments are COSMO-EU operational analyses. A description of the COSMO-EU setup can be found here, but unfortunately only in German : https://www.dwd.de/SharedDocs/downloads/DE/modelldokumentationen/nwv/cosmo_eu/cosmo_eu_dbbeschr_201406.pdf?__blob=publicationFile&v=3

The COSMO-EU model runs with a horizontal grid spacing of 7 km and parameterized convection. It uses a nudging data assimulation system.

The data are stored in the SKY database at DWD and are only accessible for registered users with a DWD account. For these experiments Dr. Christian Keil (LMU) provided the data using the following database command:

`read db=roma d=2016052700 cat=c3e_main_fc_rout p=CAPE_ML,TOT_PREC enum=1/to/ s[h]=0/to/24/by/1 f=/e/gtmp/dlrkeil/2016052700.grib2`

On our local system, the data are stored as GRIB2 files here:

`/project/meteo/w2w/Unwetter2016/icbc`

## Interpolation to COSMO-DE grid
To use the COSMO-EU analyses as initial and boundary conditions for our COSMO-DE runs, they first had to be interpolated using the interpolation program **int2lm**. The official user guide can be found here: http://cosmo-model.org/content/model/documentation/core/int2lm_2.03.pdf 

The version used for this project was version 2.02 cloned from our Waves2Weather internal repository (git commit b3894777). The actual executable used is stored here: 

`/project/meteo/w2w/A6/S.Rasp/archive/convective_variability/exe/tstint2lm`

The **namelists** for the int2lm runs are stored in this current directory. The path structure is as follows:

`< Date in format yyyymmddhh >/eu2de/run_int2lm_ceu_de`

The namelists are configured to write the interpolated data to 

`/home/scratch/users/stephan.rasp/< Date in format yyyymmddhh >/dein_ceu/`

but the data was later copied to 

`/project/meteo/scratch/users/stephan.rasp/convective_variability_data/raw_data/< Date in format yyyymmddhh >/dein_ceu/`

An external file containing constant fields (e.g. land usage) is required for the interpolation. The file used is saved here: 

`/project/meteo/w2w/A6/S.Rasp/archive/convective_variability/external/lm_d0_02800_1605x1605.mol.g1`

## COSMO-DE runs
With the interpolated data we can run the actual COSMO-DE simulations which are then used for our final analysis. A user guide can be found here: http://cosmo-model.org/content/model/documentation/core/cosmo_userguide_5.04.pdf

The base version used here is 5.04a, but with the PSP-scheme added. The actual version can be found in the Waves2Weather internal repository here (git commit 213f63d8):

`https://gitlab.wavestoweather.de/S.Rasp/cosmo`

If you are interested in using the PSP-enhanced version of COSMO, please contact me. The actual executable is stored here and in a private figshape directory:

`/project/meteo/w2w/A6/S.Rasp/archive/convective_variability/exe/lmparbin`

The **namelists** are provided in the current directory using the path:

`< Date in format yyyymmddhh >/cde_ceu_pspens/`

Then there are three different scripts for each day. Note that these are specific to our slurm cluster:

1. `var_array_1-50` is the namelist for the ensemble simulations. 
2. `var_array_det` is the namelist for the deterministic reference run.
3. `var_array_uv_1-5` is the namelist, which reruns the first 5 ensemble members with full 3D output of the horizontal wind field for the computation of the kinetic energy spectra.

The output is saved in

`/home/scratch/users/stephan.rasp/< Date in format yyyymmddhh >/deout_ceu_pspens/< ensemble member or det >/

but was later copied to

`/project/meteo/scratch/users/stephan.rasp/convective_variability_data/raw_data/< Date in format yyyymmddhh >/deout_ceu_pspens/< ensemble member or det >/`

The actual output fields are in another sub-directory `OUTPUT`. The COSMO `lfff` are split into different groups:

1. `*.nc_30` contains `W, QC, QI, QS` every 30 minutes
2. `*.nc_30_buoy` contains `QV, QR, T, P, TTEND_DIAB, TTENS_MPHY` every 30 minutes
3. `*.nc_30_surf` contains `TOT_PREC, CAPE_MK, TOT_PR, HPBL` every 30 minutes
4. `*.nc_30_uv` contains `U, V` every 60 minutes (despite its name)

Additionally, for 2016060400 the output for the PSP scheme (`TTENSSTO, WTENSSTO, QVTENSSTO, RND_BLPERT`) were written in the `*.nc_30` file for ensemble member 1 using the namelist `var_array_varout`. This is for potential illustration purposes.
