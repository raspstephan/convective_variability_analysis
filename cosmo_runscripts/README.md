# Description of INT2LM and COSMO runs
**Disclaimer:** Because Deutscher Wetterdienst (DWD) does not provide their models and data as open source, we are not able to provide any actual data or model code. With this limitation in mind, we nonetheless will try to provide a complete description of how our final data was produced. The hope is that this makes it easy for ther researchers with access to DWD data to set up similar experiments. Furthermore, it should provide internal (LMU) users with detailed instructions on the data structure and model setup. for this reason, all the specific paths are provided.


## Initial and Boundary data
The initial and boundary data for these experiments are COSMO-EU operational analyses. A description of the COSMO-EU setup can be found here, but unfortunately only in German : https://www.dwd.de/SharedDocs/downloads/DE/modelldokumentationen/nwv/cosmo_eu/cosmo_eu_dbbeschr_201406.pdf?__blob=publicationFile&v=3

The COSMO-EU model runs with a horizontal grid spacing of 7 km and parameterized convection. It uses a nudging data assimulation system.

The data are stored in the SKY database at DWD and are only accessible for registered users with a DWD account. For these experiments Dr. Christian Keil (LMU) provided the data using the following database command ( **but this is not really the one for my data, right???** ):

'read db=roma d=2016052700 cat=c3e_main_fc_rout p=CAPE_ML,TOT_PREC enum=1/to/ s[h]=0/to/24/by/1 f=/e/gtmp/dlrkeil/2016052700.grib2'

The data are stored as GRIB2 files here:

'/project/meteo/w2w/Unwetter2016/icbc'

## Interpolation to COSMO-DE grid
To use the COSMO-EU analyses as initial and boundary conditions for our COSMO-DE runs, they first had to be interpolated using the interpolation program **int2lm**. The version used for this project was version 2.02 cloned from our Waves2Weather internal repository (git commit b3894777). The actual executable used is stored here: 

'/project/meteo/w2w/A6/S.Rasp/archive/convective_variability/exe/tstint2lm'

The **namelists** for the int2lm runs are stored in this current directory. The path structure is as follows:

'< Date in format yyyymmddhh >/eu2de/run_int2lm_ceu_de'

The namelists are configured to write the interpolated data to 

'/home/scratch/users/stephan.rasp/< Date in format yyyymmddhh >/dein_ceu/'

An external file containing constant fields (e.g. land usage) is required for the interpolation. The file used is saved here ( **Where did I get that from?? ** ): 

'/project/meteo/w2w/A6/S.Rasp/archive/convective_variability/external/lm_d0_02800_1605x1605.mol.g1'

## COSMO-DE runs
With the interpolated data we can run the actual COSMO-DE simulations which are then used for our final analysis.

