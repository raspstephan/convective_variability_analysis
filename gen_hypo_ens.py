"""
This script create a hypothetical ensemble which can be analyzed with the 
compute script.

"""

# Imports
import os
import argparse
import numpy as np
from datetime import timedelta
from netCDF4 import Dataset
from cosmo_utils.helpers import make_timelist, ddhhmmss
from SubcriticalModel_Stephan import CloudPercolation_cw_i

# Setup
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--nens', metavar = 'nens', type=int, default = 20)
parser.add_argument('--type', metavar = 'type', type=str)
args = parser.parse_args()

# Convert times to timedelta objects
tstart = timedelta(hours = 0)   # Cannot be 0 because of tau_c calculation!
tend = timedelta(hours = 2)  
tinc = timedelta(minutes = 60)  # temporal resolution for analysis
timelist = make_timelist(tstart, tend, tinc)

sx = 357
sy = 357
dalpha = 0.025

#simulation set-up

N_ensemble = 2                                                                 #Number of calculated fields

#domain set-up 

phys_L_km=357*2.8                                                               #Domain size in km 
DL=2.8                                                                          #Resolution in km 
L=int(phys_L_km/DL)                                                             #Number of grid-cells in one direction

r_co_km=1.4                                                                     #Minimal radius of discs in km

#model properties

#coverage_fraction=0.05                                                          
N_clouds= 585                                                                   #Number of clouds

#Cloud properties
r_m = 2.5                                                                     #Mean radius of discs in km
cs  = 30                                                                       #Factor by which probability is increased within the rings around the clouds
fac_cw = 3                                                                     #Prob. increased within the ring r_disc < r < r_disc*fac_cw
 

ensdir = '/home/scratch/users/stephan.rasp/hypo_' + args.type + '/deout_ceu_pspens/'

# Loop over time
for it, t in enumerate(timelist):
    print 'time: ', t
    ncdffn = 'lfff' + ddhhmmss(t)
    # Loop over ensemble members
    for imem in range(1,args.nens+1):
        print 'Member:', imem
        ensstr = str(imem) + '/OUTPUT/'
        if not os.path.exists(ensdir + ensstr): os.makedirs(ensdir + ensstr)
        ensstr += ncdffn
        # Allocate NetCDF file
        rootgrp = Dataset(ensdir + ensstr + '.nc', 'w', 
                        format='NETCDF4')
        # create dimensions
        timedim = rootgrp.createDimension('time', 1)
        leveldim = rootgrp.createDimension('level', 1)
        rlondim = rootgrp.createDimension('rlon', sx)
        rlatdim = rootgrp.createDimension('rlat', sy)
        
        
        # Create variables
        m = rootgrp.createVariable('m', 'f8', ('time', 'level', 'rlat', 'rlon'))
        m.long_name = 'hypo m'
        rp = rootgrp.createVariable('rotated_pole', 'S1', ())
        rp.grid_north_pole_latitude = 0
        rp.grid_north_pole_longitude = 0
        rlon = rootgrp.createVariable('rlon', 'f8', ('rlon'))
        rlon[:] = np.linspace(0,1,sx)
        rlat = rootgrp.createVariable('rlat', 'f8', ('rlat'))
        rlat[:] = np.linspace(0,1,sy)
        
        cloud_field, cloud_centers, rA_km =CloudPercolation_cw_i(phys_L_km,DL,N_clouds,r_co_km,r_m,fac_cw,cs,plot_field=True) 
        print cloud_field.shape
        # Create field
        m[0,0,:,:] = cloud_field
        
        ## 1. Draw N from Poisson distribution
        #N = np.random.poisson(meanN, 1)
        ## 2. randomly pick location
        #randx = np.random.randint(0, sx, N)
        #randy = np.random.randint(0, sy, N)
        ## 3. draw m from exponential distribution
        #randm = np.random.exponential(meanm, N)
        #for rx, ry, rm in zip(randx, randy, randm):
            #m[0,0,rx,ry] = rm
            
        
        rootgrp.close()
