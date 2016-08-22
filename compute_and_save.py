"""
This scripts reads the original COSMO data, computes statistics and 
saves the results in a netCDF file.
"""

# Imports
import argparse
from netCDF4 import Dataset
import numpy as np
from datetime import timedelta
from cosmo_utils.pyncdf import getfobj_ncdf_ens, getfobj_ncdf
from cosmo_utils.helpers import make_timelist, ddhhmmss

# Setup
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--ana', metavar = 'ana', type=str)
parser.add_argument('--date', metavar = 'date', type=str)
parser.add_argument('--water', metavar = 'water', type=bool, default = True)
parser.add_argument('--collapse', metavar = 'collapse', type=bool, 
                    default = True)
parser.add_argument('--nens', metavar = 'nens', type=int, default = 20)
parser.add_argument('--tstart', metavar = 'tstart', type=int, default = 1)
parser.add_argument('--tend', metavar = 'tend', type=int, default = 24)
parser.add_argument('--tinc', metavar = 'tinc', type=int, default = 60)
args = parser.parse_args()

# Convert times to timedelta objects
tstart = timedelta(hours= args.tstart)   # Cannot be 0 because of tau_c calculation!
tend = timedelta(hours = args.tend)  
tinc = timedelta(minutes = args.tinc)  # temporal resolution for analysis

# Make lists for loops, dimensions
timelist = make_timelist(tstart, tend, tinc)
nlist = [256, 128, 64, 32, 16, 8, 4]
#hlist = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000,
            #8000, 9000, 10000]
hlist = [3000]

ensdir = '/home/scratch/users/stephan.rasp/' + args.date + '/deout_ceu_pspens/'

# Analysis-specific setup
if args.ana == 'm':
    HH = getfobj_ncdf(ensdir + '/1/OUTPUT/lfff00000000c.nc_30m', 'HHL')
    aslcol = HH.data[:, -1, -1]   # ATTENTION Hard coded colums above sea level 
    levlist = []
    realhlist = []
    for h in hlist:
        levlist.append(np.argmin(np.abs(aslcol-h)))   # find closest level
        realhlist.append(aslcol[np.argmin(np.abs(aslcol-h))])
    sufx = '.nc_30m'
elif args.ana == 'p':
    levlist = [None]
    realhlist = ['surf']
else:
    raise Exception, 'wrong analysis'

################################################################################
# Allocate NetCDF file
rootgrp = Dataset('./tmp', 'w', format='NETCDF4')
# Create dimensions
tdim = rootgrp.createDimension('time', len(timelist))
ndim = rootgrp.createDimension('n', len(nlist))
levdim = rootgrp.createDimension('levs', len(levlist))
xdim = rootgrp.createDimension('x', nlist[0])
ydim = rootgrp.createDimension('y', nlist[0])
xndim = rootgrp.createDimension('x_n', None)
yndim = rootgrp.createDimension('y_n', None)
nclddim = rootgrp.createDimension('N_cld', None)
drdim = rootgrp.createDimension('dr', None)   #TODO: size of dim

# Create variables
time     = rootgrp.createVariable('time', 'f8', ('time',))
n        = rootgrp.createVariable('n', 'i4', ('n',))
levs     = rootgrp.createVariable('levs', 'i4', ('levs',))
cld_size = rootgrp.createVariable('cld_size', 'f8', ('time','levs','N_cld'))
cld_sum  = rootgrp.createVariable('cld_sum', 'f8', ('time','levs','N_cld'))
rdf      = rootgrp.createVariable('rdf', 'f8', ('time','levs','dr'))
varM     = rootgrp.createVariable('varM', 'f8', ('time','levs','n','x_n','y_n'))
varN     = rootgrp.createVariable('varN', 'f8', ('time','levs','n','x_n','y_n'))
N        = rootgrp.createVariable('N', 'f8', ('time','levs','n','x_n','y_n'))
M        = rootgrp.createVariable('M', 'f8', ('time','levs','n','x_n','y_n'))
m        = rootgrp.createVariable('m', 'f8', ('time','levs','n','x_n','y_n'))
# End allocation
################################################################################


###################
## Time loop      #
###################
for it, t in enumerate(timelist):
    print 'time: ', t
    
    ############################################################################
    # Load COSMO data
    ncdffn = 'lfff' + ddhhmmss(t) + sufx
    fieldlist = getfobj_ncdf_ens(ensdir, 'sub', nens, ncdffn, 
                                 dir_suffix='/OUTPUT/', fieldn = fieldn, 
                                 nfill=1, levs = levlist, return_arrays = True)
    if ana == 'm':   # Load some extra variables 
        qclist = getfobj_ncdf_ens(ensdir, 'sub', nens, ncdffn, 
                                  dir_suffix='/OUTPUT/', fieldn = 'QC', 
                                  nfill=1, levs = levlist, return_arrays = True)
        ncdffn_rho = ncdffn + '_buoy'
        rholist = getfobj_ncdf_ens(ensdir, 'sub', nens, ncdffn_rho, 
                                   dir_suffix='/OUTPUT/', fieldn = 'RHO', 
                                   nfill=1, levs = levlist, return_arrays = True)
    else:   # Fill lists with None
        qclist = [None]*len(fieldlist)
        rholist = [None]*len(fieldlist)
    # Crop all fields to analysis domain
    sxo, syo = fobjlist[0].data.shape  # Original field shape
    lx1 = (sxo-256-1)/2 # ATTENTION first dimension is actually y
    lx2 = -(lx1+1) # Number of grid pts to exclude at border
    ly1 = (syo-256-1)/2
    ly2 = -(ly1+1)
    for field in fieldlist + qclist + rholist:
        if not field == None:
            field = field[:,lx1:lx2, ly1:ly2]
    # End loading data
    ############################################################################

    ####################
    ## lev loop        #
    ####################
    for iz, lev in enumerate(levlist):
        
        ########################################################################
        # Calculate cloud statistics
        
        # Member loop
        sizelist = []
        sumlist = []
        rdflist = []
        for field, qc, rho in zip(fieldlist, qclist, rholist):
            # Identify clouds
            if ana == 'm':
                tmp = identify_clouds(field[iz], thresh, qc[iz],
                                      opt_thresh = 0., water = water,
                                      rho = rho[iz])
                labels, cld_size_mem, cld_sum_mem = tmp
                cld_sum_mem *= dx*dx  # Rho is now already included
            else:
                tmp = identify_clouds(field[iz], thresh)
                labels, cld_size_mem, cld_sum_mem = tmp
            sizelist.append(cld_size_mem)
            sumlist.append(cld_sum_mem)
            
            # Calculate RDF
            g, r = rdf(labels, field[iz], normalize = True, rmax = 30, dr = 2)
            rdflist.append(g)
        
        # Save lists and mean rdf
        cld_size[it, iz, :] = np.array(sizelist)
        cld_sum[it, iz, :] = np.array(sumlist)
        rdf[it, iz, :] = np.mean(glist, axis = 0)
        
        ################
        ## n loop      #
        ################
        #for n in nlist:
        
        
        
        
        
        
        
