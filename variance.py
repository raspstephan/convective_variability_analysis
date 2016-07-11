"""
This is a script to analyze the convetive variance from 
COSMO ensembles
"""
# Imports
import sys
from cosmo_utils.pyncdf import getfobj_ncdf_ens
from cosmo_utils.diag import mean_spread_fieldobjlist
from cosmo_utils.plot import fig_contourf_1sp
from cosmo_utils.helpers import yyyymmddhh_strtotime, make_timelist, ddhhmmss, yyyymmddhh
from datetime import timedelta
import os
import numpy as np
import matplotlib.pyplot as plt

# Setup
ana = sys.argv[1]  # 'm' or 'p'
date = sys.argv[2]
ensdir = '/home/scratch/users/stephan.rasp/' + date + '/deout_onlypsp/'
nens = 2
nens = [1,4]
tstart = timedelta(hours=2)
tend = timedelta(hours = 4)
tinc = timedelta(hours = 1)
hl = 25   # Number of grid pts to exclude at border
plotdir = '/home/s/S.Rasp/Dropbox/figures/PhD/variance/' + date

# Specific setup for type of analysis
if ana == 'm':
    fieldn = 'W'
    thresh = 1.
    sufx = 'z.nc_1h'
    plotdir += '/m/'

if ana == 'p':
    fieldn = 'inst. prec'
    thresh = 0.001
    sufx = '.nc_5m'
    plotdir += '/p/'
    
# Create plotdir if not exist
if not os.path.exists(plotdir):
    os.makedirs(plotdir)


# Make the timelist
timelist = make_timelist(tstart, tend, tinc)
# Time loop 
for t in timelist:
    
    
    ncdffn = 'lfff' + ddhhmmss(t) + sufx
    
    # Load ensembles 
    fobjlist = getfobj_ncdf_ens(ensdir, 'sub', nens, ncdffn, 
                                dir_suffix='/OUTPUT/',
                                fieldn = fieldn, nfill=2)
    if ana == 'm':   # Additional positive QC filter
        pass
    
    
    # Calculate what needs to be calculated
    # 1. Tau_c
    taucfn = 'lfff' + ddhhmmss(t) + '.nc_5m'
    tauclist = getfobj_ncdf_ens(ensdir, 'sub', nens, taucfn, 
                                dir_suffix='/OUTPUT/',
                                fieldn = 'TAU_C', nfill=2)
    # Get mean tauc object
    meantauc, tmp = mean_spread_fieldobjlist(tauclist)
    dimeantauc = np.nanmean(meantauc.data[hl:-hl, hl:-hl])
    # 2. Cloud size and m/p distribution
    
    
    
    
    # Plot what needs to be plotted every time step
    # 1. Ensemble mean tau_c over Germany
    title_sufx = 'mean tau_c, di_mean: ' + str(dimeantauc) + 'h + ' + ddhhmmss(t)
    fig = fig_contourf_1sp(meantauc, pllevels = np.arange(0, 21, 1),
                           extend = 'max', sp_title = title_sufx)
    plt.tight_layout()
    fig.savefig(plotdir + 'test_tauc_' + ddhhmmss(t))
    # 2. Plot cloud size and m/p distribution
    
    
    
    
    
# Plot summary plots
