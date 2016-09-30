"""
This scripts reads the original COSMO data, computes statistics and 
saves the results in a netCDF file.
"""

# Imports
import warnings
warnings.filterwarnings("ignore")   # ATTENTION To suppress future warning for None
import argparse
from netCDF4 import Dataset, date2num
import numpy as np
from datetime import timedelta
from cosmo_utils.pyncdf import getfobj_ncdf_ens, getfobj_ncdf, \
                               getfobj_ncdf_timeseries
from cosmo_utils.helpers import make_timelist, ddhhmmss, yymmddhhmm, \
                                yyyymmddhh_strtotime
from cosmo_utils.diag import identify_clouds, calc_rdf, crosscor, int_rad_2d,\
                             get_totmask,powspec_2d_hor,powspec_2d_hor_alter
from scipy.ndimage.measurements import center_of_mass
from scipy.signal import correlate
import os


# Setup - Input arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--ana', metavar = 'ana', type=str)
parser.add_argument('--date', metavar = 'date', type=str)
parser.add_argument('--height', metavar = 'height', type=float, default =3000)
parser.add_argument('--water', metavar = 'water', type=str, default = 'True')
parser.add_argument('--nens', metavar = 'nens', type=int, default = 20)
parser.add_argument('--tstart', metavar = 'tstart', type=int, default = 1)
parser.add_argument('--tend', metavar = 'tend', type=int, default = 24)
parser.add_argument('--tinc', metavar = 'tinc', type=int, default = 60)
parser.add_argument('--minmem', metavar = 'minmem', type=int, default = 5)
parser.add_argument('--dr', metavar = 'dr', type=int, default = 2)
args = parser.parse_args()
# Convert water to bool 
if args.water == 'True':
    args.water = True
elif args.water == 'False':
    args.water = False
else:
    raise Exception

# Create file str
savedir = '/home/scratch/users/stephan.rasp/results/'

heightstr = str(int(args.height))
savestr = (args.date + '_ana-' + args.ana + '_wat-' + str(args.water) + 
           '_height-' + heightstr +
           '_nens-' + str(args.nens) + '_tstart-' + str(args.tstart) + 
           '_tend-' + str(args.tend) + '_tinc-' + str(args.tinc) + 
           '_minmem-' + str(args.minmem) + '_dr-' + str(args.dr) + '.nc')
print savestr
# Convert times to timedelta objects
tstart = timedelta(hours = args.tstart)   # Cannot be 0 because of tau_c calculation!
tend = timedelta(hours = args.tend)  
tinc = timedelta(minutes = args.tinc)  # temporal resolution for analysis


# Setup - Create or define lists
timelist = make_timelist(tstart, tend, tinc)
nlist = [256, 128, 64, 32, 16, 8, 4]
histbinedges = [0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 1000]

# Setup - Directories for input data
ensdir = '/home/scratch/users/stephan.rasp/' + args.date + '/deout_ceu_pspens/'
radarpref = '/project/meteo/w2w/A6/radolan/netcdf_cosmo_de/raa01-rw_10000-'
radarsufx = '-dwd---bin.nc'

# Analysis-specific setup
dx = 2800.   # Grid Spacing ATTENTION: Hard coded and wrong as of now

# Get analysis level
if not args.ana == 'hypo':
    HH = getfobj_ncdf(ensdir + '/1/OUTPUT/lfff00000000c.nc_30m', 'HHL')
    aslcol = HH.data[:, -1, -1]   # ATTENTION Hard coded colums above sea level 
    lev = np.argmin(np.abs(aslcol-args.height))   # find closest level
    realh = aslcol[np.argmin(np.abs(aslcol-args.height))]
    print 'level', lev
    sxo, syo = HH.data.shape[1:]  # Original field shape
    sufx = '.nc_30m'
    fieldn = 'W'
    thresh = 1.
else:
    thresh = 0.
    sxo, syo = (357, 357)
    sufx = '.nc'

rmax_rdf = 36
dr_rdf = args.dr

# Determine analysis domain
lx1 = (sxo-256-1)/2 # ATTENTION first dimension is actually y
lx2 = -(lx1+1) # Number of grid pts to exclude at border
ly1 = (syo-256-1)/2
ly2 = -(ly1+1)

if args.ana == 'vert':
    levlist = range(10, 50, 2)
    heightlist = []
    for lev in levlist:
        heightlist.append(aslcol[lev])
    HHcropped = HH.data[-1,lx1:lx2, ly1:ly2]
    HH50tot = np.mean(HHcropped[:,:])
    HH50south = np.mean(HHcropped[:256/2, :])
    HH50north = np.mean(HHcropped[256/2:, :])
    print 'tot', HH50tot, 'south', HH50south, 'north', HH50north

#elif args.ana == 'hypo':
    #levlist = [0]
    #sufx = '.nc'
    #fieldn = 'm'
    #ensdir = ('/home/scratch/users/stephan.rasp/hypo_' + args.date + 
              #'/deout_ceu_pspens/')
    #thresh = 0.


################################################################################
os.system('rm ' + savedir + savestr)
# Allocate NetCDF file
rootgrp = Dataset(savedir + savestr, 'w', format='NETCDF4')
# Create dimensions
tdim = rootgrp.createDimension('time', len(timelist))
ndim = rootgrp.createDimension('n', len(nlist))
xdim = rootgrp.createDimension('x', nlist[0])
ydim = rootgrp.createDimension('y', nlist[0])
nclddim = rootgrp.createDimension('N_cld', 1e6)
drdim = rootgrp.createDimension('dr', rmax_rdf/dr_rdf+1) # For RDF
drcorrdim = rootgrp.createDimension('drcorr', nlist[0]) # For 2D ACF
binsdim = rootgrp.createDimension('bins', len(histbinedges)-1)
specdim = rootgrp.createDimension('spec', 128)

# Create variables and add attributes 
time     = rootgrp.createVariable('time', 'f8', ('time',))
time[:]  = [td.total_seconds() for td in timelist]
n        = rootgrp.createVariable('n', 'i4', ('n',))
n[:]     = nlist
dr       = rootgrp.createVariable('dr', 'f4', ('dr',))


# weather
if args.ana == 'weather':
    ditauc   = rootgrp.createVariable('ditauc', 'f8', ('time'))
    dicape   = rootgrp.createVariable('dicape', 'f8', ('time'))
    diprec   = rootgrp.createVariable('diprec', 'f8', ('time'))
    dihpbl   = rootgrp.createVariable('dihpbl', 'f8', ('time'))

# prec 
if args.ana == 'prec':
    rdf_prec_model = rootgrp.createVariable('rdf_prec_model', 'f8', ('time','dr'))
    rdf_prec_obs   = rootgrp.createVariable('rdf_prec_obs', 'f8', ('time','dr'))
    hist_model   = rootgrp.createVariable('hist_model', 'f8', ('bins'))
    hist_obs   = rootgrp.createVariable('hist_obs', 'f8', ('bins'))

# spectra
if args.ana == 'spectra':
    bgkespec = rootgrp.createVariable('bgkespec', 'f8', ('time','spec'))
    dkespec  = rootgrp.createVariable('dkespec', 'f8', ('time','spec'))
    bgprecspec = rootgrp.createVariable('bgprecspec', 'f8', ('time','spec'))
    dprecspec  = rootgrp.createVariable('dprecspec', 'f8', ('time','spec'))
    speck    = rootgrp.createVariable('speck', 'f8', ('spec'))
    speclam  = rootgrp.createVariable('speclam', 'f8', ('spec'))

# clouds
if args.ana == 'clouds':
    cld_size = rootgrp.createVariable('cld_size', 'f8', ('time','N_cld'))
    cld_sum  = rootgrp.createVariable('cld_sum', 'f8', ('time','N_cld'))
    rdf      = rootgrp.createVariable('rdf', 'f8', ('time','dr'))
    exw      = rootgrp.createVariable('exw', 'f8', ('time', 'x', 'y'))
    exq      = rootgrp.createVariable('exq', 'f8', ('time', 'x', 'y'))
    #exbin    = rootgrp.createVariable('exbin', 'f8', ('time', 'levs', 'x', 'y'))
    excld    = rootgrp.createVariable('excld', 'f8', ('time', 'x', 'y'))
    exwater  = rootgrp.createVariable('exwater', 'f8', ('time', 'x', 'y'))
    totN     = rootgrp.createVariable('totN', 'f8', ('time',))

# coarse
if args.ana == 'coarse':
    varM     = rootgrp.createVariable('varM', 'f8', ('time','n','x','y'))
    varN     = rootgrp.createVariable('varN', 'f8', ('time','n','x','y'))
    varm     = rootgrp.createVariable('varm', 'f8', ('time','n','x','y'))
    meanN    = rootgrp.createVariable('meanN', 'f8', ('time','n','x','y'))
    meanM    = rootgrp.createVariable('meanM', 'f8', ('time','n','x','y'))
    meanm    = rootgrp.createVariable('meanm', 'f8', ('time','n','x','y'))
    varQmp   = rootgrp.createVariable('varQmp', 'f8', ('time','n','x','y'))
    meanQmp  = rootgrp.createVariable('meanQmp', 'f8', ('time','n','x','y'))
    #varQtot  = rootgrp.createVariable('varQtot', 'f8', ('time','n','x','y'))
    #meanQtot = rootgrp.createVariable('meanQtot', 'f8', ('time','n','x','y'))
    
if args.ana == 'vert':
    levdim = rootgrp.createDimension('levs', len(levlist))
    levs     = rootgrp.createVariable('levs', 'i4', ('levs',))
    levs[:]  = levlist
    height   = rootgrp.createVariable('height', 'i4', ('levs',))
    height[:]= heightlist
    Mtot     = rootgrp.createVariable('Mtot', 'f8', ('time', 'levs'))
    Msouth   = rootgrp.createVariable('Msouth', 'f8', ('time', 'levs'))
    Mnorth   = rootgrp.createVariable('Mnorth', 'f8', ('time', 'levs'))
    
if args.ana == 'hypo':
    cld_size = rootgrp.createVariable('cld_size', 'f8', ('time','N_cld'))
    rdf      = rootgrp.createVariable('rdf', 'f8', ('time','dr'))
    varM     = rootgrp.createVariable('varM', 'f8', ('time','n','x','y'))
    varN     = rootgrp.createVariable('varN', 'f8', ('time','n','x','y'))
    varm     = rootgrp.createVariable('varm', 'f8', ('time','n','x','y'))
    meanN    = rootgrp.createVariable('meanN', 'f8', ('time','n','x','y'))
    meanM    = rootgrp.createVariable('meanM', 'f8', ('time','n','x','y'))
    meanm    = rootgrp.createVariable('meanm', 'f8', ('time','n','x','y'))

# currently not used
#hpbl     = rootgrp.createVariable('hpbl', 'f8', ('time','levs','n','x','y'))
#enstauc  = rootgrp.createVariable('enstauc', 'f8', ('time', 'x', 'y'))
#acf2d    = rootgrp.createVariable('acf2d', 'f8', ('time','levs','n','drcorr'))
#Mmem1    = rootgrp.createVariable('Mmem1', 'f8', ('time','levs','n','x','y'))

# End allocation
################################################################################


# Preliminaries before time loop
if args.ana == 'prec':
    # Load radar data for all times
    dateobj = yyyymmddhh_strtotime(args.date)
    dtradar = timedelta(minutes = 10)
    radarts = getfobj_ncdf_timeseries(radarpref, dateobj+tstart-dtradar, 
                                    dateobj+tend-dtradar, tinc, 
                                        reftime=dateobj, ncdffn_sufx=radarsufx, 
                                        fieldn = 'pr', abs_datestr='yymmddhhmm',
                                        dwdradar = True)
    # Get mask
    radarmask = get_totmask(radarts)
    
    # Crop data
    radarmask = radarmask[lx1+62:lx2-42, ly1+22:ly2-42]
    for i in range(len(radarts)):
        radarts[i].data = radarts[i].data[lx1+62:lx2-42, ly1+22:ly2-42]
        
    # Allocate lists to get time average histograms
    tothist_model = []
    tothist_obs = []


###################
## Time loop      #
###################
for it, t in enumerate(timelist):
    print 'time: ', t
    ############################################################################
    # Load COSMO data
    ncdffn = 'lfff' + ddhhmmss(t) + sufx
    
    if args.ana in ['clouds', 'coarse']:
        # Load fields required for mass flux analysis
        fieldlist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                    dir_suffix='/OUTPUT/', fieldn = fieldn, 
                                    nfill=1, levs = lev, return_arrays=True)
        
        # Crop all fields to analysis domain
        for i in range(args.nens):
            fieldlist[i] = fieldlist[i][lx1:lx2, ly1:ly2]
            
        qclist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                  dir_suffix='/OUTPUT/', fieldn = 'QC', 
                                  nfill=1, levs = lev, return_arrays = True)
        # Add QI and QS
        qilist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                  dir_suffix='/OUTPUT/', fieldn = 'QI', 
                                  nfill=1, levs = lev, return_arrays = True)
        qslist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                  dir_suffix='/OUTPUT/', fieldn = 'QS', 
                                  nfill=1, levs = lev, return_arrays = True)
        for i in range(args.nens):
            qclist[i] = (qclist[i][lx1:lx2, ly1:ly2] + 
                         qilist[i][lx1:lx2, ly1:ly2] + 
                         qslist[i][lx1:lx2, ly1:ly2])

        del qilist
        del qslist
        ncdffn_buoy = ncdffn + '_buoy'
        rholist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_buoy, 
                                   dir_suffix='/OUTPUT/', fieldn = 'RHO', 
                                   nfill=1, levs=lev, return_arrays=True)
        
        for i in range(args.nens):
            rholist[i] = rholist[i][lx1:lx2, ly1:ly2]
            
    if args.ana == 'vert':
        # Load fields required for mass flux analysis
        fieldlist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                    dir_suffix='/OUTPUT/', fieldn = fieldn, 
                                    nfill=1, levs = levlist, return_arrays=True)
        
        # Crop all fields to analysis domain
        for i in range(args.nens):
            fieldlist[i] = fieldlist[i][:,lx1:lx2, ly1:ly2]
            
        qclist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                  dir_suffix='/OUTPUT/', fieldn = 'QC', 
                                  nfill=1, levs = levlist, return_arrays = True)
        # Add QI and QS
        qilist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                  dir_suffix='/OUTPUT/', fieldn = 'QI', 
                                  nfill=1, levs = levlist, return_arrays = True)
        qslist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                  dir_suffix='/OUTPUT/', fieldn = 'QS', 
                                  nfill=1, levs = levlist, return_arrays = True)
        for i in range(args.nens):
            qclist[i] = (qclist[i][:,lx1:lx2, ly1:ly2] + 
                         qilist[i][:,lx1:lx2, ly1:ly2] + 
                         qslist[i][:,lx1:lx2, ly1:ly2])

        del qilist
        del qslist
        ncdffn_buoy = ncdffn + '_buoy'
        rholist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_buoy, 
                                   dir_suffix='/OUTPUT/', fieldn = 'RHO', 
                                   nfill=1, levs=levlist, return_arrays=True)
        
        for i in range(args.nens):
            rholist[i] = rholist[i][:,lx1:lx2, ly1:ly2]
            
    if args.ana == 'coarse':
        # Get vertically integrated Q    
        Qmplist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_buoy, 
                                  dir_suffix='/OUTPUT/', fieldn = 'TTENS_MPHY', 
                                  nfill=1, return_arrays = True)
        #Qtotlist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_buoy, 
                                  #dir_suffix='/OUTPUT/', fieldn = 'TTENS_DIAB', 
                                  #nfill=1, return_arrays = True)
        for i in range(args.nens):
            Qmplist[i] = np.mean(Qmplist[i][:, lx1:lx2, ly1:ly2], axis = 0)
            #Qtotlist[i] = np.mean(Qtotlist[i][:, lx1:lx2, ly1:ly2], axis = 0)

    if args.ana in ['weather', 'prec', 'spectra']:
        # Load precipitation data
        ncdffn_surf = ncdffn + '_surf'
        preclist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_surf, 
                                    dir_suffix='/OUTPUT/', fieldn='PREC_ACCUM', 
                                    nfill=1, return_arrays=True)
        for i in range(args.nens):
            preclist[i] = preclist[i][lx1:lx2, ly1:ly2]
        
    if args.ana == 'weather':
        # Load weather info data
        tauclist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_surf, 
                                    dir_suffix='/OUTPUT/', fieldn = 'TAU_C', 
                                    nfill=1, return_arrays=True)
        hpbllist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_surf, 
                                    dir_suffix='/OUTPUT/', fieldn = 'HPBL', 
                                    nfill=1, return_arrays=True)
        capelist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_surf, 
                                    dir_suffix='/OUTPUT/', fieldn = 'CAPE_ML', 
                                    nfill=1, return_arrays=True)
        for i in range(args.nens):
            tauclist[i] = tauclist[i][lx1:lx2, ly1:ly2]
            hpbllist[i] = hpbllist[i][lx1:lx2, ly1:ly2]
            capelist[i] = capelist[i][lx1:lx2, ly1:ly2]
    
    if args.ana == 'spectra':
        # Load U and V
        ncdffn_uv = ncdffn + '_uv'
        ulist = getfobj_ncdf_ens(ensdir, 'sub', 5, ncdffn_uv, 
                                    dir_suffix='/OUTPUT/', fieldn = 'U', 
                                    nfill=1, return_arrays = True)
        vlist = getfobj_ncdf_ens(ensdir, 'sub', 5, ncdffn_uv, 
                                    dir_suffix='/OUTPUT/', fieldn = 'V', 
                                    nfill=1, return_arrays = True)
        for i in range(5):
            ulist[i] = ulist[i][:,lx1:lx2, ly1:ly2]
            vlist[i] = vlist[i][:,lx1:lx2, ly1:ly2]
    
    if args.ana == 'hypo':
        fieldlist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                    dir_suffix='/OUTPUT/', fieldn = 'm', 
                                    nfill=1, return_arrays = True)
        for i in range(args.nens):
            fieldlist[i] = fieldlist[i][0,lx1:lx2, ly1:ly2]
    
    # End loading data
    ############################################################################
    
    ############################################################################
    # Do analysis
    
    if args.ana == 'weather':
        # Calculate ensemble and domain mean and save for timestep
        ditauc[it] = np.nanmean(tauclist)
        dihpbl[it] = np.nanmean(hpbllist)
        dicape[it] = np.nanmean(capelist)
        diprec[it] = np.nanmean(preclist)
        # weather is DONE!!!
        
    if args.ana == 'prec':
        # 1. histograms
        hist_tmp = []
        for i in range(args.nens):
            hist_tmp.append(np.histogram(preclist[i][~radarmask], 
                                         histbinedges)[0])
        # This is now the model hist for one time step
        tothist_model.append(np.mean(hist_tmp, axis = 0)) 
        radarfield = radarts[it].data
        # This is now the obs hist for one time step
        tothist_obs.append(np.histogram(radarfield[~radarmask], 
                                         histbinedges)[0])
        
        # 2. rdf
        rdf_prec_modellist = []
        for field in preclist:
            # Identify clouds
            tmpfield = field
            tmpfield[~radarmask] == 0.
            tmp = identify_clouds(tmpfield, 1., water = args.water)
            labels, cld_size_mem, cld_sum_mem = tmp
            g, r = calc_rdf(labels, tmpfield, normalize = True, rmax = rmax_rdf, 
                            dr = dr_rdf)
            rdf_prec_modellist.append(g)
        rdf_prec_model[it, :] = np.mean(rdf_prec_modellist, axis = 0)
        
        # Now for the observation field
        tmpfield = radarfield
        tmpfield[~radarmask] == 0.
        tmp = identify_clouds(tmpfield, 1., water = args.water)
        labels, cld_size_mem, cld_sum_mem = tmp
        g, r = calc_rdf(labels, tmpfield, normalize = True, rmax = rmax_rdf, 
                        dr = dr_rdf)
        rdf_prec_obs[it, :] = g
        dr[:] = r   # km
        # rdf is DONE!!!
        
    if args.ana == 'spectra':
        # 1. Calculate DKE spectra
        vertlim = 15
        # a. Get ensemble average backgroud KE spectrum
        kelist = []
        for u, v in zip(ulist, vlist):
            vertlist = []
            for k in range(vertlim, v.shape[0]):
                p, kspec, s = powspec_2d_hor(u[k,:,:], v[k,:,:], dx, dx)
                vertlist.append(p)
            kelist.append(np.mean(vertlist, axis = 0))
        bgkespec[it,:] = np.mean(kelist, axis = 0)
        
        # b. Get ensemble mean difference KE spectrum
        dkelist = []
        for i in range(len(ulist)-1):
            for j in range(i+1, len(ulist)):
                du = ulist[i] - ulist[j]
                dv = vlist[i] - vlist[j]
                vertlist = []
                for k in range(vertlim, v.shape[0]):
                    p, kspec, s = powspec_2d_hor(du[k,:,:], dv[k,:,:], dx, dx)
                    vertlist.append(p)
                dkelist.append(np.mean(vertlist, axis = 0))
        dkespec[it,:] = np.mean(dkelist, axis = 0)
                
        # 2. Calculate Precipitation spectra
        # a. Get ensemble average backgroud KE spectrum
        kelist = []
        for prec in preclist[:5]:
            p, kspec, s = powspec_2d_hor_alter(prec[:,:], dx, dx)
            kelist.append(p)
        bgprecspec[it,:] = np.mean(kelist, axis = 0)
        
        # b. Get ensemble mean difference KE spectrum
        dkelist = []
        for i in range(len(preclist[:5])-1):
            for j in range(i+1, len(preclist[:5])):
                dprec = preclist[i] - preclist[j]
                p, kspec, s = powspec_2d_hor_alter(dprec, dx, dx)
                dkelist.append(p)
        dprecspec[it,:] = np.mean(dkelist, axis = 0)
        
        speck[:] = kspec
        speclam[:] = s
        
        
    if args.ana in ['clouds', 'coarse']:
        # Member loop
        sizelist = []
        sumlist = []
        rdflist = []
        labelslist = []   # Save for use later
        comlist = []      # Save for use later
        for field, qc, rho, imem in zip(fieldlist, qclist, rholist, 
                                  range(len(fieldlist))):
            # Identify clouds
            
            if imem == 0 and args.ana == 'clouds':
                exw[it,:,:] = field
                exq[it,:,:] = qc
                tmp = identify_clouds(field, thresh, qc,
                                    opt_thresh = 0., water = False,
                                    rho = rho)
                excld[it,:,:] = tmp[0]
                
            tmp = identify_clouds(field, thresh, qc,
                                    opt_thresh = 0., water = args.water,
                                    rho = rho)
            labels, cld_size_mem, cld_sum_mem = tmp
            if imem == 0 and args.ana == 'clouds':
                exwater[it,:,:] = labels
            cld_sum_mem *= dx*dx  # Rho is now already included
            sizelist.append(cld_size_mem)
            sumlist.append(cld_sum_mem)
            
            labelslist.append(labels)
            # Calculate centers of mass
            num = np.unique(labels).shape[0]   # Number of clouds
            com = np.array(center_of_mass(field, labels, range(1,num)))
            if com.shape[0] == 0:   # Accout for empty arrays
                com = np.empty((0,2))
            comlist.append(com)
            
            if args.ana == 'clouds':
                # Calculate RDF
                g, r = calc_rdf(labels, field, normalize = True, rmax = rmax_rdf, 
                                dr = dr_rdf)
                rdflist.append(g)
                dr[:] = r   # km
        
        if args.ana == 'clouds':
            # Save lists and mean rdf
            ntmp = len([i for sl in sumlist for i in sl])
            print ntmp
            cld_size[it, :ntmp] = [i for sl in sizelist for i in sl]  # Flatten
            cld_sum[it, :ntmp] = [i for sl in sumlist for i in sl]
            rdf[it, :] = np.mean(rdflist, axis = 0)
            totN[it] = ntmp/float(args.nens)
    
        if args.ana == 'coarse':
            ########################################################################
            # Calculate cloud statistics
            
            ################
            ## n loop      #
            ################
            for i_n, n in enumerate(nlist):
                print 'n: ', n
                ####################################################################
                # Calculate coarse variances and means
                # Determine size of coarse arrays
                nx = int(np.floor(256/n))
                ny = int(np.floor(256/n))

                # Loop over coarse grid boxes
                for ico  in range(nx):
                    for jco in range(ny):
                        # Get limits for each N box
                        xmin = ico*n
                        xmax = (ico+1)*n
                        ymin = jco*n
                        ymax = (jco+1)*n
                        
                        # These are the ensemble lists for each box
                        tmp_cldlist = []
                        tmp_Mlist = []
                        tmp_Nlist = []
                        tmp_Qmplist = []
                        # Loop over members
                        for field, labels, com, cld_sum_mem, imem, Qmpfield in \
                                                    zip(fieldlist, 
                                                        labelslist,
                                                        comlist, 
                                                        sumlist,
                                                        range(args.nens),
                                                        Qmplist):
                                                        
                            # Get the collapsed clouds for each box
                            bool_arr = ((com[:,0]>=xmin)&(com[:,0]<xmax)&
                                        (com[:,1]>=ymin)&(com[:,1]<ymax))
                            
                            # This is the array with all clouds for this box and member
                            box_cld_sum = cld_sum_mem[bool_arr]
                            
                            # This lists then contains all clouds for all members in a box
                            tmp_cldlist += list(box_cld_sum)
                            
                            # If the array is empty set M to zero
                            if len(box_cld_sum) > 0:
                                tmp_Mlist.append(np.sum(box_cld_sum))
                            else:
                                tmp_Mlist.append(0.)
                            
                            # This is the number of clouds
                            tmp_Nlist.append(box_cld_sum.shape[0])
                            
                            # This is the MEAN heating rate
                            tmp_Qmplist.append(np.mean(Qmpfield[ico*n:(ico+1)*n, 
                                                        jco*n:(jco+1)*n]))
                            # End member loop #############
                        
                        # Now convert the list with all clouds for this box
                        tmp_cldlist = np.array(tmp_cldlist)
                        
                        # Calculate statistics and save them in ncdf file
                        # Check if x number of members have clouds in them
                        
                        if np.sum(np.array(tmp_Nlist)>0) >= args.minmem:
                            varM[it,i_n,ico,jco] = np.var(tmp_Mlist, ddof = 1)
                            varN[it,i_n,ico,jco] = np.var(tmp_Nlist, ddof = 1)
                            varm[it,i_n,ico,jco] = np.var(tmp_cldlist, ddof = 1)
                            meanM[it,i_n,ico,jco] = np.mean(tmp_Mlist)
                            meanm[it,i_n,ico,jco] = np.mean(tmp_cldlist)
                            meanN[it,i_n,ico,jco] = np.mean(tmp_Nlist)
                        else:
                            varM[it,i_n,ico,jco] = np.nan
                            varN[it,i_n,ico,jco] = np.nan
                            varm[it,i_n,ico,jco] = np.nan
                            meanM[it,i_n,ico,jco] = np.nan
                            meanm[it,i_n,ico,jco] = np.nan
                            meanN[it,i_n,ico,jco] = np.nan
                        # This means NaNs only appear when minmem criterion is not met    
                        
                        varQmp[it,i_n,ico,jco] = np.var(tmp_Qmplist, ddof = 1)
                        meanQmp[it,i_n,ico,jco] = np.mean(tmp_Qmplist)

    if args.ana == 'hypo':
        # Member loop
        sizelist = []
        rdflist = []
        labelslist = []   # Save for use later
        comlist = []      # Save for use later
        for imem, field in enumerate(fieldlist):
            # Identify clouds
            tmp = identify_clouds(field, thresh, water = args.water)
            labels, cld_size_mem, cld_sum_mem = tmp

            sizelist.append(cld_size_mem)
            
            labelslist.append(labels)
            # Calculate centers of mass
            num = np.unique(labels).shape[0]   # Number of clouds
            com = np.array(center_of_mass(field, labels, range(1,num)))
            if com.shape[0] == 0:   # Accout for empty arrays
                com = np.empty((0,2))
            comlist.append(com)

            g, r = calc_rdf(labels, field, normalize = True, rmax = rmax_rdf, 
                            dr = dr_rdf)
            rdflist.append(g)
            dr[:] = r   # km
        
        
        ntmp = len([i for sl in sizelist for i in sl])
        cld_size[it, :ntmp] = [i for sl in sizelist for i in sl]  # Flatten
        rdf[it, :] = np.mean(rdflist, axis = 0)

        
        ########################################################################
        # Calculate cloud statistics
        
        ################
        ## n loop      #
        ################
        for i_n, n in enumerate(nlist):
            print 'n: ', n
            ####################################################################
            # Calculate coarse variances and means
            # Determine size of coarse arrays
            nx = int(np.floor(256/n))
            ny = int(np.floor(256/n))

            # Loop over coarse grid boxes
            for ico  in range(nx):
                for jco in range(ny):
                    # Get limits for each N box
                    xmin = ico*n
                    xmax = (ico+1)*n
                    ymin = jco*n
                    ymax = (jco+1)*n
                    
                    # These are the ensemble lists for each box
                    tmp_cldlist = []
                    tmp_Mlist = []
                    tmp_Nlist = []
                    # Loop over members
                    for field, labels, com, cld_size_mem, imem,  in \
                                                zip(fieldlist, 
                                                    labelslist,
                                                    comlist, 
                                                    sizelist,
                                                    range(args.nens),
                                                    ):
                                                    
                        # Get the collapsed clouds for each box
                        bool_arr = ((com[:,0]>=xmin)&(com[:,0]<xmax)&
                                    (com[:,1]>=ymin)&(com[:,1]<ymax))
                        
                        # This is the array with all clouds for this box and member
                        box_cld_sum = cld_size_mem[bool_arr]
                        
                        # This lists then contains all clouds for all members in a box
                        tmp_cldlist += list(box_cld_sum)
                        
                        # If the array is empty set M to zero
                        if len(box_cld_sum) > 0:
                            tmp_Mlist.append(np.sum(box_cld_sum))
                        else:
                            tmp_Mlist.append(0.)
                        
                        # This is the number of clouds
                        tmp_Nlist.append(box_cld_sum.shape[0])

                        # End member loop #############
                    
                    # Now convert the list with all clouds for this box
                    tmp_cldlist = np.array(tmp_cldlist)
                    
                    # Calculate statistics and save them in ncdf file
                    # Check if x number of members have clouds in them
                    
                    if np.sum(np.array(tmp_Nlist)>0) >= args.minmem:
                        varM[it,i_n,ico,jco] = np.var(tmp_Mlist, ddof = 1)
                        varN[it,i_n,ico,jco] = np.var(tmp_Nlist, ddof = 1)
                        varm[it,i_n,ico,jco] = np.var(tmp_cldlist, ddof = 1)
                        meanM[it,i_n,ico,jco] = np.mean(tmp_Mlist)
                        meanm[it,i_n,ico,jco] = np.mean(tmp_cldlist)
                        meanN[it,i_n,ico,jco] = np.mean(tmp_Nlist)
                    else:
                        varM[it,i_n,ico,jco] = np.nan
                        varN[it,i_n,ico,jco] = np.nan
                        varm[it,i_n,ico,jco] = np.nan
                        meanM[it,i_n,ico,jco] = np.nan
                        meanm[it,i_n,ico,jco] = np.nan
                        meanN[it,i_n,ico,jco] = np.nan
                    # This means NaNs only appear when minmem criterion is not met    


        # End coarse upscaled variances and means
        ####################################################################
    
    if args.ana == 'vert':
        for iz, lev in enumerate(levlist):
            print 'lev', lev
            # Member loop
            sumlist = []
            sumlist_N = []
            sumlist_S = []
            for field, qc, rho, imem in zip(fieldlist, qclist, rholist, 
                                    range(len(fieldlist))):
                # Identify clouds
                tmp = identify_clouds(field[iz], thresh, qc[iz],
                                        opt_thresh = 0., water = args.water,
                                        rho = rho[iz])
                labels, cld_size_mem, cld_sum_mem = tmp
                cld_sum_mem *= dx*dx  # Rho is now already included
                
                # Calculate centers of mass
                num = np.unique(labels).shape[0]   # Number of clouds
                com = np.array(center_of_mass(field[iz], labels, range(1,num)))
                if com.shape[0] == 0:   # Accout for empty arrays
                    com = np.empty((0,2))
                
                sumlist += list(cld_sum_mem)
                bool_N = com[:,0]>=256/2
                sumlist_N += list(cld_sum_mem[bool_N])
                sumlist_S += list(cld_sum_mem[~bool_N])
            Mtot[it,iz] = np.sum(sumlist) / args.nens
            Mnorth[it,iz] = np.sum(sumlist_N) / args.nens
            Msouth[it,iz] = np.sum(sumlist_S) / args.nens
        
    # End do analysis
    ############################################################################

# Get statistics over time
if args.ana == 'prec':
    tothist_model = np.mean(tothist_model, axis = 0)
    hist_model[:] = tothist_model
    tothist_obs = np.mean(tothist_obs, axis = 0)
    hist_obs[:] = tothist_obs
# Close ncdf file
rootgrp.close()
            

        
        
        
        
