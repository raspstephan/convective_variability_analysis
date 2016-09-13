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
from cosmo_utils.pyncdf import getfobj_ncdf_ens, getfobj_ncdf
from cosmo_utils.helpers import make_timelist, ddhhmmss
from cosmo_utils.diag import identify_clouds, calc_rdf, crosscor, int_rad_2d
from scipy.ndimage.measurements import center_of_mass
from scipy.signal import correlate


# Setup
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--ana', metavar = 'ana', type=str)
parser.add_argument('--date', metavar = 'date', type=str)
parser.add_argument('--height', metavar = 'height', type=float, nargs = '+',
                    default = [3000])
parser.add_argument('--water', metavar = 'water', type=bool, default = True)
parser.add_argument('--nens', metavar = 'nens', type=int, default = 20)
parser.add_argument('--tstart', metavar = 'tstart', type=int, default = 1)
parser.add_argument('--tend', metavar = 'tend', type=int, default = 24)
parser.add_argument('--tinc', metavar = 'tinc', type=int, default = 60)
args = parser.parse_args()


# Functions
def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile


# Create file str
savedir = '/home/scratch/users/stephan.rasp/results/'
heightstr = ''
for h in args.height:
    heightstr += str(int(h)) + '_'
savestr = (args.date + '_ana-' + args.ana + '_wat-' + str(args.water) + 
           '_height-' + heightstr +
           'nens-' + str(args.nens) + '_tstart-' + str(args.tstart) + 
           '_tend-' + str(args.tend) + '_tinc-' + str(args.tinc) + '.nc')

# Convert times to timedelta objects
tstart = timedelta(hours = args.tstart)   # Cannot be 0 because of tau_c calculation!
tend = timedelta(hours = args.tend)  
tinc = timedelta(minutes = args.tinc)  # temporal resolution for analysis


# Make lists for loops, dimensions
timelist = make_timelist(tstart, tend, tinc)
nlist = [256, 128, 64, 32, 16, 8, 4]
#hlist = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000,
            #8000, 9000, 10000]
hlist = args.height

ensdir = '/home/scratch/users/stephan.rasp/' + args.date + '/deout_ceu_pspens/'


# Analysis-specific setup
dx = 2800.
if args.ana == 'm':
    HH = getfobj_ncdf(ensdir + '/1/OUTPUT/lfff00000000c.nc_30m', 'HHL')
    aslcol = HH.data[:, -1, -1]   # ATTENTION Hard coded colums above sea level 
    levlist = []
    realhlist = []
    for h in hlist:
        levlist.append(np.argmin(np.abs(aslcol-h)))   # find closest level
        realhlist.append(aslcol[np.argmin(np.abs(aslcol-h))])
    sufx = '.nc_30m'
    fieldn = 'W'
    thresh = 1.
    HHcropped = HH.data[-1,50:-51, 50:-51]
    HH50tot = np.mean(HHcropped[:,:])
    HH50south = np.mean(HHcropped[:256/2, :])
    HH50north = np.mean(HHcropped[256/2:, :])
    print 'tot', HH50tot, 'south', HH50south, 'north', HH50north
elif args.ana == 'hypo':
    levlist = [0]
    sufx = '.nc'
    fieldn = 'm'
    ensdir = ('/home/scratch/users/stephan.rasp/hypo_' + args.date + 
              '/deout_ceu_pspens/')
    thresh = 0.
elif args.ana == 'p':
    levlist = [None]
    realhlist = ['surf']
    raise Exception, 'p not implemented yet'
else:
    raise Exception, 'wrong analysis'

################################################################################
# Allocate NetCDF file
rootgrp = Dataset(savedir + savestr, 'w', format='NETCDF4')
# Create dimensions
tdim = rootgrp.createDimension('time', len(timelist))
ndim = rootgrp.createDimension('n', len(nlist))
levdim = rootgrp.createDimension('levs', len(levlist))
xdim = rootgrp.createDimension('x', nlist[0])
ydim = rootgrp.createDimension('y', nlist[0])
nclddim = rootgrp.createDimension('N_cld', 1e6)
drdim = rootgrp.createDimension('dr', 30/2+1) # For RDF
drcorrdim = rootgrp.createDimension('drcorr', nlist[0]) # For 2D ACF


# Create variables and add attributes 
time     = rootgrp.createVariable('time', 'f8', ('time',))
time[:]  = [td.total_seconds() for td in timelist]

n        = rootgrp.createVariable('n', 'i4', ('n',))
n[:]     = nlist

levs     = rootgrp.createVariable('levs', 'i4', ('levs',))
levs[:]  = levlist

dr       = rootgrp.createVariable('dr', 'f4', ('dr',))

ditauc   = rootgrp.createVariable('ditauc', 'f8', ('time'))
enstauc  = rootgrp.createVariable('enstauc', 'f8', ('time', 'x', 'y'))
cld_size = rootgrp.createVariable('cld_size', 'f8', ('time','levs','N_cld'))
cld_sum  = rootgrp.createVariable('cld_sum', 'f8', ('time','levs','N_cld'))
rdf      = rootgrp.createVariable('rdf', 'f8', ('time','levs','dr'))
acf2d    = rootgrp.createVariable('acf2d', 'f8', ('time','levs','n','drcorr'))
varM     = rootgrp.createVariable('varM', 'f8', ('time','levs','n','x','y'))
varN     = rootgrp.createVariable('varN', 'f8', ('time','levs','n','x','y'))
varm     = rootgrp.createVariable('varm', 'f8', ('time','levs','n','x','y'))
meanN    = rootgrp.createVariable('meanN', 'f8', ('time','levs','n','x','y'))
meanM    = rootgrp.createVariable('meanM', 'f8', ('time','levs','n','x','y'))
meanm    = rootgrp.createVariable('meanm', 'f8', ('time','levs','n','x','y'))
varQmp   = rootgrp.createVariable('varQmp', 'f8', ('time','levs','n','x','y'))
meanQmp  = rootgrp.createVariable('meanQmp', 'f8', ('time','levs','n','x','y'))
varQtot  = rootgrp.createVariable('varQtot', 'f8', ('time','levs','n','x','y'))
meanQtot = rootgrp.createVariable('meanQtot', 'f8', ('time','levs','n','x','y'))
Mtot     = rootgrp.createVariable('Mtot', 'f8', ('time','levs'))
Msouth   = rootgrp.createVariable('Msouth', 'f8', ('time','levs'))
Mnorth   = rootgrp.createVariable('Mnorth', 'f8', ('time','levs'))

Mmem1    = rootgrp.createVariable('Mmem1', 'f8', ('time','levs','n','x','y'))

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
    fieldlist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                 dir_suffix='/OUTPUT/', fieldn = fieldn, 
                                 nfill=1, levs = levlist, return_arrays = True)
    
    # Crop all fields to analysis domain
    sxo, syo = fieldlist[0][0].shape  # Original field shape
    lx1 = (sxo-256-1)/2 # ATTENTION first dimension is actually y
    lx2 = -(lx1+1) # Number of grid pts to exclude at border
    ly1 = (syo-256-1)/2
    ly2 = -(ly1+1)
    for i in range(len(fieldlist)):
        fieldlist[i] = fieldlist[i][:, lx1:lx2, ly1:ly2]
        
    if args.ana == 'm':   # Load some extra variables 
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
        for i in range(len(qclist)):
            qclist[i] = (qclist[i][:, lx1:lx2, ly1:ly2] + 
                         qilist[i][:, lx1:lx2, ly1:ly2] + 
                         qslist[i][:, lx1:lx2, ly1:ly2])
        del qilist
        del qslist
        ncdffn_buoy = ncdffn + '_buoy'
        rholist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_buoy, 
                                   dir_suffix='/OUTPUT/', fieldn = 'RHO', 
                                   nfill=1, levs = levlist, return_arrays = True)
        
        for i in range(len(rholist)):
            rholist[i] = rholist[i][:, lx1:lx2, ly1:ly2]
            
        # Get vertically integrated Q    
        Qmplist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_buoy, 
                                  dir_suffix='/OUTPUT/', fieldn = 'TTENS_MPHY', 
                                  nfill=1, return_arrays = True)
        Qtotlist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_buoy, 
                                  dir_suffix='/OUTPUT/', fieldn = 'TTENS_DIAB', 
                                  nfill=1, return_arrays = True)
        for i in range(len(Qmplist)):
            Qmplist[i] = np.mean(Qmplist[i][:, lx1:lx2, ly1:ly2], axis = 0)
            Qtotlist[i] = np.mean(Qtotlist[i][:, lx1:lx2, ly1:ly2], axis = 0)
        
    else:   # Fill lists with None
        qclist = [None]*len(fieldlist)
        rholist = [None]*len(fieldlist)
        Qmplist = [None]*len(fieldlist)
        Qtotlist = [None]*len(fieldlist)
        
    # Load tau_c data
    if not args.ana == 'hypo':
        ncdffn_surf = ncdffn + '_surf'
        tauclist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn_surf, 
                                    dir_suffix='/OUTPUT/', fieldn = 'TAU_C', 
                                    nfill=1, levs = levlist, return_arrays = True)
        for i in range(len(tauclist)):
            tauclist[i] = tauclist[i][lx1:lx2, ly1:ly2]
    
        

    # End loading data
    ############################################################################
    
    ############################################################################
    # Calculate mean tau_c and save
    if not args.ana == 'hypo':
        ditauc[it] = np.nanmean(tauclist)
        enstauc[it] = np.nanmean(tauclist, axis = 0)
    # End calculate mean tau_c and save
    ############################################################################

    ####################
    ## lev loop        #
    ####################
    for iz, lev in enumerate(levlist):
        print 'lev: ', lev
        ########################################################################
        # Calculate cloud statistics
        
        # Member loop
        sizelist = []
        sumlist = []
        rdflist = []
        labelslist = []   # Save for use later
        comlist = []      # Save for use later
        for field, qc, rho in zip(fieldlist, qclist, rholist):
            # Identify clouds
            if args.ana == 'm':
                tmp = identify_clouds(field[iz], thresh, qc[iz],
                                      opt_thresh = 0., water = args.water,
                                      rho = rho[iz])
                labels, cld_size_mem, cld_sum_mem = tmp
                cld_sum_mem *= dx*dx  # Rho is now already included
            else:
                tmp = identify_clouds(field[iz], thresh, water = args.water)
                labels, cld_size_mem, cld_sum_mem = tmp
            sizelist.append(cld_size_mem)
            sumlist.append(cld_sum_mem)
            
            labelslist.append(labels)
            # Calculste centers of mass
            num = np.unique(labels).shape[0]   # Number of clouds
            com = np.array(center_of_mass(field[iz], labels, range(1,num)))
            if com.shape[0] == 0:   # Accout for empty arrays
                com = np.empty((0,2))
            comlist.append(com)
            
            # Calculate RDF
            g, r = calc_rdf(labels, field[iz], normalize = True, rmax = 30, 
                            dr = 2)
            rdflist.append(g)
            dr[:] = r * 2.8   # km
        
        # Save lists and mean rdf
        ntmp = len([i for sl in sumlist for i in sl])
        cld_size[it, iz, :ntmp] = [i for sl in sizelist for i in sl]  # Flatten
        cld_sum[it, iz, :ntmp] = [i for sl in sumlist for i in sl]
        rdf[it, iz, :] = np.mean(rdflist, axis = 0)
        # End calculate cloud statistics
        ########################################################################
        
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
            
            # Member loop
            varmlist = []
            mlist = []
            Mlist = []
            Nlist = []
            
            # NOTE I need all m's for every coarse box, then I can calculate M, m and var(m) and N
            
            # Loop over coarse grid boxes
            # Allocate coarse arrays
            nmem = len(fieldlist)
            Mmem_coarse = np.empty((nmem, nx, ny))# These are for the ACF correlation without the nan filter
            for ico  in range(nx):
                for jco in range(ny):
                    # Get limits for each N box
                    xmin = ico*n
                    xmax = (ico+1)*n
                    ymin = jco*n
                    ymax = (jco+1)*n
                    
                    tmp_cldlist = []
                    tmp_Mlist = []
                    tmp_Nlist = []
                    tmp_Qmplist = []
                    tmp_Qtotlist = []
                    # Loop over members
                    for field, labels, com, cld_sum_mem, imem, Qmpfield, Qtotfield in zip(fieldlist, 
                                                               labelslist,
                                                               comlist, 
                                                               sumlist,
                                                               range(nmem),
                                                               Qmplist,
                                                               Qtotlist):
                        # Get the collapsed clouds for each box
                        bool_arr = ((com[:,0]>=xmin)&(com[:,0]<xmax)&
                                    (com[:,1]>=ymin)&(com[:,1]<ymax))
                        # This lists then contains all clouds for all members in a box
                        box_cld_sum = cld_sum_mem[bool_arr]
                        tmp_cldlist += list(box_cld_sum)
                        if len(box_cld_sum) > 0:
                            tmp_Mlist.append(np.sum(box_cld_sum))
                            Mmem_coarse[imem, ico, jco] = np.sum(box_cld_sum)
                        else:
                            tmp_Mlist.append(0.)
                            Mmem_coarse[imem, ico, jco] = 0.
                        tmp_Nlist.append(box_cld_sum.shape[0])
                        tmp_Qmplist.append(np.mean(Qmpfield[ico*n:(ico+1)*n, 
                                                       jco*n:(jco+1)*n]))
                        tmp_Qtotlist.append(np.mean(Qtotfield[ico*n:(ico+1)*n, 
                                                       jco*n:(jco+1)*n]))
                        # End member loop
                    
                    tmp_cldlist = np.array(tmp_cldlist)
                    # Calculate statistics and save them in ncdf file
                    # Check if x number of members have clouds in them
                    min_mem = 5
                    if np.sum(np.array(tmp_Nlist)>0) >= min_mem:
                        varM[it,iz,i_n,ico,jco] = np.var(tmp_Mlist, ddof = 1)
                        varN[it,iz,i_n,ico,jco] = np.var(tmp_Nlist, ddof = 1)
                        varm[it,iz,i_n,ico,jco] = np.var(tmp_cldlist, ddof = 1)
                        meanM[it,iz,i_n,ico,jco] = np.mean(tmp_Mlist)
                        meanm[it,iz,i_n,ico,jco] = np.mean(tmp_cldlist)
                        meanN[it,iz,i_n,ico,jco] = np.mean(tmp_Nlist)
                    else:
                        varM[it,iz,i_n,ico,jco] = np.nan
                        varN[it,iz,i_n,ico,jco] = np.nan
                        varm[it,iz,i_n,ico,jco] = np.nan
                        meanM[it,iz,i_n,ico,jco] = np.nan
                        meanm[it,iz,i_n,ico,jco] = np.nan
                        meanN[it,iz,i_n,ico,jco] = np.nan
                    varQmp[it,iz,i_n,ico,jco] = np.var(tmp_Qmplist, ddof = 1)
                    meanQmp[it,iz,i_n,ico,jco] = np.mean(tmp_Qmplist)
                    varQtot[it,iz,i_n,ico,jco] = np.var(tmp_Qtotlist, ddof = 1)
                    meanQtot[it,iz,i_n,ico,jco] = np.mean(tmp_Qtotlist)
            
            Mmem1[it,iz,i_n,:nx,:ny] = Mmem_coarse[0]
            if n == 4:
                Mmem_mean = np.mean(Mmem_coarse, axis = 0)[:nx,:ny]
                Mmem_south = np.mean(Mmem_coarse, axis = 0)[:nx/2,:ny]
                Mtot[it,iz] = np.sum(Mmem_mean)
                Msouth[it,iz] = np.sum(Mmem_south)
                Mnorth[it,iz] = np.sum(Mmem_mean) - np.sum(Mmem_south)
                
            
            ## Calculate 2dACF 
            #if n < nlist[1]:
                #tmp_acflist = []
                #for imem in range(nmem):
                    #Mdiff = ((Mmem_coarse[imem] - np.mean(Mmem_coarse, axis = 0))/
                            #1)
                            ##np.mean(Mmem_coarse, axis = 0))
                    #Mdiff[np.isnan(Mdiff)] = 0.
                    #C = crosscor(Mdiff, Mdiff, minusmean = False)
                    #tmp_acflist.append(C)
                ##print Mdiff
                #Cmean =  np.nanmean(tmp_acflist, axis = 0)
                #tmp_acf2d = radial_profile(Cmean, (nx/2,ny/2))
                #print Mdiff[nx/2-2:nx/2+2,nx/2-2:nx/2+2]
                #print Cmean[nx/2-2:nx/2+2,nx/2-2:nx/2+2]
                #print tmp_acf2d
                #acf2d[it,iz,i_n,:tmp_acf2d.shape[0]] = tmp_acf2d
            #else:
                #acf2d[it,iz,i_n,0] = np.nan
            
                    
  
                
            
            # End coarse upscaled variances and means
            ####################################################################
            
# Close ncdf file
rootgrp.close()
            

        
        
        
        
