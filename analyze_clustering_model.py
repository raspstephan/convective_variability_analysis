"""
This script create a hypothetical ensemble which can be analyzed with the 
compute script.

"""

# Imports
import os
import argparse
import numpy as np
from SubcriticalModel_Stephan import CloudPercolation_cw_i
from cosmo_utils.diag import identify_clouds, calc_rdf
from scipy.ndimage.measurements import center_of_mass
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

# Setup
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--nens', metavar = 'nens', type=int, default = 20)
parser.add_argument('--water', metavar = 'water', type=str, default = 'True')
args = parser.parse_args()
# Convert water to bool 
if args.water == 'True':
    args.water = True
elif args.water == 'False':
    args.water = False
else:
    raise Exception

# 1. Create the hypothetical ensemble
sx = 256
nlist = [256, 128, 64, 32, 16, 8, 4]

#domain set-up 

phys_L_km=sx*2.8  #Domain size in km 
DL=2.8  #Resolution in km 
L=int(phys_L_km/DL)  #Number of grid-cells in one direction

r_co_km=0.5   #Minimal radius of discs in km

#model properties

#coverage_fraction=0.05                                                          
N_clouds= 585  #Number of clouds

#Cloud properties
r_m = 2.5  #Mean radius of discs in km
cs  = 0   #Factor by which probability is increased within the rings around the clouds
fac_cw = 4 #Prob. increased within the ring r_disc < r < r_disc*fac_cw
 

#ensdir = '/home/scratch/users/stephan.rasp/hypo_' + args.type + '/deout_ceu_pspens/'



# Loop over ensemble members
sizelist = []
comlist = []
glist = []
for imem in range(1,args.nens+1):
    print 'Member:', imem
    
    # Draw actual N cloud from Poisson distribution
    N_rand = np.random.poisson(N_clouds)
    
    cloud_field, cloud_centers, rA_km =CloudPercolation_cw_i(phys_L_km,
                                                             DL,N_rand,
                                                             r_co_km,r_m,fac_cw,
                                                             cs,plot_field=False) 
    # cloud field is the cloud field
    # cloud centers is a N_cld (-1???) by (x,y) array with the coordinates of the centers
    # rA_km is ???
    
    tmp = identify_clouds(cloud_field, thresh = 0, water = args.water)
    
    labels, cld_size_mem, cld_sum_mem = tmp
    num = np.unique(labels).shape[0]   # Number of clouds
    com = np.array(center_of_mass(cloud_field, labels, range(1,num)))
    
    g, r = calc_rdf(labels, cloud_field, normalize = True, rmax = 36, 
                            dr = 1)
    glist.append(g)
    comlist.append(com)
    sizelist.append(cld_size_mem)

# Loop over n
varMlist = []
varmlist = []
varNlist = []
meanMlist = []
meanmlist = []
meanNlist = []
for i_n, n in enumerate(nlist):
    print 'n: ', n
    ####################################################################
    # Calculate coarse variances and means
    # Determine size of coarse arrays
    nx = int(np.floor(sx/n))
    ny = int(np.floor(sx/n))
    
    varMlist.append(np.empty((nx,ny)))
    varmlist.append(np.empty((nx,ny)))
    varNlist.append(np.empty((nx,ny)))
    meanMlist.append(np.empty((nx,ny)))
    meanmlist.append(np.empty((nx,ny)))
    meanNlist.append(np.empty((nx,ny)))

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
            for com, cld_sum_mem in zip(comlist, sizelist):
                                            
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

                # End member loop #############
            
            # Now convert the list with all clouds for this box
            tmp_cldlist = np.array(tmp_cldlist)
            
            # Calculate statistics and save them in ncdf file
            
            varMlist[-1][ico,jco] = np.var(tmp_Mlist, ddof = 1)
            varNlist[-1][ico,jco] = np.var(tmp_Nlist, ddof = 1)
            varmlist[-1][ico,jco] = np.var(tmp_cldlist, ddof = 1)
            meanMlist[-1][ico,jco] = np.mean(tmp_Mlist)
            meanmlist[-1][ico,jco] = np.mean(tmp_cldlist)
            meanNlist[-1][ico,jco] = np.mean(tmp_Nlist)

            # This means NaNs only appear when minmem criterion is not met

# Analzing and plotting

# Cloud hist
totlist = []
for l in sizelist:
    totlist += list(l)
print 'mean cloud size =', (np.mean(totlist)/1e6), 'km^2'
print 'mean cloud radius =', np.sqrt((np.mean(totlist)/1e6)/np.pi), 'km'
print 'var cld size = ', np.var(totlist, ddof = 1)
print 'beta =', np.var(totlist, ddof = 1)/(np.mean(totlist)**2)
fig, ax = plt.subplots(1,1)
ax.hist(totlist)
ax.set_yscale('log')
plt.show()
plt.close('all')


# RDF
fig, ax = plt.subplots(1,1)
meang = np.mean(glist, axis = 0)
ax.plot(r/1.e3, g)
plt.show()



# Std_v_mean and parameters
clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
         "#ff00ff")
#fig, axarr = plt.subplots(1,2)
for i_n, n in enumerate(nlist):
    print n
    varM = np.ravel(varMlist[i_n])
    meanM = np.ravel(meanMlist[i_n])
    varm = np.ravel(varmlist[i_n])
    meanm = np.ravel(meanmlist[i_n])
    varN = np.ravel(varNlist[i_n])
    meanN = np.ravel(meanNlist[i_n])
    
    print 'mean m =', np.nanmean(meanm)
    print 'var m =', np.nanmean(varm)
    
    alpha = varN/meanN 
    print 'alpha =', np.nanmean(alpha)
    beta = varm/(meanm**2) 
    print 'beta =', np.nanmean(beta)
    predict = 2 * meanm * meanM
    frac = varM/predict 
    print 'frac =', np.nanmean(frac)
    
    
    #axarr[0].scatter(x, np.sqrt(y), c = clist[i_n])
    
    #predict = 2 * np.ravel(meanMlist[i_n]) * np.ravel(meanmlist[i_n])
    #frac = y/predict
    
    #axarr[1].scatter(x, frac, c = clist[i_n])
    
#axarr[0].set_xscale('log')
#axarr[0].set_yscale('log')
#axarr[0].set_xlim(1e6, 1e11)
#axarr[0].set_ylim(4e6, 2e9)

#axarr[1].set_xscale('log')
#axarr[1].set_yscale('log')

##plt.show()
#plt.close('all')








    
    
