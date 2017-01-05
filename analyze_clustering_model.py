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
parser.add_argument('--nens', metavar = 'nens', type=int, default = 4)
parser.add_argument('--water', metavar = 'water', type=str, default = 'True')
parser.add_argument('--str', metavar = 'str', type=str, default = '')
parser.add_argument('--cs', metavar = 'cs', type=float)
parser.add_argument('--cw', metavar = 'cw', type=float)
parser.add_argument('--N', metavar = 'N', type=int)
args = parser.parse_args()
# Convert water to bool 
if args.water == 'True':
    args.water = True
elif args.water == 'False':
    args.water = False
else:
    raise Exception

plotdir = '/home/s/S.Rasp/Dropbox/figures/PhD/variance/perc_model/' + args.str


# 1. Create the hypothetical ensemble
sx = 256
nlist = [256, 128, 64, 32, 16, 8, 4]

#domain set-up 

phys_L_km=sx*2.8  #Domain size in km 
DL=2.8  #Resolution in km 
L=int(phys_L_km/DL)  #Number of grid-cells in one direction

r_co_km=0.   #Minimal radius of discs in km

#model properties

#coverage_fraction=0.05                                                          
N_clouds= args.N  #Number of clouds

#Cloud properties
r_m = 2.4  #Mean radius of discs in km
cs  = args.cs  #Factor by which probability is increased within the rings around the clouds
fac_cw = args.cw #Prob. increased within the ring r_disc < r < r_disc*fac_cw
 

#ensdir = '/home/scratch/users/stephan.rasp/hypo_' + args.type + '/deout_ceu_pspens/'




# Draw one big field
Nx = args.nens
Nens = Nx*Nx
phys_L_km *= Nx
L *= Nx
print N_clouds * Nens
cloud_field_big, cloud_centers, rA_km =CloudPercolation_cw_i(phys_L_km,
                                                             DL,N_clouds*Nens,
                                                             r_co_km,r_m,fac_cw,
                                                             cs,
                                                             plot_field=False)


# Loop over ensemble members
sizelist = []
comlist = []
glist = []
glist_ns = []
numlist = []
for i in range(Nx):
    for j in range(Nx):
        xmin = i*sx
        xmax = (i+1)*sx
        ymin = j*sx
        ymax = (j+1)*sx
        cloud_field = cloud_field_big[xmin:xmax,ymin:ymax]
        
        #plt.imshow(cloud_field, cmap='jet', interpolation='nearest')
        #plt.colorbar()
        #plt.show()
        
        tmp = identify_clouds(cloud_field, thresh = 0, water = args.water)
        
        labels, cld_size_mem, cld_sum_mem = tmp
        num = np.unique(labels).shape[0]   # Number of clouds
        com = np.array(center_of_mass(cloud_field, labels, range(1,num)))
        numlist.append(num)
        g, r = calc_rdf(labels, cloud_field, normalize = True, rmax = 36, 
                                dr = 1)
        glist.append(g)
        g, r = calc_rdf(labels, cloud_field, normalize = False, rmax = 36, 
                                dr = 1)
        glist_ns.append(g)
        comlist.append(com)
        sizelist.append(cld_size_mem)

print 'Actual N mean', np.mean(numlist)
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
print 'N clouds', len(totlist)/Nens
print 'mean cloud size =', (np.mean(totlist)/1e6), 'km^2'
print 'mean cloud radius =', np.sqrt((np.mean(totlist)/1e6)/np.pi), 'km'
print 'var cld size = ', np.var(totlist, ddof = 1)
print 'beta =', np.var(totlist, ddof = 1)/(np.mean(totlist)**2)
fig, ax = plt.subplots(1,1)
dx2 = 2.8e3**2
ax.hist(totlist, bins = 30, range = [0,dx2*30])
ax.set_yscale('log')
ax.set_xlabel('Cloud size [m^2]')
ax.set_ylabel('Number')
plt.savefig(plotdir + 'hist')
plt.close('all')


# RDF
fig, axarr = plt.subplots(1,2)
meang = np.mean(glist, axis = 0)
meang_ns = np.mean(glist_ns, axis = 0)
axarr[0].plot(r/1.e3, meang)
axarr[0].set_xlabel('Distance [km]')
axarr[0].set_ylabel('Normalized RDF')
axarr[1].plot(r/1.e3, meang_ns)
axarr[1].set_xlabel('Distance [km]')
axarr[1].set_ylabel('Non-Normalized RDF')
plt.tight_layout()
plt.savefig(plotdir + 'rdf')
plt.close('all')


# Std_v_mean and parameters
clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
         "#ff00ff")
alphalist = []
betalist = []
fraclist = []
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
    alphalist.append(np.nanmean(alpha))
    print 'alpha =', np.nanmean(alpha)
    beta = varm/(meanm**2) 
    betalist.append(np.nanmean(beta))
    print 'beta =', np.nanmean(beta)
    #predict = (alpha + beta) * meanm * meanM
    predict = 2 * meanm * meanM
    frac = varM/predict 
    fraclist.append(np.nanmean(frac))
    print 'frac =', np.nanmean(frac)

fig, axarr = plt.subplots(1, 3, figsize = (12, 4))

axarr[0].plot(range(len(nlist)), alphalist, linewidth = 2)
axarr[0].set_ylabel('alpha')
axarr[0].set_xlabel('scale [km]')
axarr[0].set_xticks(range(len(nlist)))
axarr[0].set_xticklabels(np.array(nlist) * 2.8)
axarr[0].invert_xaxis()
#axarr[0].set_yscale('log')
axarr[0].set_ylim(0.1, 5)
axarr[0].plot([0, len(nlist)-1],[1,1], c = 'gray', zorder = 0.5)

axarr[1].plot(range(len(nlist)), betalist, linewidth = 2) 
axarr[1].set_ylabel('beta')
axarr[1].set_xlabel('scale [km]')
axarr[1].set_xticks(range(len(nlist)))
axarr[1].set_xticklabels(np.array(nlist) * 2.8)
axarr[1].invert_xaxis()
axarr[1].set_ylim(0.1, 2)
axarr[1].plot([0, len(nlist)-1],[1,1], c = 'gray', zorder = 0.5)


axarr[2].plot(range(len(nlist)), fraclist, linewidth = 2) 
axarr[2].set_ylabel('variance ratio')
axarr[2].set_xlabel('scale [km]')
axarr[2].set_xticks(range(len(nlist)))
axarr[2].set_xticklabels(np.array(nlist) * 2.8)
axarr[2].invert_xaxis()
axarr[2].set_ylim(0.1, 5)
axarr[2].plot([0, len(nlist)-1],[1,1], c = 'gray', zorder = 0.5)


print 'here'


plt.tight_layout()
fig.savefig(plotdir + 'prediction')








    
    
