# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 09:52:26 2016

Run SubcriticalModel on local computer.

Version for Stephan.

f2py -c -m CloudPercolation_cwr_stephan CloudPercolation_cwr_stephan.f95

"""

import matplotlib
matplotlib.use("TKAgg")
import numpy as np
import CloudPercolation_cwr_stephan
import matplotlib.pyplot as plt
#import ipdb

print CloudPercolation_cwr_stephan.m_percolations_float_stephan.perc.__doc__

#ipdb.set_trace()

#f2py -c -m CloudPercolation_float CloudPercolation_float.f95

def CloudPercolation_cw_i(phys_L,DL,N_clouds,r_co_km,r_m,fac_cw,cs, plot_field=False):

    #Underlying distribution: exponential disc size with shifted mean

    area, area_co=r_m*r_m*np.pi, 0.#r_co_km*r_co_km*np.pi           

    effective_mean_area=area-area_co                                            #adapt mean to account for cut-off radius    
    percentage_small=(1-np.exp(-area_co/effective_mean_area))                   #fraction of clouds expected with r<r_co
    estimate_fraction=1-np.minimum(2.0*percentage_small,0.99)
    
    a=np.random.exponential(effective_mean_area,int(2*N_clouds/estimate_fraction)) #draw massflux from exponential distribution     
    
    print '1st beta =', np.var(a, ddof = 1)/(np.mean(a)**2)
    
    rA_km=np.sqrt(a/np.pi)                                                         #radius corresponding to mf [km]
    rA_km=rA_km[np.where(rA_km>r_co_km)][:N_clouds]                                #choose N clouds larger than cut-off radius and sort according to size
    rA_km=np.sort(rA_km)[::-1]  
    
    a2 = np.pi*rA_km**2
    print '2nd beta =', np.var(a2, ddof = 1)/(np.mean(a2)**2)
    
    rA,rA_ring=rA_km/DL,(rA_km*fac_cw)/DL                                         #radii in units of grid cells
    
    L=int(phys_L/DL)                                                            #number of grid-cells in one direction
    
    cloud_field_binary, cloud_centers, cloud_field =CloudPercolation_cwr_stephan.m_percolations_float_stephan.perc(cs,L,len(rA),rA,rA_ring)

    if plot_field:
        plt.imshow(cloud_field, cmap='jet', interpolation='nearest')
        plt.colorbar()
        plt.show()
    
    cloud_centers=np.array(cloud_centers)-1                                     #Account for Fortran arrays starting with 1
    

    return cloud_field, cloud_centers, rA_km


#simulation set-up

N_ensemble = 2                                                                 #Number of calculated fields

#domain set-up 

phys_L_km=256*2.8                                                               #Domain size in km 
DL=2.8                                                                          #Resolution in km 
L=int(phys_L_km/DL)                                                             #Number of grid-cells in one direction

r_co_km=2.0                                                                     #Minimal radius of discs in km

#model properties

#coverage_fraction=0.05                                                          
N_clouds= 200                                                                   #Number of clouds

#Cloud properties
r_m = 2.8                                                                       #Mean radius of discs in km
cs  = 20                                                                         #Factor by which probability is increased within the rings around the clouds
fac_cw = 3                                                                      #Prob. increased within the ring r_disc < r < r_disc*fac_cw
       
estimated_coverage_fraction= 1.0-np.exp(-N_clouds*(r_m*r_m*np.pi)/(phys_L_km*phys_L_km))  
print('estimated coverage fraction %.3f'%estimated_coverage_fraction)                                                                                                                                                       

cloud_field_stack=np.zeros((N_ensemble, L, L))
cloud_centers_array=[]
radii_array=[]

for i in range(N_ensemble):
    
    cloud_field, cloud_centers, rA_km =CloudPercolation_cw_i(phys_L_km,DL,N_clouds,r_co_km,r_m,fac_cw,cs,plot_field=False) 
    
    cloud_field_stack[i,:,:]= cloud_field
    cloud_centers_array.append(cloud_centers) 
    radii_array.append(rA_km)

#ipdb.set_trace()

#for i in range(N_ensemble):
    #cloud_field = cloud_field_stack[i,:,:]
    #plt.imshow(cloud_field, cmap='jet', interpolation = 'nearest', origin = 'lower')#, interpolation='nearest')
    #plt.show()
