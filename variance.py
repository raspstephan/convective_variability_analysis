"""
This is a script to analyze the convetive variance from 
COSMO ensembles
"""
# Imports
import sys
from cosmo_utils.pyncdf import getfobj_ncdf_ens
from cosmo_utils.diag import mean_spread_fieldobjlist, identify_clouds, rdf
from cosmo_utils.plot import fig_contourf_1sp, ax_contourf
from cosmo_utils.helpers import yyyymmddhh_strtotime, make_timelist, ddhhmmss, yyyymmddhh
from datetime import timedelta
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.ndimage import measurements
import copy
from scipy.optimize import leastsq
np.seterr(divide='ignore', invalid='ignore')   # Suppress invalide divide warnings

def origin_residual(p, x, y):
    slope = p
    err = y - (slope * x)
    return err


# Plot settings
mpl.rcParams['font.size'] = 10
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['legend.fontsize'] = 8


cmTAUC = ("#FFFFC1","#FDFEB3","#FCF6A5","#FCEF97","#FCE88A","#FCE07E",
          "#FCD873","#FCD068","#FCC860","#FCC058","#FCB853","#FBB04F",
          "#FAA74E","#F89F4E","#F7964F","#F58E52","#F28655","#F07D59",
          "#ED745D","#E96C61","#E66365","#E25A6A","#DE516D","#D94771",
          "#D53D74","#D03178","#CB247B","#C6117D","#C10080","#BC0083",
         )
levelsTAUC = np.arange(0, 31, 1)

cmNVAR = ("#0030C4","#3F5BB6","#7380C0","#A0A7CE","#CCCEDC","#DDCACD",
          "#D09AA4","#BF6A7D","#AA3656","#920031")

levelsNVAR = [1, 1.14, 1.33, 1.6, 1.78, 2, 2.25, 2.5, 3, 3.5, 4]

cmMP = ("#DBA9B4","#DDA6A0","#D8A786","#CAAA68","#B1AF4C","#8DB340",
          "#50B750","#00B76E","#00B28C","#007678")
levelsMP = np.linspace(0,1+1/11,11)



# Setup
water = False
collapse = True
ana = sys.argv[1]  # 'm' or 'p'
date = sys.argv[2]
ensdir = '/home/cosmo/stephan.rasp/' + date + '/deout_onlypsp/'
try:
    nens = int(sys.argv[3])
except:    
    nens = 20
tstart = timedelta(hours=1)
tend = timedelta(hours = 24)
tinc = timedelta(hours = 1)
lx1 = 204/2 # ATTENTION first dimension is actually y
lx2 = -(lx1+1) # Number of grid pts to exclude at border
ly1 = 164/2
ly2 = -(ly1+1)
plotdir = '/home/s/S.Rasp/Dropbox/figures/PhD/variance/' + date
if not collapse:
    plotdir += '/noncollapse/'
    date = 'noncoll_' + date
if not nens == 20:
    plotdir += '/n' + str(nens).zfill(2) + '/'
    date = 'n' + str(nens).zfill(2) + '_' + date
if not water:
    plotdir += '/nowater/'
    date = 'nowater_' + date
dx = 2800.


nlist = [256, 128, 64, 32, 16, 8, 4]

# Specific setup for type of analysis
if ana == 'm':
    fieldn = 'W'
    thresh = 1.
    sufx = 'z.nc_1h'
    plotdir += '/m/'
    ana_unit = 'kg/s'
    cmWP = ("#7C0607","#903334","#A45657","#BA7B7C","#FFFFFF",
        "#8688BA","#6567AA","#46499F","#1F28A2")
    levelsWP = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
    sizemax = 1.e8
    summax = 5e8

if ana == 'p':
    fieldn = 'TOT_PR'
    thresh = 0.001
    sufx = '.nc_5m'
    plotdir += '/p/'
    ana_unit = 'mm/h'
    cmWP = ((1    , 1     , 1    ), 
            (0    , 0.627 , 1    ),
            (0.137, 0.235 , 0.98 ),
            (0.392, 0     , 0.627),
            (0.784, 0     , 0.627),
            (1    , 0.3   , 0.9  ) )
            #(0.1  , 0.1   , 0.784),
    levelsWP = [0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
    sizemax = 8.e8
    summax = 0.35
    
# Create plotdir if not exist
if not os.path.exists(plotdir): os.makedirs(plotdir)


# Make the timelist
timelist = make_timelist(tstart, tend, tinc)
# Time loop 
# Initialize time list
r_cluster_list = []
sizemmean_list = []
summean_list = []
total_list = []
dimeantauc_list = []
for t in timelist:
    print t
    
    ncdffn = 'lfff' + ddhhmmss(t) + sufx
    
    # Load ensembles 
    fobjlist = getfobj_ncdf_ens(ensdir, 'sub', nens, ncdffn, 
                                dir_suffix='/OUTPUT/',
                                fieldn = fieldn, nfill=2)
    if ana == 'm':   # Additional positive QC filter
        qcobjlist = getfobj_ncdf_ens(ensdir, 'sub', nens, ncdffn, 
                                     dir_suffix='/OUTPUT/',
                                     fieldn = 'QC', nfill=2)
    else:
        qcobjlist = [None]*len(fobjlist)
    
    
    # Calculate what needs to be calculated
    # 1. Tau_c
    taucfn = 'lfff' + ddhhmmss(t) + '.nc_5m'
    tauclist = getfobj_ncdf_ens(ensdir, 'sub', nens, taucfn, 
                                dir_suffix='/OUTPUT/',
                                fieldn = 'TAU_C', nfill=2)
    # Get mean tauc object
    meantauc, tmp = mean_spread_fieldobjlist(tauclist)
    dimeantauc = np.nanmean(meantauc.data[lx1:lx2, ly1:ly2])
    dimeantauc_list.append(dimeantauc)
    
    # Loop over members 
    sizelist = []
    sumlist = []
    glist = []
    
    # For n loop
    comlist = []
    labelslist = []
    for fobj, qcobj in zip(fobjlist, qcobjlist):
        # 2. Cloud size and m/p distribution
        if ana == 'm':
            field = fobj.data[0,lx1:lx2, ly1:ly2]
            labels, cld_size, cld_sum = identify_clouds(field,
                                                        thresh , 
                                                        qcobj.data[0,lx1:lx2, ly1:ly2],
                                                        opt_thresh = 0.,
                                                        water = True)
            cld_sum *= dx*dx*0.9575  # Here 1 stands for rho ATTENTION
        else:
            field = fobj.data[lx1:lx2, ly1:ly2]
            labels, cld_size, cld_sum = identify_clouds(field,
                                                        thresh = thresh)
        sizelist.append(cld_size)
        sumlist.append(cld_sum) 
                
        # 3. Calculate RDF
        g, r = rdf(labels, field, normalize = True)
        glist.append(g)
        
        labelslist.append(labels)
        
        # 4. Calculate Variance (Save for n loop below)
        # TODO: Put this into a function later on
        num = np.unique(labels).shape[0]   # Number of clouds
        # Get center of mass for each cluster
        com = np.array(measurements.center_of_mass(field, labels, range(1,num)))
        if com.shape[0] == 0:   # Accout for empty arrays
            com = np.empty((0,2))
        comlist.append(com)
        sx, sy = field.shape
        
        
    # cont 2. Calculate histograms
    sizelist_flat = [i for sl in sizelist for i in sl]
    sumlist_flat = [i for sl in sumlist for i in sl]
    sizehist, sizeedges = np.histogram(sizelist_flat, 
                                       bins = 10, range = [0., sizemax])
    sumhist, sumedges = np.histogram(sumlist_flat, 
                                     bins = 20, range = [0., summax])
    sizemean = np.mean(sizelist_flat)
    summean = np.mean(sumlist_flat)
    total = np.sum(sumlist_flat)
    sizemmean_list.append(sizemean)
    summean_list.append(summean)
    total_list.append(total)
    
    
    # cont 3. Get means after member loop
    #print glist
    g = np.mean(glist, axis = 0)
    # Get clustering radius
    gthresh = 1.1
    if not np.isnan(np.mean(g)):
        # This is the index where g drops below one after the first peak
        tmpwhere = np.where(g < gthresh)[0]
        gind = tmpwhere[np.where(np.diff(tmpwhere) > 1)[0][0]+1]
        r_cluster = r[gind]/1000.   # In km
    else:
        r_cluster = np.nan
    r_cluster_list.append(r_cluster)
    
    
    # cont 4. Variance
    # Loop over n
    varres1list = []
    varres2list = []
    
    for n in nlist:
        print 'n', n
        # Determine size of coarse arrays
        nx = int(np.floor(sx/n))
        ny = int(np.floor(sy/n))
        
        # Allocate ens lists
        mplist = []
        MPlist = []
        Nlist = []
        
        # Loop over members
        for fobj, labels, com, cld_sum in zip(fobjlist, labelslist, comlist, 
                                              sumlist):
            if ana == 'm':
                field = fobj.data[0,lx1:lx2, ly1:ly2]
            else: 
                field = fobj.data[lx1:lx2, ly1:ly2]
        
            # Allocate array for saving
            mp_field = np.empty((nx, ny))
            MP_field = np.empty((nx, ny))
            N_field = np.empty((nx, ny))
            # Loop over "coarse" array
            for i in range(nx):
                for j in range(ny):
                    # Get limits for each N box
                    xmin = i*n
                    xmax = (i+1)*n
                    ymin = j*n
                    ymax = (j+1)*n
                    
                    if collapse:
                        # 1. The collapsed version
                        # Create bool_arr
                        bool_arr = ((com[:,0]>=xmin)&(com[:,0]<xmax)&
                                    (com[:,1]>=ymin)&(com[:,1]<ymax))
                        # Get cld_size for subdomain
                        sub_cld_sum = cld_sum[bool_arr]
                        
                    
                    else:
                        # 2. The "normal version"
                        subfield = field[i*n:(i+1)*n, j*n:(j+1)*n]
                        sublabels = labels[i*n:(i+1)*n, j*n:(j+1)*n]
                        lrange = np.unique(sublabels)
                        sub_cld_sum = measurements.sum(subfield, sublabels,
                                                       lrange[1:])*dx*dx*0.9575
                        
                    # Get important values
                    mp_field[i, j] = np.mean(sub_cld_sum)
                    MP_field[i, j] = np.sum(sub_cld_sum)
                    N_field[i, j] = sub_cld_sum.shape[0]
                        
            
            # Write into ensemble list
            mplist.append(mp_field)
            MPlist.append(MP_field)
            Nlist.append(N_field)
            
        # Calculate ensemble means and variances, fields
        var = np.var(MPlist, axis = 0, ddof = 1)
        MPmean = np.mean(MPlist, axis = 0)
        mpmean = np.nanmean(mplist, axis = 0)
        Nmean = np.mean(Nlist, axis = 0)
        
        # Plot these fields now
        # First, have to upscale them again to the model grid
        var_fullfield = np.ones((fobj.ny, fobj.nx)) * np.nan
        MP_fullfield = np.ones((fobj.ny, fobj.nx)) * np.nan
        mp_fullfield = np.ones((fobj.ny, fobj.nx)) * np.nan
        for i in range(nx):
            for j in range(ny):
                # Get limits for each N box
                xmin = i*n+lx1
                xmax = (i+1)*n+lx1
                ymin = j*n+ly1
                ymax = (j+1)*n+ly1
                
                var_fullfield[xmin:xmax, ymin:ymax] = var[i,j]
                MP_fullfield[xmin:xmax, ymin:ymax] = MPmean[i,j]
                mp_fullfield[xmin:xmax, ymin:ymax] = mpmean[i,j]
        
        # Create new fobj
        var_fobj = copy.deepcopy(fobj)
        var_fobj.data = var_fullfield
        var_fobj.dims = 2
        var_fobj.fieldn = 'Var'
        var_fobj.unit = '(' + ana_unit + ')^2'
        
        MP_fobj = copy.deepcopy(fobj)
        MP_fobj.data = MP_fullfield
        MP_fobj.data[MP_fobj.data == 0] = np.nan
        MP_fobj.dims = 2
        MP_fobj.unit = ana_unit
        
        mp_fobj = copy.deepcopy(fobj)
        mp_fobj.data = mp_fullfield
        mp_fobj.data[mp_fobj.data == 0] = np.nan
        mp_fobj.dims = 2
        mp_fobj.unit = ana_unit
        
        nvar_fobj = copy.deepcopy(fobj)
        nvar_fobj.data = var_fullfield /MP_fullfield/mp_fullfield
        nvar_fobj.dims = 2
        nvar_fobj.unit = ''
        
        if ana == 'm':
            MP_fobj.fieldn = 'M'
            mp_fobj.fieldn = 'm'
            nvar_fobj.fieldn = 'Normalized Var = Var / M / m'
            mpmax = 2e8
        else:
            MP_fobj.fieldn = 'P'
            mp_fobj.fieldn = 'p'
            nvar_fobj.fieldn = 'Normalized Var = Var / P / p'
        
        wpobj = copy.deepcopy(fobj)
        if wpobj.dims == 3:
            wpobj.data = fobj.data[0]
            wpobj.dims = 2
        
        
        # Plotting
        fobjlist_plot = [wpobj, MP_fobj, nvar_fobj, mp_fobj]
        cmlist = [(cmWP, levelsWP), (cmMP, levelsMP*np.nanmax(MP_fullfield)), 
                  (cmNVAR, levelsNVAR), (cmMP, levelsMP*np.nanmax(mp_fullfield))]
        
        fig, axarr = plt.subplots(2, 2, figsize = (9, 8))
        for ax, fob, cm, i in zip(list(np.ravel(axarr)), fobjlist_plot, cmlist, 
                                                range(4)):
            plt.sca(ax)   # This is necessary for some reason...
            cf, tmp = ax_contourf(ax, fob, colors=cm[0], 
                                pllevels=cm[1], sp_title=fob.fieldn,
                                Basemap_drawrivers = False,
                                npars = 0, nmers = 0, ji0=(50, 50),
                                ji1=(411, 381), extend = 'both')
            cb = fig.colorbar(cf)
            cb.set_label(fob.unit)
        figtitle = date + '+' + ddhhmmss(t) + ' n = ' + str(n).zfill(3)
        fig.suptitle(figtitle, fontsize='x-large')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # To account for suptitle
        plotdirnew = plotdir + '/var/'
        if not os.path.exists(plotdirnew): os.makedirs(plotdirnew)
        fig.savefig(plotdirnew + 'var_' + ddhhmmss(t) + '_' + str(n).zfill(3),
                    dpi = 300)
        plt.close('all')
        
        # Now actually calculate the variance constant/factor
        # Method 1: Simply take the mean of nvar 
        varres1list.append([var/(MPmean**2), Nmean])
        
        # Method 2: Do the slope approach
        
        #xfit = np.ravel(MPmean**2)
        #yfit = np.ravel(var)
        #mask = yfit > 0

        
        #slope = leastsq(origin_residual, 1, args = (xfit[mask], yfit[mask]))[0][0]
        #varres2list.append([np.sqrt(slope), np.sqrt(1/np.nanmean(Nmean))])
            
    varres1list = np.array(varres1list)
    #varres2list = np.array(varres2list)
    

    
    # Plot what needs to be plotted every time step
    # 1. Ensemble mean tau_c over Germany
    title_sufx = 'mean tau_c, di_mean: ' + str(dimeantauc) + 'h + ' + ddhhmmss(t)
    fig = fig_contourf_1sp(meantauc, pllevels = np.arange(0, 21, 1),
                           extend = 'max', sp_title = title_sufx)
    plt.tight_layout()
    plotdirnew = plotdir + '/tauc/'
    if not os.path.exists(plotdirnew): os.makedirs(plotdirnew)
    fig.savefig(plotdirnew + 'tauc_' + ddhhmmss(t), dpi = 300)
    
    # 2. Plot cloud size and m/p distribution, plus RDF
    fig, axarr = plt.subplots(1, 3, figsize = (95./25.4*2.5, 3.2))
    axarr[0].bar(sizeedges[:-1], sizehist, width = np.diff(sizeedges)[0])
    axarr[0].plot([sizemean, sizemean], [0.1, 1e4], c = 'red', 
                  alpha = 0.5)
    axarr[0].set_xlabel('Cloud size [m^2]')
    axarr[0].set_ylabel('Number of clouds')
    axarr[0].set_title('Cloud size distribution')
    axarr[0].set_xlim([0., sizemax])
    axarr[0].set_ylim([0.1, 1e4])
    axarr[0].set_yscale('log')
    
    axarr[1].bar(sumedges[:-1], sumhist, width = np.diff(sumedges)[0])
    axarr[1].plot([summean, summean], [0.1, 1e4], c = 'red', 
                  alpha = 0.5)
    axarr[1].set_ylabel('Number of clouds')
    axarr[1].set_xlim([0., summax])
    axarr[1].set_ylim([0.1, 1e4])
    axarr[1].set_yscale('log')
    
    # 3. Plot RDF
    axarr[2].plot(r/1000., g)
    axarr[2].plot([r_cluster, r_cluster], [0, 4], c = 'red', alpha = 0.5)
    axarr[2].plot([0, np.max(r)/1000.], [1, 1], c = 'gray', alpha = 0.5)
    axarr[2].set_xlabel('Distance [km]')
    axarr[2].set_ylabel('Normalized RDF')
    axarr[2].set_title('Radial distribution function')
    axarr[2].set_ylim(0, 4)
    axarr[2].set_xlim(0, np.max(r)/1000.)
    
    # Analysis specific stuff
    if ana == 'm':
        axarr[1].set_xlabel('Cloud mass flux [kg/s]')
        axarr[1].set_title('Cloud mass flux distribution')
    
    else: 
        axarr[1].set_xlabel('Cloud precipitation [mm/h]')
        axarr[1].set_title('Cloud precipitation distribution')
        
        
    plt.tight_layout()
    plotdirnew = plotdir + '/cloud_stats/'
    if not os.path.exists(plotdirnew): os.makedirs(plotdirnew)
    fig.savefig(plotdirnew + 'stats_' + ddhhmmss(t), dpi = 300)
    
    
    # 4. Plot CC06 Fig4 
    fig, axarr = plt.subplots(1, 2, figsize = (95./25.4*2, 3.2))
    clist = ("#ff0000", "#ff8000", "#ffff00","#40ff00","#00ffff","#0040ff","#ff00ff")

    for x, y, n, c in zip(list(varres1list[:,1]),
                          list(varres1list[:,0]),
                          nlist, clist):
        x[x==0] = np.nan   # Set no clouds to nan for correct mean

        xcloud = np.sqrt(1/x)
        ycloud = np.sqrt(y)

        #ax.scatter(x, y, marker = 'D', c = c, label = str(n*2.8)+'km', s = 20,
                   #linewidth = 0.2)
        axarr[0].scatter(xcloud, ycloud, marker = 'o', c = c, 
                   s = 4, zorder = 0.2, linewidth = 0, alpha = 0.8,
                   label = str(n*2.8)+'km')

        ypercent = y/(2./x)*100.
        axarr[1].scatter(xcloud, ypercent, marker = 'o', c = c, 
                   s = 4, zorder = 0.2, linewidth = 0, alpha = 0.8)
    
    #for x, y, n, c in zip(list(varres2list[:,1]),
                          #list(varres2list[:,0]),
                          #nlist, clist): 
        #ax.scatter(x, y, marker = '*', c = c, label = str(n*2.8)+'km', s = 20,
                   #linewidth = 0.2)
    
    axarr[0].legend(loc =4, ncol = 2, prop={'size':6})
    tmp = np.array([0,5])
    axarr[0].plot(tmp,tmp*np.sqrt(2), c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[0].set_xlim(0,4)
    axarr[0].set_ylim(0,5)
    axarr[0].set_xlabel('Square root (1/N)')
    if ana == 'm':
        axarr[0].set_ylabel('Square root (Var(M)/M^2)')
    else:
        axarr[0].set_ylabel('Square root (Var(P)/P^2)')
    axarr[1].set_yscale('log')
    axarr[1].plot([0.,4],[100,100], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[1].set_xlim(0,4)
    axarr[1].set_ylim(1, 1000)
    axarr[1].set_xlabel('Square root (1/N)')
    axarr[1].set_ylabel('Percent of theoretical value')
    fig.suptitle(date + '+' + ddhhmmss(t), fontsize='x-large')
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    plotdirnew = plotdir + '/var_scatter/'
    if not os.path.exists(plotdirnew): os.makedirs(plotdirnew)
    fig.savefig(plotdirnew + 'scatter_' + ddhhmmss(t), dpi = 300)
    plt.close('all')
    
    
    
# Plot summary plots
timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]

fig, axarr = plt.subplots(2, 3, figsize = (95./25.4*3, 6.))

axarr[0,0].plot(timelist_plot, total_list)
axarr[0,0].set_xlabel('time [h/UTC]')
axarr[0,0].set_xlim(timelist_plot[0], timelist_plot[-1])

axarr[0,1].plot(timelist_plot, sizemmean_list)
axarr[0,1].set_xlabel('time [h/UTC]')
axarr[0,1].set_ylabel('Mean cloud size [m^2]')
axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])

axarr[0,2].plot(timelist_plot, dimeantauc_list)
axarr[0,2].set_xlabel('time [h/UTC]')
axarr[0,2].set_ylabel('Domain mean tau_c [h]')
axarr[0,2].set_xlim(timelist_plot[0], timelist_plot[-1])

axarr[1,0].plot(timelist_plot, summean_list)
axarr[1,0].set_xlabel('time [h/UTC]')
axarr[1,0].set_xlim(timelist_plot[0], timelist_plot[-1])

axarr[1,1].plot(timelist_plot, r_cluster_list)
axarr[1,1].set_xlabel('time [h/UTC]')
axarr[1,1].set_ylabel('Clustering length [km]')
axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])

if ana == 'm':
    axarr[0,0].set_ylabel('Domain total mass flux [kg/s]')
    axarr[1,0].set_ylabel('Mean cloud mass flux [kg/s]')
    
else:
    axarr[0,0].set_ylabel('Domain total precipitation rate [mm/h]')
    axarr[1,0].set_ylabel('Mean cloud precipitation rate [mm/h]')

fig.suptitle(date + ' ' + ana, fontsize='x-large')
plt.tight_layout(rect=[0, 0.0, 1, 0.95])
fig.savefig(plotdir + 'timeseries')




