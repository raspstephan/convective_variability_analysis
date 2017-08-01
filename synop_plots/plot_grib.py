"""
Script to plot ECMWF synoptic charts.
Main script written by Matthias Schindler, LMU
Adapted by Stephan Rasp, LMU
"""


import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from gribtools import grbdat, latlon
import matplotlib.colors as mcolors
import numpy as np
import datetime
import math as m
from scipy.ndimage.filters import gaussian_filter
import sys
from subprocess import check_output
from git import Repo



basedir = '/project/meteo/scratch/S.Rasp/convective_variability_data/ecmwf_data/'
result_dir = '/home/s/S.Rasp/repositories/convective_variability_analysis/figures/synop_plots/'
exp_od = '0001'
timelist = ['20160530','20160605']


latlon = latlon(basedir + '/0001_20160530_1200UTC_an_0_129_1000.grib')
lat = latlon[0]
print lat.shape
inds = (0, lat.shape[0]/2, 0, -1)
lat = lat[inds[0]:inds[1],inds[2]:inds[3]]      ### filter global data for northern hemisphere to save some memory...
lon = latlon[1]
lon = lon[inds[0]:inds[1],inds[2]:inds[3]]

    
read_data                       = True     
plot_data                       = True
###################################################################################################
###################################################################################################
if read_data:

    print 'reading data: analysis (Geopotential) -- ' + exp_od
    datadict_od_geopot_500      = [{} for i in range(len(timelist))]
    for itim, tim in enumerate(timelist):
        dirname     = basedir + '/'
        filename    = exp_od + '_' + tim + '_' + '1200UTC_an_0_129_500.grib'
        path        = dirname + filename
        data        = grbdat(path, 'Geopotential')
        data        = data[inds[0]:inds[1],inds[2]:inds[3]]
        datadict_od_geopot_500[itim].update({tim : data})
        print path


    print 'reading data: analysis (Geopotential) -- ' + exp_od
    datadict_od_geopot_1000      = [{} for i in range(len(timelist))]
    for itim, tim in enumerate(timelist):
        dirname     = basedir + '/'
        filename    = exp_od + '_' + tim + '_' + '1200UTC_an_0_129_1000.grib'
        path        = dirname + filename
        data        = grbdat(path, 'Geopotential')
        data        = data[inds[0]:inds[1],inds[2]:inds[3]]
        datadict_od_geopot_1000[itim].update({tim : data})
        print path

    
    print 'reading data: analysis (Mean sea level pressure) -- ' + exp_od
    datadict_od_slp              = [{} for i in range(len(timelist))]
    for itim, tim in enumerate(timelist):
        dirname     = basedir + '/'
        filename    = exp_od + '_' + tim + '_' + '1200UTC_an_0_151_sfc.grib'
        path        = dirname + filename
        data        = grbdat(path, 'Mean sea level pressure')
        data        = gaussian_filter(data[inds[0]:inds[1],inds[2]:inds[3]], 10,
                                      mode='wrap')
        datadict_od_slp[itim].update({tim : data})
        print path
            
    print 'relative topography'
    
    print 'read_data complete!'
###################################################################################################
###################################################################################################
if plot_data: # projection globe, grid, sp, contourlines in grey
    # Good colormap from hclwizard.org
    colors = (
    "#F5A5FF", "#C4B7FF", "#78CAFF", "#00D9FF", "#00E3FF", "#00E9E8", "#00EAC4",
    "#00E89A", "#26E46A", "#87DC25", "#B9D300", "#DEC900", "#FCBC00", "#FFAF58",
    "#FFA28D", "#FF97BA")
    print 'plt_diff_fc_abs_error'
    fig, axarr = plt.subplots(1, 2, figsize=(7.87, 3))
    for itim, tim in enumerate(timelist):
        print itim, tim
        dateobj = datetime.datetime.strptime(tim, '%Y%m%d')
        ax = axarr[itim]
        plt.sca(ax)

        m           = Basemap(projection='cyl', llcrnrlon=-40, \
                   urcrnrlon=30.,llcrnrlat=25,urcrnrlat=80, \
                   resolution='i', area_thresh=10000.)
        m.drawcoastlines(False)
        m.drawcountries()
        m.drawcoastlines(linewidth=0.5)
        x,y         = m(lon,lat)
        mylevels    = np.arange(520,590,5)
        # cs0         = m.contourf(x,y,((datadict_od_geopot_500[itim][timelist[itim]]/9.80665)/10.)-((datadict_od_geopot_1000[itim][timelist[itim]]/9.80665)/10.),
        #                          shading='flat',latlon=True, levels=mylevels,
        #                          cmap=plt.cm.gist_rainbow_r, extend = 'both',
        #                          alpha=0.8)
        # cs1         = m.contour(x,y,(datadict_od_geopot_500[itim][timelist[itim]]/9.80665)/10.,
        #                         latlon = True, levels = np.arange(520,581,5),
        #                         colors='k', linewidths=2,
        #                         zorder=3)
        cs1         = m.contourf(x,y,(datadict_od_geopot_500[itim][timelist[itim]]/9.80665)/10.,
                                 latlon=True, levels=np.arange(520,601,5),
                                 colors=colors, extend = 'both')
        cs2         = m.contour(x,y,datadict_od_slp[itim][timelist[itim]]/100.,
                                latlon = True, levels = np.arange(980,1030,5),
                                colors='w', linewidths=1, zorder=10)

        m.drawmapboundary()
        m.drawparallels(np.arange(-90.,91.,20.),  labels=[1,0,0,0], linewidth = 0.5, fontsize = 8 )
        m.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1], linewidth = 0.5, fontsize = 8 )

        # plt.clabel(cs1, fontsize=9, inline=1, fmt='%d')
        plt.clabel(cs2, fontsize=9, inline=1, fmt='%d')
        let = 'a) ' if itim == 0 else 'b) '
        titlestr = let + dateobj.strftime('%d %b')
        ax.set_title(titlestr)

    cax = fig.add_axes([0.92, 0.09, 0.015, 0.81])

    cb = fig.colorbar(cs1, cax=cax, orientation='vertical')
    cb.set_label(r'500 hPa $\Phi$ [dam]', fontsize=8)
    cb.ax.tick_params(labelsize=8)
    plt.subplots_adjust(left=0.05, bottom=0.05, top=0.94, right=0.92)
    fig.savefig(result_dir + 'synop.pdf')
    fig.clf()
    print 'plt_diff_fc_abs_error'

# Create log string; this is copied mostly from the create_log_str function
# in helpers.py
time_stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
conda_info = check_output(['conda', 'info'])
conda_env = check_output(['conda', 'env', 'export'])
pwd = check_output(['pwd'])
git_dir = pwd.rsplit('/', 1)[0]
git_hash = Repo(git_dir).heads[0].commit
exe_str = ' '.join(sys.argv)


log_str = ("""
#####################################################
Time: %s\n
Executed command\n
----------------\n
python %s\n
\n
in directory: %s\n
\n
Git hash: %s\n
\n
Anaconda install details\n
------------------------\n
%s\n
\n
Anaconda environment yml file\n
-----------------------------
%s\n\n""" % (time_stamp, exe_str, pwd, str(git_hash), conda_info, conda_env))

logfn = result_dir + 'synop.log'
logf = open(logfn, 'w+')
logf.write(log_str)
logf.close()

