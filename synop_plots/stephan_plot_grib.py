import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from gribtools import grbdat, latlon
import matplotlib.colors as mcolors
import numpy as np
import datetime
import math as m
from scipy.ndimage.filters import gaussian_filter



basedir         = '/home/s/S.Rasp/repositories/convective_variability_analysis/synop_plots'
result_dir      = '/home/s/S.Rasp/repositories/convective_variability_analysis/figures/synop_plots/'
exp_od          = '0001' 
timelist        = ['20160530','20160605']


latlon = latlon(basedir + '/0001_20160530_1200UTC_an_0_129_1000.grib')
lat = latlon[0]
print lat.shape
inds = (0, lat.shape[0]/2, 0, -1)
lat = lat[inds[0]:inds[1],inds[2]:inds[3]]                                                           ### filter global data for northern hemisphere to save some memory...
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
    
    print 'plt_diff_fc_abs_error'
    for itim, tim in enumerate(timelist):
        print itim, tim
        fig         = plt.figure(figsize=(6,5))
        ax          = fig.add_subplot(111)
        m           = Basemap(projection='cyl', llcrnrlon=-40, \
                   urcrnrlon=30.,llcrnrlat=25,urcrnrlat=80, \
                   resolution='i', area_thresh=10000.)
        m.drawcoastlines(False)
        m.drawcountries()
        x,y         = m(lon,lat)
        mylevels    = np.arange(520,590,5)
        cs0         = m.contourf(x,y,((datadict_od_geopot_500[itim][timelist[itim]]/9.80665)/10.)-((datadict_od_geopot_1000[itim][timelist[itim]]/9.80665)/10.), shading='flat',latlon=True, levels=mylevels, cmap=plt.cm.gist_rainbow_r, extend = 'both')
        cs1         = m.contour(x,y,(datadict_od_geopot_500[itim][timelist[itim]]/9.80665)/10., latlon = True, levels = np.arange(520,581,5), colors='k', linewidths=2,
                                zorder=3)
        cs2         = m.contour(x,y,datadict_od_slp[itim][timelist[itim]]/100., latlon = True, levels = np.arange(980,1030,5), colors='w', linewidths=1, zorder=2)
        m.drawcoastlines(linewidth=0.5)
        m.drawmapboundary()
        m.drawparallels(np.arange(-90.,91.,20.),  labels=[1,0,0,0], linewidth = 0.5, fontsize = 8 )
        m.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1], linewidth = 0.5, fontsize = 8 )
        fig.tight_layout()
        cb = fig.colorbar(cs0,orientation='vertical', shrink=0.8, pad=0.0)
        cb.set_label('gpdam', fontsize = 8)
        plt.clabel(cs1, fontsize=9, inline=1, fmt='%d')
        plt.clabel(cs2, fontsize=9, inline=1, fmt='%d')
        plt.title('500 hPa Geopotential [gpdam], Surface pressure [hPa]\nrelative topography [gpdam], ' + str(tim)+'00', fontsize = 10)
        fig.savefig(result_dir + str(tim) + '_srasp_test.png', dpi = 500, bbox_inches='tight')
        fig.clf()
    print 'plt_diff_fc_abs_error'
