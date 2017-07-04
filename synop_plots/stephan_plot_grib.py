import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from gribtools import grbdat, latlon
import matplotlib.colors as mcolors
import numpy as np
import datetime
import math as m



basedir         = '/project/meteo/w2w/A3/M.Schindler/PhD/experiments/data'
result_dir      = '/project/meteo/w2w/A3/M.Schindler/PhD/experiments/evaluation/'
exp_od          = '0001' 
timelist        = ['20160530','20160605']


latlon = latlon(basedir + '/0070/0070_20160925_1200UTC_an_0_129_500.grib')
lat = latlon[0]
lat = lat[:(lat.shape[0]/2),:]                                                           ### filter global data for northern hemisphere to save some memory...
lon = latlon[1]
lon = lon[:(lon.shape[0]/2),:]
    
read_data                       = True     
plot_data                       = True
###################################################################################################
###################################################################################################
if read_data:

    print 'reading data: analysis (Geopotential) -- ' + exp_od
    datadict_od_geopot_500      = [{} for i in range(len(timelist))]
    for itim, tim in enumerate(timelist):
        dirname     = basedir + '/' + 'Stephan_Rasp' + '/'
        filename    = exp_od + '_' + tim + '_' + '1200UTC_an_0_129_500.grib'
        path        = dirname + filename
        data        = grbdat(path, 'Geopotential')
        data        = data[:(lat.shape[0]),:]
        datadict_od_geopot_500[itim].update({tim : data})
        print path


    print 'reading data: analysis (Geopotential) -- ' + exp_od
    datadict_od_geopot_1000      = [{} for i in range(len(timelist))]
    for itim, tim in enumerate(timelist):
        dirname     = basedir + '/' + 'Stephan_Rasp' + '/'
        filename    = exp_od + '_' + tim + '_' + '1200UTC_an_0_129_1000.grib'
        path        = dirname + filename
        data        = grbdat(path, 'Geopotential')
        data        = data[:(lat.shape[0]),:]
        datadict_od_geopot_1000[itim].update({tim : data})
        print path

    
    print 'reading data: analysis (Mean sea level pressure) -- ' + exp_od
    datadict_od_slp              = [{} for i in range(len(timelist))]
    for itim, tim in enumerate(timelist):
        dirname     = basedir + '/' + 'Stephan_Rasp' + '/'
        filename    = exp_od + '_' + tim + '_' + '1200UTC_an_0_151_sfc.grib'
        path        = dirname + filename
        data        = grbdat(path, 'Mean sea level pressure')
        data        = data[:(lat.shape[0]),:]
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
        fig         = plt.figure(figsize=(12,6))
        ax          = fig.add_subplot(111)
        m           = Basemap(projection='cyl', llcrnrlon=-120, \
                   urcrnrlon=40.,llcrnrlat=0,urcrnrlat=lat.max(), \
                   resolution='c')
        x,y         = m(lon,lat)
        mylevels    = np.arange(480.,604.,4.)
        cs0         = m.contourf(x,y,((datadict_od_geopot_500[itim][timelist[itim]]/9.80665)/10.)-((datadict_od_geopot_1000[itim][timelist[itim]]/9.80665)/10.), shading='flat',latlon=True, levels=mylevels, cmap=plt.cm.gist_rainbow_r, extend = 'both')
        cs1         = m.contour(x,y,(datadict_od_geopot_500[itim][timelist[itim]]/9.80665)/10., latlon = True, levels = np.arange(520,581,10), colors='k', linewidths=1.4)
        cs2         = m.contour(x,y,datadict_od_slp[itim][timelist[itim]]/100., latlon = True, levels = np.arange(980,1030,10), colors='w', linewidths=1.4)
        m.drawcoastlines(linewidth=0.5)
        m.drawmapboundary()
        m.drawparallels(np.arange(-90.,91.,30.),  labels=[1,0,0,0], linewidth = 0., fontsize = 8 )
        m.drawmeridians(np.arange(-180.,181.,30.),labels=[0,0,0,1], linewidth = 0., fontsize = 8 )
        fig.tight_layout()
        cb = fig.colorbar(cs0,orientation='vertical', shrink=0.9)
        cb.set_label('gpdam', fontsize = 8)
        plt.clabel(cs1, fontsize=9, inline=1, fmt='%d')
        plt.clabel(cs2, fontsize=9, inline=1, fmt='%d')
        plt.title('500 hPa Geopotential [gpdam], Surface pressure [hPa] and relative topography [gpdam], ' + str(tim), fontsize = 10)
        fig.savefig(result_dir + str(tim) + '_srasp_test.png', dpi = 300, bbox_inches='tight')
        fig.clf()
    print 'plt_diff_fc_abs_error'
