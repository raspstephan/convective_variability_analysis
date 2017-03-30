from cosmo_utils.pyncdf import getfobj_ncdf_timeseries
radarpref = '/project/meteo/w2w/A6/radolan/netcdf_cosmo_de/raa01-rw_10000-'
radarsufx = '-dwd---bin.nc'
from datetime import timedelta
from cosmo_utils.helpers import yyyymmddhh_strtotime
tstart = '2016052801'
tend = '2016060900'

dtradar = timedelta(minutes = 10)

radarts = getfobj_ncdf_timeseries(radarpref, yyyymmddhh_strtotime(tstart)-dtradar, 
                                  yyyymmddhh_strtotime(tend)-dtradar, 
                                  timedelta(minutes=60), 
                                  reftime=yyyymmddhh_strtotime(tstart), ncdffn_sufx=radarsufx, 
                                  fieldn = 'pr', abs_datestr='yymmddhhmm',dwdradar = True)

from cosmo_utils.diag import get_totmask
radarmask = get_totmask(radarts)

sxo, syo = (357, 357)
# Determine analysis domain
lx1 = (sxo-256-1)/2 # ATTENTION first dimension is actually y
lx2 = -(lx1+1) # Number of grid pts to exclude at border
ly1 = (syo-256-1)/2
ly2 = -(ly1+1)
# Crop data
radarmask = radarmask[lx1+62:lx2-42, ly1+22:ly2-42]
import numpy as np
np.save('./radar_tot_mask.npy',radarmask)
history

