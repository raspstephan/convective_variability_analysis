"""
Saves example fields for cloud_identification_and_rdf notebook.
"""

import sys
sys.path.append('../python_scripts')   # So that we can later use the functions from helpers.py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from cosmo_utils.pyncdf import getfobj_ncdf   # This is a module developed specifically for COSMO output

savedir = '/home/s/S.Rasp/tmp/saved_fields/'

# Pick 6 examples
alims = (50, -51, 50, -51)   # Limits for analysis area
dates = ['2016052800', '2016052800', '2016060300', 
         '2016060300', '2016053000', '2016052800']
times = [15, 20, 13, 21, 6, 9]
tags = ['a) Scattered', 'b) Clustered', 'c) Scattered', 
        'd) Clustered', 'e) Single large cloud', 'f) Hardly any clouds']

# Plot the fields
fig, axmat = plt.subplots(2, 3, figsize = (13, 10))
cm = plt.cm.cool
cm.set_under(color = 'white')
field_list = []
for d, t, tag, ax in zip(dates, times, tags, np.ravel(axmat)):
    
    fn = ('/project/meteo/scratch-old/users/stephan.rasp/convective_variability_data/raw_data/' + 
            d + '/deout_ceu_pspens/det/OUTPUT/lfff00' + str(t).zfill(2) + '0000.nc_30m_surf')
    field = getfobj_ncdf(fn, 'PREC_ACCUM').data[alims[0]:alims[1], alims[2]:alims[3]]
    np.save(savedir + d + '_00' + str(t).zfill(2) + '0000_prec', field)
    
    fn = ('/project/meteo/scratch-old/users/stephan.rasp/convective_variability_data/raw_data/' + 
            d + '/deout_ceu_pspens/det/OUTPUT/lfff00' + str(t).zfill(2) + '0000.nc_30m')
    field = getfobj_ncdf(fn, 'W').data[30, alims[0]:alims[1], alims[2]:alims[3]]   # At model level 30
    np.save(savedir + d + '_00' + str(t).zfill(2) + '0000_w', field)