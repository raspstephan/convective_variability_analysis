"""
Filename:     helpers.py
Author:       Stephan Rasp, s.rasp@lmu.de

This script contains helper functions.

"""


# import modules
import os
import sys
import yaml
import numpy as np
from netCDF4 import Dataset, date2num
from datetime import datetime, timedelta
from subprocess import check_output
from git import Repo
try:
    from cosmo_utils.pyncdf import getfobj_ncdf_timeseries
    from cosmo_utils.helpers import yyyymmddhh_strtotime
except:
    print 'No cosmo_utils detected, some things might not work.'
from numpy.ma import masked_array
from scipy.ndimage import measurements
from scipy.ndimage.morphology import binary_erosion
from scipy.ndimage.filters import maximum_filter
from skimage import morphology
from scipy.optimize import leastsq
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


# Define functions
def create_log_str(inargs, step):
    """
    Function to create a log file tracking all steps from initial call to
    figure.
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    step : str 
      'Preprocessing' or 'Plotting'
    
    Returns
    -------
    log_str : str
      String with all relevant information
      
    """
    assert step in ['Preprocessing', 'Plotting'], \
        'Step must be Preprocessing or Plotting'

    time_stamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    conda_info = check_output(['conda', 'info'])
    conda_env = check_output(['conda', 'env', 'export'])
    pwd = check_output(['pwd'])
    git_dir = pwd.rsplit('/', 1)[0]
    git_hash = Repo(git_dir).heads[0].commit
    exe_str = ' '.join(sys.argv)
    config_str = open('../config/' + inargs.config_file, 'r').read()
    param_str = ''
    for key, value in vars(inargs).items():
        param_str += '--' + key + ' ' + str(value) + '\n'

    log_str = ("""
%s log\n
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
Full argparse parameters\n
------------------------\n
%s\n
\n
Config yml file\n
---------------\n
%s\n
\n
Anaconda install details\n
------------------------\n
%s\n
\n
Anaconda environment yml file\n
-----------------------------
%s\n\n""" % (step, time_stamp, exe_str, pwd, str(git_hash), param_str,
             config_str, conda_info, conda_env))

    return log_str


def get_config(inargs, top_key, bottom_key):
    """
    Reads the config JSON file

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    top_key : str
      Name of top category 
    bottom_key : key
      Name of key in top category

    Returns
    -------
    value : value for specified key pair
    """
    config = yaml.safe_load(open('../config/' + inargs.config_file))
    return config[top_key][bottom_key]


def make_datelist(inargs, out_format='yyyymmddhh'):
    """
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    out_format : str
      Output format. 'yyyymmddhh'[default] or 'netcdf'
 
    Returns
    -------
    datelist_yyyymmddhh : list
      List containing dates as strings. If date_start == date_end, returns list
      with one element.
    """
    dateob_start = datetime(int(inargs.date_start[0:4]),
                            int(inargs.date_start[4:6]),
                            int(inargs.date_start[6:8]),
                            int(inargs.date_start[8:10]))
    dateob_end = datetime(int(inargs.date_end[0:4]),
                          int(inargs.date_end[4:6]),
                          int(inargs.date_end[6:8]),
                          int(inargs.date_start[8:10]))
    date_inc = timedelta(days=1)

    datelist = []
    date = dateob_start
    while date <= dateob_end:
        if out_format == 'yyyymmddhh':
            date_str = (str(date.year) + str(date.month).zfill(2) +
                        str(date.day).zfill(2) + str(date.hour).zfill(2))
            datelist.append(date_str)
        elif out_format == 'netcdf':
            datelist.append((date - datetime(1, 1, 1)).total_seconds())
        else:
            raise Exception('Wrong format.')
        date += date_inc
    return datelist


def get_domain_limits(inargs):
    """
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------
    l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad : int
      Start and end indices for model and radar domain
    """
    ie = get_config(inargs, 'domain', 'ie')
    je = get_config(inargs, 'domain', 'je')
    ana_irange = get_config(inargs, 'domain', 'ana_irange')
    ana_jrange = get_config(inargs, 'domain', 'ana_jrange')

    # Calculate analysis limits
    assert (ie % 2 == 1) and (je % 2 == 1), 'Function assumes odd ie and je.'
    l11 = (ie - ana_irange - 1) / 2
    l12 = -(l11 + 1)
    l21 = (je - ana_jrange - 1) / 2
    l22 = -(l21 + 1)

    # Get radar limits. Hard coded as of now
    l11_rad = get_config(inargs, 'domain', 'radar_istart')
    l12_rad = get_config(inargs, 'domain', 'radar_istop')
    l21_rad = get_config(inargs, 'domain', 'radar_jstart')
    l22_rad = get_config(inargs, 'domain', 'radar_jstop')

    return l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad


def get_radar_mask(inargs):
    """
    Returns a radar_mask for all analysis days
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    Returns
    -------
    radar_mask :  numpy array
      Has dimensions [date, time, i, j]
      Mask where True values are masked (invalid)
    """
    # Check if radar_mask exists for days
    radar_mask_fn = ('../aux_files/radar_tot_mask_' + inargs.date_start +
                     '_' + inargs.date_end + '_' + str(inargs.time_start) +
                     '_' + str(inargs.time_end) + '.nc')
    if (os.path.isfile(radar_mask_fn)) and (inargs.recompute is False):
        print('Found radar mask: ' + radar_mask_fn)
        mask_rg = Dataset(radar_mask_fn, 'r')
        radar_mask = mask_rg.variables['radar_mask'][:]
        mask_rg.close()
    else:
        print('Compute radar mask: ' + radar_mask_fn)
        l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad = \
            get_domain_limits(inargs)

        # Allocate netCDF file
        mask_rg = Dataset(radar_mask_fn, 'w', format='NETCDF4')
        datearray = np.array(make_datelist(inargs, out_format='netcdf'))
        timearray = np.arange(inargs.time_start,
                              inargs.time_end + inargs.time_inc,
                              inargs.time_inc)
        mask_rg.createDimension('date', datearray.shape[0])
        mask_rg.createDimension('time', timearray.shape[0])
        mask_rg.createDimension('i', get_config(inargs, 'domain', 'ana_irange'))
        mask_rg.createDimension('j', get_config(inargs, 'domain', 'ana_jrange'))

        radar_mask = mask_rg.createVariable('radar_mask', 'i1', ('date', 'time',
                                                                 'i', 'j'))

        for idate, date in enumerate(make_datelist(inargs)):
            datalist = get_datalist_radar(inargs, date)
            for itime, data in enumerate(datalist):
                # Crop data
                radar_mask[idate, itime, :, :] = data > 100.
                # TODO This has to go in the paper!
        radar_mask = radar_mask[:]
        mask_rg.close()

    return radar_mask


def get_pp_fn(inargs, sufx='.nc', pure_fn=False, only_value=True):
    """
    Creates a filename for the pre-processed NetCDF file
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    sufx : str
      Str to attach at the end. Default = '.nc'
    pure_fn : bool  
      If true, path is not included.
    only_value : bool  
      If True [default] only values are plotted in file name
    Returns
    -------
    pp_fn : str
      Filename with path of pre-processed NetCDF file

    """
    if pure_fn:
        pp_fn = ''
    else:
        pp_fn = (get_config(inargs, 'paths', 'preproc_data') + inargs.sub_dir +
                 '/')
        if os.path.exists(pp_fn) is False:
            os.makedirs(pp_fn)
    if inargs.plot_name is not '':
        pp_fn += inargs.pp_name + sufx
    else:
        for key, value in vars(inargs).items():
            if key not in ['recompute', 'plot_name']:
                if only_value:
                    pp_fn += str(value) + '_'
                else:
                    pp_fn += key + '-' + str(value) + '_'
        pp_fn = pp_fn[:-1] + sufx  # remove last '_'
        # Replace ', ' with '_'
        pp_fn = pp_fn.replace(', ', '_')
        # Remove brackets
        pp_fn = pp_fn.replace('[', '')
        pp_fn = pp_fn.replace(']', '')
    assert len(pp_fn) <= 255, ('File name too long with ' + str(len(pp_fn)) +
                               ' ' + pp_fn)
    return pp_fn


def pp_exists(inargs):
    """
    Check whether  preprocessed file exists.
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------
    pp_exists : bool
      True if pre-processed file exitst
    """
    return os.path.isfile(get_pp_fn(inargs))


def load_raw_data(inargs, var, group, lvl=None, radar_mask_type=False):
    """
    This function loads the required COSMO fields and returns a netcdf object 
    which has dimensions [date, time, ens_no, x, y]
    
    Parameters
    ----------
    inargs
    var : str or list
      COSMO variable to be loaded. Options are [PREC_ACCUM, W, QTOT, TTENS_MPHY]
    group : str
      Which dataset. Options are [ens, det, obs]
    lvl : int
      Vertical level for 3D fields
    radar_mask_type : bool or str
      If given, appropriate radar mask is applied. Options are 
      [total, day, hour]

    Returns
    -------
    rootgroup

    """

    if var is not 'PREC_ACCUM' and group is 'obs':
        raise Exception('obs only valid for PREC_ACCUM!')

    # Define file name
    var_str = ''
    if type(var) is not list:
        var = [var]
    for v in var:
        var_str += v + '_'

    fn = (get_config(inargs, 'paths', 'preproc_data') + 'preloaded_fields/' +
          var_str + group + '_' + inargs.date_start + '_' + inargs.date_end +
          '_' + str(inargs.time_start) + '_' + str(inargs.time_end) + '_' +
          str(inargs.time_inc))
    if hasattr(inargs, 'radar_mask'):
        fn += '_' + inargs.radar_mask
    if group == 'ens':
        fn += '_' + str(inargs.nens)
    fn += '.nc'

    # check if preloaded raw data exists
    if os.path.isfile(fn) and inargs.recompute is False:
        print('Found preloaded file: ' + fn)
        rootgroup = Dataset(fn)
        return rootgroup

    # If not preload the raw data, duh
    print('Preload raw data in ' + fn)

    # Create NetCDF file
    rootgroup = Dataset(fn, 'w', format='NETCDF4')

    # Set options for group
    if group in ['det', 'obs']:
        nens = 1
    elif group is 'ens':
        nens = inargs.nens

    # Create dimensions (Partly copied from prec_stats.py)
    dimensions = {
        'date': np.array(make_datelist(inargs, out_format='netcdf')),
        'time': np.arange(inargs.time_start, inargs.time_end + inargs.time_inc,
                          inargs.time_inc),
        'ens_no': np.arange(1, nens + 1),
        'x': np.arange(get_config(inargs, 'domain', 'ana_irange')),
        'y': np.arange(get_config(inargs, 'domain', 'ana_jrange')),
    }
    for dim_name, dim_val in dimensions.items():
        rootgroup.createDimension(dim_name, dim_val.shape[0])
        tmp_var = rootgroup.createVariable(dim_name, 'f8', dim_name)
        tmp_var[:] = dim_val

    # Create variable and arrays to temporarily save data
    var_list = []
    for v in var:
        tmp_var = rootgroup.createVariable(v, 'f8', ['date', 'time', 'ens_no',
                                                     'x', 'y'])
        var_list.append(np.empty(tmp_var[:].shape))

    # If required load radar_mask
    if radar_mask_type is not False:
        radar_mask = get_radar_mask(inargs)
        mask_var = rootgroup.createVariable('mask', 'i4', ['date', 'time', 'x',
                                                           'y'])
        # Determine radar mask
        if radar_mask_type == 'total':
            mask_var[:] = np.any(radar_mask, axis=(0, 1))
        elif radar_mask_type == 'day':
            mask_var[:] = np.any(radar_mask, axis=1)
        elif radar_mask_type == 'hour':
            mask_var[:] = radar_mask

    # Load the data, process and save it in NetCDF file
    for idate, date in enumerate(make_datelist(inargs)):
        print('Loading raw data for: ' + date)

        # Loop over ensemble members and load fields for entire day
        for ie in range(nens):
            for iv, v in enumerate(var):
                if group in ['det', 'ens']:
                    if group == 'det':
                        ens_no = 'det'
                    else:
                        ens_no = ie + 1
                    datalist = get_datalist_model(inargs, date, ens_no, v,
                                                  lvl=lvl)
                elif group == 'obs':
                    datalist = get_datalist_radar(inargs, date)
                else:
                    raise Exception('Wrong group.')

                # Save list in NetCDF file
                # The loop is needed to preserve the mask!
                for it in range(rootgroup.variables['time'].size):
                    var_list[iv][idate, it, ie, :, :] = datalist[it]


    # Now write to file all at once!
    for iv, v in enumerate(var):
        rootgroup.variables[v][:] = var_list[iv]

    return rootgroup


def get_datalist_radar(inargs, date):
    """
    Get data time series for radar observation.
    Parameters
    ----------
    inargs : : argparse object
      Argparse object with all input arguments
    date : str
      Date in format yyyymmddhh. If 'all', radar data for all days are returned
    radar_mask : 2D or 3D numpy array
      Radar mask to create masked arrays.

    Returns
    -------
    datalist : list
      List of 2D masked arrays
    """
    # Get file name
    radarpref = (get_config(inargs, 'paths', 'radar_data') +
                 get_config(inargs, 'paths', 'radar_prefx'))
    radarsufx = get_config(inargs, 'paths', 'radar_sufix')
    dtradar = timedelta(minutes=10)
    if date is 'all':
        date_start = (yyyymmddhh_strtotime(inargs.date_start) +
                      timedelta(hours=inargs.time_start) - dtradar)
        date_end = (yyyymmddhh_strtotime(inargs.date_end) +
                    timedelta(hours=inargs.time_end)-dtradar)
    else:
        dateobj = yyyymmddhh_strtotime(date)
        date_start = dateobj + timedelta(hours=inargs.time_start) - dtradar
        date_end = dateobj + timedelta(hours=inargs.time_end) - dtradar
    datalist = getfobj_ncdf_timeseries(radarpref,
                                       date_start,
                                       date_end,
                                       timedelta(hours=inargs.time_inc),
                                       reftime=date_start,
                                       ncdffn_sufx=radarsufx,
                                       fieldn='pr',
                                       abs_datestr='yymmddhhmm',
                                       dwdradar=True,
                                       return_arrays=True)
    # Crop data
    l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad = \
        get_domain_limits(inargs)
    for i, data in enumerate(datalist):
        datalist[i] = data[l11_rad:l12_rad, l21_rad:l22_rad]

    return datalist


def get_datalist_model(inargs, date, ens_no, var, radar_mask=False, lvl=None):
    """
    Get data time series for model output.
    Parameters
    ----------
    inargs : : argparse object
      Argparse object with all input arguments
    date : str
      Date in format yyyymmddhh
    ens_no : int or str 
      Ensemble number or str in case of det
    var : str 
      Variable
    radar_mask : 2D or 3D numpy array
      Radar mask to create masked arrays
    lvl : int
      Vertical level for 3D data

    Returns
    -------
    datalist : list
      List of 2D masked arrays
    """
    # Get file name
    ncdffn_pref = (get_config(inargs, 'paths', 'raw_data') +
                   date + '/deout_ceu_pspens/' + str(ens_no) +
                   '/OUTPUT/lfff')
    sufx_dict = {
        'PREC_ACCUM': '.nc_30m_surf',
        'CAPE_ML': '.nc_30m_surf',
        'TAU_C': '.nc_30m_surf',
        'HPBL': '.nc_30m_surf',
        'W': '.nc_30m',
        'QC': '.nc_30m',
        'QI': '.nc_30m',
        'QS': '.nc_30m',
        'RHO': '.nc_30m_buoy',
        'TTENS_MPHY': '.nc_30m_buoy',
    }

    datalist = getfobj_ncdf_timeseries(ncdffn_pref,
                                       timedelta(hours=inargs.time_start),
                                       timedelta(hours=inargs.time_end),
                                       timedelta(hours=inargs.time_inc),
                                       ncdffn_sufx=sufx_dict[var],
                                       return_arrays=True,
                                       fieldn=var)

    # Crop data
    l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad = \
        get_domain_limits(inargs)

    # Loop over individual time steps and apply mask
    for i, data in enumerate(datalist):

        if data.ndim == 3:
            data = data[lvl]

        datalist[i] = data[l11:l12, l21:l22]
    return datalist


def read_netcdf_dataset(inargs):
    """
    Open NetCDF file and return rootgroup object.

    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments

    Returns
    -------
    rootgroup : NetCDF object
      NetCDF object
    """
    pp_fn = get_pp_fn(inargs)
    rootgroup = Dataset(pp_fn)
    return rootgroup


def save_fig_and_log(fig, rootgroup, inargs, plot_type='', datestr=None,
                     tight=False):
    """
    Saves given figure and log file with same name.
    
    Parameters
    ----------
    fig : Figure object
    rootgroup : netCDF object
    inargs : argparse object
      Argparse object with all input arguments

    plot_type : str
      str to be added in front of the figure file and log file name
    datestr : str 
      If given, datestr is attached to plot str
    tight : bool
     If True, set bbox_inches=True in figsave

    """
    # Save figure
    plotdir = get_config(inargs, 'paths', 'figures') + inargs.sub_dir + '/'
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    if inargs.plot_name == '':
        plotfn = plotdir + plot_type + get_pp_fn(inargs, sufx='', pure_fn=True)
    else:
        plotfn = plotdir + plot_type + '_' + inargs.plot_name
    if datestr is not None:
        plotfn += '_' + datestr
    plotfn += '.' + inargs.plot_format
    if inargs.plot_format == 'pdf':
        dpi = None
    else:
        dpi = 600
    print('Saving figure: ' + plotfn)
    if tight:
        fig.savefig(plotfn, bbox_inches='tight', dpi=dpi)
    else:
        fig.savefig(plotfn, dpi=dpi)
    plt.close('all')

    # Save log file
    if inargs.plot_name == '':
        logfn = plotdir + plot_type + get_pp_fn(inargs, sufx='.log',
                                                pure_fn=True)
    else:
        logfn = plotdir + plot_type + '_' + inargs.plot_name + '.log'
    logf = open(logfn, 'w+')
    if rootgroup is not None:
        netcdf_log = rootgroup.log + '\n'
        rootgroup.close()
    else:
        netcdf_log = ''
    logf.write(create_log_str(inargs, 'Plotting') + netcdf_log)
    logf.close()


def get_composite_str(inargs, rootgroup):
    """
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    rootgroup : NetCDF object

    Returns
    -------
    date_str : str
      String for composite range
    
    """

    dateobj_start = (timedelta(seconds=int(rootgroup.variables['date'][0])) +
                     datetime(1, 1, 1))
    datestr_start = dateobj_start.strftime(get_config(inargs, 'plotting',
                                                      'date_fmt'))
    dateobj_end = (timedelta(seconds=int(rootgroup.variables['date'][-1])) +
                   datetime(1, 1, 1))
    datestr_end = dateobj_end.strftime(get_config(inargs, 'plotting',
                                                  'date_fmt'))
    return datestr_start + ' - ' + datestr_end


################################################################################
# Functions to put into ensemble_tools
################################################################################
def detect_peaks(image, neighborhood = [[0,1,0],[1,1,1],[0,1,0]]):
    """
    This function is used in identify clouds and is not documented properly!!!
    Takes an image and detect the peaks usingthe local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """

    # define an 8-connected neighborhood
    #neighborhood = generate_binary_structure(2,2)

    #neighborhood = np.ones((5, 5))

    #apply the local maximum filter; all pixel of maximal value
    #in their neighborhood are set to 1
    local_max = maximum_filter(image, footprint=neighborhood)==image
    #local_max is a mask that contains the peaks we are
    #looking for, but also the background.
    #In order to isolate the peaks we must remove the background from the mask.

    #we create the mask of the background
    background = (image==0)

    #a little technicality: we must erode the background in order to
    #successfully subtract it form local_max, otherwise a line will
    #appear along the background border (artifact of the local maximum filter)
    eroded_background = binary_erosion(background, structure=neighborhood,
                                       border_value=1)

    #we obtain the final mask, containing only peaks,
    #by removing the background from the local_max mask
    detected_peaks = local_max - eroded_background

    return detected_peaks


def identify_clouds(field, thresh, opt_field = None, opt_thresh = None,
                    water = False, dx = 2800., rho = None,
                    neighborhood=[[0, 1, 0], [1, 1, 1], [0, 1, 0]],
                    return_com=False):
    """
    Parameters
    ----------
    field : numpy.ndarray
      Field from which clouds are
    thresh : float
      Threshold for field
    opt_field : numpy.ndarray, optional
      Optional field used for creating a binary mask
    opt_thresh : float, optional
      Threshold for opt_field
    water : bool, optional
      If true, watershed algorithm is applied to identify clouds
    dx : float, optional
      Grid spacing [m]
    neighborhood : int or 2D numpy array
      Defines the search perimeter for cloud separation. Only valid of water is
      True
    return_com : bool
     If true, also returns list of centers of mass


    Returns
    -------
    labels : list
      List of labels
    cld_size : list
      List of cloud sizes
    cld_sum : list
      List of summed value of field for each cloud
    cof : np.array
      2D array with centers of mass, if return_com is True

    """

    # Get binary field, 1s where there are clouds
    binfield = field > thresh
    if opt_field is not None:
        binfield *= opt_field > opt_thresh

    if water: # Apply watershed algorithm
        if type(neighborhood) is int:   # Convert integer to matrix
            neighborhood = np.ones((neighborhood, neighborhood))

        # Get local maxima
        lmax = detect_peaks(field*binfield, neighborhood=neighborhood)
        # Do the watershed segmentation
        # Get individual labels for local maxima
        lmax_labels, ncld = measurements.label(lmax)
        labels = morphology.watershed(-field, lmax_labels, mask = binfield)

    else:  # Regular algorithm
        # Find objects
        structure = [[0,1,0],[1,1,1],[0,1,0]]
        labels, ncld = measurements.label(binfield, structure = structure)
        print('Regular ncld = %i' % (ncld))

    # Get sizes and sums
    cld_size = measurements.sum(binfield, labels, range(1, ncld+1))
    if rho is not None:
        field *= rho
    cld_sum = measurements.sum(field, labels, range(1, ncld+1))

    if return_com is not True:
        return labels, cld_size * dx * dx, cld_sum

    else:
        num = np.unique(labels).shape[0]  # Number of identified objects
        # Get centers of mass for each object
        cof = measurements.center_of_mass(field, labels, range(1, num))
        cof = np.array(cof)
        return labels, cld_size * dx * dx, cld_sum, cof


def calc_rdf(labels, field, normalize=True, dx=2800., r_max=30, dr=1, mask=None):
    """
    Computes radial distribution function
    Original credit : Julia Windmiller (MPI)
    
    Parameters
    ----------
    labels : numpy.ndarray
      Array with labels
    field : numpy.ndarray
      Original field Corresponding with labels field
    normalize : bool, optional
      If True normalize RDF
    dx : float, optional
      Grid spacing [m], used for r
    r_max : int, optional
      Maximum search radius for RDF algorithm (in grid pts)
    dr : int, optional
      Search step (in grid pts)
      
    Returns
    -------
    g: numpy.ndarray
      (Normalized) RDF
    r : numpy.ndarray
      Distance
    """

    num = np.unique(labels).shape[0]   # Number of identified objects
    # Get centers of mass for each object
    cof = measurements.center_of_mass(field, labels, range(1,num))
    cof = np.array(cof)

    # If no centers of mass are found, an enpty array is passed
    if cof.shape[0] == 0:   # Accout for empty arrays
        cof = np.empty((0,2))

    g, r, tmp = pair_correlation_2d(cof[:, 0], cof[:, 1],
                                    [field.shape[0], field.shape[1]],
                                    r_max, dr, normalize=normalize, mask=mask)

    return g, r*dx


def pair_correlation_2d(x, y, S, r_max, dr, normalize=True, mask=None):
    """
    Need new doc string
    
    https://github.com/cfinch/colloid/blob/master/adsorption/analysis.py
    
    Compute the two-dimensional pair correlation function, also known
    as the radial distribution function, for a set of circular particles
    contained in a square region of a plane.  This simple function finds
    reference particles such that a circle of radius r_max drawn around the
    particle will fit entirely within the square, eliminating the need to
    compensate for edge effects.  If no such particles exist, an error is
    returned. Try a smaller r_max...or write some code to handle edge effects! ;)
    
    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        S               length of each side of the square region of the plane
        r_max            outer diameter of largest annulus
        dr              increment for increasing radius of annulus
    Returns a tuple: (g, radii, interior_indices)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        annuli used to compute g(r)
        reference_indices   indices of reference particles
    """

    # Number of particles in ring/area of ring/number of reference
    # particles/number density
    # area of ring = pi*(r_outer**2 - r_inner**2)

    # Extract domain size
    (Sx,Sy) = S if len(S) == 2 else (S, S)

    # Find particles which are close enough to the box center that a circle of radius
    # r_max will not cross any edge of the box

    # Find indices within boundaries
    if mask is None:
        bools1 = x > r_max
        bools2 = x < (Sx - r_max)
        bools3 = y > r_max
        bools4 = y < (Sy - r_max)
        interior_indices, = np.where(bools1 * bools2 * bools3 * bools4)
    else:
        # Get closes indices for parcels in a pretty non-pythonic way
        # and check whether it is inside convolved mask
        x_round = np.round(x)
        y_round = np.round(y)
        interior_indices = []
        for i in range(x_round.shape[0]):
            if mask[int(x_round[i]), int(y_round[i])] == 1:
                interior_indices.append(i)

    num_interior_particles = len(interior_indices)

    edges = np.arange(0., r_max + dr, dr)   # Was originally 1.1?
    num_increments = len(edges) - 1
    g = np.zeros([num_interior_particles, num_increments])
    radii = np.zeros(num_increments)
    number_density = float(len(x)) / float(Sx*Sy)

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = np.sqrt((x[index] - x)**2 + (y[index] - y)**2)
        d[index] = 2 * r_max   # Because sqrt(0)

        result, bins = np.histogram(d, bins=edges, normed=False)
        if normalize:
            result = result/number_density
        g[p, :] = result

    # Average g(r) for all interior particles and compute radii
    g_average = np.zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = np.mean(g[:, i]) / (np.pi * (rOuter**2 - rInner**2))

    return g_average, radii, interior_indices


def residual_b_sqrt(p, y, x):
    b = p
    err = np.abs(y - np.sqrt(b * x))

    return err


def residual_b(p, y, x):
    b = p
    err = np.abs(y - (b * x))

    return err


def residual_linear(p, y, x):
    a, b = p
    err = y - (a + b * x)

    return err


def residual_pow(p, y, x):
    a, b = p
    err = np.log(y) - (a - b * np.log(x))

    return err


def residual_exp(p, y, x):
    a, b = p
    err = np.log(y) - (a - b * x)

    return err


def fit_curve(x, y, fit_type='sqrt'):
    """
    
    Parameters
    ----------
    x
    y

    Returns
    -------

    """

    x = np.array(x)
    y = np.array(y)
    mask = np.isfinite(y)
    if fit_type == 'sqrt':
        mask = mask * x > 0
        result = leastsq(residual_b_sqrt, [1], args=(y[mask], x[mask]))
    elif fit_type == 'linear':
        result = leastsq(residual_b, [1], args=(y[mask], x[mask]))
    elif fit_type == 'exp':
        mask = mask * y > 0
        result = leastsq(residual_exp, [10, 1], args=(y[mask], x[mask]))
    elif fit_type == 'pow':
        mask = mask * y > 0
        result = leastsq(residual_pow, [10, 1], args=(y[mask], x[mask]))
    else:
        raise Exception('Wrong fit type!')
    return result[0]


def plot_stamp(inargs, fobj, colors, levels, ax, var):
    """"""
    jpl0, jpl1, ipl0, ipl1 = (50 + inargs.zoom_lat1,
                              357 - 51 + inargs.zoom_lat2,
                              50 + inargs.zoom_lon1,
                              357 - 51 + inargs.zoom_lon2)

    data = fobj.data
    lats = fobj.lats[jpl0:jpl1, ipl0:ipl1]
    lons = fobj.lons[jpl0:jpl1, ipl0:ipl1]
    polelat = fobj.polelat
    polelon = fobj.polelon
    lllat = lats[0, 0]
    lllon = lons[0, 0]
    urlat = lats[-1, -1]
    urlon = lons[-1, -1]
    Basemap_kwargs = { \
        "projection": "lcc",
        # the projection "lcc" lets the domain appear as rectangular on the 2D plot
        "lon_0": polelon,
        "lat_0": 90. + polelat,
        "llcrnrlat": lllat,
        "urcrnrlat": urlat,
        "llcrnrlon": lllon,
        "urcrnrlon": urlon,
        "resolution": 'h',
        "area_thresh": 10000.,
    }
    m = Basemap(**Basemap_kwargs)
    x, y = m(lons, lats)

    if var == 'PREC_ACCUM':
        cfplot = m.contourf(x, y, data[jpl0:jpl1, ipl0:ipl1], levels=levels,
                            colors=colors, ax=ax)
    else:
        cm_prism = plt.cm.prism
        cm_prism.set_under(color='white')
        cfplot = m.imshow(data[jpl0:jpl1, ipl0:ipl1], cmap=cm_prism,
                          origin='lower', vmin=1)

    m.drawcoastlines()  # linewidth=0.1, antialiased=0)
    m.drawcountries(linewidth=0.1, antialiased=0)
    if inargs.ind_scale:
        m.drawmapscale(inargs.ind_scale_pos[0], inargs.ind_scale_pos[1],
                       inargs.ind_scale_pos[0], inargs.ind_scale_pos[1],
                       inargs.ind_scale_len, barstyle='fancy', fontsize=7)

    if inargs.ind_box:
        jpl0_box, jpl1_box, ipl0_box, ipl1_box = (
            50 + inargs.ind_box_lat1,
            357 - 51 + inargs.ind_box_lat2,
            50 + inargs.ind_box_lon1,
            357 - 51 + inargs.ind_box_lon2)

        lon_bl = fobj.lons[jpl0_box, ipl0_box]
        lon_br = fobj.lons[jpl0_box, ipl1_box]
        lon_tl = fobj.lons[jpl1_box, ipl0_box]
        lon_tr = fobj.lons[jpl1_box, ipl1_box]

        lat_bl = fobj.lats[jpl0_box, ipl0_box]
        lat_br = fobj.lats[jpl0_box, ipl1_box]
        lat_tl = fobj.lats[jpl1_box, ipl0_box]
        lat_tr = fobj.lats[jpl1_box, ipl1_box]

        m.plot((lon_bl, lon_br), (lat_bl, lat_br), latlon=True, color='k')
        m.plot((lon_br, lon_tr), (lat_br, lat_tr), latlon=True, color='k')
        m.plot((lon_tr, lon_tl), (lat_tr, lat_tl), latlon=True, color='k')
        m.plot((lon_tl, lon_bl), (lat_tl, lat_bl), latlon=True, color='k')

    return cfplot
