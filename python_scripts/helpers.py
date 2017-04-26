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
from cosmo_utils.helpers import yyyymmddhh_strtotime, yymmddhhmm
from numpy.ma import masked_array
from cosmo_utils.pyncdf import getfobj_ncdf_timeseries, getfobj_ncdf
from scipy.ndimage import measurements


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
    conda_list = check_output(['conda', 'list'])
    # TODO: Get base dir automatically
    git_dir = '~/repositories/convective_variability_analysis'
    git_hash = Repo(git_dir).heads[0].commit
    pwd = check_output(['pwd'])
    exe_str = ' '.join(sys.argv)
    config_str = open('../config/' + inargs.config_file, 'r').read()
    param_str = ''
    for key, value in vars(inargs).items():
        param_str += '--' + key + ' ' + str(value) + '\n'

    log_str = ("""
    %s log\n
    -----------------\n
    %s\n
    %s\n
    Full argparse parameters\n
    %s
    %s\n
    %s\n
    Git hash: %s\n
    In directory: %s\n
    %s\n
    """ % (step, time_stamp, config_str, param_str, conda_info, conda_list,
           str(git_hash)[0:7], pwd, exe_str))
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


def get_datalist_radar(inargs, date, radar_mask=False):
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
        if radar_mask is False:
            tmp_mask = np.zeros(data[l11_rad:l12_rad, l21_rad:l22_rad].shape)
        elif radar_mask.ndim == 3:
            tmp_mask = radar_mask[i, :, :]
        elif radar_mask.ndim == 2:
            tmp_mask = radar_mask
        else:
            raise Exception('Wrong dimensions for radar mask.')

        datalist[i] = masked_array(data[l11_rad:l12_rad, l21_rad:l22_rad],
                                   mask=tmp_mask)

    return datalist


def get_datalist_model(inargs, date, ens_no, var, radar_mask):
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

    Returns
    -------
    datalist : list
      List of 2D masked arrays
    """
    # Get file name
    ncdffn_pref = (get_config(inargs, 'paths', 'raw_data') +
                   date + '/deout_ceu_pspens/' + str(ens_no) +
                   '/OUTPUT/lfff')
    datalist = getfobj_ncdf_timeseries(ncdffn_pref,
                                       timedelta(hours=inargs.time_start),
                                       timedelta(hours=inargs.time_end),
                                       timedelta(hours=inargs.time_inc),
                                       ncdffn_sufx='.nc_30m_surf',
                                       return_arrays=True,
                                       fieldn=var)

    # Crop data
    l11, l12, l21, l22, l11_rad, l12_rad, l21_rad, l22_rad = \
        get_domain_limits(inargs)

    # Loop over individual time steps and apply mask
    for i, data in enumerate(datalist):
        if radar_mask.ndim == 3:
            tmp_mask = radar_mask[i, :, :]
        elif radar_mask.ndim == 2:
            tmp_mask = radar_mask
        else:
            raise Exception('Wrong dimensions for radar mask.')

        datalist[i] = masked_array(data[l11:l12, l21:l22],
                                   mask=tmp_mask)
        if var == 'TAU_C':   # Additionally mask out nans
            datalist[i] = masked_array(data[l11:l12, l21:l22],
                                       mask=tmp_mask + np.isnan(
                                           data[l11:l12, l21:l22]))
    return datalist


def get_and_crop_radar_fobj(inargs, date, time):
    """
    Returns radar fobj, cropped to my specific model domain.
    
    Parameters
    ----------
    inargs : argparse object
      Argparse object with all input arguments
    date : str 
      Date in yyyymmddhh format
    time : timedelta object
      Forecast time

    Returns
    -------
    fobj : cosmo_utils field object
      Cropped radar object
    """
    radarpref = (get_config(inargs, 'paths', 'radar_data') +
                 get_config(inargs, 'paths', 'radar_prefx'))
    radarsufx = get_config(inargs, 'paths', 'radar_sufix')
    dtradar = timedelta(minutes=10)
    t_rad = yyyymmddhh_strtotime(date) + time - dtradar
    fobj = getfobj_ncdf(radarpref + yymmddhhmm(t_rad) +
                        radarsufx, 'pr', dwdradar=True)
    # Crop radar field
    fobj.data = fobj.data[62:62 + 357, 22:22 + 357]
    fobj.lats = fobj.lats[62:62 + 357, 22:22 + 357]
    fobj.lons = fobj.lons[62:62 + 357, 22:22 + 357]
    fobj.rlats = fobj.rlats[62:62 + 357, 22:22 + 357]
    fobj.rlons = fobj.rlons[62:62 + 357, 22:22 + 357]
    fobj.nx = 357
    fobj.ny = 357

    return fobj


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


def save_fig_and_log(fig, rootgroup, inargs, plot_type='', date=None,
                     time=None):
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
    date : str 
      If given, date is attached to plot str
    time : str 
      If given, time is attached to plot str

    """
    # Save figure
    plotdir = get_config(inargs, 'paths', 'figures') + inargs.sub_dir + '/'
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    if inargs.plot_name == '':
        plotfn = plotdir + plot_type + get_pp_fn(inargs, sufx='', pure_fn=True)
    else:
        plotfn = plotdir + plot_type + '_' + inargs.plot_name
    if date is not None and time is not None:
        plotfn += '_' + str(date) + '_' + str(time)
    plotfn += '.pdf'
    print('Saving figure: ' + plotfn)
    fig.savefig(plotfn)

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
    logf.write(netcdf_log + create_log_str(inargs, 'Plotting'))
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
