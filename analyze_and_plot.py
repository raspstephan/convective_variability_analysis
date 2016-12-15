"""
This script reads the saved ncdf files, does some easy analysis and produces 
plots.
2016052800 2016052900 2016053000 2016053100 2016060100 2016060200 2016060300 2016060400 2016060500 2016060600 2016060700 2016060800
"""

# Imports
import warnings
warnings.filterwarnings("ignore")   # ATTENTION For bool index warning
import os
import argparse
from netCDF4 import Dataset
import numpy as np
from datetime import timedelta
from cosmo_utils.helpers import ddhhmmss, yyyymmddhh_strtotime, yymmddhhmm
from cosmo_utils.pyncdf import getfobj_ncdf, getfobj_ncdf_ens
from cosmo_utils.pywgrib import fieldobj
from cosmo_utils.plot import ax_contourf, fig_contourf_1sp
from cosmo_utils.diag import mean_spread_fieldobjlist
import matplotlib as mpl
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
] 
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import cPickle

# Some plotting setup
pdfwidth = 7.87   # inches for AMetSoc
mpl.rcParams['font.size'] = 10
mpl.rcParams['font.family'] = 'sans-serif'


# Define functions
def residual_bc(p, y, x):
    
    b,c = p
    err = y - (b*x**c)
    
    return err

def residual_ab(p, y, x):
    a,b = p
    err = np.abs(y - (a + b*x))
    
    return err

def residual_ab_exp(p, y, x):
    a, b = p
    err = np.log(np.abs(y - (a*np.exp(b*x))))
    
    return err

def residual_abc(p, y, x):
    a,b = p
    err = np.abs(y - (a + (b*x**c)))
    
    return err

def residual_b(p, y, x):
    b = p
    err = np.abs(y - (b*x))
    
    return err

def residual_b_sqrt(p, y, x):
    b = p
    err = np.abs(y - np.sqrt(b*x))
    
    return err

def residual_exp(p, y, x):
    a,b = p
    err = np.log(y)-(a-b*x)
    
    return err

def residual_pow(p, y, x):
    a,b = p
    err = np.log(y)-(a-b*np.log(x))
    
    return err

# Setup: this needs to match the compute script!
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--ana', metavar = 'ana', type=str, default = 'm')
parser.add_argument('--date', metavar = 'date', type=str, nargs = '+')
parser.add_argument('--height', metavar = 'height', type=float, nargs = '+',
                    default = [3000])
parser.add_argument('--water', metavar = 'water', type=str, default = 'True')
parser.add_argument('--nens', metavar = 'nens', type=int, default = 20)
parser.add_argument('--tstart', metavar = 'tstart', type=int, default = 1)
parser.add_argument('--tend', metavar = 'tend', type=int, default = 24)
parser.add_argument('--tinc', metavar = 'tinc', type=int, default = 60)
parser.add_argument('--plot', metavar = 'plot', type=str, nargs = '+')
parser.add_argument('--tplot', metavar = 'tplot', type=float, nargs = '+',
                    default = [9,12])
parser.add_argument('--minmem', metavar = 'minmem', type=int, default = 5)
parser.add_argument('--dr', metavar = 'dr', type=int, default = 2)
parser.add_argument('--hypo', metavar = 'hypo', type=str, default = 'False')
parser.add_argument('--det', metavar = 'det', type=str, default = 'False')
args = parser.parse_args()

if args.hypo == 'True':
    args.hypo = True
    print 'Hypo'
else:
    args.hypo = False

if 'all' in args.plot:
    args.plot = ['cloud_stats', 'rdf', 'scatter', 'summary_stats', 'summary_var',
                 'stamps_var', 'stamps_w', 'height_var', 'std_v_mean', 'prec_stamps',
                 'prec_hist']

################################################################################
# Load datatsets
datasetlist = []
heightstr = ''
alldatestr = ''
nlist = [256, 128, 64, 32, 16, 8, 4]
for h in args.height:
    heightstr += str(int(h))
for d in args.date:
    alldatestr += d + '_'
    print 'Loading date: ', d
    # Create file str
    savedir = '/home/scratch/users/stephan.rasp/results/'
    anastr = ('_ana-' + args.ana + '_wat-' + str(args.water) + 
                '_height-' + heightstr +
                '_nens-' + str(args.nens) + '_tstart-' + str(args.tstart) + 
                '_tend-' + str(args.tend) + '_tinc-' + str(args.tinc) + 
                '_minmem-' + str(args.minmem) + '_dr-' + str(args.dr))
    anastr_det = ('_ana-' + args.ana + '_wat-' + str(args.water) + 
                '_height-' + heightstr +
                '_nens-' + str(1) + '_tstart-' + str(args.tstart) + 
                '_tend-' + str(args.tend) + '_tinc-' + str(args.tinc) + 
                '_minmem-' + str(args.minmem) + '_dr-' + str(args.dr))
    savesuf = anastr + '.nc'
    if args.det == 'True':
        savesuf += '_det'
    # Load file 
    #datasetlist.append(Dataset(savedir + savestr, 'r'))
alldatestr = alldatestr[:-1]   # revome final underscore
if alldatestr == '2016052800_2016052900_2016053000_2016053100_2016060100_2016060200_2016060300_2016060400_2016060500_2016060600_2016060700_2016060800':
    alldatestr = 'composite'
    print 'Composite!'
if args.det == 'True':
    alldatestr += '_det'
# End load datasets
################################################################################


# Some preliminaries
print savedir + args.date[0] + savesuf
tmp_dataset = Dataset(savedir + args.date[0] + savesuf, 'r')
timelist = [timedelta(seconds=ts) for ts in tmp_dataset.variables['time']]
timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
plotdir = ('/home/s/S.Rasp/Dropbox/figures/PhD/variance/' + alldatestr + 
           '/')

ensdir = '/home/scratch/users/stephan.rasp/' + args.date[0] + '/deout_ceu_pspens/'


# Define helper function
def create1Dlist(s_i):
    l = []
    for i in range(s_i):
        l.append([])
    return l

def create2Dlist(s_i, s_j):
    l = []
    for i in range(s_i):
        l.append([])
        for j in range(s_j):
            l[-1].append([])
    return l

def create3Dlist(s_i, s_j, s_k):
    l = []
    for i in range(s_i):
        l.append([])
        for j in range(s_j):
            l[-1].append([])
            for k in range(s_k):
                l[-1][-1].append([])
    return l
    


# Now comes the plotting
################################################################################
if 'cloud_stats' in args.plot:
    print 'Plotting cloud_stats'
    dataset = Dataset(savedir + args.date[0] + savesuf, 'r')
    plotdirsub = plotdir +  '/cloud_stats/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    # Setup
    dx2 = 2.8e3**2
    print 'dx2', dx2
    nbins = 60.
    sizemax = dx2*nbins
    summax = 10e8
    if args.water == 'False':
        sizemax *= 3.5
        summax *=3.5

    savename = savedir + 'cloud_stats_' + alldatestr + anastr
    if os.path.exists(savename):
        print 'Loading pre-saved file', savename
        savefile = open(savename, 'r')
        totlist1 = cPickle.load(savefile)
        totlist2 = cPickle.load(savefile)
        savefile.close()
        
    else:
        
        totlist1 = []
        totlist2 = []
        for d in args.date:
            print d
            dataset = Dataset(savedir + d + savesuf, 'r')
            cld_size_tmp = dataset.variables['cld_size'][:, :]
            #print cld_size_tmp
            cld_size_tmp = cld_size_tmp[~cld_size_tmp.mask].data
            totlist1 += list(cld_size_tmp)
            
            if not args.hypo:
                cld_sum_tmp = dataset.variables['cld_sum'][:, :]
                cld_sum_tmp = cld_sum_tmp[~cld_sum_tmp.mask].data
                totlist2 += list(cld_sum_tmp)
        savefile = open(savename, 'w')
        cPickle.dump(totlist1, savefile, -1)
        cPickle.dump(totlist2, savefile, -1)

        savefile.close()
        
        
    sizehist, sizeedges = np.histogram(totlist1, 
                                        bins = nbins, range = [0., sizemax])
    sizemean = np.mean(totlist1)
    sizevar = np.var(totlist1)
    print 'size beta', sizevar/sizemean**2, sizemean
    print 'size mean', sizemean
    
    if not args.hypo:
        sumhist, sumedges = np.histogram(totlist2, 
                                            bins = nbins, range = [0., summax])
        summean = np.mean(totlist2)
        sumvar = np.var(totlist2)
        print 'beta', sumvar/summean**2, summean
    
    # Plot the histograms
    fig, axarr = plt.subplots(1, 3, figsize = (pdfwidth, 3))
    axarr[0].bar(sizeedges[:-1], sizehist, width = np.diff(sizeedges)[0],
                 color = 'darkgray')
    axarr[0].plot([sizemean, sizemean], [1, 1e7], c = 'red', 
                alpha = 1, label = 'mean')
    ## Fit line, 1st exp
    p0 = [10, 1]
    mask = sizehist>0
    xfit = (sizeedges[:-1] + sizeedges[1:]) / 2.
    result = leastsq(residual_exp, p0, args = (sizehist[mask], xfit[mask]))
    a,b = result[0]
    axarr[0].plot(xfit,np.exp(a-b*xfit), c = 'orange', label = 'exponential')
    
    ## Fit line, 2nd power law
    p0 = [10, 1]
    mask = sizehist>0
    xfit = (sizeedges[:-1] + sizeedges[1:]) / 2.
    result = leastsq(residual_pow, p0, args = (sizehist[mask], xfit[mask]))
    a,b = result[0]
    axarr[0].plot(xfit,np.exp(a-b*np.log(xfit)), c = 'blue', label = 'power law')

    axarr[0].set_xlabel(r'Cloud size [m$^2$]')
    axarr[0].set_ylabel('Number of clouds')
    axarr[0].set_title('Cloud size', fontsize = 10)
    axarr[0].set_xlim([0., sizemax])
    axarr[0].set_ylim([1, 1e7])
    axarr[0].set_yscale('log')
    axarr[0].legend(prop={'size':8}, loc = 1)
    axarr[0].text(0.05, 0.9, '(a)', transform = axarr[0].transAxes, 
                    fontsize = 10)
    
    if not args.hypo:
        axarr[1].bar(sumedges[:-1], sumhist, width = np.diff(sumedges)[0],
                    color = 'darkgray')
        axarr[1].plot([summean, summean], [1, 1e7], c = 'red', 
                    alpha = 1)
        axarr[1].set_ylabel('Number of clouds')
        axarr[1].set_xlim([0., summax])
        axarr[1].set_ylim([1, 1e7])
        axarr[1].set_yscale('log')
        axarr[1].set_xlabel(r'm [kg s$^{-1}$]')
        axarr[1].set_title('Mass flux per cloud', fontsize = 10)
        axarr[1].text(0.05, 0.9, '(b)', transform = axarr[1].transAxes, 
                    fontsize = 10)
        
        ## Fit line, 1st exp
        p0 = [10, 1]
        mask = sumhist>0
        xfit = (sumedges[:-1] + sumedges[1:]) / 2.
        result = leastsq(residual_exp, p0, args = (sumhist[mask], xfit[mask]))
        a,b = result[0]
        axarr[1].plot(xfit,np.exp(a-b*xfit), c = 'orange', label = 'exponential')
        
        ## Fit line, 2nd power law
        p0 = [10, 1]

        result = leastsq(residual_pow, p0, args = (sumhist[mask], xfit[mask]))
        a,b = result[0]

        axarr[1].plot(xfit,np.exp(a-b*np.log(xfit)), c = 'blue', label = 'power law')
        
        Rcorr = np.corrcoef(totlist1, totlist2)[1,0]
        print 'corr', np.corrcoef(totlist1, totlist2)[1,0]
        print 'n_cld', len(totlist1)
        #axarr[2].scatter(totlist1, totlist2, c = 'grey', linewidth = 0.1, s=4)
        #tmp = np.logspace(6,10,100)
        #slope = summean/sizemean

        #axarr[2].plot(tmp, tmp*slope, c = 'k')
        #if args.water == 'False':
            #axarr[2].set_xlim([5e6, 3e9])
            #axarr[2].set_ylim([5e6, 5e9])
        #else:
            #axarr[2].set_xlim([5e6, 5e8])
            #axarr[2].set_ylim([5e6, 2e9])
        #axarr[2].set_yscale('log')
        #axarr[2].set_xscale('log')
        #axarr[2].set_xlabel('Cloud size [m$^2$]')
        #axarr[2].set_ylabel(r'm [kg s$^{-1}$]')
        #axarr[2].text(0.05, 0.9, '(c)', transform = axarr[2].transAxes, 
                    #fontsize = 10)
        #axarr[2].text(0.65, 0.05, 'R={:.2f}'.format(Rcorr), transform = axarr[2].transAxes, 
                    #fontsize = 10)
        xfit = (sizeedges[:-1] + sizeedges[1:]) / 2.
        axarr[2].scatter(xfit, sizehist)
        axarr[2].set_yscale('log')
        axarr[2].set_xscale('log')
        axarr[2].set_ylim([1, 1e7])
        axarr[2].set_xlabel(r'Cloud size [m$^2$]')
        axarr[2].set_ylabel('Number of clouds')
        axarr[2].set_title('Cloud size on log-log', fontsize = 10)
        
    
    
    # # Now the LOG hists
    # logbins = np.logspace(6, 10, 20)
    # sizehist, sizeedges = np.histogram(totlist1, bins = logbins)
    # if not args.hypo:
    #     sumhist, sumedges = np.histogram(totlist2, bins = logbins)
        
    # axarr[1,0].bar(sizeedges[:-1], sizehist, width = np.diff(sizeedges),
    #              color = 'darkgray')
    # axarr[1,0].plot([sizemean, sizemean], [1, 1e6], c = 'red', 
    #             alpha = 0.5)
    # ## Fit line
    # #p0 = [1e4, 1e-7]
    # #result = leastsq(residual_ab_exp, p0, args = (sizehist, sizeedges[:-1]))
    # #a,b = result[0]
    # #print a,b
    # ##a = 1e4
    # ##b = -0.4e-7
    # #axarr[0].plot(sizeedges[:-1],a*np.exp(b*sizeedges[:-1]), c = 'orange')
    # #print sizeedges[:-1], a*np.exp(b*sizeedges[:-1])
    # axarr[1,0].set_xlabel('Cloud size [m^2]')
    # axarr[1,0].set_ylabel('Number of clouds')
    # axarr[1,0].set_title('Cloud size', fontsize = 10)
    # #axarr[1,0].set_xlim([0., sizemax])
    # axarr[1,0].set_ylim([1, 1e6])
    # axarr[1,0].set_yscale('log')
    # axarr[1,0].set_xscale('log')
    
    # if not args.hypo:
    #     axarr[1,1].bar(sumedges[:-1], sumhist, width = np.diff(sumedges),
    #                 color = 'darkgray')
    #     axarr[1,1].plot([summean, summean], [1, 1e6], c = 'red', 
    #                 alpha = 0.5)
    #     axarr[1,1].set_ylabel('Number of clouds')
    #     #axarr[1,1].set_xlim([5e, summax])
    #     axarr[1,1].set_ylim([1, 1e6])
    #     axarr[1,1].set_yscale('log')
    #     axarr[1,1].set_xscale('log')
    #     axarr[1,1].set_xlabel('Cloud mass flux [kg/s]')
    #     axarr[1,1].set_title('Mass flux per cloud', fontsize = 10)
    
    # axarr[1,2].set_axis_off()
        
    titlestr = (alldatestr  + ', ' + args.ana + 
                ', water=' + str(args.water) +  
                ', nens=' + str(args.nens))
    fig.suptitle('Cloud size and mass flux per cloud distributions', fontsize=12)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    
    plotsavestr = ('cloud_stats_' + alldatestr + anastr)
    print 'Save as', plotdirsub + plotsavestr
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')



################################################################################
if 'rdf' in args.plot:
    print 'Plotting rdf'

    plotdirsub = plotdir +  '/rdf/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    rdf_list = []
    rdf_ns_list = []
    for d in args.date:
        dataset = Dataset(savedir + d + savesuf, 'r')
        rdf_tmp = dataset.variables['rdf'][:]
        rdf_tmp[rdf_tmp > 1e20] = np.nan
        rdf_list.append(rdf_tmp)
        rdf_ns_tmp = dataset.variables['rdf_nonscaled'][:]
        rdf_ns_tmp[rdf_ns_tmp > 1e20] = np.nan
        rdf_ns_list.append(rdf_ns_tmp)
    rdf = np.nanmean(rdf_list, axis = 0)
    rdf_ns = np.nanmean(rdf_ns_list, axis = 0)
    # Setup
    ymax = 7
    
    # Get 3 hr averages
    rdf_3hr = []
    rdf_ns_3hr = []
    tlist_3hr = []
    dt = int(3. / (args.tinc/60.))
    for i in range(len(timelist)/dt):
        rdf_3hr.append(np.nanmean(rdf[i*dt:(i+1)*dt], axis = 0))
        rdf_ns_3hr.append(np.nanmean(rdf_ns[i*dt:(i+1)*dt], axis = 0))
        tlist_3hr.append(timelist_plot[i*dt+2])
    rdf_3hr = np.array(rdf_3hr)
    rdf_ns_3hr = np.array(rdf_ns_3hr)
    cyc = [plt.cm.jet(i) for i in np.linspace(0, 1, len(tlist_3hr))]
    cyc = ("#FFE698","#BEE15B","#00D38B","#00B3C2","#0074D6","#9600AA")
    
    # Get the data
    r =   dataset.variables['dr'][:]
    
    fig, axarr = plt.subplots(1, 2, figsize = (pdfwidth, 3.5))
    
    ############# Time loop ##############
    for it, t in enumerate(tlist_3hr):
        #print 'time: ', t
        axarr[0].plot(r/1000., rdf_3hr[it, :], c = cyc[it], 
                label = str(t+1) + 'UTC pm 1h', linewidth = 1.5)
        axarr[1].plot(r/1000., rdf_ns_3hr[it, :], c = cyc[it], 
                label = str(t+1) + 'UTC pm 1h', linewidth = 1.5)
    
    axarr[0].legend(loc = 1, ncol = 1, prop={'size':10})
    axarr[0].plot([0, np.max(r)/1000.], [1, 1], c = 'gray', alpha = 0.5)
    axarr[0].set_xlabel('Distance [km]')
    axarr[0].set_ylabel('Normalized RDF')
    axarr[0].set_title('Radial distribution function')
    axarr[0].set_ylim(0, ymax)
    axarr[0].set_xlim(0, np.max(r)/1000.)
    
    axarr[1].set_xlabel('Distance [km]')
    axarr[1].set_ylabel('Non-normalized RDF')
    axarr[1].set_title('Radial distribution function (not scaled)')
    #axarr[1].set_ylim(0, ymax)
    axarr[1].set_xlim(0, np.max(r)/1000.)
    
    titlestr = (alldatestr + 
                ', nens=' + str(args.nens))
    #fig.suptitle(titlestr, fontsize='x-large')
    plt.tight_layout()
    
    plotsavestr = ('rdf_' + alldatestr + anastr)
    print 'Save as', plotdirsub + plotsavestr
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')
        

################################################################################
if 'prec_rdf' in args.plot:
    print 'Plotting rdf'

    plotdirsub = plotdir +  '/prec_rdf/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    rdf_model_list = []
    rdf_obs_list = []
    rdf_det_list = []
    for d in args.date:
        dataset = Dataset(savedir + d + savesuf, 'r')
        rdf_tmp = dataset.variables['rdf_prec_model'][:]
        rdf_tmp[rdf_tmp > 1e20] = np.nan
        rdf_model_list.append(rdf_tmp)
        rdf_tmp = dataset.variables['rdf_prec_obs'][:]
        rdf_tmp[rdf_tmp > 1e20] = np.nan
        rdf_obs_list.append(rdf_tmp)
        dataset_det = Dataset(savedir + d + anastr_det+ '.nc_det', 'r')
        rdf_tmp = dataset_det.variables['rdf_prec_model'][:]
        rdf_tmp[rdf_tmp > 1e20] = np.nan
        rdf_det_list.append(rdf_tmp)
    rdf_model = np.nanmean(rdf_model_list, axis = 0)
    rdf_obs = np.nanmean(rdf_obs_list, axis = 0)
    rdf_det = np.nanmean(rdf_det_list, axis = 0)
    r =   dataset.variables['dr'][:]
    
    # Setup
    ymin = 0.5
    ymax = 4
    
    
    # Get 3 hr averages
    rdf_3hr_model = []
    rdf_3hr_obs = []
    rdf_3hr_det = []
    tlist_3hr = []
    dt = int(3 * args.tinc/60.)
    for i in range(len(timelist)/3):
        print i
        rdf_3hr_model.append(np.nanmean(rdf_model[i*dt:(i+1)*dt], axis = 0))
        rdf_3hr_obs.append(np.nanmean(rdf_obs[i*dt:(i+1)*dt], axis = 0))
        rdf_3hr_det.append(np.nanmean(rdf_det[i*dt:(i+1)*dt], axis = 0))
        tlist_3hr.append(timelist_plot[i*dt+1])
    cyc = [plt.cm.jet(i) for i in np.linspace(0, 1, len(tlist_3hr))]
    #cyc = ("#FFE698","#BEE15B","#00D38B","#00B3C2","#0074D6","#9600AA")
    
    rdf_max_model = []
    rdf_max_obs = []
    rdf_max_det = []
    for i in range(len(timelist)):
        rdf_max_model.append(np.max(rdf_model[i]))
        rdf_max_obs.append(np.max(rdf_obs[i]))
        rdf_max_det.append(np.max(rdf_det[i]))
    
    fig, axarr = plt.subplots(2, 3, figsize = (pdfwidth, 7))
    ############# Time loop ##############
    for it, t in enumerate(tlist_3hr):
        # Normalized
        rdf_a = rdf_3hr_model[it] - 1
        max_a = np.max(rdf_a) 
        rdf_a = rdf_a / max_a
        axarr[1,1].plot(r/1000., rdf_a, c = cyc[it], 
                        linestyle = '-' , label = str(int(t+1)).zfill(2) + 'UTC pm 1h',
                        linewidth = 1.5)
        
        
        axarr[0,0].plot(r/1000., rdf_3hr_model[it], c = cyc[it], 
                        linestyle = '-' , label = str(t+1) + 'UTC pm 1h',
                        linewidth = 1.5)
        axarr[0,1].plot(r/1000., rdf_3hr_det[it], c = cyc[it], 
                        linestyle = '-' , label = str(int(t+1)).zfill(2) + 'UTC pm 1h',
                        linewidth = 1.5)
        axarr[0,2].plot(r/1000., rdf_3hr_obs[it], c = cyc[it], 
                        linestyle = '-' , label = str(int(t+1)).zfill(2) + 'UTC pm 1h',
                        linewidth = 1.5)
    
    axarr[0,2].legend(loc = 1, ncol = 1, prop={'size':10})
    axarr[0,0].plot([0, np.max(r)/1000.], [1, 1], c = 'gray', alpha = 0.5)
    axarr[0,0].set_xlabel('Distance [km]')
    axarr[0,0].set_ylabel('Normalized RDF')
    axarr[0,0].set_title('Ensemble with PSP', fontsize = 10)
    axarr[0,0].set_ylim(ymin, ymax)
    axarr[0,0].set_xlim(0, np.max(r)/1000.)
    
    axarr[0,1].plot([0, np.max(r)/1000.], [1, 1], c = 'gray', alpha = 0.5)
    axarr[0,1].set_xlabel('Distance [km]')
    #axarr[1].set_ylabel('Normalized RDF')
    axarr[0,1].set_title('Deterministic', fontsize = 10)
    axarr[0,1].set_ylim(ymin, ymax)
    axarr[0,1].set_xlim(0, np.max(r)/1000.)
    
    axarr[0,2].plot([0, np.max(r)/1000.], [1, 1], c = 'gray', alpha = 0.5)
    axarr[0,2].set_xlabel('Distance [km]')
    #axarr[1].set_ylabel('Normalized RDF')
    axarr[0,2].set_title('Observation', fontsize = 10)
    axarr[0,2].set_ylim(ymin, ymax)
    axarr[0,2].set_xlim(0, np.max(r)/1000.)
    
    axarr[1,0].plot(timelist_plot, rdf_max_model, label = 'Ens with PSP')
    axarr[1,0].plot(timelist_plot, rdf_max_det, label = 'Deterministic')
    axarr[1,0].plot(timelist_plot, rdf_max_obs, label = 'Obs')
    axarr[1,0].legend(prop={'size':8}, loc = 1)
    
    axarr[1,1].plot([0, np.max(r)/1000.], [0, 0], c = 'gray', alpha = 0.5)
    axarr[1,1].set_xlabel('Distance [km]')
    axarr[1,1].set_ylabel('Normalized RDF scaled by max')
    axarr[1,1].set_title('Ensemble with PSP', fontsize = 10)
    #axarr[1,1].set_ylim(ymin, ymax)
    axarr[1,1].set_xlim(0, np.max(r)/1000.)
    
    titlestr = (alldatestr + 
                ', nens=' + str(args.nens))
    fig.suptitle('Radial distribution function of precipitation fields', 
                 fontsize=12)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    
    plotsavestr = ('prec_rdf_' + alldatestr + anastr)
    print 'Plotting as:' , plotdirsub + plotsavestr
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')



################################################################################
if 'prec_hist' in args.plot:
    print 'Plotting prec_hist'

    plotdirsub = plotdir +  '/prec_hist/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    hist_model = []
    hist_obs = []
    hist_det = []
    for d in args.date:
        dataset = Dataset(savedir + d + savesuf, 'r')
        hist_model.append(dataset.variables['hist_model'][:])
        hist_obs.append(dataset.variables['hist_obs'][:])
        dataset_det = Dataset(savedir + d + anastr_det+ '.nc_det', 'r')
        hist_det.append(dataset_det.variables['hist_model'][:])
    hist_model = np.mean(hist_model, axis = 0)
    hist_obs = np.mean(hist_obs, axis = 0)
    hist_det = np.mean(hist_det, axis = 0)
    hist_model[0] = hist_model[0]/10.
    hist_obs[0] = hist_obs[0]/10.
    hist_det[0] = hist_det[0]/10.
    
    # Setup
    histbinedges = [0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 1000]
    x = np.arange(len(histbinedges)-1)

    fig, ax = plt.subplots(1, 1, figsize = (pdfwidth/2., 3.5))
    ax.bar(x[1:]+0.2, hist_model[1:], width = 0.2, color = '#d9d9d9', 
           label = 'Ensemble with PSP')
    ax.bar(x[1:]+0.4, hist_det[1:], width = 0.2, color = 'gray', 
           label = 'Deterministic')
    ax.bar(x[1:]+0.6, hist_obs[1:], width = 0.2, color = 'k', 
           label = 'Observations')
    
    ax.text(0.03,0.7, 'x10', transform=ax.transAxes, fontsize=10)
    
    ax.legend(loc = 1, prop={'size':10})
    ax.set_xlabel('Hourly accumulation [mm/h]')
    ax.set_ylabel('Number of grid points')
    ax.set_title('Precipitation histogram', fontsize=12)
    plt.xticks(x[1:], histbinedges[1:-1])
    
    titlestr = (alldatestr +
                ', nens=' + str(args.nens))
    #fig.suptitle(titlestr, fontsize='x-large')
    plt.tight_layout()
    
    plotsavestr = ('prec_hist_' + alldatestr + anastr)
    print 'Plotting as:', plotdirsub + plotsavestr
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')
        


################################################################################
if 'spectra' in args.plot:
    print 'Plotting spectra'

    plotdirsub = plotdir +  '/spectra/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    dke_spec_list = []
    bgke_spec_list = []
    dprec_spec_list = []
    bgprec_spec_list = []
    for d in args.date:
        dataset = Dataset(savedir + d + savesuf, 'r')
        dke_spec_list.append(dataset.variables['dkespec'][:])
        bgke_spec_list.append(dataset.variables['bgkespec'][:])
        dprec_spec_list.append(dataset.variables['dprecspec'][:])
        bgprec_spec_list.append(dataset.variables['bgprecspec'][:])
    dke_spec = np.nanmean(dke_spec_list, axis = 0)
    bgke_spec = np.nanmean(bgke_spec_list, axis = 0)
    dprec_spec = np.nanmean(dprec_spec_list, axis = 0)
    bgprec_spec = np.nanmean(bgprec_spec_list, axis = 0)

    
    cyc = [plt.cm.jet(i) for i in np.linspace(0, 1, len(timelist))]
    cyc = ("#E7A7FF","#FF84DB","#EF8974","#AF9300","#529324","#008768",
          "#006C88","#2D3184")

    speclam = dataset.variables['speclam'][:]
    
    fig, ax = plt.subplots(1, 2, figsize = (pdfwidth, 3.5))
    
    ############# Time loop ##############
    for it, t in enumerate(timelist):
        print 'time: ', t
        # Get ratio
        ratio = dke_spec[it]/bgke_spec[it]/2.
        ax[0].plot(speclam/1000., ratio, c = cyc[it], 
                label = str(int(timelist_plot[it])).zfill(2) + ' UTC',
                linewidth = 1.5)
        
        ratio = dprec_spec[it]/bgprec_spec[it]/2.
        ax[1].plot(speclam/1000., ratio, c = cyc[it], 
                label = str(int(timelist_plot[it])).zfill(2) + ' UTC',
                linewidth = 1.5)
    
    ax[0].legend(loc = 3, ncol = 2, prop={'size':8})
    ax[0].plot([5, 1000.], [1, 1], c = 'gray', alpha = 0.5)
    ax[0].set_xlabel('Wavelength [km]')
    ax[0].set_ylabel('Saturation ratio')
    ax[0].set_title("Saturation of KE spectrum")
    ax[0].set_ylim(1e-2, 2)
    ax[0].set_xlim(5, 1000.)
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    
    ax[1].plot([5, 1000.], [1, 1], c = 'gray', alpha = 0.5)
    ax[1].set_xlabel('Wavelength [km]')
    ax[1].set_title("Saturation of precipitation spectrum")
    ax[1].set_ylim(1e-2, 2)
    ax[1].set_xlim(5, 1000.)
    ax[1].set_yscale('log')
    ax[1].set_xscale('log')
    
    titlestr = (alldatestr + 
                ', nens=' + str(args.nens))
    #fig.suptitle(titlestr, fontsize='x-large')
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    
    plotsavestr = ('spectra_' + alldatestr + anastr)
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')

################################################################################
if 'prec_spec' in args.plot:
    print 'Plotting prec_spec'

    plotdirsub = plotdir +  '/prec_spec/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    dprec_spec_list = []
    bgprec_spec_list = []
    for d in args.date:
        dataset = Dataset(savedir + d + savesuf, 'r')
        dprec_spec_list.append(dataset.variables['dprecspec'][:])
        bgprec_spec_list.append(dataset.variables['bgprecspec'][:])
    dprec_spec = np.nanmean(dprec_spec_list, axis = 0)
    bgprec_spec = np.nanmean(bgprec_spec_list, axis = 0)

    cyc = [plt.cm.jet(i) for i in np.linspace(0, 1, len(timelist))]

    speclam = dataset.variables['speclam'][:]
    
    fig, ax = plt.subplots(1, 1, figsize = (pdfwidth/2., 3.5))
    
    ############# Time loop ##############
    for it, t in enumerate(timelist):
        print 'time: ', t
        # Get ratio
        ratio = dprec_spec[it]/bgprec_spec[it]/2.
        ax.plot(speclam/1000., ratio, c = cyc[it], 
                label = str(int(timelist_plot[it])).zfill(2) + ' UTC',
                linewidth = 1.5)
    
    ax.legend(loc = 3, ncol = 2, prop={'size':8})
    ax.plot([5, 1000.], [1, 1], c = 'gray', alpha = 0.5)
    ax.set_xlabel('Wavelength [km]')
    ax.set_ylabel('Saturation ratio')
    ax.set_title("Saturation of precipitation spectrum")
    ax.set_ylim(1e-2, 2)
    ax.set_xlim(5, 1000.)
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    titlestr = (alldatestr + 
                ', nens=' + str(args.nens))
    #fig.suptitle(titlestr, fontsize='x-large')
    plt.tight_layout()
    
    plotsavestr = ('prec_spec_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water)  +
                    '_nens-' + str(args.nens) + '_tstart-' + 
                    str(args.tstart) + '_tend-' + str(args.tend) + 
                    '_tinc-' + str(args.tinc))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')

            
################################################################################
if 'scatter' in args.plot:
    print 'Plotting scatter'
    plotdirsub = plotdir +  '/scatter/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)

    # Setup
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    
    
    # Load the data 
    # These lists have 4 dimensions [time, lev, n, N_box]
    varM_list = create2Dlist(len(nlist), len(timelist))
    meanM_list = create2Dlist(len(nlist), len(timelist))
    varN_list = create2Dlist(len(nlist), len(timelist))
    meanN_list = create2Dlist(len(nlist), len(timelist))
    varm_list = create2Dlist(len(nlist), len(timelist))
    meanm_list = create2Dlist(len(nlist), len(timelist))
    
    # loop over dates 
    for d in args.date:
        # Load dataset 
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        
        # Load the required data and put it in list
        for it, time in enumerate(timelist):
            for i_n, n in enumerate(dataset.variables['n']):
                nmax = 265/n
                
                meanM = dataset.variables['meanM'][it,i_n,:nmax,:nmax]
                varM = dataset.variables['varM'][it,i_n,:nmax,:nmax]
                meanN = dataset.variables['meanN'][it,i_n,:nmax,:nmax]
                varN = dataset.variables['varN'][it,i_n,:nmax,:nmax]
                meanm = dataset.variables['meanm'][it,i_n,:nmax,:nmax]
                varm = dataset.variables['varm'][it,i_n,:nmax,:nmax]
                
                meanM_list[i_n][it] += list(np.ravel(meanM))
                varM_list[i_n][it] += list(np.ravel(varM))
                meanN_list[i_n][it] += list(np.ravel(meanN))
                varN_list[i_n][it] += list(np.ravel(varN))
                meanm_list[i_n][it] += list(np.ravel(meanm))
                varm_list[i_n][it] += list(np.ravel(varm))

                
    # now I have the lists I want in the scatter plot 
        
    # Set up the figure 
    fig, axarr = plt.subplots(5, 3, figsize = (pdfwidth, 12))

    rmselist1 = []
    rmselist2 = []
    rmselist3 = []
    rmselist4 = []
    rmselist5 = []
    ####### n loop #######################
    for i_n, n in enumerate(dataset.variables['n']):
        print 'n: ', n
        z_n = 0.1 + (n/256.)*0.1   # for z_order
        
        # Extract the arrays
        M = np.array(meanM_list[i_n])
        m = np.array(meanm_list[i_n])
        N = np.array(meanN_list[i_n])
        varM = np.array(varM_list[i_n])
        varm = np.array(varm_list[i_n])
        varN = np.array(varN_list[i_n])
        alpha = varN/N
        beta = varm/(m**2)
        # These now have dimensions [time, date]        
        
        # 1. Raw PC08 prediction
        m_c = 4e7
        predict = 2*m_c*M
        
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist1.append(nrmse)
        
        axarr[0,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[0,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[0,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[0,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        
        axarr[0,2].scatter([0], [0], marker = 'o', c = clist[i_n], 
                            s = 8, label = str(n*2.8)+'km',
                            linewidth = 0)  # Ghost plot for labels
        
        # 2. CC06 prediction
        print 'mean m', np.nanmean(m)
        predict = 2*m*M
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist2.append(nrmse)
        
        axarr[1,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[1,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[1,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[1,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        var_time = np.nanmean(m, axis = 1)
        axarr[1,2].plot(timelist_plot, var_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        
        # 3. alpha adjusted
        predict = (1+alpha)*m*M
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist3.append(nrmse)
        
        axarr[2,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[2,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[2,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[2,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        var_time = np.nanmean(alpha, axis = 1)
        axarr[2,2].plot(timelist_plot, var_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        
        # 4. beta adjusted
        predict = (1+beta)*m*M
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist4.append(nrmse)
        
        axarr[3,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[3,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[3,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[3,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        var_time = np.nanmean(beta, axis = 1)
        axarr[3,2].plot(timelist_plot, var_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        
        # 5. all
        predict = (alpha+beta)*m*M
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist5.append(nrmse)
        
        axarr[4,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[4,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[4,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[4,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        #var_time = np.nanmean(m, axis = 1)
        #axarr[1,2].plot(timelist_plot, var_time, c = clist[i_n], 
                        #label = str(n*2.8)+'km', linewidth = 1.5)
        
        

    
    # Complete the figure
    #axarr[0,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[0,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[0,0].set_xlim(1e6,1e11)
    axarr[0,0].set_ylim(0, 2.5)
    axarr[0,0].set_xscale('log')
    #axarr[0,0].set_yscale('log')
    axarr[0,0].set_xlabel(r'$M$')
    axarr[0,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / 2m_c\langle M \rangle$')
    axarr[0,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[0,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[0,1].set_ylim(0, 2.5)
    axarr[0,1].set_xlabel('time [h/UTC]')
    axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    axarr[0,2].legend(loc =3, ncol = 2, prop={'size':8})
    #axarr[0,2].set_xlabel('time [h/UTC]')
    axarr[0,2].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    #axarr[1,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[1,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[1,0].set_xlim(1e6,1e11)
    axarr[1,0].set_ylim(0, 2.5)
    axarr[1,0].set_xscale('log')
    axarr[1,0].set_xlabel(r'$M$')
    axarr[1,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / 2\langle m \rangle \langle M \rangle$')
    axarr[1,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[1,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[1,1].set_ylim(0, 2.5)
    axarr[1,1].set_xlabel('time [h/UTC]')
    axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    axarr[1,2].set_ylabel(r'$\langle m \rangle$ [kg/s]')
    axarr[1,2].set_xlabel('time [h/UTC]')
    axarr[1,2].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    #axarr[2,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[2,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[2,0].set_xlim(1e6,1e11)
    axarr[2,0].set_ylim(0, 2.5)
    axarr[2,0].set_xscale('log')
    axarr[2,0].set_xlabel(r'$M$')
    axarr[2,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / (1+\alpha)\langle m \rangle \langle M \rangle$')
    axarr[2,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[2,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[2,1].set_ylim(0, 2.5)
    axarr[2,1].set_xlabel('time [h/UTC]')
    axarr[2,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    
    axarr[2,2].set_ylim(0, 2.5)
    axarr[2,2].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[2,2].set_ylabel(r'$\alpha$')
    axarr[2,2].set_xlabel('time [h/UTC]')
    axarr[2,2].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    #axarr[3,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[3,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[3,0].set_xlim(1e6,1e11)
    axarr[3,0].set_ylim(0, 2.5)
    axarr[3,0].set_xscale('log')
    axarr[3,0].set_xlabel(r'$M$')
    axarr[3,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / (1+\beta)\langle m \rangle \langle M \rangle$')
    axarr[3,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[3,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[3,1].set_ylim(0, 2.5)
    axarr[3,1].set_xlabel('time [h/UTC]')
    axarr[3,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    axarr[3,2].set_ylim(0, 2.5)
    axarr[3,2].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[3,2].set_ylabel(r'$\beta$')
    axarr[3,2].set_xlabel('time [h/UTC]')
    axarr[3,2].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    #axarr[4,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[4,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[4,0].set_xlim(1e6,1e11)
    axarr[4,0].set_ylim(0, 2.5)
    axarr[4,0].set_xscale('log')
    axarr[4,0].set_xlabel(r'$M$')
    axarr[4,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / (\alpha+\beta)\langle m \rangle \langle M \rangle$')
    axarr[4,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[4,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[4,1].set_ylim(0, 2.5)
    axarr[4,1].set_xlabel('time [h/UTC]')
    axarr[4,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    for ax, let in zip(list(np.ravel(axarr)),
                       ['a','b','c','d','e','f','g','h','i','j','l','m','n','o','p']):
        ax.text(0.05, 0.9, '('+let+')', transform = ax.transAxes, 
                    fontsize = 10)
    

    #titlestr = (alldatestr + '\n' + args.ana + 
                #', water=' + str(args.water) + ', lev= ' + str(int(args.height[0])) + 
                #', nens=' + str(args.nens))
    #fig.suptitle(titlestr, fontsize='x-large')
    fig.suptitle('Comparison of simulation results and predictions \nof increasing complexity', fontsize=12)
    plt.tight_layout(rect=[0, 0.0, 1, 0.93])
    
    plotsavestr = ('scatter_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                    '_nens-' + str(args.nens))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)

    plt.close('all')
    
    
    ##3
    fig, ax = plt.subplots(1, 1, figsize = (pdfwidth/2., 4))
    x = np.arange(5)
    width = 0.1
    
    for i_n, n in enumerate(nlist):
        plt.bar(1+i_n*width, rmselist1[i_n], width, color = clist[i_n])
        diff2 = rmselist1[i_n] - rmselist2[i_n]
        diff3 = rmselist1[i_n] - rmselist3[i_n]
        diff4 = rmselist1[i_n] - rmselist4[i_n]
        diff5 = rmselist1[i_n] - rmselist5[i_n]
        plt.bar(2+i_n*width, -diff2, width, color = clist[i_n])
        plt.bar(3+i_n*width, -diff3, width, color = clist[i_n])
        plt.bar(4+i_n*width, -diff4, width, color = clist[i_n])
        plt.bar(5+i_n*width, -diff5, width, color = clist[i_n])
        plt.bar(2+i_n*width, rmselist2[i_n], width, color = clist[i_n])
        plt.bar(3+i_n*width, rmselist3[i_n], width, color = clist[i_n])
        plt.bar(4+i_n*width, rmselist4[i_n], width, color = clist[i_n])
        plt.bar(5+i_n*width, rmselist5[i_n], width, color = clist[i_n])
    #plt.scatter(np.log2(nlist)-0.15, rmselist1, c = clist)
    #plt.scatter(np.log2(nlist)-0.05, rmselist2, c = clist)
    #plt.scatter(np.log2(nlist)+0.05, rmselist3, c = clist)
    #plt.scatter(np.log2(nlist)+0.15, rmselist4, c = clist)
    plotsavestr = ('rmse_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                    '_nens-' + str(args.nens))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)

    plt.close('all')
        

################################################################################
if 'diurnal' in args.plot:
    print 'Plotting diurnal'

    plotdirsub = plotdir +  '/diurnal/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    dx = 2.8e3
    # Setup
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    
    savename = savedir + 'diurnal_' + alldatestr + anastr
    if os.path.exists(savename):
        print 'Loading pre-saved file', savename
        # numpy load
        savefile = open(savename, 'r')
        meanM_list = cPickle.load(savefile)
        varM_list = cPickle.load(savefile)
        meanN_list = cPickle.load(savefile)
        varN_list = cPickle.load(savefile)
        meanm_list = cPickle.load(savefile)
        varm_list = cPickle.load(savefile)
        savefile.close()
            
    else:
        print 'No pre-saved file found', savename
    
        # Load the data 
        # These lists have 4 dimensions [time, lev, n, N_box]
        varM_list = create2Dlist(len(nlist), len(timelist))
        meanM_list = create2Dlist(len(nlist), len(timelist))
        varN_list = create2Dlist(len(nlist), len(timelist))
        meanN_list = create2Dlist(len(nlist), len(timelist))
        varm_list = create2Dlist(len(nlist), len(timelist))
        meanm_list = create2Dlist(len(nlist), len(timelist))
        
        # loop over dates 
        for d in args.date:
            # Load dataset 
            print 'Loading date: ', d
            dataset = Dataset(savedir + d + savesuf, 'r')
            
            # Load the required data and put it in list
            for it, time in enumerate(timelist):
                for i_n, n in enumerate(dataset.variables['n']):
                    nmax = 256/n
                    
                    meanM = dataset.variables['meanM'][it,i_n,:nmax,:nmax]
                    varM = dataset.variables['varM'][it,i_n,:nmax,:nmax]
                    meanN = dataset.variables['meanN'][it,i_n,:nmax,:nmax]
                    varN = dataset.variables['varN'][it,i_n,:nmax,:nmax]
                    meanm = dataset.variables['meanm'][it,i_n,:nmax,:nmax]
                    varm = dataset.variables['varm'][it,i_n,:nmax,:nmax]
                    
                    meanM_list[i_n][it] += list(np.ravel(meanM))
                    varM_list[i_n][it] += list(np.ravel(varM))
                    meanN_list[i_n][it] += list(np.ravel(meanN))
                    varN_list[i_n][it] += list(np.ravel(varN))
                    meanm_list[i_n][it] += list(np.ravel(meanm))
                    varm_list[i_n][it] += list(np.ravel(varm))
                    
        # now I have the lists I want in the scatter plot 
        savefile = open(savename, 'w')
        cPickle.dump(meanM_list, savefile, -1)
        cPickle.dump(varM_list, savefile, -1)
        cPickle.dump(meanN_list, savefile, -1)
        cPickle.dump(varN_list, savefile, -1)
        cPickle.dump(meanm_list, savefile, -1)
        cPickle.dump(varm_list, savefile, -1)
        savefile.close()
        
    timelist3h = [1, 4, 7, 10, 13, 16]
    
    pname = ['CC06','alpha','CC06alpha', 'beta']
    plabel = ['Variance: simulation / prediction',r'$\alpha$',
              'Variance: simulation / prediction', r'$\beta$','']
    ptitle = ['CC06 prediction',r'Organization parameter $\alpha$',
              r'$\alpha$-adjusted CC06 prediction',
              r'$m$ distribution parameter $\beta$','']
    pylim = [(0.45, 2.1),(0.7,3),(0.45, 2.1),(0.45, 2.1),(0,1)]
    
    for ip in range(len(pname)):
        print 'ip', ip
        # Set up the figure
        fig, ax = plt.subplots(1, 1, figsize = (pdfwidth/2., 3.5))
        clist = ['#3366ff', '#009933', '#ff3300']
        labellist = ['SS: 11.2 km', 'MS: 89.6 km', 'LS: 717 km']
        for i, i_n in enumerate([6,3,0]):
            print 'n', i_n
            # Extract the arrays
            M = np.array(meanM_list[i_n])
            m = np.array(meanm_list[i_n])
            N = np.array(meanN_list[i_n])
            varM = np.array(varM_list[i_n])
            varm = np.array(varm_list[i_n])
            varN = np.array(varN_list[i_n])
            alpha = varN/N
            beta = varm/(m**2)
                        
            
            if ip == 0:
                # CC06 unadjusted
                predict = 2*m*M
                frac = varM/predict
            if ip == 1:
                # Alpha
                frac = alpha
            if ip ==2:
                predict = (1+alpha)*m*M
                frac = varM/predict
            if ip == 3:
                frac = beta
                print frac.shape
            #if ip == 4:
                #tmp = []
                #for it in range(N.shape[0]):
                    #Ntmp = N[it]
                    #mtmp = m[it]
                    #mask = np.isfinite(Ntmp) * np.isfinite(mtmp)
                    
                    #tmp.append([np.corrcoef(Ntmp[mask], mtmp[mask])[1,0]])
                    
                #frac = np.array(tmp)
                #print frac.shape
                #print frac
            # Get 3 hourly means and percentiles
            mean3h = []
            per5 = []
            per25 = []
            per75 = []
            per95 = []
            for it in timelist3h:
                start = (it-1)*2; stop = (it+2)*2
                data = np.ravel(frac[start:stop])
                data = data[np.isfinite(data)]
                #print data
                mean3h.append(np.mean(data))
                per5.append(np.percentile(data, 5))
                per25.append(np.percentile(data, 25))
                per75.append(np.percentile(data, 75))
                per95.append(np.percentile(data, 95))
            
            left = np.array(timelist3h)+5.75+i/2.
            mid = left + 0.25
            height = np.array(per75)-np.array(per25)
            #ax.bar(left, height, 0.5, per25, color = clist[i], linewidth = 0,
                #label = labellist[i])
            ax.plot(mid, mean3h, c = clist[i], label = labellist[i], 
                    linewidth = 2)
            #ax.
            ax.scatter(mid, mean3h, c = clist[i], zorder = 10, s = 20)
            for j in range(len(mid)):
                ax.plot([mid[j],mid[j]],[per25[j],per75[j]], c = clist[i], 
                        zorder = 3, linewidth = 0.7)
                ax.plot([mid[j]-0.1, mid[j]+0.1],[per25[j],per25[j]], c = clist[i],
                        linewidth = 0.7, zorder = 3)
                ax.plot([mid[j]-0.1, mid[j]+0.1],[per75[j],per75[j]], c = clist[i],
                        linewidth = 0.7, zorder = 3)
            
        ax.plot([6,24], [1,1], c = 'gray', zorder = 0.5, label = 'perfect')
        ax.legend(loc = 2, fontsize = 8)
        ax.set_xticks([6,9,12,15,18,21,24])
        ax.set_xticklabels([6,9,12,15,18,21,24])
        ax.set_xlabel('Time [h/UTC]')
        ax.set_ylabel(plabel[ip])
        ax.set_xlim(6,24)
        
        ax.set_yscale('log')
        ax.set_yticks([0.5,1,2])
        ax.set_yticklabels(['.5', '1', '2'])
        ax.set_ylim(pylim[ip])
        ax.set_title(ptitle[ip])
        plt.tight_layout()
        plotsavestr = (pname[ip] + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                        '_nens-' + str(args.nens))
        fig.savefig(plotdirsub + plotsavestr, dpi = 300)

        plt.close('all')

################################################################################
if 'std_v_mean' in args.plot:
    print 'Plotting std_v_mean'

    plotdirsub = plotdir +  '/std_v_mean/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    dx = 2.8e3
    # Setup
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    
    savename = savedir + 'std_v_mean_' + alldatestr + anastr
    if os.path.exists(savename):
        print 'Loading pre-saved file', savename
        # numpy load
        savefile = open(savename, 'r')
        stdM_list = cPickle.load(savefile)
        M_list = cPickle.load(savefile)
        stdQmp_list = cPickle.load(savefile)
        Qmp_list = cPickle.load(savefile)
        savefile.close()
            
    else:
        print 'No pre-saved file found', savename
    
        # Load the data: I want to plot mu_2, N, alpha
        # These list have 3 dims [lev, n, N_box]

        stdM_list = create1Dlist(len(nlist))
        M_list = create1Dlist(len(nlist))
        stdQmp_list = create1Dlist(len(nlist))
        Qmp_list = create1Dlist(len(nlist))
        
        # Loop over dates 
        for d in args.date:
            # Load dataset 
            print 'Loading date: ', d
            dataset = Dataset(savedir + d + savesuf, 'r')
            
            # Load the required data and put it in list
            for i_n, n in enumerate(dataset.variables['n']):
                nmax = 265/n
                meanM = dataset.variables['meanM'][:,i_n,:nmax,:nmax]
                varM = dataset.variables['varM'][:,i_n,:nmax,:nmax]
                meanQmp = dataset.variables['meanQmp'][:,i_n,:nmax,:nmax]
                varQmp = dataset.variables['varQmp'][:,i_n,:nmax,:nmax]
                
                stdM_list[i_n] += list(np.ravel(np.sqrt(varM)))
                M_list[i_n] += list(np.ravel(meanM))
                stdQmp_list[i_n] += list(np.ravel(np.sqrt(varQmp)))
                Qmp_list[i_n] += list(np.ravel(meanQmp))
                    
        # now I have the lists I want in the scatter plot 
        savefile = open(savename, 'w')
        cPickle.dump(stdM_list, savefile, -1)
        cPickle.dump(M_list, savefile, -1)
        cPickle.dump(stdQmp_list, savefile, -1)
        cPickle.dump(Qmp_list, savefile, -1)
        savefile.close()

    
    # Set up the figure 
    tmp = np.logspace(5,12, 1000)
    fig, ax = plt.subplots(1, 1, figsize = (pdfwidth/2., 3.5))
    
    blist = []
    allstdM = []
    allM = []
    ####### n loop #######################
    for i_n, n in enumerate(nlist):
        print 'n: ', n
        z_n = 0.1 + (n/256.)*0.1   # for z_order
        
        # Get the data
        stdM = np.array(stdM_list[i_n])
        M = np.array(M_list[i_n])
        p0 = [1]
        y = np.array(stdM)
        x = np.array(M)
        mask = np.isfinite(y)
        # Fit the line, CC06
        result = leastsq(residual_b_sqrt, p0, args = (y[mask], x[mask]))
        b = result[0]
        blist.append(b)
        
        allstdM += list(stdM)
        allM += list(M)

        
    allstdM = np.array(allstdM)
    allM = np.array(allM)
    
    p0 = [1]
    y = np.array(allstdM)
    x = np.array(allM)
    mask = np.isfinite(y)
    # Fit the line, CC06
    result = leastsq(residual_b_sqrt, p0, args = (y[mask], x[mask]))
    b = result[0]
    tmp2 = np.logspace(6, 11, 1000)
    ax.plot(tmp2,np.sqrt(b*tmp2), alpha = 1, c = 'orangered',
                    linestyle = '-', zorder = 0.2, label = r'CC06: $y=\sqrt{bx}$')
    # Fit the line, SPPT
    result = leastsq(residual_b, p0, args = (y[mask], x[mask]))
    b = result[0]
    ax.plot(tmp2,b*tmp2, alpha = 1, color = 'cornflowerblue',
                    linestyle = '--', zorder = 0.2, label = r'SPPT: $y = bx$')
    
    
    # Bin the data
    nbins = 10
    binedges = np.logspace(6, 11, nbins+1)
    bininds = np.digitize(allM, binedges)
    binmeans = []
    bin5 = []
    bin25 = []
    bin75 = []
    bin95 = []
    binnum = []
    for i in range(1,nbins+1):
        num = (allstdM[bininds==i]).shape[0]
        if num == 0:
            binmeans.append(np.nan)
            bin25.append(np.nan)
            bin75.append(np.nan)
            bin5.append(np.nan)
            bin95.append(np.nan)
        else:
            binmeans.append(np.average(allstdM[bininds==i]))
            bin25.append(np.percentile(allstdM[bininds==i],25))
            bin75.append(np.percentile(allstdM[bininds==i],75))
            bin5.append(np.percentile(allstdM[bininds==i],5))
            bin95.append(np.percentile(allstdM[bininds==i],95))
        binnum.append(num / float((allstdM[np.isfinite(allstdM)]).shape[0]) * 50)
    xmean = (binedges[:-1] + binedges[1:])/2.
    logmean = np.exp((np.log(binedges[:-1]) + np.log(binedges[1:])) / 2.)
    logleft = np.exp(np.log(binedges[:-1])+0.2)
    logright = np.exp(np.log(binedges[1:])-0.2)
    height = np.array(bin75)-np.array(bin25)
    width = logright - logleft
    ax.bar(logleft, height, width, bin25, linewidth = 0, color = 'gray')
    for i, binmean in enumerate(binmeans):
        ax.plot([logleft[i], logright[i]], [binmean, binmean],
                color = 'black', zorder = 2)
        ax.plot([logmean[i], logmean[i]], [bin5[i], bin95[i]], 
                color = 'black', zorder = 0.5)
        
        
        
        
        

    ax.set_xlim(1e6,1e11)
    ax.set_ylim(4e6,2e9)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'$\langle M \rangle$ [kg/s]')
    ax.set_ylabel(r'$\langle (\delta M)^2 \rangle^{1/2}$ [kg/s]')

    
    ax.set_title('Scaling of standard deviation with mean')

    ax.legend(loc = 2, ncol = 1, prop={'size':8})
    
    ax2 = plt.axes([.65, .3, .2, .2], axisbg='lightgray')
    blist = np.array(blist) / 1.e8
    x = np.array(nlist)*2.8
    ax2.plot(x, blist, c = 'orangered')
    ax2.scatter(x, blist, c = 'orangered')
    ax2.set_title('Slope of CC06 fit', fontsize = 8)
    ax2.set_ylabel(r'$b\times 10^8$', fontsize = 8, labelpad = 0.05)
    ax2.set_xlabel('n [km]', fontsize = 8, labelpad = 0.07)
    ax2.set_xscale('log')
    print x
    ax2.set_xlim(5, 1000)
    ax2.set_xticks([5, 50, 500])
    ax2.set_xticklabels([5, 50, 500], fontsize = 8)
    ax2.set_yticks([0.5, 1.5])
    ax2.set_yticklabels([0.5, 1.5], fontsize = 8)

    plt.tight_layout()
    
    plotsavestr = ('std_v_mean_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                    '_nens-' + str(args.nens))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')

    #### Q ##########################3333333
    # Set up the figure 
    tmp = np.logspace(5,12, 1000)
    fig, ax = plt.subplots(1, 1, figsize = (pdfwidth/2., 3.5))
    
    blist = []
    allstdM = []
    allM = []
    ####### n loop #######################
    for i_n, n in enumerate(nlist):
        print 'n: ', n
        z_n = 0.1 + (n/256.)*0.1   # for z_order
        
        # Get the data
        stdM = np.array(stdQmp_list[i_n])*((n*2.8e3)**2)
        M = np.array(Qmp_list[i_n])*((n*2.8e3)**2)
        #print M
        p0 = [1]
        y = np.array(stdM)
        x = np.array(M)
        mask = np.isfinite(y)
        # Fit the line, CC06
        result = leastsq(residual_b_sqrt, p0, args = (y[mask], x[mask]))
        b = result[0]
        blist.append(b)
        
        allstdM += list(stdM)
        allM += list(M)

        
    allstdM = np.array(allstdM)
    allM = np.array(allM)
    
    p0 = [1]
    y = np.array(allstdM)
    x = np.array(allM)
    mask = np.isfinite(y) * x>0
    # Fit the line, CC06
    result = leastsq(residual_b_sqrt, p0, args = (y[mask], x[mask]))
    b = result[0]
    print b
    tmp2 = np.logspace(2, 7, 1000)
    ax.plot(tmp2,np.sqrt(b*tmp2), alpha = 1, c = 'orangered',
                    linestyle = '-', zorder = 0.2, label = r'CC06: $y=\sqrt{bx}$')
    # Fit the line, SPPT
    result = leastsq(residual_b, p0, args = (y[mask], x[mask]))
    b = result[0]
    ax.plot(tmp2,b*tmp2, alpha = 1, color = 'cornflowerblue',
                    linestyle = '--', zorder = 0.2, label = r'SPPT: $y = bx$')
    
    
    # Bin the data
    nbins = 10
    binedges = np.logspace(2, 7, nbins+1)
    bininds = np.digitize(allM, binedges)
    binmeans = []
    bin5 = []
    bin25 = []
    bin75 = []
    bin95 = []
    binnum = []
    for i in range(1,nbins+1):
        num = (allstdM[bininds==i]).shape[0]
        if num == 0:
            binmeans.append(np.nan)
            bin25.append(np.nan)
            bin75.append(np.nan)
            bin5.append(np.nan)
            bin95.append(np.nan)
        else:
            binmeans.append(np.average(allstdM[bininds==i]))
            bin25.append(np.percentile(allstdM[bininds==i],25))
            bin75.append(np.percentile(allstdM[bininds==i],75))
            bin5.append(np.percentile(allstdM[bininds==i],5))
            bin95.append(np.percentile(allstdM[bininds==i],95))
        binnum.append(num / float((allstdM[np.isfinite(allstdM)]).shape[0]) * 50)
    xmean = (binedges[:-1] + binedges[1:])/2.
    logmean = np.exp((np.log(binedges[:-1]) + np.log(binedges[1:])) / 2.)
    logleft = np.exp(np.log(binedges[:-1])+0.2)
    logright = np.exp(np.log(binedges[1:])-0.2)
    height = np.array(bin75)-np.array(bin25)
    width = logright - logleft
    ax.bar(logleft, height, width, bin25, linewidth = 0, color = 'gray')
    for i, binmean in enumerate(binmeans):
        ax.plot([logleft[i], logright[i]], [binmean, binmean],
                color = 'black', zorder = 2)
        ax.plot([logmean[i], logmean[i]], [bin5[i], bin95[i]], 
                color = 'black', zorder = 0.5)
        
        
        
        
        

    ax.set_xlim(1e2,1e7)
    ax.set_ylim(1e2,5e6)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'$\langle Q \rangle * A$ [kg/s * m^2]')
    ax.set_ylabel(r'$\langle (\delta Q)^2 \rangle^{1/2} * A$ [kg/s * m^2]')

    
    ax.set_title('Scaling of standard deviation with mean')

    ax.legend(loc = 2, ncol = 1, prop={'size':8})


    plt.tight_layout()
    
    plotsavestr = ('std_v_mean_Q_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                    '_nens-' + str(args.nens))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')



################################################################################
if 'frac_v_tauc' in args.plot:
    print 'Plotting frac_v_tauc'

    plotdirsub = plotdir +  '/frac_v_tauc/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    dx = 2.8e3
    # Setup
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    
    savename = savedir + 'frac_v_tauc_' + alldatestr + anastr
    if os.path.exists(savename):
        print 'Loading pre-saved file', savename
        # numpy load
        savefile = open(savename, 'r')
        varM_list = cPickle.load(savefile)
        M_list = cPickle.load(savefile)
        m_list = cPickle.load(savefile)
        tauc_list = cPickle.load(savefile)
        savefile.close()
            
    else:
        print 'No pre-saved file found', savename
    
        # Load the data:
        # These list have 2 dims [n, N_box]

        varM_list = create1Dlist(len(nlist))
        M_list = create1Dlist(len(nlist))
        m_list = create1Dlist(len(nlist))
        tauc_list = create1Dlist(len(nlist))

        
        # Loop over dates 
        for d in args.date:
            # Load dataset 
            print 'Loading date: ', d
            dataset = Dataset(savedir + d + savesuf, 'r')
            
            # Load the required data and put it in list
            for i_n, n in enumerate(dataset.variables['n']):
                nmax = 256/n
                meanM = dataset.variables['meanM'][:,i_n,:nmax,:nmax]
                varM = dataset.variables['varM'][:,i_n,:nmax,:nmax]
                meanm = dataset.variables['meanm'][:,i_n,:nmax,:nmax]
                tauc = dataset.variables['meantauc'][:,i_n,:nmax,:nmax]
                
                varM_list[i_n] += list(np.ravel(varM))
                M_list[i_n] += list(np.ravel(meanM))
                m_list[i_n] += list(np.ravel(meanm))
                tauc_list[i_n] += list(np.ravel(tauc))

                    
        # now I have the lists I want in the scatter plot 
        savefile = open(savename, 'w')
        cPickle.dump(varM_list, savefile, -1)
        cPickle.dump(M_list, savefile, -1)
        cPickle.dump(m_list, savefile, -1)
        cPickle.dump(tauc_list, savefile, -1)
        savefile.close()

    
    # Set up the figure 
    tmp = np.logspace(5,12, 1000)
    fig, ax = plt.subplots(1, 1, figsize = (pdfwidth/2., 3.5))
    
    blist = []
    allstdM = []
    allM = []
    ####### n loop #######################
    for i_n, n in enumerate(nlist):
        print 'n: ', n
        z_n = 0.1 + (n/256.)*0.1   # for z_order
        
        # Get the data
        varM = np.array(varM_list[i_n])
        M = np.array(M_list[i_n])
        m = np.array(m_list[i_n])
        tauc = np.array(tauc_list[i_n])
        
        frac = varM/(2*m*M)
        mask2 = M > 1.1 * m / 50.
        mask = np.isfinite(tauc) * np.isfinite(frac)
        mask = mask & mask2
        Rcorr = np.corrcoef(tauc[mask], frac[mask])[1,0]
        print Rcorr
        ax.scatter(tauc[mask], frac[mask], c = clist[i_n], zorder = z_n, 
                   linewidth = 0, s = z_n**2*100,
                   label = str(n) + ' corr: {:.2f}'.format(Rcorr))
        


    ax.set_xlim(0.01,1e3)
    ax.set_ylim(0.1,10)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('tau_c')
    ax.set_ylabel('variance ratio')

    
    ax.set_title('Timescale versus variance ratio')

    ax.legend(loc = 3, ncol = 3, prop={'size':6})
    

    plt.tight_layout()
    
    plotsavestr = ('frac_v_tauc_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                    '_nens-' + str(args.nens))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')



################################################################################
if 'geographical' in args.plot:
    print 'Plotting diurnal'

    plotdirsub = plotdir +  '/geographical/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    dx = 2.8e3
    # Setup

    
    savename = savedir + 'geographical_' + alldatestr + anastr
    if os.path.exists(savename):
        print 'Loading pre-saved file', savename
        # numpy load
        savefile = open(savename, 'r')
        array2Dlist = cPickle.load(savefile)
        savefile.close()
            
    else:
        print 'No pre-saved file found', savename
    
        outerlist = []
        
        for i_n, n in enumerate(nlist):
            outerlist.append([])
            nmax = 256/n
            for d in args.date:
                # Load dataset 
                print 'Loading date: ', d
                dataset = Dataset(savedir + d + savesuf, 'r')
            
                for it, time in enumerate(timelist):
                    
                    meanM = dataset.variables['meanM'][it,i_n,:nmax,:nmax]
                    varM = dataset.variables['varM'][it,i_n,:nmax,:nmax]
                    meanm = dataset.variables['meanm'][it,i_n,:nmax,:nmax]
                    
                    frac = varM/(2*meanm*meanM)
                    outerlist[-1].append(frac)
                    
        array2Dlist = []
        for innerlist in outerlist:
            array2Dlist.append(np.nanmean(innerlist, axis = 0))
        
        print array2Dlist

                    
        # now I have the lists I want in the scatter plot 
        savefile = open(savename, 'w')
        cPickle.dump(array2Dlist, savefile, -1)
        savefile.close()
        
    # Load a dummy fobj   
    import copy
    fn = '/home/scratch/users/stephan.rasp/2016052800/deout_ceu_pspens/1/OUTPUT/lfff00000000.nc_30m_surf'
    fobj = getfobj_ncdf(fn, fieldn = 'CAPE_ML')
    lx1 = (357-256-1)/2
    ly1 = (357-256-1)/2
    
    pllevels = np.exp(np.linspace(np.log(0.75), np.log(1.5), 20))
    
    
    for n, array2D in zip(nlist, array2Dlist):
        # Upscale field
        upfield = np.ones((fobj.ny, fobj.nx)) * np.nan
        nx, ny = array2D.shape
        for i in range(nx):
            for j in range(ny):
                # Get limits for each N box
                xmin = i*n+lx1
                xmax = (i+1)*n+lx1
                ymin = j*n+ly1
                ymax = (j+1)*n+ly1
                
                upfield[xmin:xmax, ymin:ymax] = array2D[i,j]
        
        newfobj = copy.deepcopy(fobj)
        newfobj.data = upfield
        newfobj.dims = 2
        newfobj.fieldn = 'Variance fraction'
        newfobj.unit = ''
        
        
        fig = fig_contourf_1sp(newfobj, pllevels = pllevels, extend = 'both')
    
        plotsavestr = ('geographical_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                        '_nens-' + str(args.nens) + '_n-' + str(n))
        fig.savefig(plotdirsub + plotsavestr, dpi = 300)
        plt.close('all')


################################################################################
if 'N_v_m' in args.plot:
    print 'Plotting N_v_m'

    plotdirsub = plotdir +  '/N_v_m/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    dx = 2.8e3
    # Setup
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    
    savename = savedir + 'N_v_m_' + alldatestr + anastr
    if os.path.exists(savename):
        print 'Loading pre-saved file', savename
        # numpy load
        savefile = open(savename, 'r')
        m_list = cPickle.load(savefile)
        N_list = cPickle.load(savefile)
        savefile.close()
            
    else:
        print 'No pre-saved file found', savename
    
        # Load the data: I want to plot mu_2, N, alpha
        # These list have 3 dims [lev, n, N_box]

        m_list = create1Dlist(len(nlist))
        N_list = create1Dlist(len(nlist))
        
        # Loop over dates 
        for d in args.date:
            # Load dataset 
            print 'Loading date: ', d
            dataset = Dataset(savedir + d + savesuf, 'r')
            
            # Load the required data and put it in list
            for i_n, n in enumerate(dataset.variables['n']):
                nmax = 265/n
                A = float((n*2.8e3)**2)
                meanm = dataset.variables['meanm'][:,i_n,:nmax,:nmax]
                meanN = dataset.variables['meanN'][:,i_n,:nmax,:nmax]
                
                m_list[i_n] += list(np.ravel(meanm))
                N_list[i_n] += list(np.ravel(meanN))
                    
        # now I have the lists I want in the scatter plot 
        savefile = open(savename, 'w')
        cPickle.dump(m_list, savefile, -1)
        cPickle.dump(N_list, savefile, -1)
        savefile.close()

    
    # Set up the figure 
    tmp = np.logspace(5,12, 1000)
    fig, ax = plt.subplots(1, 1, figsize = (pdfwidth/2., 3.5))
    
    blist = []
    allstdM = []
    allM = []
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    ####### n loop #######################
    for i_n, n in enumerate(nlist):
        print 'n: ', n
        z_n = 0.1 + (n/256.)*0.1   # for z_order
        
        # Get the data
        m = np.array(m_list[i_n])
        N = np.array(N_list[i_n])
        mask = np.isfinite(m) * np.isfinite(N)
        Rcorr = np.corrcoef(N[mask], m[mask])[1,0]
        print Rcorr
        ax.scatter(N, m, c = clist[i_n], zorder = z_n)

    
    #ax.set_xlim(1e-12, 1e-8)
    ax.set_xscale('log')
    
    plotsavestr = ('N_v_m_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                    '_nens-' + str(args.nens))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)

    plt.close('all')
    
    


################################################################################
if 'summary_stats' in args.plot:
    print 'Plotting summary_stats'
    plotdirsub = plotdir +  '/summary_stats/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    # Setup
    cyc = [plt.cm.bone(i) for i in np.linspace(0.1, 0.9, len(args.date))]

    ####### date loop #######################
    compM_list = []
    compm_list = []
    compsize_list = []
    compm_list_tot = []
    compsize_list_tot = []
    compN_list = []
    #compQ_list = []
    for d in args.date:
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        tmp1 = dataset.variables['cld_size'][:,:]
        tmp1[tmp1 > 1e20] = np.nan
        compsize_list.append(np.nanmean(tmp1, axis = 1))
        compsize_list_tot.append(tmp1)
        
        tmp2 = dataset.variables['cld_sum'][:,:]
        tmp2[tmp2 > 1e20] = np.nan
        compm_list.append(np.nanmean(tmp2, axis = 1))
        compM_list.append(np.nansum(tmp2, axis = 1)/float(args.nens))
        compm_list_tot.append(tmp2)

        
        #tmp3 = dataset.variables['meanQmp'][:,0,0,0]
        #compQ_list.append(tmp3)
        
        tmp4 = dataset.variables['totN'][:]
        tmp4[tmp4 > 1e20] = np.nan
        compN_list.append(tmp4)

    
    # Get the composite means
    #compsize = np.nanmean(np.array(compsize_list), axis = 0)
    #compm = np.nanmean(np.array(compm_list), axis = 0)
    #compstdm = np.nanstd(np.array(compm_list), axis = 0)
    #print compm, compstdm
    compM = np.nanmean(np.array(compM_list), axis = 0)
    #compQ = np.nanmean(np.array(compQ_list), axis = 0)
    compN = np.nanmean(np.array(compN_list), axis = 0)
    
    compm = []
    compstdm = []
    compsize = []
    for it in range(compsize_list_tot[0].shape[0]):
        tmp1 = []
        tmp2 = []
        for arr1, arr2 in zip(compsize_list_tot, compm_list_tot):
            tmp1.append(arr1[it])
            tmp2.append(arr2[it])
        compsize.append(np.nanmean(tmp1))
        compm.append(np.nanmean(tmp2))
        compstdm.append(np.nanstd(tmp2))
        
    
    
    timelist = [timedelta(seconds=ts) for ts in dataset.variables['time']]
    timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
    # Create the figure
    fig, axarr = plt.subplots(2, 2, figsize = (pdfwidth, 7))
    
    axarr[0,0].plot(timelist_plot, compm, c = 'orangered', linewidth = 2)
    #axarr[0,0].errorbar(timelist_plot, compm, yerr = compstdm, zorder = 0.1,
                        #)
    for ic, yplot in enumerate(compm_list):
        axarr[0,0].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    axarr[0,0].set_xlabel('time [h/UTC]')
    axarr[0,0].set_xlim(timelist_plot[0], timelist_plot[-1])
    axarr[0,0].set_ylabel('Mean cloud mass flux [kg/s]')
    
    axarr[0,1].plot(timelist_plot, compsize, c = 'orangered', linewidth = 2)
    for ic, yplot in enumerate(compsize_list):
        axarr[0,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    axarr[0,1].set_xlabel('time [h/UTC]')
    axarr[0,1].set_ylabel('Mean cloud size [m^2]')
    axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    axarr[1,0].plot(timelist_plot, compM, c = 'orangered', linewidth = 2)
    for ic, yplot in enumerate(compM_list):
        axarr[1,0].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    axarr[1,0].set_xlabel('time [h/UTC]')
    axarr[1,0].set_xlim(timelist_plot[0], timelist_plot[-1])
    axarr[1,0].set_ylabel('Domain total mass flux [kg/s]')
    
    
    axarr[1,1].plot(timelist_plot, compN, c = 'orangered', linewidth = 2)
    for ic, yplot in enumerate(compN_list):
        axarr[1,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    axarr[1,1].set_xlabel('time [h/UTC]')
    axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    axarr[1,1].set_ylabel('Mean N')
    
    #axarr[1,1].plot(timelist_plot, compQ, c = 'orangered', linewidth = 2)
    #for ic, yplot in enumerate(compQ_list):
        #axarr[1,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    #axarr[1,1].set_xlabel('time [h/UTC]')
    #axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    #axarr[1,1].set_ylabel('Mean Q [h]')
    
    
    fig.suptitle('Temporal evolution of domain averaged quantities', fontsize = 12.)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    
    plotsavestr = ('summary_stats_' + alldatestr + anastr)
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')
            
            
################################################################################
if 'summary_var' in args.plot:
    print 'Plotting summary_var'
    plotdirsub = plotdir +  '/summary_var/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    # Setup
    timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    
    
    # Load the data 
    # These lists have 4 dimensions [time, lev, n, N_box]
    varM_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    meanM_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    varN_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    meanN_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    varm_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    meanm_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    
    # loop over dates 
    for d in args.date:
        # Load dataset 
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        
        # Load the required data and put it in list
        for it, time in enumerate(timelist):
            for iz, lev in enumerate(dataset.variables['levs']):
                for i_n, n in enumerate(dataset.variables['n']):
                    nmax = 265/n
                    
                    meanM = dataset.variables['meanM'][it,iz,i_n,:nmax,:nmax]
                    varM = dataset.variables['varM'][it,iz,i_n,:nmax,:nmax]
                    meanN = dataset.variables['meanN'][it,iz,i_n,:nmax,:nmax]
                    varN = dataset.variables['varN'][it,iz,i_n,:nmax,:nmax]
                    meanm = dataset.variables['meanm'][it,iz,i_n,:nmax,:nmax]
                    varm = dataset.variables['varm'][it,iz,i_n,:nmax,:nmax]
                    
                    meanM_list[it][iz][i_n] += list(np.ravel(meanM))
                    varM_list[it][iz][i_n] += list(np.ravel(varM))
                    meanN_list[it][iz][i_n] += list(np.ravel(meanN))
                    varN_list[it][iz][i_n] += list(np.ravel(varN))
                    meanm_list[it][iz][i_n] += list(np.ravel(meanm))
                    varm_list[it][iz][i_n] += list(np.ravel(varm))
    

    mu2N_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    mu2Nalpha_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    mu2Nbeta_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    mu2Nab_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    alpha_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    beta_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    
    # Get means and convert to plotting arrays
    for it, time in enumerate(timelist):
        for iz, lev in enumerate(dataset.variables['levs']):
            for i_n, n in enumerate(dataset.variables['n']):
                # Extract the arrays
                varM = np.array(varM_list[it][iz][i_n])
                meanM = np.array(meanM_list[it][iz][i_n])
                varN = np.array(varN_list[it][iz][i_n])
                meanN = np.array(meanN_list[it][iz][i_n])
                varm = np.array(varm_list[it][iz][i_n])
                meanm = np.array(meanm_list[it][iz][i_n])
                
                mu2 = varM/(meanM**2)
                alpha = varN/meanN
                beta = varm/(meanm**2)
                
                mu2N_list3d[it][iz][i_n] = np.nanmean(mu2*meanN/2)
                mu2Nalpha_list3d[it][iz][i_n] = np.nanmean(mu2*meanN/(1+alpha))
                mu2Nbeta_list3d[it][iz][i_n] = np.nanmean(mu2*meanN/(1+beta))
                mu2Nab_list3d[it][iz][i_n] = np.nanmean(mu2*meanN/(beta+alpha))
                alpha_list3d[it][iz][i_n] = np.nanmean(alpha)
                beta_list3d[it][iz][i_n] = np.nanmean(beta)
    
    mu2N_list3d =  np.array(mu2N_list3d)
    mu2Nalpha_list3d =  np.array(mu2Nalpha_list3d)
    mu2Nbeta_list3d =  np.array(mu2Nbeta_list3d)
    mu2Nab_list3d =  np.array(mu2Nab_list3d)
    alpha_list3d =  np.array(alpha_list3d)
    beta_list3d =  np.array(beta_list3d)
    
    ######## Lev loop #####################
    for iz, lev in enumerate(dataset.variables['levs']):
        print 'lev: ', lev
        
        # Create the figure
        fig, axarr = plt.subplots(3, 2, figsize = (95./25.4*2, 10))
        
        for i_n, n in enumerate(nlist):
            axarr[0,0].plot(timelist_plot, mu2N_list3d[:,iz,i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km', linewidth = 1.5)
            axarr[0,1].plot(timelist_plot, mu2Nab_list3d[:,iz,i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km', linewidth = 1.5)
            axarr[1,0].plot(timelist_plot, mu2Nalpha_list3d[:,iz,i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km', linewidth = 1.5)
            axarr[1,1].plot(timelist_plot, mu2Nbeta_list3d[:,iz,i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km', linewidth = 1.5)
            axarr[2,0].plot(timelist_plot, alpha_list3d[:,iz,i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km', linewidth = 1.5)
            axarr[2,1].plot(timelist_plot, beta_list3d[:,iz,i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km', linewidth = 1.5)
            
        axarr[0,0].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[0,0].set_ylim(0, 2)
        #axarr[0,0].set_yscale('log')
        axarr[0,0].set_xlabel('time [h/UTC]')
        axarr[0,0].set_ylabel(r'$\mu_2 \langle N \rangle/2$')
        axarr[0,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[0,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[0,1].set_ylim(0, 2)
        #axarr[0,1].set_yscale('log')
        axarr[0,1].set_xlabel('time [h/UTC]')
        axarr[0,1].set_ylabel(r'$\mu_2 \langle N \rangle/(\alpha +\beta)$')
        axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[1,0].plot(timelist_plot, [1]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[1,0].set_ylim(0, 2)
        #axarr[1,0].set_yscale('log')
        axarr[1,0].set_xlabel('time [h/UTC]')
        axarr[1,0].set_ylabel(r'$\mu_2 \langle N \rangle/(1+\alpha)$')
        axarr[1,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[1,1].plot(timelist_plot, [1]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[1,1].set_ylim(0, 2)
        #axarr[1,1].set_yscale('log')
        axarr[1,1].set_xlabel('time [h/UTC]')
        axarr[1,1].set_ylabel(r'$\mu_2 \langle N \rangle/(1+\beta)$')
        axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[2,0].plot(timelist_plot, [1]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[2,0].set_ylim(0, 2)
        #axarr[2,0].set_yscale('log')
        axarr[2,0].set_xlabel('time [h/UTC]')
        axarr[2,0].set_ylabel(r'$\alpha$')
        axarr[2,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[2,1].plot(timelist_plot, [1]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[2,1].set_ylim(0, 2)
        #axarr[2,1].set_yscale('log')
        axarr[2,1].set_xlabel('time [h/UTC]')
        axarr[2,1].set_ylabel(r'$\beta$')
        axarr[2,1].set_xlim(timelist_plot[0], timelist_plot[-1])
            
        axarr[2,1].legend(loc =3, ncol = 2, prop={'size':6})
              
              
        titlestr = (alldatestr + '\n' + args.ana + 
                    ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                    ', nens=' + str(args.nens))
        fig.suptitle(titlestr, fontsize='x-large')
        plt.tight_layout(rect=[0, 0.0, 1, 0.93])
        
        plotsavestr = ('summary_var_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
                        '_nens-' + str(args.nens)+ '_tstart-' + 
                        str(args.tstart) + '_tend-' + str(args.tend) + 
                        '_tinc-' + str(args.tinc))
        fig.savefig(plotdirsub + plotsavestr, dpi = 300)
        plt.close('all')


################################################################################
if 'summary_weather' in args.plot:
    print 'Plotting summary_weather'
    plotdirsub = plotdir +  '/summary_weather/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    # Setup
    cyc = [plt.cm.bone(i) for i in np.linspace(0.1, 0.9, len(args.date))]


    ####### date loop #######################
    comphpbl_list = []
    compcape_list = []
    compprec_list = []
    comptauc_list = []
    for d in args.date:
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        
        hpbl_tmp = dataset.variables['dihpbl'][:]
        cape_tmp = dataset.variables['dicape'][:]
        prec_tmp = dataset.variables['diprec'][:]
        hpbl_tmp[hpbl_tmp > 1e10] = np.nan
        cape_tmp[cape_tmp > 1e10] = np.nan
        prec_tmp[prec_tmp > 1e10] = np.nan
        comphpbl_list.append(hpbl_tmp)
        compcape_list.append(cape_tmp)
        compprec_list.append(prec_tmp)
        
        tauc_tmp = dataset.variables['ditauc'][:]
        tauc_tmp[tauc_tmp > 1e10] = np.nan
        comptauc_list.append(tauc_tmp)
    
    # Get the composite means
    comphpbl = np.nanmean(np.array(comphpbl_list), axis = 0)
    compcape = np.nanmean(np.array(compcape_list), axis = 0)
    compprec = np.nanmean(np.array(compprec_list), axis = 0)
    comptauc = np.nanmean(np.array(comptauc_list), axis = 0)
    
    timelist = [timedelta(seconds=ts) for ts in dataset.variables['time']]
    timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
    
    # Create the figure
    fig, axarr = plt.subplots(2, 2, figsize = (pdfwidth, 7.))
    
    axarr[0,0].plot(timelist_plot, compprec, c = 'orangered', linewidth = 2)
    for ic, yplot in enumerate(compprec_list):
        axarr[0,0].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    axarr[0,0].set_xlabel('time [h/UTC]')
    axarr[0,0].set_xlim(timelist_plot[0], timelist_plot[-1])
    axarr[0,0].set_ylabel('Precipitation [mm/h]')
    
    axarr[0,1].plot(timelist_plot, compcape, c = 'orangered', linewidth = 2)
    for ic, yplot in enumerate(compcape_list):
        axarr[0,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    axarr[0,1].set_xlabel('time [h/UTC]')
    axarr[0,1].set_ylabel('CAPE [J/kg]')
    axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    axarr[1,0].plot(timelist_plot, comptauc, c = 'orangered', linewidth = 2)
    for ic, yplot in enumerate(comptauc_list):
        axarr[1,0].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    axarr[1,0].set_xlabel('time [h/UTC]')
    axarr[1,0].set_xlim(timelist_plot[0], timelist_plot[-1])
    axarr[1,0].set_ylabel('Convective timescale [h]')
    
    axarr[1,1].plot(timelist_plot, comphpbl, c = 'orangered', linewidth = 2)
    for ic, yplot in enumerate(comphpbl_list):
        axarr[1,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    axarr[1,1].set_xlabel('time [h/UTC]')
    axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    axarr[1,1].set_ylabel('PBL height [m]')
    
    titlestr = (alldatestr + ', ' + args.ana + 
                ', water=' + str(args.water) + ', height= ' + heightstr + 
                ', nens=' + str(args.nens))
    #fig.suptitle(titlestr, fontsize='x-large')
    fig.suptitle('Temporal evolution of domain averaged quantities', fontsize = 12.)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    
    plotsavestr = ('summary_weather_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) +
                    '_nens-' + str(args.nens) + '_tstart-' + 
                    str(args.tstart) + '_tend-' + str(args.tend) + 
                    '_tinc-' + str(args.tinc))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')
            


################################################################################
if 'prec_comp' in args.plot:
    print 'Plotting prec_comp'

    plotdirsub = plotdir +  '/prec_comp/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    preclist_obs = []
    preclist_det = []
    preclist_psp = []
    for d in args.date:
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        # ugly filename fix becuase det alsways has nens = 1
        dataset_det = Dataset(savedir + d + anastr_det+ '.nc_det', 'r')
        
        preclist_obs.append(dataset.variables['prec_mask_obs'][:])
        preclist_psp.append(dataset.variables['prec_mask_model'][:])
        preclist_det.append(dataset_det.variables['prec_mask_model'][:])

    obs_mean = np.nanmean(np.array(preclist_obs), axis = 0)
    det_mean = np.nanmean(np.array(preclist_det), axis = 0)
    psp_mean = np.nanmean(np.array(preclist_psp), axis = 0)

    timelist = [timedelta(seconds=ts) for ts in dataset.variables['time']]
    timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
    
    # Create the figure
    fig, axarr = plt.subplots(1, 3, figsize = (pdfwidth, 3))
    cyc = [plt.cm.bone(i) for i in np.linspace(0.1, 0.9, len(args.date))]
    axarr[0].plot(timelist_plot, psp_mean, c = 'orangered', linewidth = 2)
    for ic, yplot in enumerate(preclist_psp):
        axarr[0].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
    axarr[0].set_xlabel('time [h/UTC]')
    axarr[0].set_xlim(timelist_plot[0], timelist_plot[-1])
    axarr[0].set_ylabel('Precipitation [mm/h]')
    axarr[0].set_title('Ensemble mean with PSP')
    
    axarr[1].plot(timelist_plot, det_mean-psp_mean, c = 'orangered', 
                  linewidth = 2)
    axarr[1].plot(timelist_plot, [0]*len(timelist_plot), c = 'k', 
                  linewidth = 0.5, zorder = 0.1)
    for ic, yplot in enumerate(preclist_det):
        yplotdiff = yplot - preclist_psp[ic]
        axarr[1].plot(timelist_plot, yplotdiff, zorder = 0.5, c = cyc[ic])
    axarr[1].set_xlabel('time [h/UTC]')
    axarr[1].set_xlim(timelist_plot[0], timelist_plot[-1])
    axarr[1].set_ylabel('Prec diff det - psp [mm/h]')
    axarr[1].set_title('Difference det - psp')
    
    axarr[2].plot(timelist_plot, obs_mean-psp_mean, c = 'orangered', 
                  linewidth = 2)
    axarr[2].plot(timelist_plot, [0]*len(timelist_plot), c = 'k', 
                  linewidth = 0.5, zorder = 0.1)
    for ic, yplot in enumerate(preclist_obs):
        yplotdiff = yplot - preclist_psp[ic]
        axarr[2].plot(timelist_plot, yplotdiff, zorder = 0.5, c = cyc[ic])
    axarr[2].set_xlabel('time [h/UTC]')
    axarr[2].set_xlim(timelist_plot[0], timelist_plot[-1])
    axarr[2].set_ylabel('Prec diff obs - psp [mm/h]')
    axarr[2].set_title('Difference obs - psp')
    
    plt.tight_layout(rect=[0, 0.0, 1, 1])
    
    plotsavestr = ('prec_comp_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) +
                    '_nens-' + str(args.nens) + '_tstart-' + 
                    str(args.tstart) + '_tend-' + str(args.tend) + 
                    '_tinc-' + str(args.tinc))
    print 'Plotting as:' , plotdirsub + plotsavestr
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')
    

################################################################################
if 'stamps_var' in args.plot:
    print 'Plotting stamps_var'

    plotdirsub = plotdir +  '/stamps_var/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    
    # Load the data
    # dataset loop
    NvarMN_list = []
    M_list = []
    Mdiff_list = []
    enstauc_list = []
    for d in args.date:
        # Load dataset 
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        varM = dataset.variables['varM'][:]
        varN = dataset.variables['varN'][:]
        varm = dataset.variables['varm'][:]
        meanM = dataset.variables['meanM'][:]
        meanN = dataset.variables['meanN'][:]
        meanm = dataset.variables['meanm'][:]
        Mmem1 = dataset.variables['Mmem1'][:]
        # Calculate the statistics
        NvarMN = varM / (meanM**2) * meanN
        varNoN = varN / meanN
        M_list.append(meanM)
        Mdiff_list.append((Mmem1-meanM)/meanM)
        NvarMN_list.append(NvarMN)
        enstauc_list.append(dataset.variables['enstauc'][:])
    
    # Get dataset mean
    NvarMN = np.nanmean(NvarMN_list, axis = 0)
    M = np.nanmean(NvarMN_list, axis = 0)
    Mdiff = np.nanmean(Mdiff_list, axis = 0)
    enstauc = np.nanmean(enstauc_list, axis = 0)

    
    # Plot setup
    cmM = ("#0030C4","#3F5BB6","#7380C0","#A0A7CE","#CCCEDC","#DDCACD",
          "#D09AA4","#BF6A7D","#AA3656","#920031")
    levelsM = np.array([1, 1.14, 1.33, 1.6, 1.78, 2, 2.25, 2.5, 3, 3.5, 4])/2.
    
    ############# Time loop ##############
    for it, t in enumerate(timelist):
        print 'time: ', t
        
        ######## Lev loop #####################
        for iz, lev in enumerate(dataset.variables['levs']):
            print 'lev: ', lev
            
            
            ####### n loop #######################
            for i_n, n in enumerate(dataset.variables['n']):
                print 'n: ', n
                
                # Do upscaling and create fobjs
                # 0. Load one fobj to get parameters
                tmpfobj = getfobj_ncdf(ensdir + '1/OUTPUT/lfff00000000c.nc_30m',
                                       fieldn = 'HSURF')
                # Crop all fields to analysis domain
                sxo, syo = tmpfobj.data.shape  # Original field shape
                lx1 = (sxo-256-1)/2 # ATTENTION first dimension is actually y
                lx2 = -(lx1+1) # Number of grid pts to exclude at border
                ly1 = (syo-256-1)/2
                ly2 = -(ly1+1)
                rlats = tmpfobj.rlats[lx1:lx2,0]
                rlons = tmpfobj.rlons[0,ly1:ly2]
                # 1. tauc 
                taucobj = fieldobj(data = enstauc[it],
                                   fieldn = 'TAU_C',
                                   rlats = rlats,
                                   rlons = rlons,
                                   polelat = tmpfobj.polelat,
                                   polelon = tmpfobj.polelon,
                                   levs_inp = 'surf',
                                   unit = 'h')
                
                # Map fields to original grid
                NvarMN_map = np.empty((taucobj.ny, taucobj.nx))
                M_map = np.empty((taucobj.ny, taucobj.nx))
                Mdiff_map = np.empty((taucobj.ny, taucobj.nx))
                for i in range(256/n):
                    for j in range(256/n):
                        # Get limits for each N box
                        xmin = i*n
                        xmax = (i+1)*n
                        ymin = j*n
                        ymax = (j+1)*n
                        
                        NvarMN_map[xmin:xmax, ymin:ymax] = NvarMN[it,iz,i_n,i,j]
                        #print NvarMN_map[xmin:xmax, ymin:ymax]
                        Mdiff_map[xmin:xmax, ymin:ymax] = Mdiff[it,iz,i_n,i,j]
                        M_map[xmin:xmax, ymin:ymax] = M[it,iz,i_n,i,j]
                
                # Set missing values to nans
                M_map[M_map == 0.] = np.nan
                Mdiff_map[Mdiff_map == 0.] = np.nan
                NvarMN_map[NvarMN_map == 0.] = np.nan
                
                
                # 2. NvarMN
                NvarMNobj = fieldobj(data = NvarMN_map/2.,
                                   fieldn = '0.5 * NvarMN',
                                   rlats = rlats,
                                   rlons = rlons,
                                   polelat = tmpfobj.polelat,
                                   polelon = tmpfobj.polelon,
                                   levs_inp = lev,
                                   unit = '')
                
                # 2. varNoN
                Mobj = fieldobj(data = M_map,
                                   fieldn = 'meanM',
                                   rlats = rlats,
                                   rlons = rlons,
                                   polelat = tmpfobj.polelat,
                                   polelon = tmpfobj.polelon,
                                   levs_inp = lev,
                                   unit = '')
                
                # 3. NvarMN_adj
                Mdiffobj = fieldobj(data = Mdiff_map,
                                   fieldn = 'Mdiff',
                                   rlats = rlats,
                                   rlons = rlons,
                                   polelat = tmpfobj.polelat,
                                   polelon = tmpfobj.polelon,
                                   levs_inp = lev,
                                   unit = '')
            
                # Set up the figure 
                fig, axarr = plt.subplots(2, 2, figsize = (9, 8))
                
                # Plot
                # 1. tauc
                plt.sca(axarr[0,0])   # This is necessary for some reason...
                cf, tmp = ax_contourf(axarr[0,0], taucobj, 
                                      pllevels=np.arange(0, 21, 1), 
                                      sp_title=taucobj.fieldn,
                                      Basemap_drawrivers = False,
                                      Basemap_parallelslabels = [0,0,0,0],
                                      Basemap_meridiansslabels = [0,0,0,0],
                                      extend = 'both')
                cb = fig.colorbar(cf)
                cb.set_label(taucobj.unit)
                
                # 2. NvarMN
                plt.sca(axarr[0,1])
                cf, tmp = ax_contourf(axarr[0,1], NvarMNobj, 
                                      pllevels=levelsM, colors = cmM,
                                      sp_title=NvarMNobj.fieldn,
                                      Basemap_drawrivers = False,
                                      Basemap_parallelslabels = [0,0,0,0],
                                      Basemap_meridiansslabels = [0,0,0,0],
                                      extend = 'both')
                cb = fig.colorbar(cf)
                cb.set_label(NvarMNobj.unit)
                
                # 3. varNoN
                plt.sca(axarr[1,0])
                cf, tmp = ax_contourf(axarr[1,0], Mobj, 
                                      pllevels=levelsM,  colors = cmM,
                                      sp_title=Mobj.fieldn,
                                      Basemap_drawrivers = False,
                                      Basemap_parallelslabels = [0,0,0,0],
                                      Basemap_meridiansslabels = [0,0,0,0],
                                      extend = 'both')
                cb = fig.colorbar(cf)
                cb.set_label(Mobj.unit)
                
                # 4. NvarMN_adj
                plt.sca(axarr[1,1])
                cf, tmp = ax_contourf(axarr[1,1], Mdiffobj, 
                                      pllevels=np.arange(-1, 1.1, 0.2),  colors = cmM,
                                      sp_title=Mdiffobj.fieldn,
                                      Basemap_drawrivers = False,
                                      Basemap_parallelslabels = [0,0,0,0],
                                      Basemap_meridiansslabels = [0,0,0,0],
                                      extend = 'both')
                cb = fig.colorbar(cf)
                cb.set_label(Mdiffobj.unit)
                
                titlestr = (alldatestr + '+' + ddhhmmss(t) + ', ' + args.ana + 
                            ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                            ', nens=' + str(args.nens) +  ',  n=' + str(n))
                fig.suptitle(titlestr, fontsize='x-large')
                plt.tight_layout(rect=[0, 0.0, 1, 0.95])
                
                plotsavestr = ('new_stamps_var_' + alldatestr + '_ana-' + args.ana + 
                            '_wat-' + str(args.water) + '_lev-' + str(lev) +
                            '_nens-' + str(args.nens) + '_time-' + ddhhmmss(t) +
                            '_n-' + str(n))
                fig.savefig(plotdirsub + plotsavestr, dpi = 300)
                plt.close('all')
            

################################################################################
if 'stamps_w' in args.plot:
    print 'Plotting stamps_w'
    if len(datasetlist) > 1:
        raise Exception, 'More than one date is not implemented'
    dataset = Dataset(savedir + d + savesuf, 'r')
    plotdirsub = plotdir +  '/stamps_w/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    # Plot setup
    cmW = ("#7C0607","#903334","#A45657","#BA7B7C","#FFFFFF",
            "#8688BA","#6567AA","#46499F","#1F28A2")
    levelsW = np.array([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
    
    cmPREC = ((1    , 1     , 1    ), 
        (0    , 0.627 , 1    ),
        (0.137, 0.235 , 0.98 ),
        (0.392, 0     , 0.627),
        (0.784, 0     , 0.627),
        (1    , 0.3   , 0.9  ) )
    levelsPREC = [0, 0.1, 0.3, 1, 3, 10, 30]
    
    ############# Time loop ##############
    for it, t in enumerate(timelist):
        print 'time: ', t
        
        ######## Lev loop #####################
        for iz, lev in enumerate(dataset.variables['levs']):
            print 'lev: ', lev
            
            
            # Load Prec ens data
            ncdffn = 'lfff' + ddhhmmss(t) + '.nc_30m_surf'
            precobjlist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                        dir_suffix='/OUTPUT/', fieldn = 'PREC_ACCUM', 
                                        nfill=1)
            # Get mean 
            meanprec, tmp = mean_spread_fieldobjlist(precobjlist)

            # Load W ens data
            ncdffn = 'lfff' + ddhhmmss(t) + '.nc_30m'
            wobjlist = getfobj_ncdf_ens(ensdir, 'sub', args.nens, ncdffn, 
                                        dir_suffix='/OUTPUT/', fieldn = 'W', 
                                        nfill=1, 
                                        levs = dataset.variables['levs'][:])
            
            # Crop all fields to analysis domain
            sxo, syo = wobjlist[0].data.shape[1:]  # Original field shape
            lx1 = (sxo-256-1)/2 # ATTENTION first dimension is actually y
            lx2 = -(lx1+1) # Number of grid pts to exclude at border
            ly1 = (syo-256-1)/2
            ly2 = -(ly1+1)

            
            for wobj in wobjlist:
                wobj.data = wobj.data[:, lx1:lx2, ly1:ly2]
                wobj.rlats = wobj.rlats[lx1:lx2, ly1:ly2]
                wobj.rlons = wobj.rlons[lx1:lx2, ly1:ly2]
                wobj.lats = wobj.lats[lx1:lx2, ly1:ly2]
                wobj.lons = wobj.lons[lx1:lx2, ly1:ly2]
            
            meanprec.data = meanprec.data[lx1:lx2, ly1:ly2]
            meanprec.rlats = meanprec.rlats[lx1:lx2, ly1:ly2]
            meanprec.rlons = meanprec.rlons[lx1:lx2, ly1:ly2]
            meanprec.lats = meanprec.lats[lx1:lx2, ly1:ly2]
            meanprec.lons = meanprec.lons[lx1:lx2, ly1:ly2]
            
            # Set up the figure 
            fig, axarr = plt.subplots(2, 2, figsize = (9, 8))
            
            # Plot
            # 1. mean W
            plt.sca(axarr[0,0])   # This is necessary for some reason...
            cf, tmp = ax_contourf(axarr[0,0], meanprec, 
                                    pllevels=levelsPREC, 
                                    colors = cmPREC,
                                    sp_title=meanprec.fieldn,
                                    Basemap_drawrivers = False,
                                    Basemap_parallelslabels = [0,0,0,0],
                                    Basemap_meridiansslabels = [0,0,0,0],
                                    extend = 'both')
            cb = fig.colorbar(cf)
            cb.set_label(meanprec.unit)
            
            for ax, fob, i, in zip(list(np.ravel(axarr)[1:]), wobjlist[:3],
                                    range(1,4)):
            
                # 2. NvarMN
                plt.sca(ax)
                cf, tmp = ax_contourf(ax, fob, 
                                    pllevels=levelsW, colors = cmW,
                                    sp_title=fob.fieldn + str(i),
                                    Basemap_drawrivers = False,
                                    Basemap_parallelslabels = [0,0,0,0],
                                    Basemap_meridiansslabels = [0,0,0,0],
                                    extend = 'both', lev = iz)
                cb = fig.colorbar(cf)
                cb.set_label(fob.unit)
            
            
            titlestr = (args.date[0] + '+' + ddhhmmss(t) + ', ' + args.ana + 
                        ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                        ', nens=' + str(args.nens))
            fig.suptitle(titlestr, fontsize='x-large')
            plt.tight_layout(rect=[0, 0.0, 1, 0.95])
            
            plotsavestr = ('stamps_w_' + args.date[0] + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
                        '_nens-' + str(args.nens) + '_time-' + ddhhmmss(t))
            fig.savefig(plotdirsub + plotsavestr, dpi = 300)
            plt.close('all')


################################################################################
if 'height_var' in args.plot:
    print 'Plotting height_var'
    plotdirsub = plotdir +  '/height_var/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    # Setup
    t1 = int(args.tplot[0]); t2 = int(args.tplot[1])
    timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
    UTCstart = timelist[t1]
    UTCstop = timelist[t2-1]
    print 'Summing from', UTCstart, 'to', UTCstop
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")

    # ATTENTION: START COPY FROM SUMMARY_VAR
    # Load the data 
    # These lists have 4 dimensions [time, lev, n, N_box]
    varM_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    meanM_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    varN_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    meanN_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    varm_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    meanm_list = create3Dlist(len(timelist), len(args.height), len(nlist))
    
    # loop over dates 
    for d in args.date:
        # Load dataset 
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        
        # Load the required data and put it in list
        for it, time in enumerate(timelist):
            for iz, lev in enumerate(dataset.variables['levs']):
                for i_n, n in enumerate(dataset.variables['n']):
                    nmax = 265/n
                    
                    meanM = dataset.variables['meanM'][it,iz,i_n,:nmax,:nmax]
                    varM = dataset.variables['varM'][it,iz,i_n,:nmax,:nmax]
                    meanN = dataset.variables['meanN'][it,iz,i_n,:nmax,:nmax]
                    varN = dataset.variables['varN'][it,iz,i_n,:nmax,:nmax]
                    meanm = dataset.variables['meanm'][it,iz,i_n,:nmax,:nmax]
                    varm = dataset.variables['varm'][it,iz,i_n,:nmax,:nmax]
                    
                    meanM_list[it][iz][i_n] += list(np.ravel(meanM))
                    varM_list[it][iz][i_n] += list(np.ravel(varM))
                    meanN_list[it][iz][i_n] += list(np.ravel(meanN))
                    varN_list[it][iz][i_n] += list(np.ravel(varN))
                    meanm_list[it][iz][i_n] += list(np.ravel(meanm))
                    varm_list[it][iz][i_n] += list(np.ravel(varm))
    

    mu2N_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    mu2Nalpha_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    mu2Nbeta_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    mu2Nab_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    alpha_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    beta_list3d = create3Dlist(len(timelist), len(args.height), len(nlist))
    
    # Get means and convert to plotting arrays
    for it, time in enumerate(timelist):
        for iz, lev in enumerate(dataset.variables['levs']):
            for i_n, n in enumerate(dataset.variables['n']):
                # Extract the arrays
                varM = np.array(varM_list[it][iz][i_n])
                meanM = np.array(meanM_list[it][iz][i_n])
                varN = np.array(varN_list[it][iz][i_n])
                meanN = np.array(meanN_list[it][iz][i_n])
                varm = np.array(varm_list[it][iz][i_n])
                meanm = np.array(meanm_list[it][iz][i_n])
                
                mu2 = varM/(meanM**2)
                alpha = varN/meanN
                beta = varm/(meanm**2)
                
                mu2N_list3d[it][iz][i_n] = np.nanmean(mu2*meanN/2)
                mu2Nalpha_list3d[it][iz][i_n] = np.nanmean(mu2*meanN/(1+alpha))
                mu2Nbeta_list3d[it][iz][i_n] = np.nanmean(mu2*meanN/(1+beta))
                mu2Nab_list3d[it][iz][i_n] = np.nanmean(mu2*meanN/(beta+alpha))
                alpha_list3d[it][iz][i_n] = np.nanmean(alpha)
                beta_list3d[it][iz][i_n] = np.nanmean(beta)
    
    mu2N_list3d =  np.array(mu2N_list3d)
    mu2Nalpha_list3d =  np.array(mu2Nalpha_list3d)
    mu2Nbeta_list3d =  np.array(mu2Nbeta_list3d)
    mu2Nab_list3d =  np.array(mu2Nab_list3d)
    alpha_list3d =  np.array(alpha_list3d)
    beta_list3d =  np.array(beta_list3d)
    # ATTENTION: STOP COPY FROM SUMMARY_VAR
    
    # Get the required time avarage
    mu2N_list3d =  np.mean(mu2N_list3d[t1:t2], axis = 0)
    mu2Nalpha_list3d =  np.mean(mu2Nalpha_list3d[t1:t2], axis = 0)
    mu2Nbeta_list3d =  np.mean(mu2Nbeta_list3d[t1:t2], axis = 0)
    mu2Nab_list3d =  np.mean(mu2Nab_list3d[t1:t2], axis = 0)
    alpha_list3d =  np.mean(alpha_list3d[t1:t2], axis = 0)
    beta_list3d =  np.mean(beta_list3d[t1:t2], axis = 0)

    levs = dataset.variables['levs'][:]

    # Create the figure
    fig, axarr = plt.subplots(3, 2, figsize = (95./25.4*2, 10))
    
    for i_n, n in enumerate(nlist):
        axarr[0,0].plot(mu2N_list3d[:,i_n], args.height, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        axarr[0,1].plot(mu2Nab_list3d[:,i_n], args.height, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        axarr[1,0].plot(mu2Nalpha_list3d[:,i_n], args.height, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        axarr[1,1].plot(mu2Nbeta_list3d[:,i_n], args.height, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        axarr[2,0].plot(alpha_list3d[:,i_n], args.height, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        axarr[2,1].plot(beta_list3d[:,i_n], args.height, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        
    axarr[0,0].plot([1.]*len(args.height), args.height, c = 'gray', 
                    zorder = 0.1)
    axarr[0,0].set_xlim(0, 2)
    #axarr[0,0].set_yscale('log')
    axarr[0,0].set_ylabel('height [m]')
    axarr[0,0].set_xlabel(r'$\mu_2 \langle N \rangle/2$')
    axarr[0,0].set_ylim(args.height[0], args.height[-1])
    
    axarr[0,1].plot([1.]*len(args.height), args.height, c = 'gray', 
                    zorder = 0.1)
    axarr[0,1].set_xlim(0, 2)
    #axarr[0,1].set_yscale('log')
    axarr[0,1].set_ylabel('height [m]')
    axarr[0,1].set_xlabel(r'$\mu_2 \langle N \rangle/(\alpha +\beta)$')
    axarr[0,1].set_ylim(args.height[0], args.height[-1])
    
    axarr[1,0].plot([1.]*len(args.height), args.height, c = 'gray', 
                    zorder = 0.1)
    axarr[1,0].set_xlim(0, 2)
    #axarr[1,0].set_yscale('log')
    axarr[1,0].set_ylabel('height [m]')
    axarr[1,0].set_xlabel(r'$\mu_2 \langle N \rangle/(1+\alpha)$')
    axarr[1,0].set_ylim(args.height[0], args.height[-1])
    
    axarr[1,1].plot([1.]*len(args.height), args.height, c = 'gray', 
                    zorder = 0.1)
    axarr[1,1].set_xlim(0, 2)
    #axarr[1,1].set_yscale('log')
    axarr[1,1].set_ylabel('height [m]')
    axarr[1,1].set_xlabel(r'$\mu_2 \langle N \rangle/(1+\beta)$')
    axarr[1,1].set_ylim(args.height[0], args.height[-1])
    
    axarr[2,0].plot([1.]*len(args.height), args.height, c = 'gray', 
                    zorder = 0.1)
    axarr[2,0].set_xlim(0, 2)
    #axarr[2,0].set_yscale('log')
    axarr[1,0].set_ylabel('height [m]')
    axarr[2,0].set_xlabel(r'$\alpha$')
    axarr[2,0].set_ylim(args.height[0], args.height[-1])
    
    axarr[2,1].plot([1.]*len(args.height), args.height, c = 'gray', 
                    zorder = 0.1)
    axarr[2,1].set_xlim(0, 2)
    #axarr[2,1].set_yscale('log')
    axarr[2,1].set_ylabel('height [m]')
    axarr[2,1].set_xlabel(r'$\beta$')
    axarr[2,1].set_ylim(args.height[0], args.height[-1])
        
    axarr[2,1].legend(loc =4, ncol = 1, prop={'size':6})
            
            
    titlestr = (alldatestr + '\n' + args.ana + 
                ', water=' + str(args.water) + 
                ', nens=' + str(args.nens) + ', from ' + str(UTCstart) + 
                ' to ' + str(UTCstop))
    fig.suptitle(titlestr, fontsize='x-large')
    plt.tight_layout(rect=[0, 0.0, 1, 0.93])
    
    plotsavestr = ('height_var_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) +
                    '_nens-' + str(args.nens)+ '_tstart-' + 
                    str(args.tstart) + '_tend-' + str(args.tend) + 
                    '_tinc-' + str(args.tinc) + '_tplot-' + str(t1) + 
                    '-' + str(t2))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')


################################################################################
if 'M_vert' in args.plot:
    print 'Plotting M_vert'
    plotdirsub = plotdir +  '/M_vert/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    # Setup
    t1 = int(args.tplot[0]); t2 = int(args.tplot[1])
    timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
    UTCstart = timelist[t1]
    UTCstop = timelist[t2-1]
    print 'Summing from', UTCstart, 'to', UTCstop
    
    # loop over dates 
    Mtotlist = []
    Msouthlist = []
    Mnorthlist = []
    for d in args.date:
        # Load dataset 
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        
        Mtotlist.append(dataset.variables['Mtot'][:])
        Msouthlist.append(dataset.variables['Msouth'][:])
        Mnorthlist.append(dataset.variables['Mnorth'][:])
    height = dataset.variables['height'][:]/1e3

    
    Mtotlist =  np.mean(Mtotlist, axis = 0)
    Msouthlist =  np.mean(Msouthlist, axis = 0)
    Mnorthlist =  np.mean(Mnorthlist, axis = 0)
    
    # Get the required time avarage
    Mtotlist =  np.mean(Mtotlist[t1:t2], axis = 0)
    Msouthlist =  np.mean(Msouthlist[t1:t2], axis = 0)
    Mnorthlist =  np.mean(Mnorthlist[t1:t2], axis = 0)

    levs = dataset.variables['levs'][:]
    
    print Mtotlist
    print levs
    # Create the figure
    fig, ax = plt.subplots(1, 2, figsize = (pdfwidth, 4))
    

    ax[0].plot(Mtotlist, levs, c = 'r', 
                    label = 'Total domain: mean height asl 247m', linewidth = 1.5)
    ax[0].plot(Msouthlist, levs, c = 'g', 
                    label = 'Southern half: mean height asl 415m', linewidth = 1.5)
    ax[0].plot(Mnorthlist, levs, c = 'b', 
                    label = 'Northern half: mean height asl 79m', linewidth = 1.5)

    #ax.set_xlim(0, 2)
    #axarr[0,0].set_yscale('log')
    ax[0].set_ylabel('Vertical model level')
    ax[0].set_xlabel(r'$M$[kg/s]')
    ax[0].invert_yaxis()
    ax[0].plot([0,1.2e10],[30,30], c = 'gray', zorder = 0.5)
    ax[0].plot([0,1.2e10],[27,27], c = 'lightgray', zorder = 0.5)
    ax[0].plot([0,1.2e10],[34,34], c = 'lightgray', zorder = 0.5)
    #ax.set_ylim(args.height[0], args.height[-1])
    

    ax[0].legend(loc =4, ncol = 1, prop={'size':8})
    
    ax[1].plot(height, levs)
    ax[1].invert_yaxis()
    ax[1].set_xlabel(r'Height asl [km]')
    ax[1].plot([0,14],[30,30], c = 'gray', zorder = 0.5)
    ax[1].plot([0,14],[27,27], c = 'lightgray', zorder = 0.5)
    ax[1].plot([0,14],[34,34], c = 'lightgray', zorder = 0.5)
    
    ax[0].set_title('Vertical mass flux profile', fontsize=10)
    ax[1].set_title('Height for a column above sea level', fontsize=10)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    plt.suptitle(UTCstart, fontsize = 12)
    
    plotsavestr = ('M_vert_' + alldatestr + anastr + '_tplot-' + str(t1) + 
                    '-' + str(t2))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')


################################################################################
if 'prec_stamps' in args.plot:
    print 'Plotting prec_stamps'
    plotdirsub = plotdir +  '/prec_stamps/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    cmPrec = ((1    , 1     , 1    ), 
            (0    , 0.627 , 1    ),
            (0.137, 0.235 , 0.98 ),
            (0.392, 0     , 0.627),
            (0.784, 0     , 0.627))
            #(0.1  , 0.1   , 0.784),
    levelsPrec = [0, 1, 3, 10, 30, 100.]
    
    
    ensdir = '/home/scratch/users/stephan.rasp/' + args.date[0] + '/deout_ceu_pspens/'
    t = timedelta(hours = args.tplot[0])
    print t
    #HHobj = getfobj_ncdf(ensdir + '/1/OUTPUT/lfff00000000c.nc_30m', 'HHL')
    TTENSobj = getfobj_ncdf(ensdir + '/1/OUTPUT/lfff' + ddhhmmss(t) + 
                            '.nc_30m', 'QVTENSSTO')
    PREC1obj = getfobj_ncdf(ensdir + '/1/OUTPUT/lfff' + ddhhmmss(t) + 
                            '.nc_30m_surf', 'PREC_ACCUM')
    PREC2obj = getfobj_ncdf(ensdir + '/2/OUTPUT/lfff' + ddhhmmss(t) + 
                            '.nc_30m_surf', 'PREC_ACCUM')
    PREC3obj = getfobj_ncdf(ensdir + '/3/OUTPUT/lfff' + ddhhmmss(t) + 
                            '.nc_30m_surf', 'PREC_ACCUM')
    
    radarpref = '/project/meteo/w2w/A6/radolan/netcdf_cosmo_de/raa01-rw_10000-'
    radarsufx = '-dwd---bin.nc'
    dateobj = yyyymmddhh_strtotime(args.date[0])
    dtradar = timedelta(minutes = 10)
    t_rad = dateobj - dtradar + t
    print t_rad
    RADARobj = getfobj_ncdf(radarpref + yymmddhhmm(t_rad) + 
                            radarsufx, 'pr', dwdradar = True)
    
    # Set up the figure 
    fig, axarr = plt.subplots(2, 4, figsize = (pdfwidth, 5))
    plt.sca(axarr[1,0])   # This is necessary for some reason...
    print TTENSobj.data.max()
    cf, tmp = ax_contourf(axarr[1,0], TTENSobj,
                          cmap = 'bwr',pllevels=np.linspace(-3e-7, 3e-7, 100),
                            ji0=(110, 160),
                            ji1=(357-161, 357-111),
                            sp_title='QV perturbation',
                            Basemap_drawrivers = False,
                            npars = 0, nmers = 0,
                            lev =45,
                            extend = 'both')
    #cb = fig.colorbar(cf)
    #cb.set_label(HHobj.unit)
    
    for ax, fob, i, in zip(list(axarr[0,1:]),
                           [PREC1obj,PREC2obj, PREC3obj], range(1,4)):
        
        # 2. NvarMN
        plt.sca(ax)
        cf, tmp = ax_contourf(ax, fob, 
                              colors = cmPrec, pllevels = levelsPrec,
                              ji0=(50, 50),
                            ji1=(357-51, 357-51),
                            sp_title='Member ' + str(i),
                            Basemap_drawrivers = False,
                            npars = 0, nmers = 0)
        #cb = fig.colorbar(cf)
        #cb.set_label(fob.unit)
    plt.sca(axarr[0,0])
    cf, tmp = ax_contourf(axarr[0,0], RADARobj, 
                              colors = cmPrec, pllevels = levelsPrec,
                              ji0=(50+62, 50+22),
                            ji1=(357-51+62, 357-51+22),
                            sp_title='RADAR',
                            Basemap_drawrivers = False,
                            npars = 0, nmers = 0)
    #cb = fig.colorbar(cf)
    #cb.set_label(fob.unit)
    
    #####################################33
    ## Set up the figure 2
    #fig, axarr = plt.subplots(2, 2, figsize = (9, 8))
    #plt.sca(axarr[0,0])   # This is necessary for some reason...
     
    #tmpfield = np.ones((51, 357, 357), dtype = bool)
    #tmpfield[:,110:357-161,160:357-111] = 0
    #HHobj.data[tmpfield] = np.nan
    
    #cf, tmp = ax_contourf(axarr[0,0], HHobj,
                          #cmap = 'gist_earth', pllevels = np.linspace(0, 1000, 100),
                            #ji0=(50, 50),
                            #ji1=(357-51, 357-51),
                            #sp_title=HHobj.fieldn,
                            #Basemap_drawrivers = False,
                            #npars = 0, nmers = 0,
                            #lev = 50,
                            #extend = 'both')
    #cb = fig.colorbar(cf)
    #cb.set_label(HHobj.unit)
    
    for ax, fob, i, in zip(list(axarr[1,1:]), 
                           [PREC1obj,PREC2obj, PREC3obj], range(1,4)):
        
        # 2. NvarMN
        plt.sca(ax)
        cf, tmp = ax_contourf(ax, fob, 
                              colors = cmPrec, pllevels = levelsPrec,
                              ji0=(110, 160),
                            ji1=(357-161, 357-111),
                            sp_title='Member ' + str(i),
                            Basemap_drawrivers = False,
                            npars = 0, nmers = 0)
        #cb = fig.colorbar(cf)
        #cb.set_label(fob.unit)
    
    titlestr = (args.date[0] + '+' + ddhhmmss(t) + ', ' + args.ana + 
                ', water=' + str(args.water) + 
                ', nens=' + str(args.nens))
    fig.suptitle('Precipitation stamps ' + str(dateobj+t), fontsize=12)
    plt.tight_layout(rect=[0, 0.0, 1, 0.93])
    
    plotsavestr = ('stamps_prec_' + args.date[0] + '_ana-' + args.ana + 
                '_wat-' + str(args.water) +
                '_nens-' + str(args.nens) + '_time-' + ddhhmmss(t))
    print plotsavestr
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')
    

if 'identification' in args.plot:
    print 'Plotting identification'
    plotdirsub = plotdir +  '/identification/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    dataset = Dataset(savedir + args.date[0] + savesuf, 'r')
    
    # ATTENTION For now this just takes one level 
    
    y1 = 10; y2 = 45; x1 = 180; x2 = 215
    exw = dataset.variables['exw'][args.tplot[0], x1:x2, y1:y2]
    exq = dataset.variables['exq'][args.tplot[0], x1:x2, y1:y2]*1000.
    excld = dataset.variables['excld'][args.tplot[0], x1:x2, y1:y2]
    exwater = dataset.variables['exwater'][args.tplot[0], x1:x2, y1:y2]
    excld[excld == 0] = np.nan
    exwater[exwater == 0] = np.nan
    exbin = (exw > 1.) * (exq > 0.)
    
    fig, axarr = plt.subplots(2, 3, figsize = (pdfwidth, 7))
    
    
    # First plot
    from matplotlib import colors
    mpl.rcParams['axes.linewidth'] = 0
    cmw = ("#0043C1","#405EBF","#6476C5","#808DCC","#99A2D4","#B0B6DC",
          "#C4C8E3","#D6D8E8","#E4E5ED","#FFFFFF","#FFFFFF","#EEE3E5",
          "#EAD3D7","#E5C0C7","#DFABB4","#D792A0","#CC7789","#C05A72",
          "#B23659","#A10040")
    cmw = colors.ListedColormap(cmw)
    boundsw = np.arange(-5,5.5,0.5)
    normw = colors.BoundaryNorm(boundsw, cmw.N)
    print exw
    Cw = axarr[0,0].imshow(exw, interpolation = 'nearest', origin = 'lower',
                  cmap = 'bwr', alpha = 1, vmin = -5, vmax = 5)
    cb = fig.colorbar(Cw, ax = axarr[0,0], orientation = 'horizontal', 
                      fraction = 0.05, pad = 0.05)
    cb.set_label('[m/s]')
    cb.set_ticks([-5, -2.5, 0, 2.5, 5])
    axarr[0,0].set_title('Vertical velocity')
    
    tmp = np.copy(exw)
    tmp[exbin == False] = np.nan
    print exbin
    print tmp
    Cw = axarr[1,0].imshow(tmp, interpolation = 'nearest', origin = 'lower',
                  cmap = 'bwr', alpha = 1, vmin = -5, vmax = 5)
    cb = fig.colorbar(Cw, ax = axarr[1,0], orientation = 'horizontal', 
                      fraction = 0.05, pad = 0.05)
    
    tmp2 = np.copy(exq)
    tmp2[exbin == False] = np.nan
    Cw = axarr[0,1].imshow(tmp2, interpolation = 'nearest', origin = 'lower',
                  cmap = 'cool', alpha = 1, vmin = 0, vmax = 3)
    cb = fig.colorbar(Cw, ax = axarr[0,1], orientation = 'horizontal', 
                      fraction = 0.05, pad = 0.05)
    cb.set_label('[g/kg]')
    cb.set_ticks([0, 1, 2, 3])
    cb.set_ticklabels(['>0', 1, 2, 3])
    axarr[0,1].set_title('Cloud water')
    
    Cw = axarr[1,1].imshow(exbin, interpolation = 'nearest', origin = 'lower',
                  cmap = 'Greys', alpha = 1, vmin = 0, vmax = 1)
    cb = fig.colorbar(Cw, ax = axarr[1,1], orientation = 'horizontal', 
                      fraction = 0.05, pad = 0.05)
    cb.set_label('')
    axarr[1,1].set_title('Binary field')
    
    Cw = axarr[0,2].imshow(excld, interpolation = 'nearest', origin = 'lower',
                  cmap = 'prism', alpha = 1, vmin = np.nanmin(excld),
                  vmax = np.nanmax(excld))
    cb = fig.colorbar(Cw, ax = axarr[0,2], orientation = 'horizontal', 
                      fraction = 0.05, pad = 0.05)
    cb.set_label('')
    axarr[0,2].set_title('Clouds identified')
    
    Cw = axarr[1,2].imshow(exwater, interpolation = 'nearest', origin = 'lower',
                  cmap = 'prism', alpha = 1, vmin = np.nanmin(exwater),
                  vmax = np.nanmax(exwater)*1.1)
    cb = fig.colorbar(Cw, ax = axarr[1,2], orientation = 'horizontal', 
                      fraction = 0.05, pad = 0.05)
    cb.set_label('')
    axarr[1,2].set_title('Clouds separated')

    for ax in list(np.ravel(axarr)):
        mpl.rcParams['axes.linewidth'] = 0
        ax.tick_params(bottom='off',
                       top='off',
                       left='off',
                       right='off',
                       labelbottom='off',
                       labeltop='off',
                       labelleft='off',
                       labelright='off')
    
    t = timedelta(hours = args.tplot[0] + args.tstart)
    titlestr = (args.date[0] + '+' + ddhhmmss(t) + ', ' + args.ana + 
                ', water=' + str(args.water) + 
                ', nens=' + str(args.nens))
    fig.suptitle('Example for cloud identification and separation', fontsize=12)
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    
    plotsavestr = ('identification_' + args.date[0] + anastr + '_' + ddhhmmss(t))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')
    
    
###################################################33
# HYPO

################################################################################
if 'hypo_std_v_mean' in args.plot:
    print 'Plotting std_v_mean'

    plotdirsub = plotdir +  '/std_v_mean/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    dx = 2.8e3
    # Setup
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    
    
    # Load the data: I want to plot mu_2, N, alpha
    # These list have 3 dims [lev, n, N_box]

    stdM_list = create1Dlist(len(nlist))
    M_list = create1Dlist(len(nlist))

    
    # Loop over dates 
    for d in args.date:
        # Load dataset 
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        
        # Load the required data and put it in list
        for i_n, n in enumerate(dataset.variables['n']):
            nmax = 265/n
            meanM = dataset.variables['meanM'][:,i_n,:nmax,:nmax]
            varM = dataset.variables['varM'][:,i_n,:nmax,:nmax]

            
            stdM_list[i_n] += list(np.ravel(np.sqrt(varM)))
            M_list[i_n] += list(np.ravel(meanM))

                
    # now I have the lists I want in the scatter plot 

    
    # Set up the figure 
    tmp = np.logspace(5,12, 1000)
    fig, axarr = plt.subplots(1, 2, figsize = (pdfwidth, 4))
    
    ####### n loop #######################
    for i_n, n in enumerate(dataset.variables['n']):
        print 'n: ', n
        z_n = 0.1 + (n/256.)*0.1   # for z_order
        
        # Get the data
        stdM = np.array(stdM_list[i_n])
        M = np.array(M_list[i_n])

        
        # M v var(M) 
        axarr[0].scatter(M, stdM, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        
        # Fit the line, SPPT
        tmp2 = np.logspace(6, 11, 1000)
        p0 = [1]
        y = np.array(stdM)
        x = np.array(M)
        mask = np.isfinite(y)
        result = leastsq(residual_b, p0, args = (y[mask], x[mask]))
        b = result[0]
        axarr[0].plot(tmp2,b*tmp2, c = clist[i_n], alpha = 1, 
                        linestyle = '--', zorder = 2, label = 'SPPT')
        # How good is the fit
        nrmse = (np.sqrt(np.mean(residual_b(b, y[mask], x[mask])**2))/
                 np.mean(y[mask]))
        print b, nrmse
        axarr[1].errorbar(np.log2(n)-0.1, b, yerr = nrmse/2., c = clist[i_n], 
                            fmt = 'o', label = 'SPPT')
        
        # Fit the line, CC06
        result = leastsq(residual_b_sqrt, p0, args = (y[mask], x[mask]))
        b = result[0]
        axarr[0].plot(tmp2,np.sqrt(b*tmp2), c = clist[i_n], alpha = 1, 
                        linestyle = '-.', zorder = 2, label = 'CC06')
        # How good is the fit
        nrmse = (np.sqrt(np.mean(residual_b_sqrt(b, y[mask], x[mask])**2))/
                 np.mean(y[mask]))
        print b, nrmse
        axarr[1].errorbar(np.log2(n)+0.1, b/1e8, yerr = nrmse/2., c = clist[i_n], 
                            fmt = 'x', label = 'CC06 / 1e8')
        
       
        

    axarr[0].set_xlim(1e6,1e11)
    axarr[0].set_ylim(1e6,1e10)
    axarr[0].set_xscale('log')
    axarr[0].set_yscale('log')
    #axarr[0,0].xaxis.grid(True)
    #ax2.invert_xaxis()
    axarr[0].set_xlabel(r'$\langle M \rangle$ [kg/s]')
    axarr[0].set_ylabel(r'$\langle (\delta M)^2 \rangle^{1/2}$ [kg/s]')
    axarr[0].set_title('Mass flux: variability against mean', fontsize = 10)
    axarr[0].legend(loc = 4, ncol = 2, prop={'size':6})
    
    
    axarr[1].set_xlabel(r'log$_2$ n')
    axarr[1].set_ylabel(r'fit parameter and NRMSE')
    axarr[1].set_ylim(0,2)
    axarr[1].set_title('How good are the CC06 and SPPT fits', fontsize = 10)
    axarr[1].legend(loc = 1, ncol = 2, prop={'size':6})
   
    
    axarr[0].text(0.05, 0.9, '(a)', transform = axarr[0].transAxes, 
                    fontsize = 10)
    axarr[1].text(0.05, 0.9, '(b)', transform = axarr[1].transAxes, 
                    fontsize = 10)

    
    titlestr = (alldatestr + '\n' + args.ana + 
                ', water=' + str(args.water) + ', lev= ' + str(int(args.height[0])) + 
                ', nens=' + str(args.nens))
    titlestr += '\nCC06 variance scaling fits data reasonably well'
    fig.suptitle(titlestr)
    plt.tight_layout(rect=[0, 0.0, 1, 0.90])
    
    plotsavestr = ('hypo_std_v_mean_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                    '_nens-' + str(args.nens))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')
    


################################################################################
if 'hypo_scatter' in args.plot:
    print 'Plotting scatter'
    plotdirsub = plotdir +  '/scatter/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)

    # Setup
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    
    
    # Load the data 
    # These lists have 4 dimensions [time, lev, n, N_box]
    varM_list = create2Dlist(len(nlist), len(timelist))
    meanM_list = create2Dlist(len(nlist), len(timelist))
    varN_list = create2Dlist(len(nlist), len(timelist))
    meanN_list = create2Dlist(len(nlist), len(timelist))
    varm_list = create2Dlist(len(nlist), len(timelist))
    meanm_list = create2Dlist(len(nlist), len(timelist))
    
    # loop over dates 
    for d in args.date:
        # Load dataset 
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        
        # Load the required data and put it in list
        for it, time in enumerate(timelist):
            for i_n, n in enumerate(dataset.variables['n']):
                nmax = 265/n
                
                meanM = dataset.variables['meanM'][it,i_n,:nmax,:nmax]
                varM = dataset.variables['varM'][it,i_n,:nmax,:nmax]
                meanN = dataset.variables['meanN'][it,i_n,:nmax,:nmax]
                varN = dataset.variables['varN'][it,i_n,:nmax,:nmax]
                meanm = dataset.variables['meanm'][it,i_n,:nmax,:nmax]
                varm = dataset.variables['varm'][it,i_n,:nmax,:nmax]
                
                meanM_list[i_n][it] += list(np.ravel(meanM))
                varM_list[i_n][it] += list(np.ravel(varM))
                meanN_list[i_n][it] += list(np.ravel(meanN))
                varN_list[i_n][it] += list(np.ravel(varN))
                meanm_list[i_n][it] += list(np.ravel(meanm))
                varm_list[i_n][it] += list(np.ravel(varm))

                
    # now I have the lists I want in the scatter plot 
        
    # Set up the figure 
    fig, axarr = plt.subplots(5, 3, figsize = (pdfwidth, 12))

    rmselist1 = []
    rmselist2 = []
    rmselist3 = []
    rmselist4 = []
    rmselist5 = []
    ####### n loop #######################
    for i_n, n in enumerate(dataset.variables['n']):
        print 'n: ', n
        z_n = 0.1 + (n/256.)*0.1   # for z_order
        
        # Extract the arrays
        M = np.array(meanM_list[i_n])
        m = np.array(meanm_list[i_n])
        N = np.array(meanN_list[i_n])
        varM = np.array(varM_list[i_n])
        varm = np.array(varm_list[i_n])
        varN = np.array(varN_list[i_n])
        alpha = varN/N
        beta = varm/(m**2)
        # These now have dimensions [time, date]        
        
        # 1. Raw PC08 prediction
        m_c = 4e7
        predict = 2*m_c*M
        
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist1.append(nrmse)
        
        axarr[0,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[0,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[0,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[0,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        
        axarr[0,2].scatter([0], [0], marker = 'o', c = clist[i_n], 
                            s = 8, label = str(n*2.8)+'km',
                            linewidth = 0)  # Ghost plot for labels
        
        # 2. CC06 prediction
        print 'mean m', np.nanmean(m)
        predict = 2*m*M
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist2.append(nrmse)
        
        axarr[1,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[1,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[1,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[1,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        var_time = np.nanmean(m, axis = 1)
        axarr[1,2].plot(timelist_plot, var_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        
        # 3. alpha adjusted
        predict = (1+alpha)*m*M
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist3.append(nrmse)
        
        axarr[2,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[2,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[2,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[2,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        var_time = np.nanmean(alpha, axis = 1)
        axarr[2,2].plot(timelist_plot, var_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        
        # 4. beta adjusted
        predict = (1+beta)*m*M
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist4.append(nrmse)
        
        axarr[3,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[3,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[3,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[3,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        var_time = np.nanmean(beta, axis = 1)
        axarr[3,2].plot(timelist_plot, var_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        
        # 5. all
        predict = (alpha+beta)*m*M
        frac = varM/predict
        frac_mean = np.nanmean(frac)
        frac_std = np.nanstd(frac)
        nrmse = np.sqrt(np.nanmean(((varM-predict)/varM)**2))
        rmselist5.append(nrmse)
        
        axarr[4,0].scatter(M, frac, marker = 'o', c = clist[i_n], 
                            s = 4, zorder = z_n, linewidth = 0, 
                            alpha = 0.8)
        axarr[4,0].scatter(np.nanmean(M), frac_mean, marker = 'o', 
                           c = clist[i_n], 
                           s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                        label = r'm:{:.2f}, s:{:.2f}'.format(frac_mean, 
                                                    frac_std))
        axarr[4,0].errorbar(np.nanmean(M), frac_mean, marker = 'o', 
                            mec = clist[i_n], 
                        ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                        yerr = frac_std, c = 'black')
        
        frac_time = np.nanmean(frac, axis = 1)
        axarr[4,1].plot(timelist_plot, frac_time, c = clist[i_n], 
                        label = str(n*2.8)+'km', linewidth = 1.5)
        #var_time = np.nanmean(m, axis = 1)
        #axarr[1,2].plot(timelist_plot, var_time, c = clist[i_n], 
                        #label = str(n*2.8)+'km', linewidth = 1.5)
        
        

    
    # Complete the figure
    #axarr[0,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[0,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[0,0].set_xlim(1e6,1e11)
    axarr[0,0].set_ylim(0, 2.5)
    axarr[0,0].set_xscale('log')
    #axarr[0,0].set_yscale('log')
    axarr[0,0].set_xlabel(r'$M$')
    axarr[0,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / 2m_c\langle M \rangle$')
    axarr[0,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[0,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[0,1].set_ylim(0, 2.5)
    axarr[0,1].set_xlabel('time [h/UTC]')
    axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    axarr[0,2].legend(loc =3, ncol = 2, prop={'size':8})
    #axarr[0,2].set_xlabel('time [h/UTC]')
    axarr[0,2].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    #axarr[1,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[1,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[1,0].set_xlim(1e6,1e11)
    axarr[1,0].set_ylim(0, 2.5)
    axarr[1,0].set_xscale('log')
    axarr[1,0].set_xlabel(r'$M$')
    axarr[1,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / 2\langle m \rangle \langle M \rangle$')
    axarr[1,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[1,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[1,1].set_ylim(0, 2.5)
    axarr[1,1].set_xlabel('time [h/UTC]')
    axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    axarr[1,2].set_ylabel(r'$\langle m \rangle$ [kg/s]')
    axarr[1,2].set_xlabel('time [h/UTC]')
    axarr[1,2].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    #axarr[2,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[2,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[2,0].set_xlim(1e6,1e11)
    axarr[2,0].set_ylim(0, 2.5)
    axarr[2,0].set_xscale('log')
    axarr[2,0].set_xlabel(r'$M$')
    axarr[2,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / (1+\alpha)\langle m \rangle \langle M \rangle$')
    axarr[2,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[2,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[2,1].set_ylim(0, 2.5)
    axarr[2,1].set_xlabel('time [h/UTC]')
    axarr[2,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    
    axarr[2,2].set_ylim(0, 2.5)
    axarr[2,2].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[2,2].set_ylabel(r'$\alpha$')
    axarr[2,2].set_xlabel('time [h/UTC]')
    axarr[2,2].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    #axarr[3,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[3,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[3,0].set_xlim(1e6,1e11)
    axarr[3,0].set_ylim(0, 2.5)
    axarr[3,0].set_xscale('log')
    axarr[3,0].set_xlabel(r'$M$')
    axarr[3,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / (1+\beta)\langle m \rangle \langle M \rangle$')
    axarr[3,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[3,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[3,1].set_ylim(0, 2.5)
    axarr[3,1].set_xlabel('time [h/UTC]')
    axarr[3,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    axarr[3,2].set_ylim(0, 2.5)
    axarr[3,2].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[3,2].set_ylabel(r'$\beta$')
    axarr[3,2].set_xlabel('time [h/UTC]')
    axarr[3,2].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    #axarr[4,0].legend(loc =3, ncol = 2, prop={'size':6})
    tmp = np.array([0,10])
    axarr[4,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    axarr[4,0].set_xlim(1e6,1e11)
    axarr[4,0].set_ylim(0, 2.5)
    axarr[4,0].set_xscale('log')
    axarr[4,0].set_xlabel(r'$M$')
    axarr[4,0].set_ylabel(r'$\langle (\delta M)^2 \rangle / (\alpha+\beta)\langle m \rangle \langle M \rangle$')
    axarr[4,0].plot([1e6,1e11],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
            zorder = 0.1)
    
    axarr[4,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
    axarr[4,1].set_ylim(0, 2.5)
    axarr[4,1].set_xlabel('time [h/UTC]')
    axarr[4,1].set_xlim(timelist_plot[0], timelist_plot[-1])
    
    for ax, let in zip(list(np.ravel(axarr)),
                       ['a','b','c','d','e','f','g','h','i','j','l','m','n','o','p']):
        ax.text(0.05, 0.9, '('+let+')', transform = ax.transAxes, 
                    fontsize = 10)
    

    #titlestr = (alldatestr + '\n' + args.ana + 
                #', water=' + str(args.water) + ', lev= ' + str(int(args.height[0])) + 
                #', nens=' + str(args.nens))
    #fig.suptitle(titlestr, fontsize='x-large')
    fig.suptitle('Comparison of simulation results and predictions \nof increasing complexity', fontsize=12)
    plt.tight_layout(rect=[0, 0.0, 1, 0.93])
    
    plotsavestr = ('hypo_scatter_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) + '_lev-' + str(int(args.height[0])) +
                    '_nens-' + str(args.nens))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)

    plt.close('all')
