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
from cosmo_utils.helpers import ddhhmmss
from cosmo_utils.pyncdf import getfobj_ncdf, getfobj_ncdf_ens
from cosmo_utils.pywgrib import fieldobj
from cosmo_utils.plot import ax_contourf
from cosmo_utils.diag import mean_spread_fieldobjlist
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

# Define functions
def residual(p, y, x):
    
    a,b,c = p
    #err = y - (a + b*x**c)
    err = y - (0 + b*x**c)
    
    return err


# Setup: this needs to match the compute script!
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--ana', metavar = 'ana', type=str, default = 'm')
parser.add_argument('--date', metavar = 'date', type=str, nargs = '+')
parser.add_argument('--height', metavar = 'height', type=float, nargs = '+',
                    default = [3000])
parser.add_argument('--water', metavar = 'water', type=bool, default = True)
parser.add_argument('--nens', metavar = 'nens', type=int, default = 20)
parser.add_argument('--tstart', metavar = 'tstart', type=int, default = 1)
parser.add_argument('--tend', metavar = 'tend', type=int, default = 24)
parser.add_argument('--tinc', metavar = 'tinc', type=int, default = 60)
parser.add_argument('--plot', metavar = 'plot', type=str, nargs = '+')
parser.add_argument('--tplot', metavar = 'tplot', type=float, nargs = '+',
                    default = [9,12])
args = parser.parse_args()

if 'all' in args.plot:
    args.plot = ['cloud_stats', 'rdf', 'scatter', 'summary_stats', 'summary_var',
                 'stamps_var', 'stamps_w', 'height_var', 'std_v_mean']

################################################################################
# Load datatsets
datasetlist = []
heightstr = ''
alldatestr = ''
nlist = [256, 128, 64, 32, 16, 8, 4]
for h in args.height:
    heightstr += str(int(h)) + '_'
for d in args.date:
    alldatestr += d + '_'
    print 'Loading date: ', d
    # Create file str
    savedir = '/home/scratch/users/stephan.rasp/results/'
    savesuf = ('_ana-' + args.ana + '_wat-' + str(args.water) + 
               '_height-' + heightstr +
               'nens-' + str(args.nens) + '_tstart-' + str(args.tstart) + 
               '_tend-' + str(args.tend) + '_tinc-' + str(args.tinc) + '.nc')
    # Load file 
    #datasetlist.append(Dataset(savedir + savestr, 'r'))
alldatestr = alldatestr[:-1]   # revome final underscore
if alldatestr == '2016052800_2016052900_2016053000_2016053100_2016060100_2016060200_2016060300_2016060400_2016060500_2016060600_2016060700_2016060800':
    alldatestr = 'composite'
    print 'Composite!'
# End load datasets
################################################################################


# Some preliminaries
tmp_dataset = Dataset(savedir + args.date[0] + savesuf, 'r')
timelist = [timedelta(seconds=ts) for ts in tmp_dataset.variables['time']]
timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
plotdir = ('/home/s/S.Rasp/Dropbox/figures/PhD/variance/' + alldatestr + 
           '/' + args.ana)

ensdir = '/home/scratch/users/stephan.rasp/' + args.date[0] + '/deout_ceu_pspens/'


# Define helper function
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
    sizemax = 4e8
    summax = 10e8
    totlist1 = []
    totlist2 = []
    for iz, lev in enumerate(dataset.variables['levs']):
        tmplist1 = []
        tmplist2 = []
        for d in args.date:
            dataset = Dataset(savedir + d + savesuf, 'r')
            cld_size_tmp = dataset.variables['cld_size'][:,iz, :]
            cld_size_tmp = cld_size_tmp[~cld_size_tmp.mask].data
            cld_sum_tmp = dataset.variables['cld_sum'][:,iz, :]
            cld_sum_tmp = cld_sum_tmp[~cld_sum_tmp.mask].data
            tmplist1 += list(cld_size_tmp)
            tmplist2 += list(cld_sum_tmp)
        totlist1.append(tmplist1)
        totlist2.append(tmplist2)
    

    ######## Lev loop #####################
    for iz, lev in enumerate(dataset.variables['levs']):
        # Get the data
        print 'lev: ', lev

        sizehist, sizeedges = np.histogram(totlist1[iz], 
                                            bins = 15, range = [0., sizemax])
        sizemean = np.mean(totlist1[iz])

        sumhist, sumedges = np.histogram(totlist2[iz], 
                                            bins = 15, range = [0., summax])
        summean = np.mean(totlist2[iz])
        
        # Plot the histograms
        fig, axarr = plt.subplots(1, 2, figsize = (95./25.4*2.5, 4.2))
        axarr[0].bar(sizeedges[:-1], sizehist, width = np.diff(sizeedges)[0])
        axarr[0].plot([sizemean, sizemean], [1, 1e6], c = 'red', 
                    alpha = 0.5)
        axarr[0].set_xlabel('Cloud size [m^2]')
        axarr[0].set_ylabel('Number of clouds')
        axarr[0].set_title('Cloud size distribution')
        axarr[0].set_xlim([0., sizemax])
        axarr[0].set_ylim([1, 1e6])
        axarr[0].set_yscale('log')
        
        axarr[1].bar(sumedges[:-1], sumhist, width = np.diff(sumedges)[0])
        axarr[1].plot([summean, summean], [1, 1e6], c = 'red', 
                    alpha = 0.5)
        axarr[1].set_ylabel('Number of clouds')
        axarr[1].set_xlim([0., summax])
        axarr[1].set_ylim([1, 1e6])
        axarr[1].set_yscale('log')
        axarr[1].set_xlabel('Cloud mass flux [kg/s]')
        axarr[1].set_title('Cloud mass flux distribution')
        
        titlestr = (alldatestr  + ', ' + args.ana + 
                    ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                    ', nens=' + str(args.nens))
        fig.suptitle(titlestr, fontsize='x-large')
        plt.tight_layout(rect=[0, 0.0, 1, 0.95])
        
        plotsavestr = ('cloud_stats_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
                        '_nens-' + str(args.nens))
        fig.savefig(plotdirsub + plotsavestr, dpi = 300)
        plt.close('all')



################################################################################
if 'rdf' in args.plot:
    print 'Plotting rdf'

    plotdirsub = plotdir +  '/rdf/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    rdf_list = []
    for d in args.date:
        dataset = Dataset(savedir + d + savesuf, 'r')
        rdf_tmp = dataset.variables['rdf'][:]
        rdf_tmp[rdf_tmp > 1e20] = np.nan
        rdf_list.append(rdf_tmp)
    rdf = np.nanmean(rdf_list, axis = 0)
    
    # Setup
    ymax = 5
    
    # Get 3 hr averages
    rdf_3hr = []
    tlist_3hr = []
    dt = 3 * args.tinc/60.
    for i in range(len(timelist)/3):
        rdf_3hr.append(np.nanmean(rdf[i*dt:(i+1)*dt], axis = 0))
        tlist_3hr.append(i*dt+1)
    rdf_3hr = np.array(rdf_3hr)
    cyc = [plt.cm.jet(i) for i in np.linspace(0, 1, len(tlist_3hr))]
    ######## Lev loop #####################
    for iz, lev in enumerate(dataset.variables['levs']):
        print 'lev: ', lev
        
        # Get the data
        r =   dataset.variables['dr'][:]
        
        fig, ax = plt.subplots(1, 1, figsize = (95./25.4*1.25, 4.5))
        
        ############# Time loop ##############
        for it, t in enumerate(tlist_3hr):
            #print 'time: ', t
            ax.plot(r/1000., rdf_3hr[it, iz, :], c = cyc[it], 
                    label = str(t) + 'UTC pm 1h')
        
        ax.legend(loc = 1, ncol = 2, prop={'size':6})
        ax.plot([0, np.max(r)/1000.], [1, 1], c = 'gray', alpha = 0.5)
        ax.set_xlabel('Distance [km]')
        ax.set_ylabel('Normalized RDF')
        ax.set_title('Radial distribution function')
        ax.set_ylim(0, ymax)
        ax.set_xlim(0, np.max(r)/1000.)
        
        titlestr = (alldatestr + '\nlev= ' + str(lev) + 
                    ', nens=' + str(args.nens))
        fig.suptitle(titlestr, fontsize='x-large')
        plt.tight_layout(rect=[0, 0.0, 1, 0.85])
        
        plotsavestr = ('rdf_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
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
    
    
    # Load the data: I want to plot mu_2, N, alpha
    # These list have 3 dims [lev, n, N_box]
    mu2_list = create2Dlist(len(args.height), len(nlist))
    N_list = create2Dlist(len(args.height), len(nlist))
    alpha_list = create2Dlist(len(args.height), len(nlist))
    beta_list = create2Dlist(len(args.height), len(nlist))

    
    # Loop over dates 
    for d in args.date:
        # Load dataset 
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        
        # Load the required data and put it in list
        for iz, lev in enumerate(dataset.variables['levs']):
            for i_n, n in enumerate(dataset.variables['n']):
                nmax = 265/n
                meanM = dataset.variables['meanM'][:,iz,i_n,:nmax,:nmax]
                varM = dataset.variables['varM'][:,iz,i_n,:nmax,:nmax]
                meanN = dataset.variables['meanN'][:,iz,i_n,:nmax,:nmax]
                varN = dataset.variables['varN'][:,iz,i_n,:nmax,:nmax]
                meanm = dataset.variables['meanm'][:,iz,i_n,:nmax,:nmax]
                varm = dataset.variables['varm'][:,iz,i_n,:nmax,:nmax]
                
                mu2 = np.ravel(varM / (meanM**2))
                alpha = np.ravel(varN / meanN)
                beta = np.ravel(varm / (meanm**2))
                mu2_list[iz][i_n] += list(mu2)
                N_list[iz][i_n] += list(np.ravel(meanN))
                alpha_list[iz][i_n] += list(alpha)
                beta_list[iz][i_n] += list(beta)

                
    # now I have the lists I want in the scatter plot 
    

    ######## Lev loop #####################
    for iz, lev in enumerate(dataset.variables['levs']):
        print 'lev: ', lev
        
        # Set up the figure 
        fig, axarr = plt.subplots(4, 2, figsize = (95./25.4*2, 13))

        
        ####### n loop #######################
        for i_n, n in enumerate(dataset.variables['n']):
            print 'n: ', n
            z_n = 0.1 + (n/256.)*0.1   # for z_order
            
            # Extract the arrays
            N = np.array(N_list[iz][i_n])
            mu2 = np.array(mu2_list[iz][i_n])
            alpha = np.array(alpha_list[iz][i_n])
            beta = np.array(beta_list[iz][i_n])

            
            
            
            # Scatterplots 1
            x = np.sqrt(2/N)
            y = np.sqrt(mu2)
            xmean = np.sqrt(2/np.nanmean(N))
            
            yfrac = mu2 * N / 2
            yfrac_mean = np.nanmean(yfrac)
            yfrac_std = np.nanstd(yfrac)
            ymean = yfrac_mean * xmean
            print 'before', yfrac_mean, yfrac_std

            axarr[0,0].scatter(x, y, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            axarr[0,0].scatter(xmean, ymean, marker = 'o', c = clist[i_n], 
                            s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                            label = str(n*2.8)+'km')
            
            axarr[0,1].scatter(x, yfrac, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            axarr[0,1].scatter(xmean, yfrac_mean, marker = 'o', c = clist[i_n], 
                            s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                            label = r'mean: {:.2f}, std: {:.2f}'.format(yfrac_mean, yfrac_std))
            axarr[0,1].errorbar(xmean, yfrac_mean, marker = 'o', mec = clist[i_n], 
                            ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                            yerr = yfrac_std, c = 'black')
            
            # Scatterplots 2
            # Middle
            axarr[1,0].scatter(alpha, yfrac, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            alpha_mean = np.nanmean(alpha)
            axarr[1,0].scatter(alpha_mean, yfrac_mean, marker = 'o', c = clist[i_n], 
                            s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1)
            
            #x = np.sqrt((1+alpha)/N)
            #xmean = np.sqrt(np.nanmean((1+alpha)/N))
            
            yfrac = mu2 * N / (1 + alpha)
            yfrac_mean = np.nanmean(yfrac)
            yfrac_std = np.nanstd(yfrac)
            ymean = yfrac_mean * xmean
            print 'after', yfrac_mean, yfrac_std
            
            axarr[1,1].scatter(x, yfrac, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            axarr[1,1].scatter(xmean, yfrac_mean, marker = 'o', c = clist[i_n], 
                            s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                            label = r'mean: {:.2f}, std: {:.2f}'.format(yfrac_mean, yfrac_std))
            axarr[1,1].errorbar(xmean, yfrac_mean, marker = 'o', mec = clist[i_n], 
                            ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                            yerr = yfrac_std, c = 'black')
            
            # Scatterplots 3
            # Bottom
            axarr[2,0].scatter(beta, yfrac, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            beta_mean = np.nanmean(beta)
            axarr[2,0].scatter(beta_mean, yfrac_mean, marker = 'o', c = clist[i_n], 
                            s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1)
            
            #x = np.sqrt((1+alpha)/N)
            #xmean = np.sqrt(np.nanmean((1+alpha)/N))
            
            yfrac = mu2 * N / (1 + beta)
            yfrac_mean = np.nanmean(yfrac)
            yfrac_std = np.nanstd(yfrac)
            ymean = yfrac_mean * xmean
            print 'after', yfrac_mean, yfrac_std
            
            axarr[2,1].scatter(x, yfrac, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            axarr[2,1].scatter(xmean, yfrac_mean, marker = 'o', c = clist[i_n], 
                            s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                            label = r'mean: {:.2f}, std: {:.2f}'.format(yfrac_mean, yfrac_std))
            axarr[2,1].errorbar(xmean, yfrac_mean, marker = 'o', mec = clist[i_n], 
                            ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                            yerr = yfrac_std, c = 'black')
            
            
            # Scatterplots 4
            # Bottom
            axarr[3,0].scatter(alpha+beta, yfrac, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            ab_mean = np.nanmean(alpha+beta)
            axarr[3,0].scatter(ab_mean, yfrac_mean, marker = 'o', c = clist[i_n], 
                            s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1)
            
            #x = np.sqrt((1+alpha)/N)
            #xmean = np.sqrt(np.nanmean((1+alpha)/N))
            
            yfrac = mu2 * N / (alpha + beta)
            yfrac_mean = np.nanmean(yfrac)
            yfrac_std = np.nanstd(yfrac)
            ymean = yfrac_mean * xmean
            print 'after', yfrac_mean, yfrac_std
            
            axarr[3,1].scatter(x, yfrac, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            axarr[3,1].scatter(xmean, yfrac_mean, marker = 'o', c = clist[i_n], 
                            s = 40, zorder = 0.5, linewidth = 0.8, alpha = 1,
                            label = r'mean: {:.2f}, std: {:.2f}'.format(yfrac_mean, yfrac_std))
            axarr[3,1].errorbar(xmean, yfrac_mean, marker = 'o', mec = clist[i_n], 
                            ms = 0, zorder = 0.4, linewidth = 1.2, alpha = 1,
                            yerr = yfrac_std, c = 'black')
            
            
        
        # Complete the figure
        axarr[0,0].legend(loc =3, ncol = 2, prop={'size':6})
        tmp = np.array([0,10])
        axarr[0,0].plot(tmp,tmp, c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[0,0].set_xlim(0.05,5)
        axarr[0,0].set_ylim(0.02,10)
        axarr[0,0].set_xscale('log')
        axarr[0,0].set_yscale('log')
        axarr[0,0].invert_xaxis()
        axarr[0,0].set_xlabel(r'$\sqrt{2/N}$')
        axarr[0,0].set_ylabel(r'$\sqrt{\mu_2}$')
        
        axarr[0,1].legend(loc =3, ncol = 2, prop={'size':6})
        axarr[0,1].plot([0.,10],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[0,1].set_xlim(0.05,5)
        axarr[0,1].set_ylim(-1, 2.5)
        axarr[0,1].set_xscale('log')
        #axarr[0,1].set_yscale('log')
        axarr[0,1].invert_xaxis()
        axarr[0,1].set_xlabel(r'$\sqrt{2/N}$')
        axarr[0,1].set_ylabel(r'$\mu_2 \langle N \rangle/2$')
        
        # Complete the figure, middle row
        axarr[1,0].plot([0,10],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[1,0].plot([1,1],[0,10], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[1,0].set_xlim(0.05,5)
        axarr[1,0].set_ylim(0.02,10)
        axarr[1,0].set_xscale('log')
        axarr[1,0].set_yscale('log')
        #axarr[1,0].invert_xaxis()
        axarr[1,0].set_xlabel(r'$\alpha$')
        axarr[1,0].set_ylabel(r'$\mu_2 \langle N \rangle/2$')
        
        axarr[1,1].legend(loc =3, ncol = 2, prop={'size':6})
        axarr[1,1].plot([0.,10],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[1,1].set_xlim(0.05,5)
        axarr[1,1].set_ylim(-1, 2.5)
        axarr[1,1].set_xscale('log')
        #axarr[1,1].set_yscale('log')
        axarr[1,1].invert_xaxis()
        axarr[1,1].set_xlabel(r'$\sqrt{2/N}$')
        axarr[1,1].set_ylabel(r'$\mu_2 \langle N \rangle/(1+\alpha)$')
        
        # Complete the figure, bottom row
        axarr[2,0].plot([0,10],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[2,0].plot([1,1],[0,10], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[2,0].set_xlim(0.05,5)
        axarr[2,0].set_ylim(0.02,10)
        axarr[2,0].set_xscale('log')
        axarr[2,0].set_yscale('log')
        #axarr[1,0].invert_xaxis()
        axarr[2,0].set_xlabel(r'$\beta$')
        axarr[2,0].set_ylabel(r'$\mu_2 \langle N \rangle/2$')
        
        axarr[2,1].legend(loc =3, ncol = 2, prop={'size':6})
        axarr[2,1].plot([0.,10],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[2,1].set_xlim(0.05,5)
        axarr[2,1].set_ylim(-1, 2.5)
        axarr[2,1].set_xscale('log')
        #axarr[2,1].set_yscale('log')
        axarr[2,1].invert_xaxis()
        axarr[2,1].set_xlabel(r'$\sqrt{2/N}$')
        axarr[2,1].set_ylabel(r'$\mu_2 \langle N \rangle/(1+\beta)$')
        
        # Complete the figure, bottom row
        axarr[3,0].plot([0,10],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[3,0].plot([2,2],[0,10], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[3,0].set_xlim(0.1,10)
        axarr[3,0].set_ylim(0.02,10)
        axarr[3,0].set_xscale('log')
        axarr[3,0].set_yscale('log')
        #axarr[1,0].invert_xaxis()
        axarr[3,0].set_xlabel(r'$\alpha + \beta$')
        axarr[3,0].set_ylabel(r'$\mu_2 \langle N \rangle/2$')
        
        axarr[3,1].legend(loc =3, ncol = 2, prop={'size':6})
        axarr[3,1].plot([0.,10],[1,1], c = 'gray', alpha = 0.5, linestyle = '--',
                zorder = 0.1)
        axarr[3,1].set_xlim(0.05,5)
        axarr[3,1].set_ylim(-1, 2.5)
        axarr[3,1].set_xscale('log')
        #axarr[2,1].set_yscale('log')
        axarr[3,1].invert_xaxis()
        axarr[3,1].set_xlabel(r'$\sqrt{2/N}$')
        axarr[3,1].set_ylabel(r'$\mu_2 \langle N \rangle/(\alpha +\beta)$')
        
        titlestr = (alldatestr + '\n' + args.ana + 
                    ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                    ', nens=' + str(args.nens))
        fig.suptitle(titlestr, fontsize='x-large')
        plt.tight_layout(rect=[0, 0.0, 1, 0.93])
        
        plotsavestr = ('scatter_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
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
    
    
    # Load the data: I want to plot mu_2, N, alpha
    # These list have 3 dims [lev, n, N_box]

    stdM_list = create2Dlist(len(args.height), len(nlist))
    M_list = create2Dlist(len(args.height), len(nlist))
    stdQmp_list = create2Dlist(len(args.height), len(nlist))
    Qmp_list = create2Dlist(len(args.height), len(nlist))
    stdQtot_list = create2Dlist(len(args.height), len(nlist))
    Qtot_list = create2Dlist(len(args.height), len(nlist))
    
    # Loop over dates 
    for d in args.date:
        # Load dataset 
        print 'Loading date: ', d
        dataset = Dataset(savedir + d + savesuf, 'r')
        
        # Load the required data and put it in list
        for iz, lev in enumerate(dataset.variables['levs']):
            for i_n, n in enumerate(dataset.variables['n']):
                nmax = 265/n
                meanM = dataset.variables['meanM'][:,iz,i_n,:nmax,:nmax]
                varM = dataset.variables['varM'][:,iz,i_n,:nmax,:nmax]
                meanQmp = dataset.variables['meanQmp'][:,iz,i_n,:nmax,:nmax]
                varQmp = dataset.variables['varQmp'][:,iz,i_n,:nmax,:nmax]
                meanQtot = dataset.variables['meanQtot'][:,iz,i_n,:nmax,:nmax]
                varQtot = dataset.variables['varQtot'][:,iz,i_n,:nmax,:nmax]
                
                stdM_list[iz][i_n] += list(np.ravel(np.sqrt(varM)))
                M_list[iz][i_n] += list(np.ravel(meanM))
                stdQmp_list[iz][i_n] += list(np.ravel(np.sqrt(varQmp)))
                Qmp_list[iz][i_n] += list(np.ravel(meanQmp))
                stdQtot_list[iz][i_n] += list(np.ravel(np.sqrt(varQtot)))
                Qtot_list[iz][i_n] += list(np.ravel(meanQtot))
                
    # now I have the lists I want in the scatter plot 
    

    ######## Lev loop #####################
    for iz, lev in enumerate(dataset.variables['levs']):
        print 'lev: ', lev
        
        # Set up the figure 

        fig, axarr = plt.subplots(1, 3, figsize = (95./25.4*3, 4))
        
        allstdM = []
        allM = []
        ####### n loop #######################
        for i_n, n in enumerate(dataset.variables['n']):
            print 'n: ', n
            z_n = 0.1 + (n/256.)*0.1   # for z_order
            

            stdM = np.array(stdM_list[iz][i_n])
            M = np.array(M_list[iz][i_n])
            allstdM += list(np.ravel(stdM))
            allM += list(np.ravel(M))
            axarr[0].scatter(M, stdM, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            
            stdQmp = np.array(stdQmp_list[iz][i_n]) * 3600. * 24.
            Qmp = np.array(Qmp_list[iz][i_n]) * 3600. * 24.
            axarr[1].scatter(Qmp, stdQmp, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            
            stdQtot = np.array(stdQtot_list[iz][i_n]) * 3600. * 24.
            Qtot = np.array(Qtot_list[iz][i_n]) * 3600. * 24.
            axarr[2].scatter(Qtot, stdQtot, marker = 'o', c = clist[i_n], 
                                s = 4, zorder = z_n, linewidth = 0, 
                                alpha = 0.8)
            
        
        # Fit the line
        p0 = [1,1,1]
        y = np.array(allstdM)
        x = np.array(allM)
        mask = np.isfinite(y)
        result = leastsq(residual, p0, args = (y[mask], x[mask]))
        a,b,c = result[0]
        print a, b, c
        
        # Complete the figure
        #ax2.legend(loc =3, ncol = 2, prop={'size':6})
        tmp = np.linspace(0,1e10, 100000)
        axarr[0].plot(tmp,tmp, c = 'gray', alpha = 1, linestyle = '--',
                zorder = 2)
        axarr[0].plot(tmp,np.sqrt(tmp*5e7), c = 'gray', alpha = 1, linestyle = '-.',
                zorder = 2)
        axarr[0].plot(tmp,a+b*tmp**c, c = 'gray', alpha = 1, linestyle = '-',
                zorder = 2)
        axarr[0].set_xlim(1e6,1e10)
        axarr[0].set_ylim(1e6,1e10)
        axarr[0].set_xscale('log')
        axarr[0].set_yscale('log')
        #ax2.invert_xaxis()
        axarr[0].set_xlabel(r'$\langle M \rangle$')
        axarr[0].set_ylabel(r'$\sqrt{\langle (\delta M)^2 \rangle}$')
        
        axarr[1].set_xlabel(r'$\langle Q_{\mathrm{mphy}} \rangle$')
        axarr[1].set_ylabel(r'$\sqrt{\langle (\delta Q_{\mathrm{mphy}})^2 \rangle}$')
        
        axarr[2].set_xlabel(r'$\langle Q_{\mathrm{tot}} \rangle$')
        axarr[2].set_ylabel(r'$\sqrt{\langle (\delta Q_{\mathrm{tot}})^2 \rangle}$')
        
        
        titlestr = (alldatestr + '\n' + args.ana + 
                    ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                    ', nens=' + str(args.nens))
        fig.suptitle(titlestr)
        plt.tight_layout(rect=[0, 0.0, 1, 0.90])
        
        plotsavestr = ('std_v_mean_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
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

    ######## Lev loop #####################
    for iz in range(len(args.height)):
        print 'lev: ', iz
        
        ####### date loop #######################
        compM_list = []
        compm_list = []
        compsize_list = []
        compN_list = []
        compQ_list = []
        for d in args.date:
            print 'Loading date: ', d
            dataset = Dataset(savedir + d + savesuf, 'r')
            tmp1 = dataset.variables['cld_size'][:,iz,:]
            compsize_list.append(np.mean(tmp1, axis = 1))
            
            tmp2 = dataset.variables['cld_sum'][:,iz,:]
            compm_list.append(np.mean(tmp2, axis = 1))
            compM_list.append(np.sum(tmp2, axis = 1))
            
            tmp3 = dataset.variables['meanQmp'][:,iz,0,0,0]
            compQ_list.append(tmp3)
            
            tmp4 = dataset.variables['meanN'][:,iz,0,0,0]
            compN_list.append(tmp4)

        
        # Get the composite means
        compsize = np.nanmean(np.array(compsize_list), axis = 0)
        compm = np.nanmean(np.array(compm_list), axis = 0)
        compM = np.nanmean(np.array(compM_list), axis = 0)
        compQ = np.nanmean(np.array(compQ_list), axis = 0)
        compN = np.nanmean(np.array(compN_list), axis = 0)
        
        lev = dataset.variables['levs'][iz]
        timelist = [timedelta(seconds=ts) for ts in dataset.variables['time']]
        timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
        # Create the figure
        fig, axarr = plt.subplots(3, 2, figsize = (95./25.4*3, 11))
        
        axarr[0,0].plot(timelist_plot, compM, c = 'orangered', linewidth = 2)
        for ic, yplot in enumerate(compM_list):
            axarr[0,0].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
        axarr[0,0].set_xlabel('time [h/UTC]')
        axarr[0,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        axarr[0,0].set_ylabel('Domain total mass flux [kg/s]')
        
        axarr[0,1].plot(timelist_plot, compsize, c = 'orangered', linewidth = 2)
        for ic, yplot in enumerate(compsize_list):
            axarr[0,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
        axarr[0,1].set_xlabel('time [h/UTC]')
        axarr[0,1].set_ylabel('Mean cloud size [m^2]')
        axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[1,0].plot(timelist_plot, compm, c = 'orangered', linewidth = 2)
        for ic, yplot in enumerate(compm_list):
            axarr[1,0].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
        axarr[1,0].set_xlabel('time [h/UTC]')
        axarr[1,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        axarr[1,0].set_ylabel('Mean cloud mass flux [kg/s]')
        
        axarr[1,1].plot(timelist_plot, compQ, c = 'orangered', linewidth = 2)
        for ic, yplot in enumerate(compQ_list):
            axarr[1,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
        axarr[1,1].set_xlabel('time [h/UTC]')
        axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])
        axarr[1,1].set_ylabel('Mean Q [h]')
        
        axarr[2,1].plot(timelist_plot, compN, c = 'orangered', linewidth = 2)
        for ic, yplot in enumerate(compN_list):
            axarr[2,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
        axarr[2,1].set_xlabel('time [h/UTC]')
        axarr[2,1].set_xlim(timelist_plot[0], timelist_plot[-1])
        axarr[2,1].set_ylabel('Mean N [h]')
        
        titlestr = (alldatestr + ', ' + args.ana + 
                    ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                    ', nens=' + str(args.nens))
        fig.suptitle(titlestr, fontsize='x-large')
        plt.tight_layout(rect=[0, 0.0, 1, 0.95])
        
        plotsavestr = ('summary_stats_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
                        '_nens-' + str(args.nens) + '_tstart-' + 
                        str(args.tstart) + '_tend-' + str(args.tend) + 
                        '_tinc-' + str(args.tinc))
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

    ######## Lev loop #####################
    for iz in range(len(args.height)):
        print 'lev: ', iz
        
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
        
        lev = dataset.variables['levs'][iz]
        timelist = [timedelta(seconds=ts) for ts in dataset.variables['time']]
        timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
        # Create the figure
        fig, axarr = plt.subplots(2, 2, figsize = (95./25.4*3, 7.))
        
        axarr[0,0].plot(timelist_plot, compprec, c = 'orangered', linewidth = 2)
        for ic, yplot in enumerate(compprec_list):
            axarr[0,0].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
        axarr[0,0].set_xlabel('time [h/UTC]')
        axarr[0,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        axarr[0,0].set_ylabel('Domain average precipitation [mm/h]')
        
        axarr[0,1].plot(timelist_plot, compcape, c = 'orangered', linewidth = 2)
        for ic, yplot in enumerate(compcape_list):
            axarr[0,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
        axarr[0,1].set_xlabel('time [h/UTC]')
        axarr[0,1].set_ylabel('Domain average CAPE [J/kg]')
        axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[1,0].plot(timelist_plot, comptauc, c = 'orangered', linewidth = 2)
        for ic, yplot in enumerate(comptauc_list):
            axarr[1,0].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
        axarr[1,0].set_xlabel('time [h/UTC]')
        axarr[1,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        axarr[1,0].set_ylabel('Domain average tau_c [h]')
        
        axarr[1,1].plot(timelist_plot, comphpbl, c = 'orangered', linewidth = 2)
        for ic, yplot in enumerate(comphpbl_list):
            axarr[1,1].plot(timelist_plot, yplot, zorder = 0.5, c = cyc[ic])
        axarr[1,1].set_xlabel('time [h/UTC]')
        axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])
        axarr[1,1].set_ylabel('Domain average PBL height [m]')
        
        titlestr = (alldatestr + ', ' + args.ana + 
                    ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                    ', nens=' + str(args.nens))
        fig.suptitle(titlestr, fontsize='x-large')
        plt.tight_layout(rect=[0, 0.0, 1, 0.95])
        
        plotsavestr = ('summary_weather_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
                        '_nens-' + str(args.nens) + '_tstart-' + 
                        str(args.tstart) + '_tend-' + str(args.tend) + 
                        '_tinc-' + str(args.tinc))
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
    

    
    Mtotlist =  np.mean(Mtotlist, axis = 0)
    Msouthlist =  np.mean(Msouthlist, axis = 0)
    Mnorthlist =  np.mean(Mnorthlist, axis = 0)
    
    # Get the required time avarage
    Mtotlist =  np.mean(Mtotlist[t1:t2], axis = 0)
    Msouthlist =  np.mean(Msouthlist[t1:t2], axis = 0)
    Mnorthlist =  np.mean(Mnorthlist[t1:t2], axis = 0)

    levs = dataset.variables['levs'][:]

    # Create the figure
    fig, ax = plt.subplots(1, 1, figsize = (95./25.4*1, 4))
    

    ax.plot(Mtotlist, args.height, c = 'r', 
                    label = 'all', linewidth = 1.5)
    ax.plot(Msouthlist, args.height, c = 'g', 
                    label = 'south', linewidth = 1.5)
    ax.plot(Mnorthlist, args.height, c = 'b', 
                    label = 'north', linewidth = 1.5)

    #ax.set_xlim(0, 2)
    #axarr[0,0].set_yscale('log')
    ax.set_ylabel('height [m]')
    ax.set_xlabel(r'$M$')
    ax.set_ylim(args.height[0], args.height[-1])
    

    ax.legend(loc =4, ncol = 1, prop={'size':6})
            
            
    titlestr = (alldatestr + '\n' + args.ana + 
                ', water=' + str(args.water) + 
                ', nens=' + str(args.nens) + ', from ' + str(UTCstart) + 
                ' to ' + str(UTCstop))
    fig.suptitle(titlestr, fontsize='x-large')
    plt.tight_layout(rect=[0, 0.0, 1, 0.93])
    
    plotsavestr = ('M_vert_' + alldatestr + '_ana-' + args.ana + 
                    '_wat-' + str(args.water) +
                    '_nens-' + str(args.nens)+ '_tstart-' + 
                    str(args.tstart) + '_tend-' + str(args.tend) + 
                    '_tinc-' + str(args.tinc) + '_tplot-' + str(t1) + 
                    '-' + str(t2))
    fig.savefig(plotdirsub + plotsavestr, dpi = 300)
    plt.close('all')

