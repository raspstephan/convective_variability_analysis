"""
This script reads the saved ncdf files, does some easy analysis and produces 
plots.
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
import matplotlib.pyplot as plt


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
args = parser.parse_args()


################################################################################
# Load datatsets
datasetlist = []
heightstr = ''
alldatestr = ''
for h in args.height:
    heightstr += str(int(h)) + '_'
for d in args.date:
    alldatestr += d + '_'
    print 'Loading date: ', d
    # Create file str
    savedir = '/home/s/S.Rasp/repositories/variance/results/'
    savestr = (d + '_ana-' + args.ana + '_wat-' + str(args.water) + 
               '_height-' + heightstr +
               'nens-' + str(args.nens) + '_tstart-' + str(args.tstart) + 
               '_tend-' + str(args.tend) + '_tinc-' + str(args.tinc) + '.nc')
    # Load file 
    datasetlist.append(Dataset(savedir + savestr, 'r'))
alldatestr = alldatestr[:-1]   # revome final underscore
# End load datasets
################################################################################


# Some preliminaries
timelist = [timedelta(seconds=ts) for ts in datasetlist[0].variables['time']]
plotdir = ('/home/s/S.Rasp/Dropbox/figures/PhD/variance/' + alldatestr + 
           '/' + args.ana)


# Now comes the plotting
################################################################################
if 'cloud_stats' in args.plot:
    if len(datasetlist) > 1:
        raise Exception, 'More than one date is not implemented'
    dataset = datasetlist[0]
    plotdirsub = plotdir +  '/cloud_stats/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    # Setup
    sizemax = 2.e8
    summax = 7.5e8
    
    ############# Time loop ##############
    for it, t in enumerate(timelist):
        print 'time: ', t
        ######## Lev loop #####################
        for iz, lev in enumerate(dataset.variables['levs']):
            # Get the data
            print 'lev: ', lev
            cld_size_tmp = dataset.variables['cld_size'][it,iz, :]
            cld_size_tmp = cld_size_tmp[~cld_size_tmp.mask].data
            sizehist, sizeedges = np.histogram(cld_size_tmp, 
                                               bins = 15, range = [0., sizemax])
            sizemean = np.mean(cld_size_tmp)
            cld_sum_tmp = dataset.variables['cld_sum'][it,iz, :]
            cld_sum_tmp = cld_size_tmp[~cld_sum_tmp.mask].data
            sumhist, sumedges = np.histogram(cld_sum_tmp, 
                                             bins = 15, range = [0., summax])
            summean = np.mean(cld_sum_tmp)
            
            # Plot the histograms
            fig, axarr = plt.subplots(1, 2, figsize = (95./25.4*2.5, 4.2))
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
            axarr[1].set_xlabel('Cloud mass flux [kg/s]')
            axarr[1].set_title('Cloud mass flux distribution')
            
            titlestr = (args.date[0] + '+' + ddhhmmss(t) + ', ' + args.ana + 
                        ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                        ', nens=' + str(args.nens))
            fig.suptitle(titlestr, fontsize='x-large')
            plt.tight_layout(rect=[0, 0.0, 1, 0.95])
            
            plotsavestr = ('cloud_stats_' + args.date[0] + '_ana-' + args.ana + 
                           '_wat-' + str(args.water) + '_lev-' + str(lev) +
                           '_nens-' + str(args.nens) + '_time-' + ddhhmmss(t))
            fig.savefig(plotdirsub + plotsavestr, dpi = 300)



################################################################################
if 'rdf' in args.plot:
    if len(datasetlist) > 1:
        raise Exception, 'More than one date is not implemented'
    dataset = datasetlist[0]
    plotdirsub = plotdir +  '/rdf/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    # Setup
    ymax = 10.
    
    ############# Time loop ##############
    for it, t in enumerate(timelist):
        print 'time: ', t
        
        ######## Lev loop #####################
        for iz, lev in enumerate(dataset.variables['levs']):
            print 'lev: ', lev
            
            # Get the data
            rdf = dataset.variables['rdf'][it, iz, :]
            r =   dataset.variables['dr'][:]
            
            fig, ax = plt.subplots(1, 1, figsize = (95./25.4*1.25, 3.7))
            
            ax.plot(r/1000., rdf, c = 'g')
            ax.plot([0, np.max(r)/1000.], [1, 1], c = 'gray', alpha = 0.5)
            ax.set_xlabel('Distance [km]')
            ax.set_ylabel('Normalized RDF')
            ax.set_title('Radial distribution function')
            ax.set_ylim(0, ymax)
            ax.set_xlim(0, np.max(r)/1000.)
            
            titlestr = (args.date[0] + '+' + ddhhmmss(t) + ', ' + args.ana + 
                        ',\nwater=' + str(args.water) + ', lev= ' + str(lev) + 
                        ', nens=' + str(args.nens))
            fig.suptitle(titlestr, fontsize='x-large')
            plt.tight_layout(rect=[0, 0.0, 1, 0.85])
            
            plotsavestr = ('rdf_' + args.date[0] + '_ana-' + args.ana + 
                           '_wat-' + str(args.water) + '_lev-' + str(lev) +
                           '_nens-' + str(args.nens) + '_time-' + ddhhmmss(t))
            fig.savefig(plotdirsub + plotsavestr, dpi = 300)
            
            
################################################################################
if 'scatter' in args.plot:
    if len(datasetlist) > 1:
        raise Exception, 'More than one date is not implemented'
    dataset = datasetlist[0]
    plotdirsub = plotdir +  '/scatter/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    
    # Setup
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")
    # Load the data
    varM = dataset.variables['varM'][:]
    meanM = dataset.variables['meanM'][:]
    meanN = dataset.variables['meanN'][:]
    
    ############# Time loop ##############
    for it, t in enumerate(timelist):
        print 'time: ', t
        
        ######## Lev loop #####################
        for iz, lev in enumerate(dataset.variables['levs']):
            print 'lev: ', lev
            
            # Set up the figure 
            fig, axarr = plt.subplots(1, 2, figsize = (95./25.4*2, 3.7))
            
            ####### n loop #######################
            for i_n, n in enumerate(dataset.variables['n']):
                print 'n: ', n
                
                # Convert the data
                varM_tmp = varM[it,iz,i_n,:,:]
                varM_tmp = np.ravel(varM_tmp[~varM_tmp.mask].data)
                meanM_tmp = meanM[it,iz,i_n,:,:]
                meanM_tmp = np.ravel(meanM_tmp[~meanM_tmp.mask].data)
                meanN_tmp = meanN[it,iz,i_n,:,:]
                meanN_tmp = np.ravel(meanN_tmp[~meanN_tmp.mask].data)
                
                # Get the x and y values and filter where there are no clouds
                cld_mask = meanN_tmp > 0
                x = (varM_tmp/(meanM_tmp**2))[cld_mask]
                y = (1/meanN_tmp)[cld_mask]
                
                # Scatterplots
                axarr[0].scatter(x, y, marker = 'o', c = clist[i_n], 
                                 s = 4, zorder = 0.2, linewidth = 0, 
                                 alpha = 0.8, label = str(n*2.8)+'km')
            
            # Complete the figure
            axarr[0].legend(loc =3, ncol = 2, prop={'size':6})
            tmp = np.array([0,10])
            axarr[0].plot(tmp,tmp*np.sqrt(2), c = 'gray', alpha = 0.5, linestyle = '--',
                    zorder = 0.1)
            axarr[0].set_xlim(0.05,10)
            axarr[0].set_ylim(0.01,100)
            axarr[0].set_xscale('log')
            axarr[0].set_yscale('log')
            axarr[0].invert_xaxis()
            axarr[0].set_xlabel('Square root (1/N)')
            axarr[0].set_ylabel('Square root (Var(M)/M^2)')
            
            titlestr = (args.date[0] + '+' + ddhhmmss(t) + ', ' + args.ana + 
                        ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                        ', nens=' + str(args.nens))
            fig.suptitle(titlestr, fontsize='x-large')
            plt.tight_layout(rect=[0, 0.0, 1, 0.95])
            
            plotsavestr = ('scatter_' + args.date[0] + '_ana-' + args.ana + 
                           '_wat-' + str(args.water) + '_lev-' + str(lev) +
                           '_nens-' + str(args.nens) + '_time-' + ddhhmmss(t))
            fig.savefig(plotdirsub + plotsavestr, dpi = 300)
            

################################################################################
if 'summary_stats' in args.plot:
    plotdirsub = plotdir +  '/summary_stats/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    # Setup
    timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]

    ######## Lev loop #####################
    for iz, lev in enumerate(datasetlist[0].variables['levs']):
        print 'lev: ', lev
        
        ####### date loop #######################
        compM = []
        compm = []
        compsize = []
        for dataset in datasetlist:
            tmp1 = dataset.variables['cld_size'][:,iz,:]
            compsize.append(np.mean(tmp1, axis = 1))
            
            tmp2 = dataset.variables['cld_sum'][:,iz,:]
            compm.append(np.mean(tmp2, axis = 1))
            compM.append(np.sum(tmp2, axis = 1))
        
        # Get the composite means
        compsize = np.mean(np.array(compsize), axis = 0)
        compm = np.mean(np.array(compm), axis = 0)
        compM = np.mean(np.array(compM), axis = 0)
        
        # Create the figure
        fig, axarr = plt.subplots(2, 2, figsize = (95./25.4*3, 7.))
        
        axarr[0,0].plot(timelist_plot, compM)
        axarr[0,0].set_xlabel('time [h/UTC]')
        axarr[0,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        axarr[0,0].set_ylabel('Domain total mass flux [kg/s]')
        
        axarr[0,1].plot(timelist_plot, compsize)
        axarr[0,1].set_xlabel('time [h/UTC]')
        axarr[0,1].set_ylabel('Mean cloud size [m^2]')
        axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[1,0].plot(timelist_plot, compm)
        axarr[1,0].set_xlabel('time [h/UTC]')
        axarr[1,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        axarr[1,0].set_ylabel('Mean cloud mass flux [kg/s]')
        
        titlestr = (alldatestr + ', ' + args.ana + 
                    ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                    ', nens=' + str(args.nens))
        fig.suptitle(titlestr, fontsize='x-large')
        plt.tight_layout(rect=[0, 0.0, 1, 0.95])
        
        plotsavestr = ('summary_stats_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
                        '_nens-' + str(args.nens))
        fig.savefig(plotdirsub + plotsavestr, dpi = 300)

            
            
################################################################################
if 'summary_var' in args.plot:
    plotdirsub = plotdir +  '/summary_var/'
    if not os.path.exists(plotdirsub): os.makedirs(plotdirsub)
    # Setup
    timelist_plot = [(dt.total_seconds()/3600) for dt in timelist]
    # Clist for n
    clist = ("#ff0000", "#ff8000", "#e6e600","#40ff00","#00ffff","#0040ff",
             "#ff00ff")

    ######## Lev loop #####################
    for iz, lev in enumerate(datasetlist[0].variables['levs']):
        print 'lev: ', lev
        
        ####### date loop #######################
        compvarNoN = []
        compNvarMN = []
        compvarmomm = []
        compNvarMN_adj = []
        for i_d, dataset in enumerate(datasetlist):
            print 'Date: ', i_d
            # Load the data
            varM = dataset.variables['varM'][:]
            varN = dataset.variables['varN'][:]
            varm = dataset.variables['varm'][:]
            meanM = dataset.variables['meanM'][:]
            meanN = dataset.variables['meanN'][:]
            meanm = dataset.variables['meanm'][:]

            ########## n loop ########
            NvarMN_nlist = []
            varNoN_nlist = []
            varmomm_nlist = []
            NvarMN_adj_nlist = []
            for i_n, n in enumerate(datasetlist[0].variables['n']):
                print 'analysis n: ', n
                # Extract the correct data
                varM_tmp = varM[:,iz,i_n,:,:]
                varN_tmp = varN[:,iz,i_n,:,:]
                varm_tmp = varm[:,iz,i_n,:,:]
                meanM_tmp = meanM[:,iz,i_n,:,:]
                meanN_tmp = meanN[:,iz,i_n,:,:]
                meanm_tmp = meanm[:,iz,i_n,:,:]
                
                NvarM = varM_tmp/(meanM_tmp**2)
                NvarMN = NvarM * meanN_tmp
                
                varNoN = varN_tmp / meanN_tmp
                
                varmomm = varm_tmp / (meanm_tmp**2)
                
                NvarMN_adj = NvarMN / (1+varNoN)
                
                # Get timeseries mean
                tmp1 = []
                tmp2 = []
                tmp3 = []
                tmp4 = []
                for it in range(NvarMN.shape[0]):
                    tmp1.append(np.mean(NvarMN[it]))
                    tmp2.append(np.mean(varNoN[it]))
                    tmp3.append(np.mean(varmomm[it]))
                    tmp4.append(np.mean(NvarMN_adj[it]))
                NvarMN_nlist.append(tmp1)
                varNoN_nlist.append(tmp2)
                varmomm_nlist.append(tmp3)
                NvarMN_adj_nlist.append(tmp4)
            
            compNvarMN.append(NvarMN_nlist)
            compvarNoN.append(varNoN_nlist)
            compvarmomm.append(varmomm_nlist)
            compNvarMN_adj.append(NvarMN_adj_nlist)
        
        # Calculate composits
        compNvarMN = np.mean(compNvarMN, axis = 0) / 2.
        compvarNoN = np.mean(compvarNoN, axis = 0)
        compvarmomm = np.mean(compvarmomm, axis = 0)
        compNvarMN_adj = np.mean(compNvarMN_adj, axis = 0)
        
        # Create the figure
        fig, axarr = plt.subplots(2, 2, figsize = (95./25.4*3, 7.))
        
        for i_n, n in enumerate(datasetlist[0].variables['n']):
            axarr[0,0].plot(timelist_plot, compNvarMN[i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km')
            axarr[0,1].plot(timelist_plot, compNvarMN_adj[i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km')
            axarr[1,0].plot(timelist_plot, compvarmomm[i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km')
            axarr[1,1].plot(timelist_plot, compvarNoN[i_n], c = clist[i_n], 
                            label = str(n*2.8)+'km')
            
        axarr[0,0].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[0,0].set_xlabel('time [h/UTC]')
        axarr[0,0].set_ylabel('NVar(M)N/2')
        axarr[0,0].set_xlim(0,2)
        axarr[0,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[0,1].plot(timelist_plot, [1.]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[0,1].set_xlabel('time [h/UTC]')
        axarr[0,1].set_ylabel('NVar(M) <N> / (1+Var(N)/N)')
        axarr[0,1].set_xlim(0,2)
        axarr[0,1].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[1,0].plot(timelist_plot, [1]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[1,0].set_xlabel('time [h/UTC]')
        axarr[1,0].set_ylabel('Var(m) / m^2')
        axarr[1,0].set_xlim(0,2)
        axarr[1,0].set_xlim(timelist_plot[0], timelist_plot[-1])
        
        axarr[1,1].plot(timelist_plot, [1]*len(timelist_plot), c = 'gray', 
                        zorder = 0.1)
        axarr[1,1].set_xlabel('time [h/UTC]')
        axarr[1,1].set_ylabel('Var(N)/N')
        axarr[1,1].set_xlim(0,2)
        axarr[1,1].set_xlim(timelist_plot[0], timelist_plot[-1])
            
        axarr[1,1].legend(loc =3, ncol = 2, prop={'size':6})
              
              
        titlestr = (alldatestr + ', ' + args.ana + 
                    ', water=' + str(args.water) + ', lev= ' + str(lev) + 
                    ', nens=' + str(args.nens))
        fig.suptitle(titlestr, fontsize='x-large')
        plt.tight_layout(rect=[0, 0.0, 1, 0.95])
        
        plotsavestr = ('summary_var_' + alldatestr + '_ana-' + args.ana + 
                        '_wat-' + str(args.water) + '_lev-' + str(lev) +
                        '_nens-' + str(args.nens))
        fig.savefig(plotdirsub + plotsavestr, dpi = 300)
                
            

            
            
            
 
            
