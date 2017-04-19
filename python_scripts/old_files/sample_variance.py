# Imports
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
pdfwidth = 7.87

def calc_beta(sample):
    return np.var(sample, ddof = 1)/np.mean(sample)**2

def calc_alpha(sample):
    return np.var(sample, ddof = 1)/np.mean(sample)

# 
mean_m = 4.7e7
n_iter = 10000
n_sample = (range(2,100,1) + range(100, 1000, 20) + range(1000, 2000, 200))

#alpha_list = []
beta_list = []

for ns in n_sample:
    print ns
    #tmplist_alpha = []
    tmplist_beta = []
    for ni in range(n_iter):
        
        #sample_alpha = np.random.poisson(ns, 50)
        #tmplist_alpha.append(calc_alpha(sample_alpha))
        
        sample_beta = np.random.exponential(mean_m, ns)
        tmplist_beta.append(calc_beta(sample_beta))
        
    #alpha_list.append(np.mean(tmplist_alpha))
    beta_list.append(np.mean(tmplist_beta))

#print alpha_list
#print beta_list

fn = './beta_lookup.npy'
np.save(fn, (n_sample, beta_list))

fig, ax = plt.subplots(1,1, figsize = (pdfwidth/2., 3))
ax.plot(n_sample, beta_list, linewidth = 2, c = 'k')
ax.set_xlabel('Sample size')
ax.set_ylabel(r'$\beta$')
ax.set_xscale('log')
ax.set_title(r'Dependency of $\beta$ on sample size')
plt.tight_layout()
fig.savefig('./sample_variance')

