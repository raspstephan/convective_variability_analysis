from netCDF4 import Dataset
import multiprocessing as mp
from joblib import Parallel, delayed

fn = '/home/s/S.Rasp/tmp/test.nc'

rootgroup = Dataset(fn, 'w', format='NETCDF4')

dimlen = 6
rootgroup.createDimension('dim1', dimlen)
rootgroup.createVariable('var1', 'f8', 'dim1')


def change_var(args):
    rg, i, x = args
    rg.variables['var1'][i] = x


#for i in range(dimlen):
#    change_var(rootgroup, i, i)

nproc = 3
#pool = mp.Pool(processes=nproc)

arg_list = []
for i in range(dimlen):
    arg_list.append((rootgroup, i, i))

#pool.map(change_var, arg_list)

Parallel(n_jobs=nproc)(delayed(change_var)(arg) for arg in arg_list)

print rootgroup.variables['var1'][:]
