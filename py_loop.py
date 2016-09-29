import os
import sys

import argparse
# Setup - Input arguments
parser = argparse.ArgumentParser(description = 'Process input')
parser.add_argument('--ana', metavar = 'ana', type=str)
parser.add_argument('--date', metavar = 'date', type=int)
parser.add_argument('--height', metavar = 'height', type=float, default =3000)
parser.add_argument('--water', metavar = 'water', type=str, default = 'True')
parser.add_argument('--nens', metavar = 'nens', type=int, default = 20)
parser.add_argument('--tstart', metavar = 'tstart', type=int, default = 1)
parser.add_argument('--tend', metavar = 'tend', type=int, default = 24)
parser.add_argument('--tinc', metavar = 'tinc', type=int, default = 60)
parser.add_argument('--minmem', metavar = 'minmem', type=int, default = 5)
parser.add_argument('--dr', metavar = 'dr', type=int, default = 2)
args = parser.parse_args()


anastr = (' --ana ' + args.ana + ' --nens ' + str(args.nens) + 
          ' --height ' + str(args.height) + ' --water ' + str(args.water) + 
          ' --tstart ' + str(args.tstart) + ' --tend ' + str(args.tend) + 
          ' --tinc ' + str(args.tinc) + ' --minmem ' + str(args.minmem) + 
          ' --dr ' + str(args.dr))

datelist = ['2016052800', '2016052900', '2016053000', '2016053100', '2016060100', '2016060200',
         '2016060300', '2016060400', '2016060500', '2016060600', '2016060700', '2016060800']
d = datelist[args.date-1]
print 'python compute_and_save.py --date ' + d + anastr
os.system('python compute_and_save.py --date ' + d + anastr)
