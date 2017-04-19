import glob
import os
import argparse

parser = argparse.ArgumentParser(description = 'Process input')
#parser.add_argument('--date', metavar = 'date', type=str)
parser.add_argument('--nens', metavar = 'nens', type=int, default = 50)
args = parser.parse_args()

dates = ['2016052800', '2016052900', '2016053000', '2016053100', '2016060100', '2016060200',
         '2016060300', '2016060400', '2016060500', '2016060600', '2016060700', '2016060800']
for d in dates:
    for n in range(1,args.nens+1):
        dir = ('/home/scratch/users/stephan.rasp/' + d + 
            '/deout_ceu_pspens/' + str(n) + '/OUTPUT/')
        # Remove all .dat files
        os.system('rm ' + dir + '*.dat')
        ## Get list of all files
        #filelist = glob.glob(dir+'*')

        #for fn in filelist:
            #print fn
            #newfn = fn + '_new'
            #os.system('mv ' + fn + ' ' + newfn)
            
            #os.system('nccopy -k 4 -d 5 ' + oldfn + ' ' + fn)
            
            #os.system('rm ' + oldfn)
            
