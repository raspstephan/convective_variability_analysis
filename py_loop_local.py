import os
datelist = ['2016052800', '2016052900', '2016053000', '2016053100', '2016060100', '2016060200',
         '2016060300', '2016060400', '2016060500', '2016060600', '2016060700', '2016060800']
#datelist = ['2016060300', '2016060400', '2016060500', '2016060600', '2016060700', '2016060800']
for date in datelist:
    os.system('python compute_and_save.py --ana clouds --date ' + date + 
              ' --nens 50 --tstart 6 --tend 24 --tinc 30 --minmem 5 --dr 1 --water False &')
              
    #os.system('python ../compute_and_save.py --ana prec --date ' + date + 
              #' --nens 1 --tstart 6 --tend 24 --tinc 60 --dr 1 --water False --det True&')
