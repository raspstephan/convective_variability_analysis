import os

dates = ['2016052800', '2016053000', '2016053100']
ana = ['p', 'm']
water = ['', 'nowater']

for d in dates:
    for a in ana:
        for w in water:
            os.system('python variance.py ' + a + ' ' + d + ' ' + w + '&')
                      #'&> /dev/null &')   # No output to console
"""            
python compute_and_save.py --ana m --date 2016052800 --nens 20 --tstart 1 --tend 24 --height 2000 3000 5000
python analyze_and_plot.py --ana m --date 2016052800 --nens 20 --tstart 1 --tend 24 --height 2000 3000 5000 --plot all 

python compute_and_save.py --ana m --date 2016053000 --nens 20 --tstart 1 --tend 24 --height 2000 3000 5000
python analyze_and_plot.py --ana m --date 2016053000 --nens 20 --tstart 1 --tend 24 --height 2000 3000 5000 --plot all 

python compute_and_save.py --ana m --date 2016053100 --nens 20 --tstart 1 --tend 24 --height 2000 3000 5000
python analyze_and_plot.py --ana m --date 2016053100 --nens 20 --tstart 1 --tend 24 --height 2000 3000 5000 --plot all 
"""
