import os

dates = ['2009070100', '2016060200']
ana = ['p', 'm']
water = ['', 'nowater']

for d in dates:
    for a in ana:
        for w in water:
            os.system('python variance.py ' + a + ' ' + d + ' ' + w + '&')
                      #'&> /dev/null &')   # No output to console
            
