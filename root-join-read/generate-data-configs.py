#!/usr/bin/env python

# 2019-01-09

import os

path = './'

setnum = 'SET1'
proc = 'V00-04-04'

#setnum = 'SET1'
#setnum = 'SET2'
#setnum = 'SET12'
#proc = 'V00-04-12'

datfile = path+setnum+'_'+proc+'.dat'

ws = '\t'
for i in range(1, 10):
    output = 'D' + ws + \
             str(i) + ws + \
             setnum + ws + \
             proc + ws + \
             datfile + ws + \
             './dummy.root\n'
    
    filename = 'c'+str(i)+'-data-'+setnum+'-'+proc
    with open(os.path.join('./temp', filename), 'w') as fconf:
        fconf.write(output)

