#!/usr/bin/env python

# 2019-09-03
# for v400 and later

import os

path = './'

setnum = 'SET3'
proc = 'V00-04-14'

ws='  '
for i in range(1, 10):
    datfile = path+'C'+str(i)+'_'+str(setnum)+'_'+str(proc)+'.dat'
    output = 'D' + ws + \
             str(i) + ws + \
             setnum + ws + \
             proc + ws + \
             datfile + ws + \
             './dummy.root\n'
    
    filename = 'c'+str(i)+'-data-'+setnum+'-'+proc
    with open(os.path.join('./temp', filename), 'w') as fconf:
        fconf.write(output)

