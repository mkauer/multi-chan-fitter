#!/usr/bin/env python

# 2018-12-27

import os

path = './'

#setnum = 'set1'
#proc = 'V00-04-04'

setnum = 'set1'
#setnum = 'set2'
#setnum = 'set12'
proc = 'V00-04-12'

datfile = path+setnum+'_'+proc+'.dat'

ws = '\t'
for i in range(1, 10):
    output = 'D' + ws + \
             'SM' + ws + \
             str(i) + ws + \
             datfile + ws + \
             setnum + ws + \
             './dummy.root\n'
    
    filename = 'c'+str(i)+'-data'+setnum+'-'+proc
    with open(os.path.join('./temp', filename), 'w') as fconf:
        fconf.write(output)

