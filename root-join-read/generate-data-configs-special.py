#!/usr/bin/env python

# 2019-06-03
# more special case for V00-04-14 comparison
# special case to test specific run numbers

import os

#for run in [1544, 1546, 1616, 1617, 1618, 1625, 1626, 1627,
#            1634, 1652, 1654, 1666, 1671, 1672, 1675, 1678,
#            1683, 1690, 1718, 1719, 1720, 1771, 1777, 1856,
#            1858]:

for run in [1544]:

    #production = "V00-04-04"
    production = "V00-04-12"
    #production = "V00-04-14"
    run = str(run)
    datfile = './'+run+'_'+production+'.dat'

    ws='  '
    for i in range(1, 10):
        output = 'D' + ws + \
                 str(i) + ws + \
                 run + ws + \
                 production + ws + \
                 datfile + ws + \
                 './dummy.root\n'
    
        filename = 'c'+str(i)+'-data-'+run+'-'+production
        with open(os.path.join('./temp', filename), 'w') as fconf:
            fconf.write(output)

