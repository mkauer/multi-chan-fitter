#!/usr/bin/env python

# 2019-06-24

# add in gamma Tl208 just because it's quicker this way
# surface stuff takes too long and needs to be separated
# add teflon-surf depths
# add other nai-surf depths
# steel and innersteel still need v3.1.1 MC
# lsveto still needs v3.1.1 MC
# add pmt U235
# add internal Te121
# add internal U235
# updated format for v300

import os

allTheThings = [
    #['nai-surf-10um','Pb210','Pb210','GRND'],
    #['nai-surf-1um','Pb210','Pb210','GRND'],
    #['nai-surf-0.01um','Pb210','Pb210','GRND'],
    #['teflon-bulk','Pb210','Pb210','GRND'],
    #['teflon-surf-10um','Pb210','Pb210','GRND'],
    #['teflon-surf-0.1um','Pb210','Pb210','GRND'],
    #['reflector','Pb210','Pb210','GRND'],
    #['coppercase','Co60','Co60','GRND'],
    ['gamma','Tl208','Tl208','GRND']
]


for thing in allTheThings:
    
    ws='  '
    for i in [1,2,3,4,6,7,9]:
        
        filename = 'c'+str(i)+'-'+thing[0]+thing[2]+'_'+thing[3]
        print filename
        
        output=''
        output += 'B'+ws
        output += str(i)+ws
        output += thing[0]+ws
        output += thing[1]+ws
        output += thing[2]+ws
        output += thing[3]+ws
        output += '1'+ws
        output += '0.1'+ws
        output += '10'+ws
        if thing[0] == 'lsveto' or 'steel' in thing[0]: output += 'v3.1.1'+ws
        else: output += 'v4.0.2'+ws
        ### pmt U235 still needs to be v4.0.1
        #if thing[1] == 'U235': output += 'v4.0.1'+ws
        #else: output += 'v3.1.1'+ws
        output += './dummy.root\n'
        
        with open(os.path.join('./temp', filename), 'w') as fconf:
            fconf.write(output)

