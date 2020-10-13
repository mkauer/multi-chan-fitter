#!/usr/bin/env python

# 2020-05-12

# generate nai and teflon surf for all! (and teflon bulk)

# testing teflon expo profiles - just C6
# reprocesses with alpha cuts
# face components are garbage
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
    ['teflon-surf-expo-1-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-2-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-3-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-4-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-5-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-6-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-7-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-8-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-9-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-10-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-20-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-30-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-40-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-50-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-60-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-70-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-80-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-90-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-100-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-200-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-300-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-400-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-500-in','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-600-in','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-700-in','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-800-in','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-900-in','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-1000-in','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-2000-in','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-3000-in','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-4000-in','Pb210','Pb210','GRND'],

    ['teflon-surf-expo-1-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-2-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-3-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-4-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-5-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-6-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-7-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-8-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-9-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-10-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-20-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-30-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-40-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-50-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-60-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-70-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-80-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-90-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-100-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-200-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-300-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-400-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-500-out','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-600-out','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-700-out','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-800-out','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-900-out','Pb210','Pb210','GRND'],
    ['teflon-surf-expo-1000-out','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-2000-out','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-3000-out','Pb210','Pb210','GRND'],
    #['teflon-surf-expo-4000-out','Pb210','Pb210','GRND'],
    

    ['nai-surf-expo-1-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-2-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-3-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-4-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-5-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-6-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-7-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-8-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-9-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-10-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-20-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-30-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-40-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-50-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-60-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-70-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-80-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-90-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-100-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-200-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-300-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-400-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-500-side','Pb210','Pb210','GRND'],
    ['nai-surf-expo-1000-side','Pb210','Pb210','GRND'],

    ['teflon-bulk','Pb210','Pb210','GRND'],
    
]

allTheThings = [
    ['coppercase','U238','U238','Th230'],
    ['coppercase','U238','Th230','Ra226'],
    ['coppercase','U238','Ra226','Rn222'],
    ['coppercase','U238','Rn222','Pb210'],
    ['coppercase','U238','Pb210','GRND'],
    ['coppercase','Th232','Th232','Ra228'],
    ['coppercase','Th232','Ra228','Th228'],
    ['coppercase','Th232','Th228','GRND'],

    ['teflon-bulk','U238','U238','Th230'],
    ['teflon-bulk','U238','Th230','Ra226'],
    ['teflon-bulk','U238','Ra226','Rn222'],
    ['teflon-bulk','U238','Rn222','Pb210'],
    ['teflon-bulk','U238','Pb210','GRND'],
    ['teflon-bulk','Th232','Th232','Ra228'],
    ['teflon-bulk','Th232','Ra228','Th228'],
    ['teflon-bulk','Th232','Th228','GRND'],

    ['nai-surf-10um','U238','U238','Th230'],
    ['nai-surf-10um','U238','Th230','Ra226'],
    ['nai-surf-10um','U238','Ra226','Rn222'],
    ['nai-surf-10um','U238','Rn222','Pb210'],
    ['nai-surf-10um','U238','Pb210','GRND'],
    ['nai-surf-10um','Th232','Th232','Ra228'],
    ['nai-surf-10um','Th232','Ra228','Th228'],
    ['nai-surf-10um','Th232','Th228','GRND'],
]
    
for thing in allTheThings:
    
    ws='  '
    #for i in [1,2,3,4,5,6,7,8,9]:
    #for i in [1,2,3,4,5,6,7,8]:
    #for i in [2,3,4,6,7]:
    for i in [4]:
        
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

