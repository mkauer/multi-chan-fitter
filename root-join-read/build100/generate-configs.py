#!/usr/bin/env python

# 2018-06-28

import os

allTheThings = [
    ['nai-surf-10um','Pb210','Pb210','GRND'],
    ['teflon-bulk','Pb210','Pb210','GRND'],
    ['copper','Co60','Co60','GRND'],
    ['internal','Sn113','Sn113','GRND'],
    ['internal','Cd109','Cd109','GRND'],
    ['internal','H3','H3','GRND'],
    ['internal','Na22','Na22','GRND'],
    ['internal','I125','I125','GRND'],
    ['internal','Te121m','Te121m','GRND'],
    ['internal','Te123m','Te123m','GRND'],
    ['internal','Te125m','Te125m','GRND'],
    ['internal','Te127m','Te127m','GRND'],
    ['internal','K40','K40','GRND'],
    ['internal','U238','U238','GRND'],
    ['internal','U238','U238','Pb210'],
    ['internal','U238','U238','Th230'],
    ['internal','U238','Th230','Ra226'],
    ['internal','U238','Ra226','Rn222'],
    ['internal','U238','Rn222','Pb210'],
    #['internal','U238','Pb210','GRND'], # must be unique names
    ['internal','Pb210','Pb210','GRND'],
    ['internal','Th232','Th232','GRND'],
    ['internal','Th232','Th232','Ra228'],
    ['internal','Th232','Ra228','Th228'],
    ['internal','Th232','Th228','GRND'],
    ['pmt','K40','K40','GRND'],
    ['pmt','U238','U238','GRND'],
    ['pmt','U238','U238','Pb210'],
    ['pmt','U238','U238','Th230'],
    ['pmt','U238','Th230','Ra226'],
    ['pmt','U238','Ra226','Rn222'],
    ['pmt','U238','Rn222','Pb210'],
    ['pmt','U238','Pb210','GRND'],
    ['pmt','Th232','Th232','GRND'],
    ['pmt','Th232','Th232','Ra228'],
    ['pmt','Th232','Ra228','Th228'],
    ['pmt','Th232','Th228','GRND'],
    ['lsveto','K40','K40','GRND'],
    ['lsveto','U238','U238','GRND'],
    ['lsveto','U238','U238','Pb210'],
    ['lsveto','U238','U238','Th230'],
    ['lsveto','U238','Th230','Ra226'],
    ['lsveto','U238','Ra226','Rn222'],
    ['lsveto','U238','Rn222','Pb210'],
    ['lsveto','U238','Pb210','GRND'],
    ['lsveto','Th232','Th232','GRND'],
    ['lsveto','Th232','Th232','Ra228'],
    ['lsveto','Th232','Ra228','Th228'],
    ['lsveto','Th232','Th228','GRND'],
    ['innersteel','U238','U238','GRND'],
    ['innersteel','U238','U238','Pb210'],
    ['innersteel','U238','U238','Th230'],
    ['innersteel','U238','Th230','Ra226'],
    ['innersteel','U238','Ra226','Rn222'],
    ['innersteel','U238','Rn222','Pb210'],
    ['innersteel','U238','Pb210','GRND'],
    ['innersteel','Th232','Th232','GRND'],
    ['innersteel','Th232','Th232','Ra228'],
    ['innersteel','Th232','Ra228','Th228'],
    ['innersteel','Th232','Th228','GRND']
]


for thing in allTheThings:
    
    filename = thing[0]+thing[2]+'_'+thing[3]
    print filename
    
    ws='  '
    output=''
    for i in range(1, 10):
        output += 'B'+ws+'SM'+ws+str(i)+ws+thing[0]+ws+thing[1]+ws+thing[2]+ws+thing[3]+ws+'1'+ws+'0.1'+ws+'10'+ws+'./dummy.root\n'
    
    with open(os.path.join('./configs', filename), 'w') as fconf:
        fconf.write(output)


