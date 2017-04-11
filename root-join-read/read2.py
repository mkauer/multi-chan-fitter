#!/usr/bin/env python

######################################################################
# Pack all needed histos into one root file
# hi/lo energy data should be calibrated
# hi/lo energy MC should be resolution smeared
#
# version: 2016-12-05
#
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ updated to new format of readROOT2 returning bkgs and sigs
# ~ switch to readROOT2()
# ~ switch to funcs2 import
# ~ everything should have unique keys now
# + use buildMC()
# ~ fixed hostname selection
# ~ might need a totaly different framework for the hist joining...
# ~ shit - each hist needs its own unique name for this to work,
#   I need to revisit my rootana code on how to do this
# + hacking away at it
#
# email me: mkauer@physics.wisc.edu
######################################################################
#
# Where is MC?
# /data/MC/KIMS-NaI/user-scratch/sim/processed/K40/set2
# /data/MC/KIMS-NaI/user-scratch/sim/processed/U238/set2
# /data/MC/KIMS-NaI/user-scratch/sim/processed/Th232/set2
#
# Where is raw data?
# /data/KIMS/COSINE/PHY_RUN
#
# My processed data is currently in
# /home/mkauer/temp
#
######################################################################

import os,sys
import numpy as np
from ROOT import *
import ROOT
#import funcs2
from funcs2 import *


def _myself_(argv):
    
    #energy = 0
    #rootfile = './joined-master.root'
    rootfile = './join2-test.root'
    data, bkgs, sigs = readROOT2(rootfile)
    
    canv = TCanvas('canv', 'canv', 800, 600)
    canv.SetLogy()

    this = bkgs
    
    init=1
    for name in this:
        #print name
        if not 'x1' in name: continue
        if not 'e1' in name: continue
        #if not 'pmt' in name: continue
        if not 'internal' in name: continue
        if not 'K40' in name: continue
        print name
        if init:
            this[name]['hist'].Draw()
            init=0
        else:
            this[name]['hist'].Draw('same')

    canv.Update()
    
    raw_input("Enter to quit...\n")

######################################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])

