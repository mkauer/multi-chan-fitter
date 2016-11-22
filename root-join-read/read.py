#!/usr/bin/env python

######################################################################
# Pack all needed histos into one root file
# hi/lo energy data should be calibrated
# hi/lo energy MC should be resolution smeared
#
# version: 2016-10-27
#
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
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
import funcs
from funcs import *


def _myself_(argv):

    energy = 0
    rootfile = './joined-master.root'
    #rootfile = './join-test.root'
    data, mc, locs, isos = readROOT(rootfile, energy)

    canv = TCanvas('canv', 'canv', 800, 600)
    canv.SetLogy()
    
    i=0
    data[i].Draw()
    for loc in locs:
        for iso in isos:
            if not mc[loc][iso]['hist']:
                print 'Warning:', loc, iso, 'not found...'
                continue
            mc[loc][iso]['hist'][i].Draw('same')
            
    
    raw_input("Enter to quit...\n")

######################################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])

