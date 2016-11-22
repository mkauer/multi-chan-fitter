#!/usr/bin/env python

######################################################################
# Pack all needed histos into one root file
# hi/lo energy data should be calibrated
# hi/lo energy MC should be resolution smeared
#
# version: 2016-11-01
#
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + added resolution/2 to funcs.py - see my toy mc
# ~ needed absolute paths here
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

local = amLocal()
gROOT.SetBatch(1)

def _myself_(argv):
    
    rfile = TFile("join-test.root","RECREATE")
    
    for eng in range(2):
        #print eng
        
        if local: path = "/home/mkauer/COSINE/CUP/mc-fitting/data/phys/"
        else:     path = "/home/mkauer/temp/"
        
        data = getData([path+"*1324*root*"], eng)
        
        locs = ['internal','pmt']
        isos = ['K40','U238','Th232']
        mc = buildMC(locs, isos, eng)
        
        rfile.Write()
        
    rfile.Close()
    
    print "\nDONE\n"
    
######################################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])

