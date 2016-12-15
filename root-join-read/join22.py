#!/usr/bin/env python

######################################################################
# Pack all needed histos into one root file
# hi/lo energy data should be calibrated
# hi/lo energy MC should be resolution smeared
# 
# version: 2016-12-15
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + implement Pushpa's calib and resol functions
# ~ use data run 1544
# ~ new paths for raw data and processed data
# ~ using the new version of building data and mc from backgrounds file
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

import os,sys
import numpy as np
from ROOT import *
import ROOT
#import funcs2
from funcs2 import *

local = amLocal()
gROOT.SetBatch(1)


def _myself_(argv):

    runNum = 1544
    mcfile = 'backgrounds2.txt'
    
    rootfile = 'join22-test-'+str(runNum)+'.root'

    rfile = TFile(rootfile, 'RECREATE')
    print 'creating rootfile',rootfile

    data = getData22(runNum)
    bkgs, sigs = buildMC2(mcfile, 2)
    
    rfile.Write()
    rfile.Close()
    
    print "\nDONE\n"
    

if __name__ == "__main__":
    _myself_(sys.argv[1:])

