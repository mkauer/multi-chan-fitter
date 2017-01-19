#!/usr/bin/env python

######################################################################
# Pack all needed histos into one root file
#
# Works with v32 and later versions
# 
# version: 2017-01-18
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + command line arg for the input file
# ~ use getData32() for new V00-02-00 location
# ~ also using tweaked resolution funtion for C8
# ~ import funcs32.py and use buildMC32() and use run 1546
# ~ import funcs3.py and use buildMC3()
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
import ROOT
from ROOT import *

from funcs32 import *


local = amLocal()
gROOT.SetBatch(1)

def _myself_(argv):

    if len(sys.argv)>1:
        mcfile = str(sys.argv[1])
    else:
        mcfile = 'backgrounds32.txt'

    runNum = 1546
    #mcfile = 'backgrounds32.txt'
    
    #rootfile = 'join32-'+str(runNum)+'-test.root'
    fname = mcfile.split('/')[-1].split('.')[0]
    rootfile = 'join32-'+fname+'-master.root'
    
    rfile = TFile(rootfile, 'RECREATE')
    print 'creating rootfile',rootfile
    
    #data = getData32(runNum, 'V00-02-00')
    bkgs, sigs = buildMC32(mcfile, 2)
    
    rfile.Write()
    rfile.Close()
    
    print "\nDONE\n"
    
    
if __name__ == "__main__":
    _myself_(sys.argv[1:])

