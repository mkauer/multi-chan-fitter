#!/usr/bin/env python

######################################################################
# Pack all needed histos into one root file
#
# Works with v40 and later versions
# 
# version: 2017-02-01
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + mega revamp for new v40 format
#
# email me: mkauer@physics.wisc.edu
######################################################################

import os,sys
import numpy as np
import ROOT
from ROOT import *

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs40 import *


local = amLocal()
gROOT.SetBatch(1)

def _myself_(argv):
    
    if len(sys.argv)>1:
        mcfile = str(sys.argv[1])
    else:
        print "ERROR: specify a backgrounds file or 'data' for data"
        return
    
    build = 'build40'
    
    if local:
        outpath="/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/"+build+'/'
    else:
        outpath="/home/mkauer/mc-fitting/root-join-read/"+build+'/'
    
    fname = mcfile.split('/')[-1].split('.')[0]
    rootfile = build+'-'+fname+'-test.root'
    rfile = TFile(outpath+rootfile, 'RECREATE')
    
    data, bkgs, sigs = build40(mcfile)
    
    """
    if 'data' in mcfile:
        runNum = 1546
        rootfile = 'join32-data-'+str(runNum)+'-test.root'
        rfile = TFile(outpath+rootfile, 'RECREATE')
        print 'creating rootfile',rootfile
        data = getData32(runNum, 'V00-02-00')
    else:
        fname = mcfile.split('/')[-1].split('.')[0]
        rootfile = 'join32-'+fname+'-test.root'
        rfile = TFile(outpath+rootfile, 'RECREATE')
        print 'creating rootfile',rootfile
        bkgs, sigs = buildMC32(mcfile, 2)
    """
    
    rfile.Write()
    rfile.Close()
    
    print "\nDONE\n"
    
    
if __name__ == "__main__":
    _myself_(sys.argv[1:])

