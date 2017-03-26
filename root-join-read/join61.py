#!/usr/bin/env python

######################################################################
# Pack all needed histos into one root file
#
# Works with v61 and later versions
# 
# version: 2017-03-26
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + now includes primPMTid cuts
# ~ changed to join61.py and build61
# ~ changed to join60.py and build60
# ~ changed to join41.py and build41
# + build everything - all/single/multi hit into one rootfile
# ~ use build40(mcfile,2) to force no reusing of data
# + mega revamp for new v40 format
#
# email me: mkauer@physics.wisc.edu
######################################################################

import os,sys
import numpy as np

import ROOT
from ROOT import *
ROOT.gROOT.Reset()
ROOT.gErrorIgnoreLevel = kWarning

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs61 import *

local = amLocal()
gROOT.SetBatch(1)

def _myself_(argv):
    
    if len(sys.argv)>1:
        mcfile = str(sys.argv[1])
    else:
        print "ERROR: specify a backgrounds file or 'data' for data"
        return
    
    build = 'build61'
    
    if local:
        outpath="/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/"+build+'/'
    else:
        outpath="/home/mkauer/mc-fitting/root-join-read/"+build+'/'
    
    fname = mcfile.split('/')[-1]

    if local: ver = 'test'
    else: ver = 'master'
    
    rootfile = build+'-'+fname+'-'+ver+'.root'
    rfile = TFile(outpath+rootfile, 'RECREATE')
    
    data1, bkgs1, sigs1 = build61(mcfile, 2, 'MSA')
    
    rfile.Write()
    rfile.Close()
    
    print "\nDONE\n"
    
    
if __name__ == "__main__":
    _myself_(sys.argv[1:])

