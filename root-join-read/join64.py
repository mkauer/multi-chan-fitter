#!/usr/bin/env python

######################################################################
# Pack all needed histos into one root file
#
# Works with v64 and later versions
# 
# version: 2017-06-06
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ main change is different default histogram parameters
# ~ see the changelog for version v64
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
from funcs64 import *

local = amLocal()
gROOT.SetBatch(1)

def _myself_(argv):
    
    if len(sys.argv)>1:
        mcfile = str(sys.argv[1])
    else:
        print "ERROR: specify a backgrounds file or 'data' for data"
        return
    
    build = 'build64'
    
    if local:
        outpath="/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/"+build+'/'
    else:
        outpath="/home/mkauer/mc-fitting/root-join-read/"+build+'/'
    
    fname = mcfile.split('/')[-1]
    
    if local: ver = 'test'
    else: ver = 'master'
    
    rootfile = build+'-'+fname+'-'+ver+'.root'
    rfile = TFile(outpath+rootfile, 'RECREATE')
    
    data1, bkgs1, sigs1 = build64(mcfile, 2, 'MSA')
    
    rfile.Write()
    rfile.Close()
    
    print "\nDONE\n"
    
    
if __name__ == "__main__":
    _myself_(sys.argv[1:])

