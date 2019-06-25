#!/usr/bin/env python

######################################################################
# Pack extracted histos into rootfiles
# Works with v300 and later versions
# 
# version: 2019-04-12
# 
# regenerate everything for G4.9
######################################################################

import os,sys
import numpy as np

import ROOT
from ROOT import *
ROOT.gROOT.Reset()
ROOT.gErrorIgnoreLevel = kWarning

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")

### REMEMBER to change the version here!
from funcs300 import *

local = amLocal()
gROOT.SetBatch(1)

def _myself_(argv):

    ### REMEMBER to change the version here!
    build = 'build304'
    
    xstal = 0
    
    if len(sys.argv)>1:
        mcfile = str(sys.argv[1])
    else:
        print "ERROR: specify a backgrounds or data file"
        return
    
    #if len(sys.argv)>2:
    #    xstal = int(sys.argv[2])
    
    ### select a crystal or comment out
    xstal = 9
    
    if local: outpath = "/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/"+build+'/'
    else:     outpath = "/home/mkauer/mc-fitting/root-join-read/"+build+'/'
    if xstal: outpath += 'c'+str(xstal)+'/'
    if not os.path.exists(outpath):
        try:    os.path.mkdir(outpath)
        except: os.mkdir(outpath)
    fname = mcfile.split('/')[-1]
    rootfile = build+'-'+fname+'.root'
    rfile = TFile(outpath+rootfile, 'RECREATE')
    
    ### REMEMBER to change the version here!
    if xstal: data, bkgs, sigs, runtime = build300(mcfile, 1, 2, 'SM', [xstal])
    else:     data, bkgs, sigs, runtime = build300(mcfile, 1, 2, 'SM')

    rfile.Write()
    rfile.Close()
    
    print "\nDONE\n"
    
    
if __name__ == "__main__":
    _myself_(sys.argv[1:])

