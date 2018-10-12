#!/usr/bin/env python

######################################################################
# Pack extracted histos into rootfiles
# Works with v101 and later versions
# 
# version: 2018-09-04
# 
# see CHANGELOG for changes
######################################################################

import os,sys
import numpy as np

import ROOT
from ROOT import *
ROOT.gROOT.Reset()
ROOT.gErrorIgnoreLevel = kWarning

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs101 import *

local = amLocal()
gROOT.SetBatch(1)

def _myself_(argv):
    xstal = 0
    
    if len(sys.argv)>1:
        mcfile = str(sys.argv[1])
    else:
        print "ERROR: specify a backgrounds or data file"
        return

    #if len(sys.argv)>2:
    #    xstal = int(sys.argv[2])
    
    # select a specific crystal
    # or comment out
    #xstal = 2
    
    build = 'build101'
    
    if local: outpath = "/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/"+build+'/'
    else:     outpath = "/home/mkauer/mc-fitting/root-join-read/"+build+'/'

    #outpath += 'noAlphaCut/'
    
    if xstal: outpath += 'c'+str(xstal)+'/'

    if not os.path.exists(outpath):
        try:    os.path.mkdir(outpath)
        except: os.mkdir(outpath)
    
    fname = mcfile.split('/')[-1]
    
    rootfile = build+'-'+fname+'.root'
    rfile = TFile(outpath+rootfile, 'RECREATE')

    ### REMEMBER to change the build version here if it changed!!!
    if xstal: data, bkgs, sigs, runtime = build101(mcfile, 1, 2, 'SM', [xstal])
    else:     data, bkgs, sigs, runtime = build101(mcfile, 1, 2, 'SM')
    
    rfile.Write()
    rfile.Close()
    
    print "\nDONE\n"
    
    
if __name__ == "__main__":
    _myself_(sys.argv[1:])

