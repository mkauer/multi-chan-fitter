#!/usr/bin/env python

######################################################################
# Pack all needed histos into one root file
#
# Works with v90 and later versions
# 
# version: 2017-11-21
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ see the changelog for version v90
#
# email: mkauer@physics.wisc.edu
######################################################################

import os,sys
import numpy as np

import ROOT
from ROOT import *
ROOT.gROOT.Reset()
ROOT.gErrorIgnoreLevel = kWarning

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs90 import *

local = amLocal()
gROOT.SetBatch(1)

def _myself_(argv):
    xstal = 0
    
    if len(sys.argv)>1:
        mcfile = str(sys.argv[1])
    else:
        print "ERROR: specify a backgrounds or data file"
        return

    # select a specific crystal
    # or comment out
    xstal = 3
    
    build = 'build90'
    
    if local:
        outpath="/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/"+build+'/'
    else:
        outpath="/home/mkauer/mc-fitting/root-join-read/"+build+'/'

    if xstal:
        outpath += 'c'+str(xstal)+'/'

    if not os.path.exists(outpath):
        try:
            os.path.mkdir(outpath)
        except:
            os.mkdir(outpath)
    
    fname = mcfile.split('/')[-1]
    
    rootfile = build+'-'+fname+'.root'
    rfile = TFile(outpath+rootfile, 'RECREATE')

    if xstal:
        data, bkgs, sigs, runtime = build90(mcfile, 2, 'MS', [xstal])
    else:
        data, bkgs, sigs, runtime = build90(mcfile, 2, 'MS')
    
    rfile.Write()
    rfile.Close()
    
    print "\nDONE\n"
    
    
if __name__ == "__main__":
    _myself_(sys.argv[1:])

