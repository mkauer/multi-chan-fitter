#!/usr/bin/env python
######################################################################
# updateBkgs61.py
# 
# Update a background file of choice to fit results file of choice
# 
# Works with v60 and later versions
# 
# version: 2017-04-03
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ moved the main code into funcs61.py as updateBkgsFile()
# + hack something together
# 
######################################################################

import os,sys

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs61 import *


def _myself_(argv):
    
    BF = 'F'
    bkgsfile = "./backgrounds61.txt"
    resultsfile = "./plots/local_1546_Nchan-fit_loEfit-20-180_hiEfit-300-2400_hiEfitRebin-10_fitRebinScale1_mcscale1_mcsumw20_datsumw20_dru1_hiEplotRebin-10_reuse1_chansMS_v61_fit-results.txt"
    newbkgs = bkgsfile[:-4]+'-updated.txt'
    
    updateBkgsFile61(bkgsfile, resultsfile, newbkgs, BF)
        
    
if __name__ == "__main__":
    _myself_(sys.argv[1:])

