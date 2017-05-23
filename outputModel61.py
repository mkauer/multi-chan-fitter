#!/usr/bin/env python
######################################################################
# outputModel61.py
# 
# Create a table of the backgrounds model and activities
# 
# Works with v60 and later versions
# 
# version: 2017-04-03
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ moved the main code into funcs61.py as updateBkgsFile()
# ~ tweak the output write statements
# + hack something together
# 
######################################################################

import os,sys

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs61 import *


def _myself_(argv):

    #modelfile = "./backgrounds61-updated.txt"
    modelfile = "./backgrounds65-broken-internal.txt"
    outtable = modelfile[:-4]+'-table.txt'

    outputModelTable61(modelfile, outtable)
        

if __name__ == "__main__":
    _myself_(sys.argv[1:])

