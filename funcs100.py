#!/usr/bin/env python
######################################################################
# funcs100.py
# 
# Adding LS-veto functionality!
# 
# version: 2018-06-26
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ import funcs93
# 
# email me: mkauer@physics.wisc.edu
######################################################################

import os,sys,re
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs93 import *

