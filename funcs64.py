#!/usr/bin/env python
######################################################################
# funcs64.py
# 
# New data format and selection criteria
# 
# Works with v64 and later versions
# 
# version: 2017-06-05
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# 
# email me: mkauer@physics.wisc.edu
######################################################################
# 
# Where is MC?
# /data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/
# 
# Where is the processed data?
# Right now I'm using:
# /data/COSINE/NTP/phys/V00-02-03_MERGED
# and run ntp_I001546*
# 
######################################################################

import os,sys
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs63 import *


