#!/usr/bin/env python

############################################################
# Matt Kauer - mkauer@icecube.wisc.edu
#-----------------------------------------------------------
# 
# 2024-01-29 - tweaked for alpha analysis
############################################################

import os,sys,re
import shutil
import socket
import copy
import math
import numpy as np
import datetime

# py2-3 compat
try:
    input = raw_input
except NameError:
    pass

# needed for ROOT-6 compat
import ctypes
from ctypes import *

script = os.path.basename(__file__)
#V = 'v'+(script.split('-')[0])
V = 'v501'
print('INFO: running script --> {0} ({1})'.format(script, V))

# check python version
import platform
py_ver = platform.python_version()
print('INFO: using python --> {0}'.format(py_ver))

# import ROOT and check version
import ROOT
from ROOT import *
try:
    root_ver = (ROOT.__version__)
except:
    print('ERROR: Requires ROOT >= 6')
    sys.exit()
print('INFO: using root --> {0}'.format(root_ver))

ROOT.gErrorIgnoreLevel = kWarning
#ROOT.gErrorIgnoreLevel = kError

# import local functions
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from funcs500 import *
from funcs_old import *
from funcs_new import *


### batch job?
#print 'HOSTNAME:', socket.gethostname()
batch = onCup()
#batch = 1

### get the total number of possible crystals (including lsveto)
numx = numX()

### print debug info
debug = 0

# draw options
dopt = 'hist same'
#dopt = 'P same'

# save histograms to rootfile?
SAVE_RAW_HISTS    = 0
SAVE_SCALED_HISTS = 0
scale_to_1mBq     = 0
SAVE_TOTAL_HISTS  = 0


### ==========  GENERAL INPUTS  ======================================
### note to add to saved plot names?
#note = 0
note = ''


#mcfile = 'backgrounds_502-2.txt' # the rn222 alpha model v20
#mcfile = 'backgrounds_502-3.txt' # the rn220 alpha model v30
#mcfile = 'backgrounds_502-4.txt' # try xpmt and innersteel
#mcfile = 'backgrounds_502-5.txt' # added pmtbase
#mcfile = 'backgrounds_502-6.txt' # added film k40
#mcfile = 'backgrounds_502-7.txt' # added vetopmt
#mcfile = 'backgrounds_502-8.txt' # testing new calib and add cu-case
#mcfile = 'backgrounds_502-9.txt' # fit xpmt without others
#mcfile = 'backgrounds_502-10.txt' # fit xpmt with others
#mcfile = 'backgrounds_502-11.txt' # fit bulk cucase
#mcfile = 'backgrounds_502-12.txt' # fix Th228 excess in LSveto
#mcfile = 'backgrounds_502-13.txt' # fix internal Ra228
#mcfile = 'backgrounds_502-14.txt' # shift I129 test (not using)
#mcfile = 'backgrounds_502-14-2.txt' # refit pmt test (not using)
#mcfile = 'backgrounds_502-12_update.txt' # stuff for Gyunho
#mcfile = 'backgrounds_502-12_test.txt' # stuff for Gyunho
#mcfile = 'backgrounds_502-data.txt' # testing new data format
#mcfile = 'backgrounds_502-15.txt' # try new Pb206
mcfile = 'backgrounds_502-15-po210.txt' # new Po210 functions for all
#mcfile = 'backgrounds_502-16.txt' # try refitting with Pb206+betas

#mcfile = 'backgrounds_502-gyunho-1.txt' # using Gyunho's data
#mcfile = 'backgrounds_502-gyunho-2.txt' # try refitting

#mcfile = 'backgrounds_502-gyunho-mlp.txt' # try refitting



#mcfile = 'testing-gyunho-data.txt'
#mcfile = 'testing-BiPo.txt'
#mcfile = 'testing-pb210-alphas.txt'
#mcfile = 'testing-neutrons.txt'
#mcfile = 'testing-c58.txt'
#mcfile = 'testing-surf-Th232.txt'
#mcfile = 'testing-pmtbase.txt'
#mcfile = 'testing-lsveto.txt'
#mcfile = 'testing-film.txt'
#mcfile = 'testing-vetopmt.txt'
#mcfile = 'testing-teflon-rn222.txt'
#mcfile = 'testing-surf-cucase.txt'
#mcfile = 'testing-internals.txt'
#mcfile = 'testing-Th228.txt'
#mcfile = 'testing-xplastic.txt'
#mcfile = 'testing-I128.txt'
#mcfile = 'testing-new-Pb206.txt'
#mcfile = 'testing-data-sets.txt'
#mcfile = 'testing-surface-alpha.txt'
#mcfile = 'testing-alphas.txt'


# get a short name for the plot dir
if mcfile.startswith('backgrounds'): shortname = mcfile[12:-4]
elif mcfile.startswith('testing'): shortname = mcfile[:-4]
else: shortname = ''

print('INFO: using backgrounds config file --> {0}'.format(mcfile))


### ==========  OPTIMIZATION OPTIONS  ================================
### test MC energy shift in bins (12bins = 1keV)
binShift = {}
###             for xstals :    1    2    3    4    5    6    7    8    9
binShift[  'teflon-Pb210'] = [  0,  -6,  -8,   0,   0, -18,  -9,   0,   0]
binShift['internal-Pb210'] = [ -7,  -5,  -5,  -8,   0,  -7,  -7,   0,   0]
binShift[ 'mystery-My007'] = [ 18,  12,  12,  18,   0,  18,  18,   0,   0]
#binShift[ 'internal-I129'] = [  0,   0,   0,   0,   0, -6,   0,   0,   0]

### shift MC for Gyunho's recalibrted data
if 'gyunho-1' in mcfile or 'gyunho-2' in mcfile:
    print('INFO: shifting Gyunho energy')
    binShift = {}
    ###             for xstals :    1    2    3    4    5    6    7    8    9
    binShift[  'teflon-Pb210'] = [  0,  -6,  -8,  -2,   0, -17,  -9,   0,   0]
    #binShift['internal-Pb210'] = [  0,   0,   0,   0,   0,  -3,   0,   0,   0]

### shift MC for Gyunho's MLP data
if 'gyunho-mlp' in mcfile:
    print('INFO: shifting MLP data energy')
    binShift = {}
    ###             for xstals :    1    2    3    4    5    6    7    8    9
    #binShift[  'teflon-Pb210'] = [  0,   0,   0,   0,   0,   0,   0,   0,   0]
    #binShift['internal-Pb210'] = [  0,   0,   0,   0,   0,   0,   0,   0,   0]



### test MC smoothing?
### smoothing window in +/- number of bins
smoothing = 0
#smoothing = 10
smoothWhat = ['film']
#smoothWhat = ['cushield', 'steel']
#smoothWhat = ['pmt', 'cushield', 'steel']

### fix MC hists nbins if things aren't consistent
fix_nbins = 1


### ==========  PLOTTING OPTIONS  ====================================
### show stuff not in a group?
showNoneGroup = 0

### plot components in groups? [0,1]
ingroups   = 1

### show the legends? [0,1]
showLeg    = 1

### show the total? [0,1]
showTotal  = 1

### plot the total in gray==0, red==1? [0,1]
redtotal   = 1

### plot the residuals? [0,1]
showResid  = 0

### show individual channel plots? [0,1]
save_indi  = 0
show_indi  = 0

### show low energy zoomed plots? [0,1]
save_zoom  = 0
show_zoom  = 0
zoomLog    = 0

### plot alpha channels? [0,1]
save_alpha = 0
show_alpha = 1
and_multi  = 0

### only show qX alpha peaks? [0,1,2]
justQ      = 0

### combine surf-expo-??? in alpha plots? [0,1]
combineExp = 1

### combine 'others' into the makePlots() plots? [0,1]
combineOth = 1

### print out chi2 of data to total? [0,1]
printChi2  = 0
gyunhoChi2 = 0

### set data label for legend? [0,'name']
dataname   = 0
#dataname   = 'Data'

### set total label for legend? [0,'name']
totname    = 0
totname    = 'Total MC'

### force the reuse of all joined rootfiles in mcfile? [0,1,2]
### very nice for debugging
### [0] default - use whatever is specified in the backgrounds file
### [1] forces reusing of all data/bkgs/sigs
### [2] forces NOT reusing any data/bkgs/sigs
reuse = 0

### helper checks
if show_indi: save_indi = 1
if show_zoom: save_zoom = 1
if show_alpha: save_alpha = 1

### ==========  FITTING OPTIONS  =====================================
### show data/hists before/after the fitting
show_prefit  = 0
show_postfit = 1
showFitSigs  = 1
showFitLeg   = 0
fitYrange = [1e-3, 50]

### skip fitting to xstals
skipfit = [5, 8]

### select channels to fit
fitchans = ['S', 'M']

### try fitting Qs separately?
### must set to 0 if not fitting alphas
fit_Q_separate = 1

### set the fitter step size
### can help to go smaller than default 1e-2
stepSize = 0
#stepSize = 1e-1
#stepSize = 1e-2
#stepSize = 1e-3
stepSize = 1e-4

### Let's finally try different fit ranges!!
fitranges = [{} for x in range(numx)]
for i in range(numx):

    ### all channels default except alphas
    # --------------------------------------------------------------
    fitranges[i]['S0'] = [2,   2,   80]  # single-hit low-energy
    fitranges[i]['S1'] = [4, 100, 3500]  # single-hit high-energy
    fitranges[i]['S2'] = [1,   0,    0]  # single-hit alphas
    fitranges[i]['M0'] = [2,   2,  100]  # multi-hit low-energy
    fitranges[i]['M1'] = [4, 100, 3500]  # multi-hit high-energy
    fitranges[i]['M2'] = [1,   0,    0]  # multi-hit alphas
    
    ### for just low E fit
    """
    #fitranges[i]['S0'] = [1, 0, 0]
    #fitranges[i]['S0'] = [1, 60, 96]
    fitranges[i]['S0'] = [1, 40, 50]
    #fitranges[i]['S0'] = [1, 2, 31]
    fitranges[i]['S1'] = [1, 0, 0]
    fitranges[i]['M0'] = [1, 0, 0]
    #fitranges[i]['M0'] = [4, 24, 40]
    fitranges[i]['M1'] = [1, 0, 0]
    # with alpha?
    #fitranges[i]['S2'] = [1, 1000, 2300]
    """
    ### for general fitting
    """
    #fitranges[i]['S0'] = [1, 0, 0]
    fitranges[i]['S0'] = [1, 65, 85]
    #fitranges[i]['S0'] = [1, 70, 85]
    #fitranges[i]['S1'] = [1, 0, 0]
    fitranges[i]['S1'] = [1, 150, 400] # Rn222, Th228, U235
    #fitranges[i]['S1'] = [1, 210, 270] # Th228 (Te123m)
    #fitranges[i]['S1'] = [1, 800, 1700] # K40, Ra228
    #fitranges[i]['S1'] = [1, 150, 4000]
    #fitranges[i]['M0'] = [1, 0, 0]
    fitranges[i]['M0'] = [1, 5, 80]
    #fitranges[i]['M0'] = [2, 40, 80]
    #fitranges[i]['M1'] = [1, 0, 0]
    fitranges[i]['M1'] = [2, 2300, 3500]
    #fitranges[i]['M1'] = [1, 2000, 4000]
    #fitranges[i]['M1'] = [2, 200, 4000]
    """
    ### 4-peak fitting Rn222, Th228, U235 pmt/cu-case/xplastic
    """
    fitranges[i]['S0'] = [1, 70, 85]
    fitranges[i]['S1'] = [1, 150, 400]
    fitranges[i]['M0'] = [1, 5, 80]
    fitranges[i]['M1'] = [2, 2500, 3500]
    """
    ### for just alpha fit
    
    fitranges[i]['S0'] = [1, 0, 0]
    #fitranges[i]['S0'] = [1, 40, 60]
    fitranges[i]['S1'] = [1, 0, 0]
    fitranges[i]['M0'] = [1, 0, 0]
    fitranges[i]['M1'] = [1, 0, 0]
    fitranges[i]['S2'] = [1, 2500, 3500]
    fitranges[i]['M2'] = [1, 0, 0]
    #fitranges[i]['M2'] = [4, 1100, 4100]
    
    ### for just lsveto fit
    """
    fitranges[i]['S0'] = [1, 0, 0]
    fitranges[i]['S1'] = [1, 0, 0]
    fitranges[i]['M0'] = [1, 0, 0]
    fitranges[i]['M1'] = [1, 0, 0]
    """
    
### special case for C1, C5, C8
#----------------------------------------
#for i in [0, 4, 7]:
#    fitranges[i]['S0'] = [1, 10, 60]
#    fitranges[i]['S1'] = [6, 0, 0]
#    fitranges[i]['M0'] = [1, 10, 90]
#    fitranges[i]['M1'] = [6, 80, 4000]

### lsveto
#----------------------------------------

fitranges[8]['S0'] = [1, 0, 0]
#fitranges[8]['S1'] = [1, 0, 0]
#fitranges[8]['S1'] = [1, 450, 4000]
fitranges[8]['S1'] = [1, 450, 5000]
#fitranges[8]['S1'] = [4, 5000, 7500]
fitranges[8]['M0'] = [1, 0, 0]
#fitranges[8]['M1'] = [1, 0, 0]
#fitranges[8]['M1'] = [1, 100, 4500]
fitranges[8]['M1'] = [1, 100, 5000]
#fitranges[8]['M1'] = [4, 5000, 7500]


# not lsveto
"""
fitranges[8]['S0'] = [1, 0, 0]
fitranges[8]['S1'] = [1, 0, 0]
fitranges[8]['M0'] = [1, 0, 0]
fitranges[8]['M1'] = [1, 0, 0]
"""


### exclude xstals from fit
#----------------------------------------
if skipfit:
    print('INFO: not fitting to {0}'.format(skipfit))
    for X in skipfit:
        i = X-1
        fitranges[i]['S0'] = [1, 0, 0]
        fitranges[i]['S1'] = [1, 0, 0]
        fitranges[i]['S2'] = [1, 0, 0]
        fitranges[i]['M0'] = [1, 0, 0]
        fitranges[i]['M1'] = [1, 0, 0]
        fitranges[i]['M2'] = [1, 0, 0]


### ==========  MORE FITTING OPTIONS  ===============================
### which MC to fit globally (to all crystals simultaneously)?
globalmc = [
    'pmt',
    'pmtbase',
    'plastic',
    'lsveto',
    'vetopmt',
    'film',
    'cushield',
    'innersteel',
    'steel',
    'gamma',
    'neutron',
]

### include bkgs from 'other' crystals? [0,1]
others = 1

### skip others for some locations
skip_locas = []
#skip_locas = ['xplastic']
#skip_locas = ['xpmt', 'xpmtbase']

### keep others for some isotopes
keep_isos = []
#keep_isos = ['Th228_GRND']

### use fit bounds from backgrounds file? [0,1,2,3]
### [0] max bounds are 0-1
### [1] use bounds specified in backgrounds file (as percent of activity)
### [2] use 'newBounds' specified below (as percent of activity)
### [3] use 'otherBnds' specified below (as a fractional scaling)
useBounds = 1

### new bounds to overwrite from file (as a percent of activity)
newBounds = [0.7, 1.3]

### else use these other bounds (as fractional scaling)
otherBnds = [1e-6, 0.9]


### ==========  MORE PLOTTING OPTIONS  =================================
### select channels to plot
pltchans = ['S', 'M']

xstal_settings = {
    'S0': {
        'xrange': [0, 100],
        'log': 1,
        'yrange': [
            [0, 20], # linear y range
            [1e-2, 100] # log y range
        ],
        'rebin': 2
    },
    'M0': {
        'xrange': [0, 100],
        'log': 1,
        'yrange': [
            [0, 5], # linear y range
            [1e-2, 50] # log y range
        ],
        'rebin': 2
    },
    'S1': {
        'xrange': [0, 4000],
        #'xrange': [3000, 9000],
        'log': 1,
        'yrange': [
            [0, 20], # linear y range
            [1e-6, 30] # log y range
            #[1e-6, 2e-3] # log y range
        ],
        'rebin': 2
        #'rebin': 20
    },
    'M1': {
        'xrange': [0, 4000],
        #'xrange': [3000, 9000],
        'log': 1,
        'yrange': [
            [0, 5], # linear y range
            [1e-6, 10] # log y range
            #[1e-6, 2e-3] # log y range
        ],
        'rebin': 2
        #'rebin': 20
    },
    'S2': {
        'xrange': [1000, 5000],
        'log': 1,
        'yrange': [
            [0, 0.8], # linear y range
            [1e-5, 2] # log y range
        ],
        'rebin': 10
    },
    'M2': {
        'xrange': [0, 6000],
        'log': 1,
        'yrange': [
            [0, 1e-3], # linear y range
            [1e-6, 1e-2] # log y range
        ],
        'rebin': 10
    }
}

ls_settings = {
    'S0': {
        'xrange': [0, 200],
        'log': 1,
        'yrange': [
            [0, 2], # linear y range
            [2e-4, 3] # log y range
        ],
        'rebin': 2
    },
    'M0': {
        'xrange': [0, 200],
        'log': 1,
        'yrange': [
            [0, 0.1], # linear y range
            [1e-5, 1] # log y range
        ],
        'rebin': 2
    },
    'S1': {
        'xrange': [0, 4500],
        'log': 1,
        'yrange': [
            [0, 2], # linear y range
            [2e-4, 3] # log y range
        ],
        'rebin': 10
    },
    'M1': {
        'xrange': [0, 4500],
        'log': 1,
        'yrange': [
            [0, 0.1], # linear y range
            [8e-7, 2e-1] # log y range
        ],
        'rebin': 10
    },
    'S2': {
        'xrange': [0, 4500],
        'log': 1,
        'yrange': [
            [0, 2], # linear y range
            [2e-4, 3] # log y range
        ],
        'rebin': 10
    },
    'M2': {
        'xrange': [0, 4500],
        'log': 1,
        'yrange': [
            [0, 0.1], # linear y range
            [8e-7, 1e-1] # log y range
        ],
        'rebin': 10
    }
}


# multi-hit test over-ride
if  0  :
    xstal_settings['S0']['xrange'] = [0, 200]
    xstal_settings['S0']['yrange'] = [[], [1e-1, 10]]
    xstal_settings['S0']['rebin']  = 6
    xstal_settings['M0']['xrange'] = [0, 200]
    xstal_settings['M0']['log']    = 0
    xstal_settings['M0']['yrange'] = [[0, 2.5], [1e-2, 3]]
    xstal_settings['M0']['rebin']  = 6
    xstal_settings['S1']['xrange'] = [0, 400]
    xstal_settings['S1']['yrange'] = [[], [1e-1, 10]]
    xstal_settings['M1']['xrange'] = [0, 400]
    xstal_settings['M1']['log']    = 0
    xstal_settings['M1']['yrange'] = [[0, 2.5], [1e-2, 3]]

# raw data in adc over-ride
if  0  :
    xstal_settings['S0']['xrange'] = False
    xstal_settings['S0']['yrange'] = False
    xstal_settings['S0']['rebin']  = 100
    xstal_settings['M0']['xrange'] = False
    xstal_settings['M0']['yrange'] = False
    xstal_settings['M0']['rebin']  = 100
    xstal_settings['S1']['xrange'] = False
    xstal_settings['S1']['yrange'] = False
    xstal_settings['S1']['rebin']  = 100
    xstal_settings['M1']['xrange'] = False
    xstal_settings['M1']['yrange'] = False
    xstal_settings['M1']['rebin']  = 100
    xstal_settings['S2']['rebin']  = 100
    xstal_settings['M2']['rebin']  = 100


### special range for zoomed in plot
zmaxE = 40

### zoomed in residual as dat-mc [0] or dat/mc [1]
zdr = 0

### use linear residual scale? [0,1]
linres = 1

### set y mean line on resid plot
lrm = 1

### set y scale on resid plot
lrs = [0, 2]
#lrs = [0.1, 10]

### main plots residual in linear scale [0,1]
liny = 0


### ==========  THINGS THAT CAN EFFECT FIT RESULTS  ==================
### scale to dru?
dru = 1

### This doesn't seem to effect the fit results at all
### set MC sumw2()? [0,1]
mcsumw2 = 0

### set data sumw2()? [0,1]
datsumw2 = 0

### set error on the total? [0,1]
toterr = 0

### chi2 test option
### ["UU", "UW", "WW", "NORM"]
chiopt = 'WU'

### set the background subtracted data errors to zero?
### also, errors are super large when scaling to DRU
### should look into this deeper at some point...
zeroFitDataError = 1


### ==========  MAIN FUNCTION  =======================================
### ==================================================================

def main(argv):
    
    #-----------------------------------------
    # set as an empty list for default action
    #-----------------------------------------
    if len(argv) > 0:
        temp = []
        for c in argv:
            temp.extend(list(c))
        temp = [int(s) for s in temp]
        xstals = list(set(temp))
        xstals.sort()
        #print 'xstals -->', xstals
    else:
        xstals = []
    #-----------------------------------------
    #-----------------------------------------
    
    gROOT.SetBatch(batch)
    gStyle.SetPalette(1)
    gStyle.SetOptStat('')
    gStyle.SetOptFit(0)

    ### where am I running?
    #if onCup(): here = '/home/mkauer/mc-fitting'
    #else: here = '/home/mkauer/COSINE/CUP/mc-fitting'
    """
    ### for saving the plots...
    plotdir = here+'/plots/c'
    for x in xstals:
        plotdir += str(x)
    if not os.path.exists(plotdir): 
        os.makedirs(plotdir)
    """
    #mcfile = os.path.join(here, mcfile)
    if not os.path.exists(os.path.join(HERE, mcfile)):
        print('ERROR: could not find backgrounds file --> {0}'.format(os.path.join(HERE, mcfile)))
        sys.exit()
    
    
    ### where everything gets loaded into dictionary
    #-----------------------------------------------------------------
    #-----------------------------------------------------------------
    #allchans = uniqString(fitchans+pltchans)
    #allchans = list(set(fitchans+pltchans))
    allchans = ['S', 'M']
    print('INFO: loading/building all the histograms')
    data, bkgs, sigs = build500(os.path.join(HERE, mcfile), others, reuse, allchans, xstals)

    ### set all hist nbins to be the same
    if fix_nbins:
        print('INFO: unifying hist bin ranges')
        data = unifyHistBins(data)
        bkgs = unifyHistBins(bkgs)
        sigs = unifyHistBins(sigs)
    
    ### only plot/fit crystals that have data
    justthese = []
    for i in range(1, numx+1):
        for key in data:
            if 'x'+str(i) in key and i not in justthese:
                justthese.append(i)
    justthese.sort()
    if len(justthese) == 0:
        print('ERROR: nothing to plot, quitting')
        sys.exit()
    print('INFO: plotting/fitting crystals --> {0}'.format(justthese))
    
    ### create dir for saving the plots...
    plotdir = HERE+'/plots/'+shortname+'_c'
    for x in justthese:
        plotdir += str(x)
    if not os.path.exists(plotdir): 
        os.makedirs(plotdir)

    ### copy the backgrounds file before fitting so it isn't overwritten later
    base_mcfile = mcfile.split('/')[-1]
    try:
        shutil.copyfile(os.path.join(HERE, mcfile), os.path.join(plotdir, base_mcfile))
    except:
        pass
    
    ### copy the actual script as well just to have it on hand
    shutil.copyfile(os.path.join(HERE, script), os.path.join(plotdir, script))

    ### save all raw histograms to a rootfile so it can be re-used
    if SAVE_RAW_HISTS:
        rootoutfile = TFile(plotdir+"/histograms-raw.root", "RECREATE")
        for hist_dict in [data, bkgs, sigs]:
            writeHists500(hist_dict)
        rootoutfile.Write()
        rootoutfile.Close()

    
    ### make of list runtimes per crystal, channel, energy
    runtimes = getRuntimes(data)
    
    # print runtimes as sanity check
    """
    for key in sortDataKeys92(data):
        X, C, E = getXCEFromKey(key)
        print('{0} = {1}'.format(key, runtimes[X][C][E]))
    return
    """
    
    ### sort keys for convenience
    datkeys = sortDataKeys92(data)
    if datsumw2:
        for key in datkeys:
            data[key]['hist'].Sumw2()
    
    ### scale into dru or not
    data = scaleData411(data, dru)
    if dru:
        if scale_to_1mBq:
            print('INFO: scaling all bkgs to 1 mBq/kg')
        bkgs = scaleBkgs411(bkgs, one_mBq=scale_to_1mBq)
        sigs = scaleBkgs411(sigs, one_mBq=scale_to_1mBq)
    else:
        bkgs = scaleBkgs411(bkgs, runtimes=runtimes)
        sigs = scaleBkgs411(sigs, runtimes=runtimes)

    ### testing MC energy scaling
    for key in binShift:
        bkgs = ScaleEnergy500(bkgs, key, binShift[key])
        sigs = ScaleEnergy500(sigs, key, binShift[key])
    
    ### make plots before combining?
    #makePlots93(bkgs, combineOth, others)
    #makePlots93(sigs, combineOth, others)
    #sys.exit()
    
    ### combine after scaling
    ### FIXME - combine but don't delete others?
    print('INFO: combining histograms')
    bkgs = combineOthers500(bkgs, globalmc, skip_locas, keep_isos)
    sigs = combineOthers500(sigs, globalmc, skip_locas, keep_isos)
    
    ### now sort and remove empty histos
    bkgs, bakkeys = sortSimKeys92(bkgs)
    sigs, sigkeys = sortSimKeys92(sigs)
    
    ### make plots after combining?
    #makePlots93(bkgs, combineOth, others)
    #makePlots93(sigs, combineOth, others)
    #sys.exit()
    
    ### do histogram smoothing?
    if smoothing:
        print('INFO: smoothing histograms --> {0}'.format(smoothWhat))
        bkgs = smooth(bkgs, smoothWhat, smoothing)
        sigs = smooth(sigs, smoothWhat, smoothing)

    ### save all scaled histograms (for other people to use?)
    if SAVE_SCALED_HISTS:
        rootoutfile = TFile(plotdir+"/histograms-scaled.root", "RECREATE")
        for hist_dict in [data, bkgs, sigs]:
            writeHists500(hist_dict)
        rootoutfile.Write()
        rootoutfile.Close()

    #-----------------------------------------------------------------
    #-----------------------------------------------------------------

    
    ### assume all data is using same runs and hist params
    try:    runtag = data[datkeys[0]]['info']['tag']
    except: runtag = 'none'
    #try:    params = globalParams(data)
    #except: params = globalParams(bkgs)

    params = globalParams()
    
    ### find unique names for color scheme?
    ### "internal-K40" for example
    uniqBkgs = []
    uniqSigs = []
    uniqAll  = []
    for key in bakkeys:
        uniqBkgs.append(key.split('-')[1]+'-'+key.split('-')[2])
        uniqAll.append(key.split('-')[1]+'-'+key.split('-')[2])
    for key in sigkeys:
        uniqSigs.append(key.split('-')[1]+'-'+key.split('-')[2])
        uniqAll.append(key.split('-')[1]+'-'+key.split('-')[2])

    uniqBkgs = sorted(list(set(uniqBkgs)))
    uniqSigs = sorted(list(set(uniqSigs)))
    uniqAll  = sorted(list(set(uniqAll)))
    #print 'INFO: Unique bkgs =',uniqBkgs
    #print 'INFO: Unique sigs =',uniqSigs
    #print 'INFO: Unique bkgs and sigs =',uniqAll
    
    ### make a string list of the globals
    globstr = ''
    if len(globalmc) == 0:
        globstr = 'none'
    else:
        for txt in globalmc:
            globstr += txt+'-'
        globstr = globstr[:-1]
    #print 'INFO: global fits to -->', globstr

    
    ### create unique colors for all bkgs
    Nc = len(uniqAll)
    print('INFO: Total number of unique bkgs and sigs = {0}'.format(Nc))
    # ROOT6 FIX
    colors, cis = rainbowSix(uniqAll)
    uniqColor = cis

    ### create unique colors for the fit hists
    sigNc = len(uniqSigs)
    sigcolors, sigcis = rainbowSix(uniqSigs)
    siguniqColor = sigcis
    
    ### colors for the groups
    gis = {
        'internal':  kBlue,
        'cosmo':     kMagenta+1,
        'surface':   kCyan+1,
        'cucase':    kYellow+1,
        #'xpmts':     kGreen,
        'pmts':      kGreen+1,
        'plastic':   kOrange,
        'lsveto':    kOrange+1,
        'vetopmt':   kGreen,
        'film':      kYellow,
        'cushield':  kOrange+4,
        'steel':     kYellow,
        #'gamma':     kOrange+2,  # now with cushield
        'neutron':   kBlue-2,
        'none':      kBlack
    }

    
    ### colors for multiple data sets
    dcs = [kBlack, kBlue, kGreen, kOrange]
    dclrs, dcs_tmp = rainbow(10)
    dcs.extend(dcs_tmp)
    
    ### legend length = MC + data + total
    Nlg = Nc+2
    lnc = 1
    if Nlg >  6: lnc = 2
    if Nlg > 12: lnc = 3
    #if Nlg > 18: lnc = 4
    xlegstop  = 0.94
    xlegstart = xlegstop-(0.2*lnc)
    #ylegstop  = 0.89
    ylegstop  = 0.91
    ylegstart = ylegstop-(0.04*6)
    
    
    ### This part puts the histos together to get ready for the fitting
    ##################################################################
    ### only do the fit if you have signals and data!
    
    resultsfile = ''
    fitting = 0
    fitStartTime = 0
    status = 0
    fitStopTime = 0
    fitTime = 0
    
    if len(sigs) > 0:
        
        print('INFO: building histograms for fit')

        fitting = 1
        
        fitsigs = {}
        fitbkgs = {}

        ### have separate global signals and backgrounds
        ### 2018-11-15
        fitglobsigs = {}
        fglobsigkeys = []
        fitglobbkgs = {}
        fglobbkgkeys = []

        #fitdata = TH1F('globData', 'globData', 1, 0, 1)

        # get/set all bin ranges first
        fitbinranges, fitbins = makeFitBinRanges(justthese, fitranges)
        #print('fit hist bins = {0}'.format(fitbins))
        #input('enter to continue')
        
        fitdata = TH1F('globData', 'globData', fitbins, 0, fitbins)
        fitdata.SetLineColor(kBlack)
        fitdata.SetMarkerColor(kBlack)
        fitdata.SetLineWidth(1)

        ftotal = TH1F('globTotal', 'globTotal', fitbins, 0, fitbins)
        ftotal.SetLineColor(kGray+1)
        ftotal.SetMarkerColor(kGray+1)
        ftotal.SetLineWidth(1)

        fresid = TH1F('globResid', 'globResid', fitbins, 0, fitbins)
        fresid.SetLineColor(kBlack)
        fresid.SetMarkerColor(kBlack)
        fresid.SetLineWidth(1)
        
        for i in [j-1 for j in justthese]:
            for C in fitchans:
                for E in ['0', '1', '2']:

                    # data
                    for dkey in datkeys:
                        if 'x'+str(i+1) in dkey and '-c'+C in dkey and '-e'+E in dkey:
                            fitdata = addToFitHist(fitdata,
                                                   cropHist(data, dkey, fitranges[i][C+E]),
                                                   fitbinranges[i][C+E])
                    
                    # signals
                    for skey in sigkeys:
                        
                        ### init histograms for global signals (pmt, lsveto, etc.)
                        bits = skey.split('-')
                        if bits[1] in globalmc:
                            fgskey = bits[1]+'-'+bits[2]
                            if fgskey not in fglobsigkeys:
                                fglobsigkeys.append(fgskey)
                                fitglobsigs[fgskey] = {}
                                fsglob = TH1F(fgskey, fgskey, fitbins, 0, fitbins)
                                fitglobsigs[fgskey]['hist'] = fsglob

                        if fit_Q_separate:
                            for Q in ['1', '2']:
                                if 'x'+str(i+1) in skey and '-c'+C in skey and '-e'+E in skey and '-q'+Q in skey:
                                    #fskey = skey.split('-c'+C+'-e'+E)[0]
                                    # need to get -q1/q2 out of main key
                                    fskey = bits[0]+'-'+bits[1]+'-'+bits[2]+'-'+bits[3]+'-'+bits[4]
                                    if fskey not in fitsigs:
                                        fitsigs[fskey] = {}
                                        stemp = TH1F(fskey, fskey, fitbins, 0, fitbins)
                                        fitsigs[fskey]['hist'] = deepcopy(stemp)
                                    if debug: print('DEBUG: adding {0} to {1}'.format(skey, fskey))
                                    fitsigs[fskey]['hist'] = addToFitHist(fitsigs[fskey]['hist'],
                                                                          cropHist(sigs, skey, fitranges[i][C+E]),
                                                                          fitbinranges[i][C+E])
                        else:
                            if 'x'+str(i+1) in skey and '-c'+C in skey and '-e'+E in skey:
                                #fskey = skey.split('-c'+C+'-e'+E)[0]
                                # need to get -q1/q2 out of main key
                                fskey = bits[0]+'-'+bits[1]+'-'+bits[2]
                                if fskey not in fitsigs:
                                    fitsigs[fskey] = {}
                                    stemp = TH1F(fskey, fskey, fitbins, 0, fitbins)
                                    fitsigs[fskey]['hist'] = deepcopy(stemp)
                                if debug: print('DEBUG: adding {0} to {1}'.format(skey, fskey))
                                fitsigs[fskey]['hist'] = addToFitHist(fitsigs[fskey]['hist'],
                                                                      cropHist(sigs, skey, fitranges[i][C+E]),
                                                                      fitbinranges[i][C+E])
                    
                    # backgrounds
                    for bkey in bakkeys:

                        ### init histograms for global backgrounds (pmt, lsveto, etc.)
                        bits = bkey.split('-')
                        if bits[1] in globalmc:
                            fgbkey = bits[1]+'-'+bits[2]
                            if fgbkey not in fglobbkgkeys:
                                fglobbkgkeys.append(fgbkey)
                                fitglobbkgs[fgbkey] = {}
                                fbglob = TH1F(fgbkey, fgbkey, fitbins, 0, fitbins)
                                fitglobbkgs[fgbkey]['hist'] = fbglob
                    
                        if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e'+E in bkey:
                            #fbkey = bkey.split('-c'+C+'-e'+E)[0]
                            # need to get -q1/q2 out of main key
                            fbkey = bits[0]+'-'+bits[1]+'-'+bits[2]
                            if fbkey not in fitbkgs:
                                fitbkgs[fbkey] = {}
                                btemp = TH1F(fbkey, fbkey, fitbins, 0, fitbins)
                                fitbkgs[fbkey]['hist'] = deepcopy(btemp)
                            fitbkgs[fbkey]['hist'] = addToFitHist(fitbkgs[fbkey]['hist'],
                                                                  cropHist(bkgs, bkey, fitranges[i][C+E]),
                                                                  fitbinranges[i][C+E])

                            


        ### ====================================================================
        ###     DEAL WITH THE GLOBAL SIGNALS AND BACKGROUNDS
        ### ====================================================================
        
        ### --------------------------------------------------------------------
        ###         GLOBAL SIGNALS 
        ### --------------------------------------------------------------------
        
        # (1) delete/sort the fsigkeys
        delete = []
        fsigkeys = []
        for fskey in fitsigs:
            if fitsigs[fskey]['hist'].Integral() > 0:
                fsigkeys.append(fskey)
            else: delete.append(fskey)
        for key in delete:
            if debug: print('DEBUG: zero events so deleting sig key [{0}]'.format(key))
            del fitsigs[key]
        fsigkeys.sort()

        # (2) sort the global sig keys
        fglobsigkeys.sort()
        
        # (3) build global sigs
        for fskey in fsigkeys:
            for fgskey in fglobsigkeys:
                #print(fskey, fgskey)
                ### make fgskey unique (pmt vs vetopmt vs xpmt) 2022-09-25
                if ('-'+fgskey) in fskey:
                    #print(fskey, ' --> ', fgskey)
                    fitglobsigs[fgskey]['hist'].Add(fitsigs[fskey]['hist'])
        
        # (4) remove global sigs from fit sigs
        L = len(fsigkeys)-1
        for k, fskey in enumerate(reversed(fsigkeys)):
            for gmckey in globalmc:
                ### make gmckey unique (steel vs innersteel) 2018-10-23
                if ('-'+gmckey+'-') in fskey:
                    if debug: print('DEBUG: remove global key from fit sigs [{0}]'.format(fsigkeys[L-k]))
                    del fsigkeys[L-k]
        
        # (5) now delete the empty sig histos
        delete = []
        fglobsigkeys = []
        for fgskey in fitglobsigs:
            if fitglobsigs[fgskey]['hist'].Integral() > 0:
                fglobsigkeys.append(fgskey)
            else: delete.append(fgskey)
        for key in delete:
            if debug: print('DEBUG: zero events so deleting global sigs key [{0}]'.format(key))
            del fitglobsigs[key]
        fglobsigkeys.sort()

        
        ### --------------------------------------------------------------------
        ###         GLOBAL BACKGROUNDS
        ### --------------------------------------------------------------------
        
        # (1) delete/sort the fbakkeys
        delete = []
        fbkgkeys = []
        for fbkey in fitbkgs:
            if fitbkgs[fbkey]['hist'].Integral() > 0:
                fbkgkeys.append(fbkey)
            else: delete.append(fbkey)
        for key in delete:
            if debug: print('DEBUG: zero events so deleting bkg key [{0}]'.format(key))
            del fitbkgs[key]
        fbkgkeys.sort()

        # (2) sort the global bkg keys
        fglobbkgkeys.sort()
        
        # (3) build global bkgs
        for fbkey in fbkgkeys:
            for fgbkey in fglobbkgkeys:
                ### make fgskey unique (pmt vs vetopmt vs xpmt) 2022-09-25
                if ('-'+fgbkey) in fbkey:
                    fitglobbkgs[fgbkey]['hist'].Add(fitbkgs[fbkey]['hist'])
        
        # (4) remove global bkgs from fit bkgs
        L = len(fbkgkeys)-1
        for k, fbkey in enumerate(reversed(fbkgkeys)):
            for gmckey in globalmc:
                ### make gmckey unique (steel vs innersteel) 2018-10-23
                if ('-'+gmckey+'-') in fbkey:
                    if debug: print('DEBUG: remove global key from fit bkgs [{0}]'.format(fbkgkeys[L-k]))
                    del fbkgkeys[L-k]

        # (5) now delete the empty bkg histos
        delete = []
        fglobbkgkeys = []
        for fgbkey in fitglobbkgs:
            if fitglobbkgs[fgbkey]['hist'].Integral() > 0:
                fglobbkgkeys.append(fgbkey)
            else: delete.append(fgbkey)
        for key in delete:
            if debug: print('DEBUG: zero events so deleting global bkgs key [{0}]'.format(key))
            del fitglobbkgs[key]
        fglobbkgkeys.sort()

        
        ### --------------------------------------------------------------------
        ###         SUBTRACT BACKGROUNDS FROM DATA
        ### --------------------------------------------------------------------

        # (1) subtract normal backgrounds
        for key in fbkgkeys:
            if debug: print('DEBUG: subtract bkg hist --> {0}'.format(key))
            fitdata.Add(fitbkgs[key]['hist'], -1)

        # (2) subtract global backgrounds
        for key in fglobbkgkeys:
            if debug: print('DEBUG: subtract global bkg hist --> {0}'.format(key))
            fitdata.Add(fitglobbkgs[key]['hist'], -1)
        
        ### ====================================================================
        ###     DONE !!!
        ### ====================================================================
        
        
        sigObj = []
        fit = []
        fitresults = {}
        fitchi2ndf = []
        wasFit = []
        bounds = []
        
        ### exclude negative bins
        print('INFO: excluding negative data bins from fit')
        for n in range(1, fitbins):
            if fitdata.GetBinContent(n) < 0:
                fitdata.SetBinContent(n, 0)
        
        ### data integral to normalize signal to
        dat_int = fitdata.Integral()
        
        ### set up the fitting object for TFractionFitter
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
            fitresults[str(i)] = []
            for fskey in fsigkeys:
                #print('signal', fskey)
                if 'x'+str(i+1) in fskey:

                    mc_int = fitsigs[fskey]['hist'].Integral()
                                        
                    ### to weight or not to weight...
                    if mcsumw2:
                        fitsigs[fskey]['hist'].Sumw2() # set stat weights
                    
                    ### normalize MC to total data
                    ### needed for TFractionFitter to work right
                    ### still don't fully understand why
                    try:
                        fitsigs[fskey]['hist'].Scale(dat_int/mc_int)
                    except:
                        print('\nERROR: No events for --> {0}'.format(fskey))
                        print('       Remove it from the fit? \n')
                        sys.exit()
                    
                    for C in allchans:
                        for E in range(3):
                            E = str(E)
                            
                            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            
                            newkeys = []
                            tmpkey = fskey+'-c'+C+'-e'+E
                            if tmpkey in sigs:
                                newkeys.append(tmpkey)
                            for Q in range(1, 3):
                                Q = str(Q)
                                tmpkey = fskey+'-q'+Q+'-c'+C+'-e'+E
                                if tmpkey in sigs:
                                    #print(tmpkey)
                                    newkeys.append(tmpkey)
                            for F in range(numx):
                                F = str(F+1)
                                tmpkey = fskey+'-f'+F+'-c'+C+'-e'+E
                                if tmpkey in sigs:
                                    newkeys.append(tmpkey)
                                for Q in range(1, 3):
                                    Q = str(Q)
                                    tmpkey = fskey+'-f'+F+'-q'+Q+'-c'+C+'-e'+E
                                    if tmpkey in sigs:
                                        newkeys.append(tmpkey)
                            """
                            try:
                                tmpkey = fskey+'-c'+C+'-e'+E
                                test = sigs[tmpkey]['hist']
                                newkeys.append(tmpkey)
                            except:
                                for F in range(numx):
                                    F = str(F+1)
                                    try:
                                        tmpkey = fskey+'-f'+F+'-c'+C+'-e'+E
                                        test = sigs[tmpkey]['hist']
                                        newkeys.append(tmpkey)
                                    except:
                                        continue
                            """
                            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            
                            for newkey in newkeys:
                                #print(newkey)
                                sigs[newkey]['hist'].Scale(dat_int/mc_int)
                                sigs[newkey]['fitscale'] = sigs[newkey]['scale'] * dat_int/mc_int
                                #print(newkey, '-->', sigs[newkey]['fitscale'])
                                
                                ### rescale the bounds to the normalized fraction
                                ### and save new values to the sigs info
                                #---------------------------------------------------------------------
                                sigs[newkey]['info']['newfbnd'] = [0,0]
                                for k in range(2):
                                    
                                    if sigs[newkey]['fitscale'] != 0:
                                        renorm = sigs[newkey]['scale'] / float(sigs[newkey]['fitscale'])
                                    
                                    if useBounds == 0:
                                        these = [0.00, 1.00]
                                        sigs[newkey]['info']['fbnd'][k] = 1./renorm * these[k]
                                        sigs[newkey]['info']['newfbnd'] = these
                                    elif useBounds == 1:
                                        sigs[newkey]['info']['newfbnd'][k] = \
                                            sigs[newkey]['info']['fbnd'][k] * renorm
                                    elif useBounds == 2:
                                        sigs[newkey]['info']['fbnd'] = newBounds
                                        sigs[newkey]['info']['newfbnd'][k] = \
                                            newBounds[k] * renorm
                                    elif useBounds == 3:
                                        sigs[newkey]['info']['fbnd'][k] = 1./renorm * otherBnds[k]
                                        sigs[newkey]['info']['newfbnd'] = otherBnds
                                    else:
                                        print('ERROR: do not know what to do with useBounds = {0}'
                                              .format(useBounds))
                                        sys.exit()
                    #print '!!!', sigs[newkeys[0]]['info']['newfbnd']
                    #bounds.append(sigs[newkeys[0]]['info']['newfbnd'])
                    # or just use last good key
                    #print('!!!!!!!! last good key --> {0}'.format(newkey))
                    bounds.append(sigs[newkey]['info']['newfbnd'])
                    #---------------------------------------------------------------------
                    
                    ### set errors to zero?
                    fitsigs[fskey]['hist'] = zeroBinError(fitsigs[fskey]['hist'])

            
            fitresults[str(i)].append('Crystal-'+str(i+1)+' fit results')
            fitresults[str(i)].append('runtime = '+str(round(runtimes[i]['S'][0]/60./60./24., 2))+' days')
            if note: fitresults[str(i)].append('note = '+str(note))
            fitresults[str(i)].append('script = '+script)
            fitresults[str(i)].append('python = '+str(py_ver))
            fitresults[str(i)].append('root = '+str(root_ver))
            fitresults[str(i)].append('selected xstals = '+str(justthese))
            fitresults[str(i)].append('channels fit = '+str(fitchans))
            fitresults[str(i)].append('globals = '+str(globalmc))
            #fitresults[str(i)].append('use others = '+str(others))
            #fitresults[str(i)].append('hist extend = '+str(extend))
            fitresults[str(i)].append('norm to dru = '+str(dru))
            #fitresults[str(i)].append('energy shift in bins = '+str(binShift1[i]))
            #fitresults[str(i)].append('shifting these = '+str(shiftWhat1))
            fitresults[str(i)].append('fitter step size = '+str(stepSize))
            for chan in ['S0', 'S1', 'S2', 'M0', 'M1', 'M2']:
                fitresults[str(i)].append('{0} fit --> rebin = [{1}] --> range = [{2} - {3}] keV'
                                          .format(chan,
                                                  fitranges[i][chan][0],
                                                  fitranges[i][chan][1],
                                                  fitranges[i][chan][2]))
        
        
        ### do the same for the global sigs
        boundskeys = []
        for fgkey in fglobsigkeys:
            #print('global', fgkey)
            mc_int = fitglobsigs[fgkey]['hist'].Integral()
            
            ### to weight or not to weight...
            if mcsumw2:
                fitglobsigs[fgkey]['hist'].Sumw2() # set stat weights
            
            ### normalize MC to total data
            ### needed for TFractionFitter to work right
            ### still don't fully understand why
            try:
                fitglobsigs[fgkey]['hist'].Scale(dat_int/mc_int)
            except:
                print('\nERROR: No events for --> {0}'.format(fgkey))
                print('       Remove it from the fit!\n')
                sys.exit()
            
            #for i in range(numx):
            for i in [j-1 for j in justthese]:
                X = str(i+1)
                for C in allchans:
                    for E in range(3):
                        E = str(E)
                        
                        newkeys=[]
                        tmpkey = 'x'+X+'-'+fgkey+'-c'+C+'-e'+E
                        if tmpkey in sigs:
                            newkeys.append(tmpkey)
                        for F in range(numx):
                            F=str(F+1)
                            tmpkey = 'x'+X+'-'+fgkey+'-f'+F+'-c'+C+'-e'+E
                            if tmpkey in sigs:
                                newkeys.append(tmpkey)
                        """
                        try:
                            tmpkey = 'x'+X+'-'+fgkey+'-c'+C+'-e'+E
                            test = sigs[tmpkey]['hist']
                            newkeys.append(tmpkey)
                        except:
                            for F in range(numx):
                                F=str(F+1)
                                try:
                                    tmpkey = 'x'+X+'-'+fgkey+'-f'+F+'-c'+C+'-e'+E
                                    test = sigs[tmpkey]['hist']
                                    newkeys.append(tmpkey)
                                except:
                                    continue
                        """
                        
                        if len(newkeys) > 0:
                            for newkey in newkeys:
                                
                                sigs[newkey]['hist'].Scale(dat_int/mc_int)
                                sigs[newkey]['fitscale'] = sigs[newkey]['scale'] * dat_int/mc_int
                                                                
                                ### rescale the bounds to the normalized fraction
                                ### and save new values to the sigs info
                                #---------------------------------------------------------------------
                                sigs[newkey]['info']['newfbnd'] = [0,0]
                                for k in range(2):
                                    renorm = sigs[newkey]['scale'] / float(sigs[newkey]['fitscale'])
                                    
                                    if useBounds == 0:
                                        these = [0.00, 1.00]
                                        sigs[newkey]['info']['fbnd'][k] = 1./renorm * these[k]
                                        sigs[newkey]['info']['newfbnd'] = these
                                    elif useBounds == 1:
                                        sigs[newkey]['info']['newfbnd'][k] = \
                                            sigs[newkey]['info']['fbnd'][k] * renorm
                                    elif useBounds == 2:
                                        sigs[newkey]['info']['fbnd'] = newBounds
                                        sigs[newkey]['info']['newfbnd'][k] = \
                                            newBounds[k] * renorm
                                    elif useBounds == 3:
                                        sigs[newkey]['info']['fbnd'][k] = 1./renorm * otherBnds[k]
                                        sigs[newkey]['info']['newfbnd'] = otherBnds
                                    else:
                                        print('ERROR: do not know what to do with useBounds = {0}'
                                              .format(useBounds))
                                        sys.exit()
                                #---------------------------------------------------------------------

                            
                            ### add bounds to unique bouds list
                            if newkeys[0] and fgkey not in boundskeys:
                                #print 'bounds for', fgkey, sigs[newkey]['info']['newfbnd']
                                boundskeys.append(fgkey)
                                bounds.append(sigs[newkeys[0]]['info']['newfbnd'])
            
            
            ### set errors to zero?
            fitglobsigs[fgkey]['hist'] = zeroBinError(fitglobsigs[fgkey]['hist'])
        
        
        ### conflict of interests going on here
        #=========================================================
        if datsumw2:
            fitdata.Sumw2()

        if zeroFitDataError:
            fitdata = zeroBinError(fitdata)
        #=========================================================

            
        ### exclude negative bins
        print('INFO: excluding negative MC bins from fit')
        for n in range(1, fitbins):
            if fitdata.GetBinContent(n) <= 0:
                for fskey in fsigkeys:
                    fitsigs[fskey]['hist'].SetBinContent(n, 0)
                for fgkey in fglobsigkeys:
                    fitglobsigs[fgkey]['hist'].SetBinContent(n, 0)
                
        totalNumFits = len(fsigkeys)+len(fglobsigkeys)
        print('INFO: total number of hists being fit = {0}'.format(totalNumFits))
        if totalNumFits < 2:
            print('ERROR: need at least 2 MC histograms')
            return
        sigObj = TObjArray(totalNumFits)
        for fskey in fsigkeys:
            #print('loading signal to fit -->', fskey)
            #print(fitsigs[fskey]['hist'].GetNbinsX())
            sigObj.append(fitsigs[fskey]['hist'])
        for fgkey in fglobsigkeys:
            #print('loading global to fit -->', fgkey)
            sigObj.append(fitglobsigs[fgkey]['hist'])
        
        # pre-fit plotting to see how it looks?
        if show_prefit:
            pfcanv = TCanvas('pfcanv', 'pfcanv', 0, 0, 1400, 900)
            fitdata.Draw()
            for fskey in fsigkeys:
                cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
                fitsigs[fskey]['hist'].SetMarkerColor(siguniqColor[cname])
                fitsigs[fskey]['hist'].SetLineColor(siguniqColor[cname])
                fitsigs[fskey]['hist'].Draw(dopt)
            for fgkey in fglobsigkeys:
                cname = fgkey
                fitglobsigs[fgkey]['hist'].SetMarkerColor(siguniqColor[cname])
                fitglobsigs[fgkey]['hist'].SetLineColor(siguniqColor[cname])
                fitglobsigs[fgkey]['hist'].Draw(dopt)
            pfcanv.SetLogy(1)
            pfcanv.Update()
            pfsave  = ''
            pfsave += str(runtag)
            pfsave += '_pre-fit'
            if note: pfsave += '_'+str(note)
            pfsave += '_'+V
            pfcanv.Print(plotdir+'/'+pfsave+'.png')
        
        
        #print 'data', fitdata.GetNbinsX()
        
        # init the fitter class
        fit = TFractionFitter(fitdata, sigObj)
        fitter = fit.GetFitter()
        #sys.exit()
        
        ### set fit bounds!!!
        for l in range(len(bounds)):
            fit.Constrain(l, bounds[l][0], bounds[l][1])
            if stepSize != 0:
                fitter.Config().ParSettings(l).SetStepSize(stepSize)
        
        ### set the fit range
        # not sure if bin 0 is underflow in this case? Yes it is!
        #fit.SetRangeX(0, fitbins-1)
        fit.SetRangeX(1, fitbins)
        #print('fitdata bin {0} = {1}'.format(0, fitdata.GetBinContent(0)))
        #print('fitdata bin {0} = {1}'.format(1, fitdata.GetBinContent(1)))
        #print('fitdata bin {0} = {1}'.format(fitbins, fitdata.GetBinContent(fitbins)))
        #print('fitdata bin {0} = {1}'.format(fitbins+1, fitdata.GetBinContent(fitbins+1)))
        

        ### if you want to get initial fit values from fitter...
        #fitter = fit.GetFitter()
        #fit_values = fitter.Config().ParamsValues()
        #fit_values = np.asarray(fit_values)
        #print fit_values

        
        #=======================================================================
        #        MACHEN SIE DAS FIT!!!
        #=======================================================================
        print('INFO: starting the fit')
        fitStartTime = datetime.datetime.now()
        status = fit.Fit()
        #status = 0  # force the plotting for testing
        status = int(status)  # because 6.14.00 returns <ROOT.TFitResultPtr object>
        fitStopTime = datetime.datetime.now()
        fitTime = int((fitStopTime-fitStartTime).total_seconds())
        #=======================================================================

        if status != 0:
            print('\n\n*******************  FIT FAILURE  *******************\n')
            print('Time to fail = {} sec \n'.format(fitTime))
            #input('\n[Enter] to quit\n')
            #sys.exit()
            return
        print('\n\n*******************  SUCCESSFUL FIT  *******************\n\n')
        
        
        chi2 = fit.GetChisquare()
        ndf  = fit.GetNDF()
        pval = fit.GetProb()

        fitresult = fit.GetPlot()
        
        #print('chi2 =', chi2)
        #print('ndf =', ndf)
        #print 'pval =', pval
        
        ### get non-zero chi2 and ndf
        NDF = 0
        CHI2 = 0
        for n in range(1, fitbins+1):
            if fitdata.GetBinContent(n) > 0 \
               and fitresult.GetBinContent(n) > 0:
                NDF += 1
                #CHI2 += abs(fitresult.GetBinContent(n) - fitdata.GetBinContent(n))**2 / fitdata.GetBinContent(n)
                CHI2 += (fitdata.GetBinContent(n) - fitresult.GetBinContent(n))**2 / fitresult.GetBinContent(n)
                
        #print('Fit Chi^2 = {0}'.format(CHI2))
        
        #ndf = NDF
        ndf = totalNumFits
        chi2 = CHI2
        try:
            fitchi2ndf = (chi2/ndf)
        except:
            fitchi2ndf = 0
            
        count = 0
        
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
            
            fitresults[str(i)].append('total number of hists being fit = '+str(totalNumFits))
            #fitresults[str(i)].append('returned fit status = '+str(status))
            fitresults[str(i)].append('chi2/ndf = %.3g/%s = %.3g'%(chi2, ndf, fitchi2ndf))
            fitresults[str(i)].append('time to fit = '+str(fitTime)+' seconds')
            
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    
                    #fscale = ROOT.Double(0.0)
                    #ferror = ROOT.Double(0.0)
                    # ROOT6 FIX
                    fscale = ctypes.c_double(0.0)
                    ferror = ctypes.c_double(0.0)
                    
                    #print 'count',count
                    fit.GetResult(count, fscale, ferror)
                    # ROOT6 FIX
                    fscale = float(fscale.value)
                    ferror = float(ferror.value)
                    
                    count += 1
                    
                    fitsigs[fskey]['hist'].Scale(fscale)

                    # this leaves out -q1/q2 [2021-07-06]
                    """
                    for C in allchans:
                        for E in range(3):
                            E = str(E)
                            #print fskey
                            newkey = fskey+'-c'+C+'-e'+E
                            try:
                                ### save the raw scaling factor from the fit
                                sigs[newkey]['hist'].Scale(fscale)
                                
                                ### save converted scaling factor
                                sigs[newkey]['fitscale'] = sigs[newkey]['fitscale'] * fscale
                            except:
                                #print '!!! could not find', newkey
                                continue

                            try:
                                ### set error as a percent of the scaling factor
                                sigs[newkey]['fiterror'] = ferror/fscale
                            except:
                                sigs[newkey]['fiterror'] = 1.
                    """
                    # do this a smarter way [2021-07-06]
                    for skey in sigs:
                        if fskey in skey:
                            if debug: print('DEBUG: scaling {0} from fit key {1}'.format(skey, fskey))
                            ### scale the hist by the fit factor
                            sigs[skey]['hist'].Scale(fscale)
                                
                            ### save converted scaling factor
                            #print('signal', skey, fskey)
                            sigs[skey]['fitscale'] = sigs[skey]['fitscale'] * fscale
                            
                            ### set error as a percent of the scaling factor
                            sigs[skey]['fiterror'] = ferror/fscale

                        
        ### do the same for the globals lsveto
        for fgkey in fglobsigkeys:

            #fscale = ROOT.Double(0.0)
            #ferror = ROOT.Double(0.0)
            # ROOT6 FIX
            fscale = ctypes.c_double(0.0)
            ferror = ctypes.c_double(0.0)

            #print 'count',count
            fit.GetResult(count, fscale, ferror)
            # ROOT6 FIX
            fscale = float(fscale.value)
            ferror = float(ferror.value)
            
            count += 1
            
            fitglobsigs[fgkey]['hist'].Scale(fscale)

            # this is really messy [2021-07-06]
            """
            for i in [j-1 for j in justthese]:
                X = str(i+1)
                for C in allchans:
                    for E in range(3):
                        E = str(E)
                        
                        newkeys=[]
                        try:
                            tmpkey = 'x'+X+'-'+fgkey+'-c'+C+'-e'+E
                            test = sigs[tmpkey]['hist']
                            newkeys.append(tmpkey)
                        except:
                            for F in range(numx):
                                F=str(F+1)
                                try:
                                    tmpkey = 'x'+X+'-'+fgkey+'-f'+F+'-c'+C+'-e'+E
                                    test = sigs[tmpkey]['hist']
                                    newkeys.append(tmpkey)
                                except:
                                    #tmpkey=0
                                    continue
                        
                        if len(newkeys) > 0:
                            for newkey in newkeys:
                                ### save the raw scaling factor from the fit
                                sigs[newkey]['hist'].Scale(fscale)

                                ### save converted scaling factor
                                sigs[newkey]['fitscale'] = sigs[newkey]['fitscale'] * fscale

                                ### set error as a percent of the scaling factor
                                try:
                                    sigs[newkey]['fiterror'] = ferror/fscale
                                except:
                                    sigs[newkey]['fiterror'] = 1.
            """
            # do this a smarter way [2021-07-06]
            for skey in sigs:
                # need to add "-" to fgkey otherwise pmt is in xpmt and xpmt is double scaled
                # fixed this 2022-06-29
                if '-'+fgkey in skey:
                    if debug: print('DEBUG: scaling {0} from fit key {1}'.format(skey, fgkey))
                    ### scale the hist by the fit factor
                    sigs[skey]['hist'].Scale(fscale)

                    ### save converted scaling factor
                    #print('global', skey, fgkey)
                    sigs[skey]['fitscale'] = sigs[skey]['fitscale'] * fscale

                    ### set error as a percent of the scaling factor
                    sigs[skey]['fiterror'] = ferror/fscale

            
        ### scale the signals to mBq/kg
        if dru:
            #print 'dru scaling'
            sigs = scaleSigs411(sigkeys, sigs)
        else:
            #print 'runtime scaling'
            sigs = scaleSigs411(sigkeys, sigs, runtimes)

        #print sigs['x6-internal-Pb210_GRND-f6-q1-cS-e2']['info']
        #print sigs['x6-internal-Pb210_GRND-f6-q2-cS-e2']['info']
        
        ### print the fit activities
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
            for fskey in fsigkeys:
                #print 'fskey =', fskey
                if 'x'+str(i+1) in fskey:
                    finit=1
                    for C in allchans:
                        for E in range(3):
                            E = str(E)
                            if finit:
                                """
                                newkey = fskey+'-c'+C+'-e'+E
                                try: test = sigs[newkey]['hist']
                                except: continue
                                """
                                newkeys = []
                                
                                # make this smarter
                                # no more try:except
                                """
                                try:
                                    tmpkey = fskey+'-c'+C+'-e'+E
                                    test = sigs[tmpkey]['hist']
                                    newkeys.append(tmpkey)
                                except:
                                    for F in range(numx):
                                        F=str(F+1)
                                        try:
                                            tmpkey = fskey+'-f'+F+'-c'+C+'-e'+E
                                            test = sigs[tmpkey]['hist']
                                            newkeys.append(tmpkey)
                                        except:
                                            #tmpkey=0
                                            continue
                                """
                                
                                tmpkey = fskey+'-c'+C+'-e'+E
                                if tmpkey in sigs:
                                    newkeys.append(tmpkey)
                                for Q in range(1, 3):
                                    Q = str(Q)
                                    tmpkey = fskey+'-q'+Q+'-c'+C+'-e'+E
                                    if tmpkey in sigs:
                                        #print(tmpkey)
                                        newkeys.append(tmpkey)
                                for F in range(numx):
                                    F = str(F+1)
                                    tmpkey = fskey+'-f'+F+'-c'+C+'-e'+E
                                    if tmpkey in sigs:
                                        newkeys.append(tmpkey)
                                    for Q in range(1, 3):
                                        Q = str(Q)
                                        tmpkey = fskey+'-f'+F+'-q'+Q+'-c'+C+'-e'+E
                                        if tmpkey in sigs:
                                            newkeys.append(tmpkey)

                                        
                                if len(newkeys) > 0:
                                    
                                    #for key in newkeys:
                                    #    print('{0} --> {1}'.format(key, sigs[key]['info']['fitacti']))
                                    
                                    newkey = newkeys[0]
                                    
                                    #print 'found', newkey
                                    #print newkey, sigs[newkey]['info']['fbnd']

                                    ### print out activity and error and bounds
                                    fitacti = sigs[newkey]['info']['fitacti']
                                    fiterro = sigs[newkey]['info']['fiterro']
                                    lobnd = sigs[newkey]['info']['acti'] * \
                                            sigs[newkey]['info']['fbnd'][0]
                                    hibnd = sigs[newkey]['info']['acti'] * \
                                            sigs[newkey]['info']['fbnd'][1]

                                    ### check if at limit
                                    limit = 0
                                    if float('%.2e'%(fitacti)) <= float('%.2e'%(lobnd)):
                                        limit = '[LOWER LIMIT]'
                                    if float('%.2e'%(fitacti)) >= float('%.2e'%(hibnd)):
                                        limit = '[UPPER LIMIT]'
                                    if limit:
                                        fitresults[str(i)].append(
                                            # show the actual error
                                            '%43s = %.2e +/- %.2e mBq  (%.2e, %.2e) %s'
                                            %('fit '+fskey, fitacti, fiterro, lobnd, hibnd, limit))
                                            # show error as %
                                            # oh this messes up my pretty plotting of fit activities
                                            #'%42s = %.2e mBq (%.1f%%) %s'
                                            #%('fit '+fskey, fitacti, 100*fiterro/fitacti, limit))

                                    else:
                                        fitresults[str(i)].append(
                                            # show the actual error
                                            '%43s = %.2e +/- %.2e mBq  (%.2e, %.2e)'
                                            %('fit '+fskey, fitacti, fiterro, lobnd, hibnd))
                                            # show error as %
                                            # oh this messes up my pretty plotting of fit activities
                                            #'%42s = %.2e mBq (%.1f%%)'
                                            #%('fit '+fskey, fitacti, 100*fiterro/fitacti))

                                    ### turn off
                                    finit = 0

                                
        ### do the same for the globals
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
            for fgkey in fglobsigkeys:
                X=str(i+1)
                finit=1
                for C in allchans:
                    for E in range(3):
                        E=str(E)
                        if finit:
                            
                            newkeys = []

                            # make this smarter
                            # no more try:except
                            """
                            try:
                                tmpkey = 'x'+X+'-'+fgkey+'-c'+C+'-e'+E
                                test = sigs[tmpkey]['hist']
                                newkeys.append(tmpkey)
                            except:
                                for F in range(numx):
                                    F=str(F+1)
                                    try:
                                        tmpkey = 'x'+X+'-'+fgkey+'-f'+F+'-c'+C+'-e'+E
                                        test = sigs[tmpkey]['hist']
                                        newkeys.append(tmpkey)
                                    except:
                                        #tmpkey=0
                                        continue
                            """
                            
                            tmpkey = 'x'+X+'-'+fgkey+'-c'+C+'-e'+E
                            if tmpkey in sigs:
                                newkeys.append(tmpkey)
                            for F in range(numx):
                                F=str(F+1)
                                tmpkey = 'x'+X+'-'+fgkey+'-f'+F+'-c'+C+'-e'+E
                                if tmpkey in sigs:
                                    newkeys.append(tmpkey)
                            
                            
                            if len(newkeys) > 0:

                                #for key in newkeys:
                                #    print('{0} --> {1}'.format(key, sigs[key]['info']['fitacti']))
                                    
                                # special case
                                #--------------------------------------
                                # print/save neutron multi-hit activity
                                found = False
                                good = False
                                for key in newkeys:
                                    if 'neutron' in key:
                                        found = True
                                        if '-cM-e1' in key:
                                            good = True
                                if found and not good: continue
                                #--------------------------------------
                                
                                newkey = newkeys[0]
                                                                
                                ### print out activity and error and bounds
                                fitacti = sigs[newkey]['info']['fitacti']
                                fiterro = sigs[newkey]['info']['fiterro']
                                lobnd = sigs[newkey]['info']['acti'] * \
                                        sigs[newkey]['info']['fbnd'][0]
                                hibnd = sigs[newkey]['info']['acti'] * \
                                        sigs[newkey]['info']['fbnd'][1]

                                ### check if at limit
                                limit = 0
                                if float('%.2e'%(fitacti)) <= float('%.2e'%(lobnd)):
                                    limit = '[LOWER LIMIT]'
                                if float('%.2e'%(fitacti)) >= float('%.2e'%(hibnd)):
                                    limit = '[UPPER LIMIT]'
                                if limit:
                                    fitresults[str(i)].append(
                                        # show error in mBq
                                        '%43s = %.2e +/- %.2e mBq  (%.2e, %.2e) %s'
                                        #%('fit '+'x'+X+'-'+fgkey, fitacti, fiterro, lobnd, hibnd, limit))
                                        # print as 'x0'
                                        %('fit '+'x'+'0'+'-'+fgkey, fitacti, fiterro, lobnd, hibnd, limit))
                                        # show error as %
                                        # oh this messes up my pretty plotting of fit activities
                                        #'%42s = %.2e mBq (%.1f%%) %s'
                                        #%('fit '+'x'+X+'-'+fgkey, fitacti, 100*fiterro/fitacti, limit))
                                
                                else:
                                    fitresults[str(i)].append(
                                        # show error in mBq
                                        '%43s = %.2e +/- %.2e mBq  (%.2e, %.2e)'
                                        #%('fit '+'x'+X+'-'+fgkey, fitacti, fiterro, lobnd, hibnd))
                                        # print as 'x0'
                                        %('fit '+'x'+'0'+'-'+fgkey, fitacti, fiterro, lobnd, hibnd))
                                        # show error as %
                                        # oh this messes up my pretty plotting of fit activities
                                        #'%42s = %.2e mBq (%.1f%%)'
                                        #%('fit '+'x'+X+'-'+fgkey, fitacti, 100*fiterro/fitacti))
                                
                                ### turn off
                                finit = 0
            
            fitresults[str(i)].append('\n')
            
            
        save = ''
        #if local: save += 'local'
        #else:     save += 'on-cup'
        save += str(runtag)
        save += '_Nchan-fit'
        #save += '_globals-'+globstr
        #save += '_loEfit-'+str(int(fLoE[0]))+'-'+str(int(fLoE[1]))
        #save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        #save += '_loEfitRebin-'+str(loEfitRebin)
        #save += '_hiEfitRebin-'+str(hiEfitRebin)
        #save += '_loEfitRebinScale'+str(loEfitRebinScale)
        #save += '_hiEfitRebinScale'+str(hiEfitRebinScale)
        #save += '_useBounds'+str(useBounds)
        #save += '_mcsumw2'+str(mcsumw2)
        #save += '_datsumw2'+str(datsumw2)
        #save += '_dru'+str(dru)
        #save += '_loEplotRebin-'+str(loEplotRebin)
        #save += '_hiEplotRebin-'+str(hiEplotRebin)
        #save += '_reuse'+str(reuse)
        #save += '_chans'+str(fitchans)
        #save += '_extend'+str(extend)
        #save += '_others'+str(others)
        save += '_fit-results'
        if note: save += '_'+str(note)
        save += '_'+V

        ### make sure the plotdir still exists
        if not os.path.exists(plotdir): 
            os.makedirs(plotdir)

        
        ### only write out files if fit is successful
        #-------------------------------------------------------------
        if status == 0:
            
            ### sort the fitresults keys
            resultskeys = []
            for rskey in fitresults:
                resultskeys.append(rskey)
            resultskeys.sort()
                    
            ### write results to file
            resultsfile = os.path.join(plotdir, save+'.txt')
            outfile = open(resultsfile, 'w')
            for key in resultskeys:
                if int(key)+1 in justthese:
                    for line in fitresults[key]:
                        outfile.write(line+'\n')
            outfile.close()

            ### save scaled histograms to rootfile?
            if SAVE_SCALED_HISTS:
                rootoutfile = TFile(plotdir+"/histograms-scaled.root", "RECREATE")
                for key in sigkeys:
                    sigs[key]['hist'].Write(key)
                for key in bakkeys:
                    bkgs[key]['hist'].Write(key)
                for key in data:
                    data[key]['hist'].Write(key)
                rootoutfile.Write()
                rootoutfile.Close()
        #-------------------------------------------------------------
        
        ### create the updated backgrounds file
        print()
        #updateBkgsFile300(justthese, os.path.join(plotdir, mcfile), resultsfile, plotdir)
        updateBkgsFile500(justthese, os.path.join(plotdir, mcfile), resultsfile, plotdir)
        
        
        #=============================================================
        #=============================================================
        #
        #   plot the full fit histogram
        #
        #=============================================================
        #=============================================================
        
        mcanv = TCanvas('mcanv', 'mcanv', 0, 0, 1400, 900)
        #mcanv.Divide(4,2)

        gStyle.SetPadTopMargin    (0.07)
        gStyle.SetPadBottomMargin (0.10)
        gStyle.SetPadLeftMargin   (0.08)
        gStyle.SetPadRightMargin  (0.03)

        font = 63
        size = 20
        yoff = 1.5
        xoff = 6

        fraction = 0.3
        mtoppad = TPad('pad1','pad1',0,fraction,1,1)
        mbotpad = TPad('pad2','pad2',0,0,1,fraction)
        mtoppad.SetBottomMargin(0.01)
        mtoppad.SetBorderMode(0)
        mtoppad.SetLogy(1)
        
        mbotpad.SetTopMargin(0.05)
        mbotpad.SetBottomMargin(0.3)
        mbotpad.SetBorderMode(0)
        mbotpad.SetLogy(1)

        mtoppad.Draw()
        mbotpad.Draw()
        mtoppad.cd()
        
        newFitTitle = str('Global Fit Histograms')
        fitdata.SetTitle(newFitTitle)
        
        fitdata.SetLineColor(kBlack)
        fitdata.SetMarkerColor(kBlack)
        fitdata.SetLineWidth(1)

        if dru: fitdata.GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
        else: fitdata.GetYaxis().SetTitle('arb. counts')
        fitdata.GetYaxis().SetTitleFont(font)
        fitdata.GetYaxis().SetTitleSize(size)
        fitdata.GetYaxis().SetTitleOffset(yoff)
        fitdata.GetYaxis().SetLabelFont(font)
        fitdata.GetYaxis().SetLabelSize(size)
        fitdata.GetYaxis().SetLabelOffset(0.01)
        
        fitdata.Draw()
        #flegs[i].AddEntry(fitdata, 'data - bkgs', legopt)
        """
        if dru:
            fitdata.SetAxisRange(2e-3, 2e3, 'y')
            #if i == 8:
            #    fitdata.SetAxisRange(2e-5, 2e1, 'y')
        """
        fitdata.SetAxisRange(fitYrange[0], fitYrange[1], 'y')
        fitdata.SetAxisRange(0, fitbins, 'x')
        
        
        Nfsigs = 0
        for fskey in fsigkeys:

            Nfsigs += 1

            # find the unique name for color and set color
            cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
            fitsigs[fskey]['hist'].SetMarkerColor(siguniqColor[cname])
            fitsigs[fskey]['hist'].SetLineColor(siguniqColor[cname])

            ### draw the sigs
            if showFitSigs: fitsigs[fskey]['hist'].Draw(dopt)

            ### add MC to total MC hist
            ftotal.Add(fitsigs[fskey]['hist'])

        for fgkey in fglobsigkeys:
            
            Nfsigs += 1

            # find the unique name for color and set color
            #cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
            cname = fgkey
            fitglobsigs[fgkey]['hist'].SetMarkerColor(siguniqColor[cname])
            fitglobsigs[fgkey]['hist'].SetLineColor(siguniqColor[cname])

            ### draw the sigs
            if showFitSigs: fitglobsigs[fgkey]['hist'].Draw(dopt)

            ### add MC to total MC hist
            ftotal.Add(fitglobsigs[fgkey]['hist'])

                
        ### select the right range
        ftotal.SetAxisRange(0, fitbins, 'x')
        if redtotal:
            ftotal.SetLineColor(kRed)
            ftotal.SetMarkerColor(kRed)
        ftotal.Draw(dopt)
        
        ### build legend after getting the numbers of signals
        Nfs = Nfsigs+2
        flnc = 1
        if Nfs > 10: flnc = 2
        if Nfs > 20: flnc = 3
        if Nfs > 30: flnc = 4
        if Nfs > 40: flnc = 5
        #fxlegstop  = 0.94
        fxlegstop  = 0.96
        #fxlegstart = fxlegstop-(0.2*flnc)
        fxlegstart = fxlegstop-(0.14*flnc)
        #fylegstop  = 0.89
        fylegstop  = 0.91
        fylegstart = fylegstop-(0.04*6)

        mleg = TLegend(fxlegstart, fylegstart, fxlegstop, fylegstop)
        #flegs.append(fleg)
        mleg.SetNColumns(flnc)
        mleg.SetFillColor(0)
        mleg.SetBorderSize(0)
        legopt = 'LPE'

        mleg.AddEntry(fitdata, 'data - bkgs', legopt)

        ### add legend entries in order
        #for name in uniqAll:
        for name in uniqSigs:
            for fskey in fsigkeys:
                # need to add a "-" here for pmt in xpmt overlap
                # fixed 2022-06-29
                if '-'+name in fskey:
                    #print('legend signal', name, fskey)
                    mleg.AddEntry(fitsigs[fskey]['hist'], fskey, legopt)
            # and for globals
            for fgkey in fglobsigkeys:
                #print(name, fgkey)
                ### make this unique (pmt vs vetopmt) 2022-09-25
                #if name in fgkey:
                if name == fgkey:
                    #print('legend global', name, fgkey)
                    mleg.AddEntry(fitglobsigs[fgkey]['hist'], fgkey, legopt)

        #ftotal.Draw('same')

        # returned from fit
        mleg.AddEntry(ftotal, 'Fit Total (chi2/ndf = '+str(round(fitchi2ndf,2))+')', legopt)
        # calc by Chi2TestX()
        #flegs[i].AddEntry(ftotal, 'Fit Total (chi2/ndf = '+str(round(fitchi2ndfv2,2))+')', legopt)
        #print 'INFO: Fit total MC from Chi2TestX chi2/ndf = '+str(round(fitchi2ndfv2,2))

        if showFitSigs and showFitLeg:
            mleg.Draw('same')


        ### plot the fit residuals
        #-------------------------------------------------------------
        mbotpad.cd()
        mleg2 = TLegend(0.72, 0.78, 0.94, 0.94)
        #flegs2.append(leg)
        mleg2.SetFillColor(0)
        mleg2.SetBorderSize(0)
        legopt = 'LPE'
        if debug: print('DEBUG: hist divide #1')
        try:
            fresid.Divide(fitdata, ftotal)
        except:
            print('WARNING: could not divide fit data by fit total')
        
        fresid.SetTitle('')
        #fresid.SetXTitle('Energy (keVee)')
        fresid.SetXTitle('fit bins')
        fresid.GetXaxis().SetTitleFont(font)
        fresid.GetXaxis().SetTitleSize(size)
        fresid.GetXaxis().SetTitleOffset(xoff)
        fresid.GetXaxis().SetLabelFont(font)
        fresid.GetXaxis().SetLabelSize(size)
        fresid.GetXaxis().SetLabelOffset(0.03)
        
        #fresid.SetYTitle('counts / keV')
        fresid.SetYTitle('data / MC')
        #fresid.SetYTitle('MC / data')
        fresid.GetYaxis().SetTitleFont(font)
        fresid.GetYaxis().SetTitleSize(size)
        fresid.GetYaxis().SetTitleOffset(yoff)
        fresid.GetYaxis().SetLabelFont(font)
        fresid.GetYaxis().SetLabelSize(size)
        fresid.GetYaxis().SetLabelOffset(0.01)
        # '5' secondary and '05' primary
        fresid.GetYaxis().SetNdivisions(505)
        
        fresid.SetAxisRange(0.1, 10, 'y')
        
        fresid.SetAxisRange(0, fitbins, 'x')
        
        if linres:
            mbotpad.SetLogy(0)
            fresid.SetAxisRange(lrs[0], lrs[1], 'y')
        fresid.Draw("P")

        ### set my reference line to lrm variable
        mzero = TLine(0, lrm, fitbins, lrm)
        mzero.SetLineColor(kRed)
        mzero.SetLineWidth(1)
        mzero.Draw("same")
        
        mleg2.AddEntry(fresid,'data / MC',legopt)
        #flegs2[i].Draw()
        #-------------------------------------------------------------
        
        mcanv.Update()

        cfsave  = ''
        cfsave += str(runtag)
        #cfsave += '_all-xstals'
        cfsave += '_post-fit'
        #cfsave += '_globals-'+globstr
        #cfsave += '_loEfit-'+str(int(fLoE[0]))+'-'+str(int(fLoE[1]))
        #cfsave += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        #cfsave += '_loEfitRebin-'+str(loEfitRebin)
        #cfsave += '_hiEfitRebin-'+str(hiEfitRebin)
        #cfsave += '_loEfitRebinScale'+str(loEfitRebinScale)
        #cfsave += '_hiEfitRebinScale'+str(hiEfitRebinScale)
        #cfsave += '_useBounds'+str(useBounds)
        #cfsave += '_mcsumw2'+str(mcsumw2)
        #cfsave += '_datsumw2'+str(datsumw2)
        #cfsave += '_dru'+str(dru)
        #cfsave += '_loEplotRebin-'+str(loEplotRebin)
        #cfsave += '_hiEplotRebin-'+str(hiEplotRebin)
        #cfsave += '_reuse'+str(reuse)
        #cfsave += '_chans'+str(fitchans)
        #cfsave += '_extend'+str(extend)
        #cfsave += '_others'+str(others)
        if note: cfsave += '_'+str(note)
        cfsave += '_'+V
        
        mcanv.Print(plotdir+'/'+cfsave+'.png')
        
        #=============================================================
        #=============================================================

        #raw_input('[Enter] to quit \n')
        
    ### end of fitting bit if you have signals

    
    
    #=================================================================
    #  plot histograms for all crystals + lsveto
    #=================================================================
    
    # number of energy ranges (low, high, alpha)
    numE = 3
    
    # number of channels (single-hit, multi-hit)
    numC = 2
    
    canvs = [[[0 for i in range(numx)] for e in range(numE)] for c in range(numC)]
        
    ### seperate memory space for the pads is key!!!!
    toppad = [[[0 for i in range(numx)] for e in range(numE)] for c in range(numC)]
    botpad = [[[0 for i in range(numx)] for e in range(numE)] for c in range(numC)]
    
    legs   = [[[0 for i in range(numx)] for e in range(numE)] for c in range(numC)]
    legs2  = [[[0 for i in range(numx)] for e in range(numE)] for c in range(numC)]
    zeros  = [[[0 for i in range(numx)] for e in range(numE)] for c in range(numC)]
    
    total  = [[0 for e in range(numE)] for c in range(numC)]
    resid  = [[0 for e in range(numE)] for c in range(numC)]
    
    gbkgs  = [[[{} for x in range(numx)] for x in range(numE)] for x in range(numC)]
    gsigs  = [[[{} for x in range(numx)] for x in range(numE)] for x in range(numC)]
    
    gStyle.SetPadTopMargin    (0.07)
    gStyle.SetPadBottomMargin (0.11)
    gStyle.SetPadRightMargin  (0.02)
    gStyle.SetPadLeftMargin   (0.10)
    
    font = 43
    size = 18
    yoff = 2
    xoff = 11
    if showResid:
        fraction = 0.25
        botmarg = 0.01
    else:
        fraction = 0
        botmarg = 0.1
        
    if not show_indi: gROOT.SetBatch(1)
    
    for C, chan in enumerate(pltchans): 
        for E in range(numE):
            
            # inits 9 hists for given chan and energy
            total[C][E] = makeTotal100(chan, E, params[E])
            resid[C][E] = makeResid100(chan, E, params[E])
            
            #for i in range(numx):
            for i in [j-1 for j in justthese]:
                canvs[C][E][i] = TCanvas('canv'+str(i+1)+chan+str(E),
                                         'canv'+str(i+1)+chan+str(E),
                                         0, 0, 1200, 700)
                
                # get plot settings
                if i < 8:
                    pkey = chan+str(E)
                    popts = xstal_settings[pkey]
                else:
                    pkey = chan+str(E)
                    popts = ls_settings[pkey]
                
                pad1 = TPad('pad1'+chan+str(E),'pad1'+chan+str(E),0,fraction,1,1)
                toppad[C][E][i] = (pad1)
                toppad[C][E][i].SetBottomMargin(botmarg)
                toppad[C][E][i].SetBorderMode(0)
                #toppad[C][E][i].SetLeftMargin(0.08)
                toppad[C][E][i].SetLogy(popts['log'])
                toppad[C][E][i].Draw()
                
                if showResid:
                    pad2 = TPad('pad2'+chan+str(E),'pad2'+chan+str(E),0,0,1,fraction)
                    botpad[C][E][i] = (pad2)
                    botpad[C][E][i].SetTopMargin(0.05)
                    botpad[C][E][i].SetBottomMargin(0.35)
                    botpad[C][E][i].SetBorderMode(0)
                    #botpad[C][E][i].SetLeftMargin(0.08)
                    botpad[C][E][i].SetLogy(1)
                    botpad[C][E][i].Draw()
                    

                toppad[C][E][i].cd()

                if ingroups:
                    #leg = TLegend(0.35, 0.75, 0.94, 0.89)
                    leg = TLegend(0.56, 0.76, 0.98, 0.92)
                    lnc = 3
                else:
                    leg = TLegend(xlegstart, ylegstart, xlegstop, ylegstop)
                    
                legs[C][E][i] = (leg)
                legs[C][E][i].SetFillColor(0)
                legs[C][E][i].SetBorderSize(0)
                legs[C][E][i].SetNColumns(lnc)
                legopt = 'LPE'

                total[C][E][i].Rebin(popts['rebin'])
                #total[C][E][i].SetLineWidth(4)
                #total[C][E][i].SetMarkerSize(4)
                if redtotal:
                    total[C][E][i].SetMarkerColor(kRed)
                    total[C][E][i].SetLineColor(kRed)

                
                dkey = 0
                datasets = 0
                for key in datkeys:
                    if 'x'+str(i+1) in key and '-c'+chan in key and '-e'+str(E) in key:

                        dkey = key
                        datasets += 1
                        
                        if dru:
                            #data[dkey]['hist'].GetYaxis().SetTitle('counts / day / kg / keV')
                            data[dkey]['hist'].GetYaxis().SetTitle('counts (dru)')
                        else:
                            data[dkey]['hist'].GetYaxis().SetTitle('arb. counts')

                        newTitle = (cnames(i)+'   '
                                    +chanNames(chan)+'   '
                                    +energyNames(E))
                        
                        data[dkey]['hist'].SetTitle(newTitle)
                        data[dkey]['hist'].SetTitleFont(font)
                        data[dkey]['hist'].SetTitleSize(size)
                        
                        data[dkey]['hist'].Rebin(popts['rebin'])
                        data[dkey]['hist'].Scale(1./float(popts['rebin']))

                        data[dkey]['hist'].SetMarkerColor(dcs[datasets-1])
                        data[dkey]['hist'].SetLineColor(dcs[datasets-1])
                        data[dkey]['hist'].SetLineWidth(1)
                        
                        data[dkey]['hist'].GetXaxis().SetTitle('Energy (keV)')
                        
                        data[dkey]['hist'].GetXaxis().SetTitleFont(font)
                        data[dkey]['hist'].GetXaxis().SetTitleSize(size)
                        data[dkey]['hist'].GetXaxis().SetTitleOffset(yoff)
                        data[dkey]['hist'].GetXaxis().SetLabelFont(font)
                        data[dkey]['hist'].GetXaxis().SetLabelSize(size)
                        
                        data[dkey]['hist'].GetYaxis().SetTitleFont(font)
                        data[dkey]['hist'].GetYaxis().SetTitleSize(size)
                        data[dkey]['hist'].GetYaxis().SetTitleOffset(yoff)
                        data[dkey]['hist'].GetYaxis().SetLabelFont(font)
                        data[dkey]['hist'].GetYaxis().SetLabelSize(size)
                        
                        if popts['xrange']:
                            data[dkey]['hist'].SetAxisRange(popts['xrange'][0], popts['xrange'][1], 'x')
                        if dru:
                            if popts['yrange']:
                                log = int(popts['log'])
                                yrange = popts['yrange'][log]
                                data[dkey]['hist'].SetAxisRange(yrange[0], yrange[1], 'y')
                        
                        data[dkey]['hist'].Draw('same')
                        days = round(data[dkey]['runtime']/86400., 2)
                        dataLegName = data[dkey]['info']['tag']\
                            +'  '+data[dkey]['info']['build']
                        if dataname:
                            dataLegName = dataname
                        legs[C][E][i].AddEntry(data[dkey]['hist'], dataLegName, legopt)
                        
                        
                tcount = 0
                bkey = 0
                for key in bakkeys:
                    if 'x'+str(i+1) in key and '-c'+chan in key and '-e'+str(E) in key:

                        bkey = key
                        
                        # find the unique name for color and set color
                        cname = key.split('-')[1]+'-'+key.split('-')[2]
                        bkgs[key]['hist'].SetMarkerColor(uniqColor[cname])
                        bkgs[key]['hist'].SetLineColor(uniqColor[cname])
                        
                        #if E and plotRebin:
                        bkgs[key]['hist'].Rebin(popts['rebin'])
                        bkgs[key]['hist'].Scale(1./float(popts['rebin']))
                        #bkgs[key]['hist'].Sumw2()
                        
                        if ingroups:
                            group = bkgs[key]['info']['group']
                            if group == 'none':
                                bits = key.split('-')
                                nkey = bits[1]+'-'+bits[2]
                                try:
                                    gbkgs[C][E][i]['none'][nkey] = deepcopy(bkgs[key]['hist'])
                                except:
                                    gbkgs[C][E][i]['none'] = {}
                                    gbkgs[C][E][i]['none'][nkey] = deepcopy(bkgs[key]['hist'])
                            else:
                                try:
                                    gbkgs[C][E][i][group].Add(bkgs[key]['hist'])
                                except:
                                    gbkgs[C][E][i][group] = deepcopy(bkgs[key]['hist'])
                                    gbkgs[C][E][i][group].SetNameTitle(group, group)
                        else:
                            bkgs[key]['hist'].Draw(dopt)

                        ### add MC to total MC hist
                        if debug: print('DEBUG: adding to total hist --> {0}'.format(key))
                        total[C][E][i].Add(bkgs[key]['hist'])
                        tcount += 1
                        
                skey = 0
                for key in sigkeys:
                    if 'x'+str(i+1) in key and '-c'+chan in key and '-e'+str(E) in key:

                        skey = key
                        
                        # find the unique name for color and set color
                        cname = key.split('-')[1]+'-'+key.split('-')[2]
                        sigs[key]['hist'].SetMarkerColor(uniqColor[cname])
                        sigs[key]['hist'].SetLineColor(uniqColor[cname])

                        sigs[key]['hist'].Rebin(popts['rebin'])
                        sigs[key]['hist'].Scale(1./float(popts['rebin']))
                        #sigs[key]['hist'].Sumw2()

                        if ingroups:
                            group = sigs[key]['info']['group']
                            if group == 'none':
                                bits = key.split('-')
                                nkey = bits[1]+'-'+bits[2]
                                try:
                                    gbkgs[C][E][i]['none'][nkey] = deepcopy(sigs[key]['hist'])
                                except:
                                    gbkgs[C][E][i]['none'] = {}
                                    gbkgs[C][E][i]['none'][nkey] = deepcopy(sigs[key]['hist'])
                            else:
                                try:
                                    gbkgs[C][E][i][group].Add(sigs[key]['hist'])
                                except:
                                    gbkgs[C][E][i][group] = deepcopy(sigs[key]['hist'])
                                    gbkgs[C][E][i][group].SetNameTitle(group, group)
                        else:
                            sigs[key]['hist'].Draw(dopt)

                        ### add MC to total MC hist
                        total[C][E][i].Add(sigs[key]['hist'])
                        tcount += 1

                if ingroups:
                    groupSort = []
                    for group in gbkgs[C][E][i]:
                        if group not in groupSort:
                            groupSort.append(group)
                    groupSort.sort()
                    for group in groupSort:
                        if group != 'none':
                            gbkgs[C][E][i][group].SetMarkerColor(gis[group])
                            gbkgs[C][E][i][group].SetLineColor(gis[group])
                            gbkgs[C][E][i][group].Draw(dopt)
                            legs[C][E][i].AddEntry(gbkgs[C][E][i][group], group, legopt)
                    if 'none' in groupSort:
                        for key in gbkgs[C][E][i]['none']:
                            gbkgs[C][E][i]['none'][key].SetMarkerColor(gis['none'])
                            gbkgs[C][E][i]['none'][key].SetLineColor(gis['none'])
                            gbkgs[C][E][i]['none'][key].Draw(dopt)
                            if showNoneGroup:
                                legs[C][E][i].AddEntry(gbkgs[C][E][i]['none'][key], key, legopt)
                else:
                    # add legend entries in order
                    for name in uniqAll:
                        for bkey in bakkeys:
                            #print bkey
                            if bkey == 'x'+str(i+1)+'-'+name+'-c'+chan+'-e'+str(E) \
                               or bkey == 'x'+str(i+1)+'-'+name+'-f'+str(i+1)+'-c'+chan+'-e'+str(E):
                                activ = '(%.2e) '%(bkgs[bkey]['info']['acti'])
                                # make legend names smaller
                                bits = bkey.split('-')
                                legname = bits[1]+' '+bits[2]
                                legs[C][E][i].AddEntry(bkgs[bkey]['hist'], legname, legopt)
                        for skey in sigkeys:
                            if skey == 'x'+str(i+1)+'-'+name+'-c'+chan+'-e'+str(E) \
                               or skey == 'x'+str(i+1)+'-'+name+'-f'+str(i+1)+'-c'+chan+'-e'+str(E):
                                activ = '(%.2e) '%(sigs[skey]['info']['fitacti'])
                                # make legend names smaller
                                bits = skey.split('-')
                                legname = bits[1]+' '+bits[2]
                                legs[C][E][i].AddEntry(sigs[skey]['hist'], legname, legopt)
                
                ### you need to scale the error by the dru scaling and/or the rebinning
                #-----------------------------------------------------------------------------
                # don't set the error on total until I understand what's going on?
                if toterr:
                    if dru:
                        for n in range(total[C][E][i].GetNbinsX()):
                            total[C][E][i].SetBinError(n+1,
                                total[C][E][i].GetBinError(n+1)*data[dkey]['druScale'])
                    else:
                        for n in range(total[C][E][i].GetNbinsX()):
                            total[C][E][i].SetBinError(n+1,
                                total[C][E][i].GetBinError(n+1)/(float(popts['rebin'])/math.sqrt(2.)))
                
                
                ### chi2 test
                ### get the chi2 of the total mc compared to data
                #=============================================================================
                #-----------------------------------------------------------------------------
                if printChi2:
                    #chi2  = ROOT.Double(-1)
                    #ndf   = ROOT.Long(1)
                    #igood = ROOT.Long(0)
                    # ROOT6 FIX
                    chi2  = ctypes.c_double(0.0)
                    ndf   = ctypes.c_int(0)
                    igood = ctypes.c_int(0)
                    if dkey and tcount:
                        if data[dkey]['hist'].GetEntries() > 0 \
                           and total[C][E][i].GetEntries() > 0:
                            pval  = data[dkey]['hist'].Chi2TestX(
                                total[C][E][i], chi2, ndf, igood, chiopt)
                            # ROOT6 FIX
                            chi2  = float(chi2.value)
                            ndf   = int(ndf.value)
                            igood = int(igood.value)
                            print('INFO: {0} total MC chi2/ndf = {1}'
                                  .format(dkey, round(chi2/ndf, 2)))

                # for Gyunho comparison
                if E == 0 and gyunhoChi2:
                    start_bin = int(round((1*12)/popts['rebin'],0))
                    stop_bin = int(round((100*12)/popts['rebin'],0))
                    #print('start/stop bins = {0}-{1}'.format(start_bin, stop_bin))
                    chi2 = 0
                    ndf = 0
                    for n in range(start_bin, stop_bin+1):
                        chi2 += ( (data[dkey]['hist'].GetBinContent(n)
                                   - total[C][E][i].GetBinContent(n))**2
                                  / total[C][E][i].GetBinContent(n) )
                        ndf += 1
                    print('{0} chi^2/ndf = {1}/{2} = {3}'
                          .format(dkey, round(chi2, 1), ndf, round(chi2/ndf, 1)))
                #-----------------------------------------------------------------------------
                #=============================================================================



                #-----------------------------------------------------------------------------
                # still draw total even if nothing so hist box is created
                if showTotal:
                    total[C][E][i].Draw(dopt)
                    if tcount:
                        total[C][E][i].SetLineWidth(1)
                        totLegName = "Total"
                        if totname: totLegName = totname
                        legs[C][E][i].AddEntry(total[C][E][i], totLegName, legopt)
                        """
                        if ingroups:
                            legs[C][E][i].AddEntry(total[C][E][i], 'Total', legopt)
                        else:
                            legs[C][E][i].AddEntry(total[C][E][i], 'Total', legopt)
                            #legs[C][E][i].AddEntry(total[C][E][i],
                            #    'Total MC (chi2/ndf = '+str(round(chi2/ndf,2))+')', legopt)
                        """
                else:
                    # draw an empty total hist to preserve plot box layout?
                    if not tcount and not dkey: total[C][E][i].Draw()
                
                ### show the legends?
                if showLeg and (dkey or bkey):
                    legs[C][E][i].Draw('same')

                ### update canvas
                toppad[C][E][i].Update()
                canvs[C][E][i].Update()

                
                ### try to get the residuals in!
                #---------------------------------------------------------
                if showResid:
                    resid[C][E][i].Rebin(popts['rebin'])
                    
                    botpad[C][E][i].cd()
                    leg2 = TLegend(0.72, 0.78, 0.94, 0.94)
                    legs2[C][E][i] = (leg2)
                    legs2[C][E][i].SetFillColor(0)
                    legs2[C][E][i].SetBorderSize(0)
                    legopt = 'LPE'

                    if tcount and dkey:
                        #print('here?')
                        total[C][E][i] = cleanSmallBins(total[C][E][i])
                        #print('done?')
                        if debug: print('DEBUG: hist divide #2')
                        try:
                            resid[C][E][i].Divide(data[dkey]['hist'], total[C][E][i])
                        except:
                            print('WARNING: could not divide', dkey, 'by total')
                            pass
                    
                    resid[C][E][i].SetTitle('')
                    resid[C][E][i].SetXTitle('Energy (keV)')
                    resid[C][E][i].SetYTitle('data / MC')
                    
                    resid[C][E][i].GetXaxis().SetTitleFont(font)
                    resid[C][E][i].GetXaxis().SetTitleSize(size)
                    resid[C][E][i].GetXaxis().SetTitleOffset(xoff-4)
                    resid[C][E][i].GetXaxis().SetLabelFont(font)
                    resid[C][E][i].GetXaxis().SetLabelSize(size)

                    resid[C][E][i].GetYaxis().SetTitleFont(font)
                    resid[C][E][i].GetYaxis().SetTitleSize(size)
                    resid[C][E][i].GetYaxis().SetTitleOffset(yoff)
                    resid[C][E][i].GetYaxis().SetLabelFont(font)
                    resid[C][E][i].GetYaxis().SetLabelSize(size)
                    
                    # '5' secondary and '05' primary
                    resid[C][E][i].GetYaxis().SetNdivisions(503)

                    if popts['xrange']:
                        resid[C][E][i].SetAxisRange(popts['xrange'][0], popts['xrange'][1], 'x')

                    if linres:
                        botpad[C][E][i].SetLogy(0)
                        resid[C][E][i].SetAxisRange(lrs[0], lrs[1], 'y')
                    else:
                        resid[C][E][i].SetAxisRange(0.1, 10, 'y')

                    ###==================================================================
                    ### set errors on resid R where R=Data/Total
                    ### sigR = R*sqrt((sigD/D**2) + (sigT/T)**2)
                    ### but sigT=0 (for now) so...
                    ### sigR = R*(sigD/D)
                    if dkey:
                        for n in range(resid[C][E][i].GetNbinsX()+1):
                            R  = resid[C][E][i].GetBinContent(n)
                            D  = data[dkey]['hist'].GetBinContent(n)
                            sD = data[dkey]['hist'].GetBinError(n)

                            try:
                                resid[C][E][i].SetBinError(n, R * (sD/D))
                            except:
                                resid[C][E][i].SetBinError(n, 0)
                            """
                            # assume sqrt(N) error on the total
                            T  = total[C][E][i].GetBinContent(n)
                            sT = np.sqrt(total[C][E][i].GetBinContent(n))
                            # assume a correlation of 1
                            c  = 1.
                            cDT = c*sD*sT
                            # this over estimates the error a lot
                            try:
                                resid[C][E][i].SetBinError(n, R * \
                                    np.sqrt(((sD/D)**2) + ((sT/T)**2) - (2*(cDT/(D*T)))))
                            except:
                                resid[C][E][i].SetBinError(n, 0)
                            """
                    else:
                        resid[C][E][i] = zeroBinError(resid[C][E][i])
                    ###==================================================================

                    resid[C][E][i].Draw()

                    # set my ref line to lrm variable
                    if popts['xrange']:
                        zero = TLine(popts['xrange'][0], lrm, popts['xrange'][1], lrm)
                    else:
                        zpars = getPars(resid[C][E][i])
                        zero = TLine(zpars[1], lrm, zpars[2], lrm)

                    zeros[C][E][i] = (zero)
                    zeros[C][E][i].SetLineColor(kRed)
                    zeros[C][E][i].SetLineWidth(1)
                    zeros[C][E][i].Draw()

                    # update canvas
                    botpad[C][E][i].Update()
                    canvs[C][E][i].Update()
    
    # saving plots
    for j in justthese:
        i = j-1
        for C, chan in enumerate(['S', 'M']):
            for E in [0, 1]: # skip alphas
                if save_indi or (i==8 and E==1):
                    if (i==8 and E!=1): # skip ls low-energy
                        continue
                    else:
                        isave = 'x{0}-c{1}-e{2}_{3}.png'.format(i+1, chan, E, V)
                        canvs[C][E][i].Print(os.path.join(plotdir, isave))

    # save the total hists to rootfile
    if SAVE_TOTAL_HISTS:
        rootoutfile = TFile(plotdir+"/histograms-totals.root", "RECREATE")
        for j in justthese:
            i = j-1
            for C in [0,1]:
                for E in [0,1]:
                    total[C][E][i].Write()
        rootoutfile.Write()
        rootoutfile.Close()

    
    gROOT.SetBatch(batch)
    # done with individual plots
    #============================================================================
    
    ### Save 4 (or 6) plots to one canvas
    #-------------------------------------------------------------
    if True:  # always show these plots
        combPlots = [0 for x in range(numx)]
        tpad = [[0 for x in range(numx)] for x in range(6)]
        bpad = [[0 for x in range(numx)] for x in range(6)]
        
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
            combPlots[i] = TCanvas('ccan-'+str(i+1), 'ccan-'+str(i+1),
                                   0, 0, 1400, 900)
            
            combPlots[i].Divide(2,2)
            Es = 2
            
            p=0
            for C, chan in enumerate(pltchans):
                for E in range(Es):
                    combPlots[i].cd(p+1)
                    tpad[p][i] = toppad[C][E][i].Clone()
                    tpad[p][i].Draw()
                    if showResid:
                        bpad[p][i] = botpad[C][E][i].Clone()
                        bpad[p][i].Draw()
                    combPlots[i].Update()
                    p += 1
            csave  = ''
            csave += 'all-channels'
            csave += '_x'+str(i+1)
            if note: csave += '_'+str(note)
            csave += '_'+str(V)

            combPlots[i].Print(plotdir+'/'+csave+'.png')

    
    ### special plots variables
    #---------------------------------------
    font = 43
    size = 24
    loff = 0.005
    xoff = 3.6
    yoff = 1.1
    if showResid:
        fraction = 0.25
        botmarg = 0.01
    else:
        fraction = 0
        botmarg = 0.1
    #---------------------------------------

    
    ### Plot the 0-40 keV region of single hit
    #------------------------------------------
    #if indi:
    if save_zoom:
        if not show_zoom: gROOT.SetBatch(1)
        
        zsh_canv = [0 for x in range(numx)]
        zsh_hist = [0 for x in range(numx)]
        zsh_data = [0 for x in range(numx)]
        zsh_total = makeTotal100('S', 0, params[0])
        zsh_resid = makeResid100('S', 0, params[0])
        zsh_line = [0 for x in range(numx)]
        zsh_top = [0 for x in range(numx)]
        zsh_bot = [0 for x in range(numx)]
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:

                zsh_canv[i] = TCanvas('zscanv-'+str(i+1),'zscanv-'+str(i+1),
                                0, 0, 1200, 700)

                zspad1 = TPad('zspad1_'+str(i),'zspad1_'+str(i),0,fraction,1,1)
                zsh_top[i] = (zspad1)
                zsh_top[i].SetBottomMargin(botmarg)
                zsh_top[i].SetBorderMode(0)
                #zsh_top[i].SetLogy(1)
                zsh_top[i].SetLogy(0)
                zsh_top[i].Draw()
                
                if showResid:
                    zspad2 = TPad('zspad2_'+str(i),'zspad2_'+str(i),0,0,1,fraction)
                    zsh_bot[i] = (zspad2)
                    zsh_bot[i].SetTopMargin(0.05)
                    zsh_bot[i].SetBottomMargin(0.3)
                    zsh_bot[i].SetBorderMode(0)
                    #zsh_bot[i].SetLogy(1)
                    zsh_bot[i].SetLogy(0)
                    zsh_bot[i].Draw()
                    

                zsh_top[i].cd()

                # rebin total here?
                zsh_total[i].Rebin(xstal_settings['S0']['rebin'])
                
                prim_list = toppad[0][0][i].GetListOfPrimitives()
                for bla in prim_list:
                    #print(bla)

                    if 'data' in bla.GetName():
                        zsh_data[i] = (toppad[0][0][i].Clone()).GetPrimitive(bla.GetName())
                        
                    if 'total' in bla.GetName():
                        #zsh_total[i] = (toppad[0][0][i].Clone()).GetPrimitive(bla.GetName())
                        continue
                    
                    try:
                        zsh_hist[i] = (toppad[0][0][i].Clone()).GetPrimitive(bla.GetName())
                        
                        zsh_hist[i].SetXTitle('energy (keV)')
                        zsh_hist[i].GetXaxis().SetTitleFont(font)
                        zsh_hist[i].GetXaxis().SetTitleSize(size)
                        #zsh_hist[i].GetXaxis().SetTitleOffset(xoff)
                        zsh_hist[i].GetXaxis().SetTitleOffset(1.1)
                        zsh_hist[i].GetXaxis().SetLabelFont(font)
                        zsh_hist[i].GetXaxis().SetLabelSize(size)
                        zsh_hist[i].GetXaxis().SetLabelOffset(loff)
                        
                        zsh_hist[i].SetYTitle('counts (dru)')
                        zsh_hist[i].GetYaxis().SetTitleFont(font)
                        zsh_hist[i].GetYaxis().SetTitleSize(size)
                        zsh_hist[i].GetYaxis().SetTitleOffset(yoff)
                        zsh_hist[i].GetYaxis().SetLabelFont(font)
                        zsh_hist[i].GetYaxis().SetLabelSize(size)
                        zsh_hist[i].GetYaxis().SetLabelOffset(loff)
                        
                        zsh_hist[i].SetAxisRange(0, zmaxE, 'x')

                        if   i+1 == 1: zsh_hist[i].SetAxisRange(0.01,  4.0, 'y')
                        elif i+1 == 2: zsh_hist[i].SetAxisRange(0.01,  4.0, 'y')
                        elif i+1 == 3: zsh_hist[i].SetAxisRange(0.01,  4.0, 'y')
                        elif i+1 == 4: zsh_hist[i].SetAxisRange(0.01,  4.0, 'y')
                        elif i+1 == 5: zsh_hist[i].SetAxisRange(0.01, 10.0, 'y')
                        elif i+1 == 6: zsh_hist[i].SetAxisRange(0.01,  3.5, 'y')
                        elif i+1 == 7: zsh_hist[i].SetAxisRange(0.01,  3.5, 'y')
                        elif i+1 == 8: zsh_hist[i].SetAxisRange(0.01, 10.0, 'y')
                        else:          zsh_hist[i].SetAxisRange(0.01, 10.0, 'y')
                                                
                        if 'data' in bla.GetName():
                            zsh_hist[i].Draw('same')
                        else:
                            zsh_hist[i].Draw(dopt)
                            zsh_total[i].Add(zsh_hist[i])
                            
                    except:
                        #print('failed on', bla)
                        continue
                    
                if showTotal:
                    if redtotal:
                        zsh_total[i].SetMarkerColor(kRed)
                        zsh_total[i].SetLineColor(kRed)
                    zsh_total[i].Draw(dopt)
                    #legs[0][0][i].AddEntry(zsh_total[i], "Total", legopt)
                    
                if showLeg:
                    try:
                        #leg = TLegend(toppad[0][0][i].GetPrimitive(bla.GetName()))
                        #leg.Draw("same")
                        legs[0][0][i].Draw('same')
                    except: continue

                # update canvas
                zsh_top[i].Update()
                
                if showResid: 
                    zsh_bot[i].cd()

                    zsh_resid[i].Rebin(xstal_settings['S0']['rebin'])

                    if zdr:
                        if debug: print('DEBUG: hist divide #3')
                        try:
                            zsh_resid[i].Divide(zsh_data[i], zsh_total[i])
                        except:
                            print('WARNING: could not divide zoomed low-energy single-hit')
                            print('         data nbins = {0} : total nbins = {1}'
                                  .format(zsh_data[i].GetXaxis().GetNbins(),
                                          zsh_total[i].GetXaxis().GetNbins()))
                            pass
                    else:
                        try:
                            zsh_resid[i] = zsh_data[i] - zsh_total[i]
                        except:
                            print('WARNING: could not subtract zoomed low-energy single-hit')
                            print('         data nbins = {0} : total nbins = {1}'
                                  .format(zsh_data[i].GetXaxis().GetNbins(),
                                          zsh_total[i].GetXaxis().GetNbins()))
                            pass

                    ### tweak the errors on the residual
                    for n in range(zsh_resid[i].GetNbinsX()+1):
                        R  = zsh_resid[i].GetBinContent(n)
                        D  = zsh_data[i].GetBinContent(n)
                        sD = zsh_data[i].GetBinError(n)

                        if zdr:
                            try:
                                zsh_resid[i].SetBinError(n, R*(sD/D))
                            except:
                                zsh_resid[i].SetBinError(n, 0)
                        else:
                            zsh_resid[i].SetBinError(n, sD)


                    ### now just formatting stuff
                    zsh_resid[i].SetTitle('')

                    zsh_resid[i].SetAxisRange(0, zmaxE, 'x')
                    if zdr: zsh_resid[i].SetAxisRange(0.8, 1.2, 'y')
                    else:   zsh_resid[i].SetAxisRange(-0.4, 0.4, 'y')

                    zsh_resid[i].SetXTitle('energy (keV)')
                    zsh_resid[i].GetXaxis().SetTitleFont(font)
                    zsh_resid[i].GetXaxis().SetTitleSize(size)
                    zsh_resid[i].GetXaxis().SetTitleOffset(xoff)
                    zsh_resid[i].GetXaxis().SetLabelFont(font)
                    zsh_resid[i].GetXaxis().SetLabelSize(size)
                    zsh_resid[i].GetXaxis().SetLabelOffset(loff)

                    if zdr: zsh_resid[i].SetYTitle('data/mc')
                    else:   zsh_resid[i].SetYTitle('data-mc (dru)')
                    zsh_resid[i].GetYaxis().SetTitleFont(font)
                    zsh_resid[i].GetYaxis().SetTitleSize(size)
                    zsh_resid[i].GetYaxis().SetTitleOffset(yoff)
                    zsh_resid[i].GetYaxis().SetLabelFont(font)
                    zsh_resid[i].GetYaxis().SetLabelSize(size)
                    zsh_resid[i].GetYaxis().SetLabelOffset(loff)

                    # '5' secondary and '05' primary
                    # means 5 divisions will be shown
                    zsh_resid[i].GetYaxis().SetNdivisions(505)

                    ### make a line at '0' or lrm variable
                    if zdr: zero = TLine(0, lrm, zmaxE, lrm)
                    else:   zero = TLine(0, 0, zmaxE, 0)
                    zsh_line[i] = zero
                    zsh_line[i].SetLineColor(kRed)
                    zsh_line[i].SetLineWidth(1)

                    zsh_resid[i].Draw()
                    zsh_line[i].Draw('same')
                    zsh_bot[i].Update()

                # update and save
                zsh_canv[i].Update()
                zsh_canv[i].Print(plotdir+'/a_singleHitZoom_c'+str(i+1)+'.png')
                # save log scale plot
                zsh_top[i].SetLogy(1)
                zsh_canv[i].Update()
                zsh_canv[i].Print(plotdir+'/a_singleHitZoom_c'+str(i+1)+'_log.png')
                # switch zoomLog setting for display
                zsh_top[i].SetLogy(zoomLog)
                zsh_canv[i].Update()
                
                
    ### Plot the 0-40 keV region of multi hit
    #-----------------------------------------
    #if indi:
    if save_zoom:
        if not show_zoom: gROOT.SetBatch(1)
        
        zmh_canv = [0 for x in range(numx)]
        zmh_hist = [0 for x in range(numx)]
        zmh_data = [0 for x in range(numx)]
        zmh_total = makeTotal100('M', 0, params[0])
        zmh_resid = makeResid100('M', 0, params[0])
        zmh_line = [0 for x in range(numx)]
        zmh_top = [0 for x in range(numx)]
        zmh_bot = [0 for x in range(numx)]
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:
                
                zmh_canv[i] = TCanvas('zmcanv-'+str(i+1),'zmcanv-'+str(i+1),
                                0, 0, 1200, 700)

                zmpad1 = TPad('zmpad1_'+str(i),'zmpad1_'+str(i),0,fraction,1,1)
                zmh_top[i] = (zmpad1)
                zmh_top[i].SetBottomMargin(botmarg)
                zmh_top[i].SetBorderMode(0)
                #zmh_top[i].SetLogy(1)
                zmh_top[i].SetLogy(0)
                zmh_top[i].Draw()
                
                if showResid:
                    zmpad2 = TPad('zmpad2_'+str(i),'zmpad2_'+str(i),0,0,1,fraction)
                    zmh_bot[i] = (zmpad2)
                    zmh_bot[i].SetTopMargin(0.05)
                    zmh_bot[i].SetBottomMargin(0.3)
                    zmh_bot[i].SetBorderMode(0)
                    #zmh_bot[i].SetLogy(1)
                    zmh_bot[i].SetLogy(0)
                    zmh_bot[i].Draw()
                
                
                zmh_top[i].cd()

                # rebin total here?
                zmh_total[i].Rebin(xstal_settings['M0']['rebin'])
                
                prim_list = toppad[1][0][i].GetListOfPrimitives()
                for bla in prim_list:
                    #print(bla.GetName(), bla)

                    if 'data' in bla.GetName():
                        zmh_data[i] = (toppad[1][0][i].Clone()).GetPrimitive(bla.GetName())
                        
                    if 'total' in bla.GetName():
                        #zmh_total[i] = (toppad[1][0][i].Clone()).GetPrimitive(bla.GetName())
                        continue
                    
                    try:
                        zmh_hist[i] = (toppad[1][0][i].Clone()).GetPrimitive(bla.GetName())
                        
                        zmh_hist[i].SetXTitle('energy (keV)')
                        zmh_hist[i].GetXaxis().SetTitleFont(font)
                        zmh_hist[i].GetXaxis().SetTitleSize(size)
                        #zmh_hist[i].GetXaxis().SetTitleOffset(xoff)
                        zmh_hist[i].GetXaxis().SetTitleOffset(1.1)
                        zmh_hist[i].GetXaxis().SetLabelFont(font)
                        zmh_hist[i].GetXaxis().SetLabelSize(size)
                        zmh_hist[i].GetXaxis().SetLabelOffset(loff)
                        
                        zmh_hist[i].SetYTitle('counts (dru)')
                        zmh_hist[i].GetYaxis().SetTitleFont(font)
                        zmh_hist[i].GetYaxis().SetTitleSize(size)
                        zmh_hist[i].GetYaxis().SetTitleOffset(yoff)
                        zmh_hist[i].GetYaxis().SetLabelFont(font)
                        zmh_hist[i].GetYaxis().SetLabelSize(size)
                        zmh_hist[i].GetYaxis().SetLabelOffset(loff)
                        
                        zmh_hist[i].SetAxisRange(0, zmaxE, 'x')

                        if   i+1 == 1: zmh_hist[i].SetAxisRange(0.01, 4.5, 'y')
                        elif i+1 == 2: zmh_hist[i].SetAxisRange(0.01, 4.5, 'y')
                        elif i+1 == 3: zmh_hist[i].SetAxisRange(0.01, 2.5, 'y')
                        elif i+1 == 4: zmh_hist[i].SetAxisRange(0.01, 2.0, 'y')
                        elif i+1 == 5: zmh_hist[i].SetAxisRange(0.01, 5.0, 'y')
                        elif i+1 == 6: zmh_hist[i].SetAxisRange(0.01, 1.2, 'y')
                        elif i+1 == 7: zmh_hist[i].SetAxisRange(0.01, 1.2, 'y')
                        elif i+1 == 8: zmh_hist[i].SetAxisRange(0.01, 5.0, 'y')
                        else:          zmh_hist[i].SetAxisRange(0.01, 5.0, 'y')

                        if 'data' in bla.GetName():
                            zmh_hist[i].Draw('same')
                        else:
                            zmh_hist[i].Draw(dopt)
                            zmh_total[i].Add(zmh_hist[i])
                            
                    except: continue

                if showTotal:
                    if redtotal:
                        zmh_total[i].SetMarkerColor(kRed)
                        zmh_total[i].SetLineColor(kRed)
                    zmh_total[i].Draw(dopt)
                    #legs[1][0][i].AddEntry(zmh_total[i], "Total", legopt)
                    
                if showLeg:
                    try:
                        #leg = TLegend(toppad[0][0][i].GetPrimitive(bla.GetName()))
                        #leg.Draw("same")
                        legs[1][0][i].Draw('same')
                    except: continue

                # update pad
                zmh_top[i].Update()
                
                if showResid:
                    zmh_bot[i].cd()

                    zmh_resid[i].Rebin(xstal_settings['M0']['rebin'])

                    if zdr:
                        if debug: print('DEBUG: hist divide #4')
                        try:
                            zmh_resid[i].Divide(zmh_data[i], zmh_total[i])
                        except:
                            print('WARNING: could not divide zoomed low-energy multi-hit')
                            print('         data nbins = {0} : total nbins = {1}'
                                  .format(zmh_data[i].GetXaxis().GetNbins(),
                                          zmh_total[i].GetXaxis().GetNbins()))
                            pass
                    else:
                        try:
                            zmh_resid[i] = zmh_data[i] - zmh_total[i]
                        except:
                            print('WARNING: could not subtract zoomed low-energy multi-hit')
                            print('         data nbins = {0} : total nbins = {1}'
                                  .format(zmh_data[i].GetXaxis().GetNbins(),
                                          zmh_total[i].GetXaxis().GetNbins()))
                            pass

                    ### tweak the errors on the residual
                    for n in range(zmh_resid[i].GetNbinsX()+1):
                        R  = zmh_resid[i].GetBinContent(n)
                        D  = zmh_data[i].GetBinContent(n)
                        sD = zmh_data[i].GetBinError(n)

                        if zdr:
                            try:
                                zmh_resid[i].SetBinError(n, R*(sD/D))
                            except:
                                zmh_resid[i].SetBinError(n, 0)
                        else:
                            zmh_resid[i].SetBinError(n, sD)

                    ### now just formatting stuff
                    zmh_resid[i].SetTitle('')

                    zmh_resid[i].SetAxisRange(0, zmaxE, 'x')
                    if zdr: zmh_resid[i].SetAxisRange(0.8, 1.2, 'y')
                    else:   zmh_resid[i].SetAxisRange(-0.4, 0.4, 'y')

                    zmh_resid[i].SetXTitle('energy (keV)')
                    zmh_resid[i].GetXaxis().SetTitleFont(font)
                    zmh_resid[i].GetXaxis().SetTitleSize(size)
                    zmh_resid[i].GetXaxis().SetTitleOffset(xoff)
                    zmh_resid[i].GetXaxis().SetLabelFont(font)
                    zmh_resid[i].GetXaxis().SetLabelSize(size)
                    zmh_resid[i].GetXaxis().SetLabelOffset(loff)

                    if zdr: zmh_resid[i].SetYTitle('data/mc')
                    else:   zmh_resid[i].SetYTitle('data-mc (dru)')
                    zmh_resid[i].GetYaxis().SetTitleFont(font)
                    zmh_resid[i].GetYaxis().SetTitleSize(size)
                    zmh_resid[i].GetYaxis().SetTitleOffset(yoff)
                    zmh_resid[i].GetYaxis().SetLabelFont(font)
                    zmh_resid[i].GetYaxis().SetLabelSize(size)
                    zmh_resid[i].GetYaxis().SetLabelOffset(loff)

                    # '5' secondary and '05' primary
                    # means 5 divisions will be shown
                    zmh_resid[i].GetYaxis().SetNdivisions(505)


                    ### make a line at '0' or lrm variable
                    if zdr: zero = TLine(0, lrm, zmaxE, lrm)
                    else:   zero = TLine(0, 0, zmaxE, 0)
                    zmh_line[i] = zero
                    zmh_line[i].SetLineColor(kRed)
                    zmh_line[i].SetLineWidth(1)

                    zmh_resid[i].Draw()
                    zmh_line[i].Draw('same')
                    zmh_bot[i].Update()
                    
                zmh_canv[i].Update()
                zmh_canv[i].Print(plotdir+'/a_multiHitZoom_c'+str(i+1)+'.png')

    gROOT.SetBatch(batch)
    
    ### Plot the single hit alpha channel
    #-------------------------------------
    #if indi:
    if save_alpha:
        if not show_alpha: gROOT.SetBatch(1)
        
        ash_canv = [0 for x in range(numx)]
        ash_hist = [{} for x in range(numx)]
        ash_data = [0 for x in range(numx)]
        ash_total = makeTotal100('S', 2, params[2])
        ash_resid = makeResid100('S', 2, params[2])
        ash_line = [0 for x in range(numx)]
        ash_top = [0 for x in range(numx)]
        ash_bot = [0 for x in range(numx)]
        ash_leg = [0 for x in range(numx)]
        
        popts = xstal_settings['S2']
        log = int(popts['log'])
        yrange = popts['yrange'][log]
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:
                
                ash_canv[i] = TCanvas('ascanv-'+str(i+1),'ascanv-'+str(i+1),
                                0, 0, 1200, 700)
                
                aspad1 = TPad('aspad1_'+str(i),'aspad1_'+str(i),0,fraction,1,1)
                ash_top[i] = (aspad1)
                ash_top[i].SetBottomMargin(botmarg)
                ash_top[i].SetBorderMode(0)
                ash_top[i].SetLogy(log)
                #ash_top[i].SetLogy(0)
                ash_top[i].Draw()

                if showResid:
                    aspad2 = TPad('aspad2_'+str(i),'aspad2_'+str(i),0,0,1,fraction)
                    ash_bot[i] = (aspad2)
                    ash_bot[i].SetTopMargin(0.05)
                    ash_bot[i].SetBottomMargin(0.3)
                    ash_bot[i].SetBorderMode(0)
                    #ash_bot[i].SetLogy(1)
                    ash_bot[i].SetLogy(0)
                    ash_bot[i].Draw()

                
                ash_top[i].cd()
                
                for key in datkeys:
                    if 'x'+str(i+1) in key and '-cS' in key and '-e2' in key:
                        ash_data[i] = deepcopy(data[key]['hist'])
                        
                        ash_data[i].SetXTitle('energy (keV)')
                        ash_data[i].GetXaxis().SetTitleFont(font)
                        ash_data[i].GetXaxis().SetTitleSize(size)
                        #ash_data[i].GetXaxis().SetTitleOffset(xoff)
                        ash_data[i].GetXaxis().SetTitleOffset(1.1)
                        ash_data[i].GetXaxis().SetLabelFont(font)
                        ash_data[i].GetXaxis().SetLabelSize(size)
                        ash_data[i].GetXaxis().SetLabelOffset(loff)
                                                
                        ash_data[i].SetYTitle('counts (dru)')
                        ash_data[i].GetYaxis().SetTitleFont(font)
                        ash_data[i].GetYaxis().SetTitleSize(size)
                        ash_data[i].GetYaxis().SetTitleOffset(yoff)
                        ash_data[i].GetYaxis().SetLabelFont(font)
                        ash_data[i].GetYaxis().SetLabelSize(size)
                        ash_data[i].GetYaxis().SetLabelOffset(loff)
                        
                        ash_data[i].SetAxisRange(popts['xrange'][0], popts['xrange'][1], 'x')
                        ash_data[i].SetAxisRange(yrange[0], yrange[1], 'y')

                        ash_data[i].Draw('same')
                        ldatname = data[key]['info']['tag']+'  '+data[key]['info']['build']
                        if dataname: ldatname = dataname
                            
                # join hists of same type together -fx -qx
                for key in bakkeys:
                    if 'x'+str(i+1) in key and '-cS' in key and '-e2' in key:
                        if justQ and '-q'+str(justQ)+'-' not in key: continue
                        if bkgs[key]['hist'].Integral() > 0:
                            histkey = re.split('-f(?:[0-9])-', key)[0]
                            if combineExp and 'expo' in histkey:
                                bits = histkey.split('-')
                                bits[1] = bits[1].split('expo')[0]
                                histkey = '-'.join(bits)
                            if histkey in ash_hist[i]:
                                ash_hist[i][histkey].Add(bkgs[key]['hist'])
                            else:
                                ash_hist[i][histkey] = (bkgs[key]['hist'])
                
                for key in sigkeys:
                    if 'x'+str(i+1) in key and '-cS' in key and '-e2' in key:
                        if justQ and '-q'+str(justQ)+'-' not in key: continue
                        if sigs[key]['hist'].Integral() > 0:
                            histkey = re.split('-f(?:[0-9])-', key)[0]
                            if combineExp and 'expo' in histkey:
                                bits = histkey.split('-')
                                bits[1] = bits[1].split('expo')[0]
                                histkey = '-'.join(bits)
                            if histkey in ash_hist[i]:
                                ash_hist[i][histkey].Add(sigs[key]['hist'])
                            else:
                                ash_hist[i][histkey] = (sigs[key]['hist'])
                
                N_hists = len(ash_hist[i])
                
                leg_ymin = 1-((N_hists+2)*0.08)
                if leg_ymin < 0.05: leg_ymin = 0.05
                leg = TLegend(0.75, leg_ymin, 0.98, 0.90)
                ash_leg[i] = (leg)
                ash_leg[i].SetFillColor(0)
                ash_leg[i].SetBorderSize(0)
                legopt = 'LPE'
                
                ash_leg[i].AddEntry(ash_data[i], ldatname, legopt)

                histkeys = sorted([key for key in ash_hist[i]])

                # fix the colors for the paper
                isoColors = {
                    'internal Pb210': kBlue    +1,
                    'internal Th228': kCyan    +1,
                    'nai-surf Pb210': kGreen   +1,
                    'teflon Pb210': kOrange  +1,
                    'teflon-surf Pb210': kMagenta +1,
                }

                
                for N, histkey in enumerate(histkeys):
                    #ash_hist[i][histkey].SetMarkerColor(alphaColors(N))
                    #ash_hist[i][histkey].SetLineColor(alphaColors(N))
                    #ash_hist[i][histkey].Draw(dopt)
                    # tweak legend labels
                    bits = histkey.split('-')
                    if combineExp and 'surf' in histkey:
                        bits[1] = bits[1].replace('surf', '-surf')
                    # use just primary isotope
                    bits[2] = bits[2].split('_')[0]
                    
                    # fix U238 group names - 2023-05-23
                    if bits[2] == 'Th230': bits[2] = 'U234'
                    if bits[2] == 'Ra226': bits[2] = 'Th230'
                    if bits[2] == 'Rn222': bits[2] = 'Ra226'
                    
                    legname = bits[1]+' '+bits[2]
                    if legname in isoColors:
                        ash_hist[i][histkey].SetMarkerColor(isoColors[legname])
                        ash_hist[i][histkey].SetLineColor(isoColors[legname])
                    else:
                        ash_hist[i][histkey].SetMarkerColor(alphaColors(N))
                        ash_hist[i][histkey].SetLineColor(alphaColors(N))
                    ash_hist[i][histkey].Draw(dopt)
                    ash_leg[i].AddEntry(ash_hist[i][histkey], legname, legopt)
                    ash_total[i].Add(ash_hist[i][histkey])

                    
                # rebin total for resid plot
                data_bins = ash_data[i].GetNbinsX()
                total_bins = ash_total[i].GetNbinsX()
                if total_bins != data_bins:
                    ash_total[i].Rebin(popts['rebin'])
                ash_total[i].SetAxisRange(popts['xrange'][0], popts['xrange'][1], 'x')
                if showTotal:
                    if redtotal:
                        ash_total[i].SetMarkerColor(kRed)
                        ash_total[i].SetLineColor(kRed)
                    ash_total[i].Draw(dopt)
                    ltotname = 'Total'
                    if totname: ltotname = totname
                    ash_leg[i].AddEntry(ash_total[i], ltotname, legopt)

                if showLeg:
                    ash_leg[i].Draw('same')
                
                if showResid:
                    # plot residuals
                    ash_bot[i].cd()
                    ash_resid[i].Rebin(popts['rebin'])

                    # remove very small bins or the Divide returns 'inf'
                    ash_total[i] = cleanSmallBins(ash_total[i])
                    if debug:
                        print('DEBUG: hist divide #5')
                    try:
                        ash_resid[i].Divide(ash_data[i], ash_total[i])
                    except:
                        print('WARNING: could not divide single-hit alpha')

                    ### tweak the errors on the residual
                    for n in range(ash_resid[i].GetNbinsX()+1):
                        R  = ash_resid[i].GetBinContent(n)
                        D  = ash_data[i].GetBinContent(n)
                        sD = ash_data[i].GetBinError(n)
                        try:
                            ash_resid[i].SetBinError(n, R*(sD/D))
                        except:
                            ash_resid[i].SetBinError(n, sD)

                    ### now just formatting stuff
                    ash_resid[i].SetTitle('')

                    ash_resid[i].SetAxisRange(popts['xrange'][0], popts['xrange'][1], 'x')
                    ash_resid[i].SetAxisRange(0.0, 2, 'y')

                    ash_resid[i].SetXTitle('energy (keV)')
                    ash_resid[i].GetXaxis().SetTitleFont(font)
                    ash_resid[i].GetXaxis().SetTitleSize(size)
                    ash_resid[i].GetXaxis().SetTitleOffset(xoff)
                    ash_resid[i].GetXaxis().SetLabelFont(font)
                    ash_resid[i].GetXaxis().SetLabelSize(size)
                    ash_resid[i].GetXaxis().SetLabelOffset(loff)

                    ash_resid[i].SetYTitle('data / MC')
                    ash_resid[i].GetYaxis().SetTitleFont(font)
                    ash_resid[i].GetYaxis().SetTitleSize(size)
                    ash_resid[i].GetYaxis().SetTitleOffset(yoff)
                    ash_resid[i].GetYaxis().SetLabelFont(font)
                    ash_resid[i].GetYaxis().SetLabelSize(size)
                    ash_resid[i].GetYaxis().SetLabelOffset(loff)

                    # '5' secondary and '05' primary
                    # means 5 divisions will be shown
                    ash_resid[i].GetYaxis().SetNdivisions(505)

                    ### make a line at lrm variable
                    zero = TLine(popts['xrange'][0], lrm, popts['xrange'][1], lrm)
                    ash_line[i] = zero
                    ash_line[i].SetLineColor(kRed)
                    ash_line[i].SetLineWidth(1)

                    ash_resid[i].Draw()
                    ash_line[i].Draw('same')
                
                ash_canv[i].Update()
                ash_canv[i].Print(plotdir+'/a_alphaSingleHit_c'+str(i+1)+'.png')

                
    ### Plot the multi hit alpha channel
    #--------------------------------------
    #if indi:
    if save_alpha and and_multi:
        if not show_alpha: gROOT.SetBatch(1)
        
        amh_canv = [0 for x in range(numx)]
        amh_hist = [{} for x in range(numx)]
        amh_data = [0 for x in range(numx)]
        amh_total = makeTotal100('M', 2, params[2])
        amh_resid = makeResid100('M', 2, params[2])
        amh_line = [0 for x in range(numx)]
        amh_top = [0 for x in range(numx)]
        amh_bot = [0 for x in range(numx)]
        amh_leg = [0 for x in range(numx)]
        
        popts = xstal_settings['M2']
        log = int(popts['log'])
        yrange = popts['yrange'][log]
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:
                
                amh_canv[i] = TCanvas('amcanv-'+str(i+1),'amcanv-'+str(i+1),
                                0, 0, 1200, 700)
                
                ampad1 = TPad('ampad1_'+str(i),'ampad1_'+str(i),0,fraction,1,1)
                amh_top[i] = (ampad1)
                amh_top[i].SetBottomMargin(botmarg)
                amh_top[i].SetBorderMode(0)
                amh_top[i].SetLogy(log)
                #amh_top[i].SetLogy(0)
                amh_top[i].Draw()

                if showResid:
                    ampad2 = TPad('ampad2_'+str(i),'ampad2_'+str(i),0,0,1,fraction)
                    amh_bot[i] = (ampad2)
                    amh_bot[i].SetTopMargin(0.05)
                    amh_bot[i].SetBottomMargin(0.3)
                    amh_bot[i].SetBorderMode(0)
                    #amh_bot[i].SetLogy(1)
                    amh_bot[i].SetLogy(0)
                    amh_bot[i].Draw()

                
                amh_top[i].cd()
                
                for key in datkeys:
                    if 'x'+str(i+1) in key and '-cM' in key and '-e2' in key:
                        amh_data[i] = deepcopy(data[key]['hist'])
                        
                        amh_data[i].SetXTitle('energy (keV)')
                        amh_data[i].GetXaxis().SetTitleFont(font)
                        amh_data[i].GetXaxis().SetTitleSize(size)
                        #amh_data[i].GetXaxis().SetTitleOffset(xoff)
                        amh_data[i].GetXaxis().SetTitleOffset(1.1)
                        amh_data[i].GetXaxis().SetLabelFont(font)
                        amh_data[i].GetXaxis().SetLabelSize(size)
                        amh_data[i].GetXaxis().SetLabelOffset(loff)
                                                
                        amh_data[i].SetYTitle('counts (dru)')
                        amh_data[i].GetYaxis().SetTitleFont(font)
                        amh_data[i].GetYaxis().SetTitleSize(size)
                        amh_data[i].GetYaxis().SetTitleOffset(yoff)
                        amh_data[i].GetYaxis().SetLabelFont(font)
                        amh_data[i].GetYaxis().SetLabelSize(size)
                        amh_data[i].GetYaxis().SetLabelOffset(loff)
                        
                        amh_data[i].SetAxisRange(popts['xrange'][0], popts['xrange'][1], 'x')
                        amh_data[i].SetAxisRange(yrange[0], yrange[1], 'y')

                        amh_data[i].Draw('same')
                        ldatname = data[key]['info']['tag']+'  '+data[key]['info']['build']
                        if dataname: ldatname = dataname
                        
                # join hists of same type together -fx -qx
                for key in bakkeys:
                    if 'x'+str(i+1) in key and '-cM' in key and '-e2' in key:
                        if justQ and '-q'+str(justQ)+'-' not in key: continue
                        if bkgs[key]['hist'].Integral() > 0:
                            histkey = re.split('-f(?:[0-9])-', key)[0]
                            if combineExp and 'expo' in histkey:
                                bits = histkey.split('-')
                                bits[1] = bits[1].split('expo')[0]
                                histkey = '-'.join(bits)
                            if histkey in amh_hist[i]:
                                amh_hist[i][histkey].Add(bkgs[key]['hist'])
                            else:
                                amh_hist[i][histkey] = (bkgs[key]['hist'])
                
                for key in sigkeys:
                    if 'x'+str(i+1) in key and '-cM' in key and '-e2' in key:
                        if justQ and '-q'+str(justQ)+'-' not in key: continue
                        if sigs[key]['hist'].Integral() > 0:
                            histkey = re.split('-f(?:[0-9])-', key)[0]
                            if combineExp and 'expo' in histkey:
                                bits = histkey.split('-')
                                bits[1] = bits[1].split('expo')[0]
                                histkey = '-'.join(bits)
                            if histkey in amh_hist[i]:
                                amh_hist[i][histkey].Add(sigs[key]['hist'])
                            else:
                                amh_hist[i][histkey] = (sigs[key]['hist'])
                
                N_hists = len(amh_hist[i])

                leg_ymin = 1-((N_hists+2)*0.08)
                if leg_ymin < 0.05: leg_ymin = 0.05
                leg = TLegend(0.75, leg_ymin, 0.98, 0.90)
                amh_leg[i] = (leg)
                amh_leg[i].SetFillColor(0)
                amh_leg[i].SetBorderSize(0)
                legopt = 'LPE'
                
                amh_leg[i].AddEntry(amh_data[i], ldatname, legopt)

                histkeys = sorted([key for key in amh_hist[i]])
                for N, histkey in enumerate(histkeys):
                    amh_hist[i][histkey].SetMarkerColor(alphaColors(N))
                    amh_hist[i][histkey].SetLineColor(alphaColors(N))
                    amh_hist[i][histkey].Draw(dopt)
                    # tweak legend labels
                    bits = histkey.split('-')
                    if combineExp and 'surf' in histkey:
                        bits[1] = bits[1].replace('surf', '-surf')
                    # use just primary isotope
                    bits[2] = bits[2].split('_')[0]
                    legname = bits[1]+' '+bits[2]
                    amh_leg[i].AddEntry(amh_hist[i][histkey], legname, legopt)
                    amh_total[i].Add(amh_hist[i][histkey])

                    
                # rebin total for resid plot
                data_bins = amh_data[i].GetNbinsX()
                total_bins = amh_total[i].GetNbinsX()
                if total_bins != data_bins:
                    amh_total[i].Rebin(popts['rebin'])
                amh_total[i].SetAxisRange(popts['xrange'][0], popts['xrange'][1], 'x')
                if showTotal:
                    if redtotal:
                        amh_total[i].SetMarkerColor(kRed)
                        amh_total[i].SetLineColor(kRed)
                    amh_total[i].Draw(dopt)
                    ltotname = 'Total'
                    if totname: ltotname = totname
                    amh_leg[i].AddEntry(amh_total[i], ltotname, legopt)

                if showLeg:
                    amh_leg[i].Draw('same')
                
                if showResid:
                    # plot residuals
                    amh_bot[i].cd()
                    amh_resid[i].Rebin(popts['rebin'])

                    # remove very small bins or the Divide returns 'inf'
                    amh_total[i] = cleanSmallBins(amh_total[i])
                    if debug: print('DEBUG: hist divide #6')
                    try:
                        amh_resid[i].Divide(amh_data[i], amh_total[i])
                    except:
                        print('WARNING: could not divide multi-hit alpha')

                    ### tweak the errors on the residual
                    for n in range(amh_resid[i].GetNbinsX()+1):
                        R  = amh_resid[i].GetBinContent(n)
                        D  = amh_data[i].GetBinContent(n)
                        sD = amh_data[i].GetBinError(n)
                        try:
                            amh_resid[i].SetBinError(n, R*(sD/D))
                        except:
                            amh_resid[i].SetBinError(n, sD)

                    ### now just formatting stuff
                    amh_resid[i].SetTitle('')

                    amh_resid[i].SetAxisRange(popts['xrange'][0], popts['xrange'][1], 'x')
                    amh_resid[i].SetAxisRange(0.0, 2, 'y')

                    amh_resid[i].SetXTitle('energy (keV)')
                    amh_resid[i].GetXaxis().SetTitleFont(font)
                    amh_resid[i].GetXaxis().SetTitleSize(size)
                    amh_resid[i].GetXaxis().SetTitleOffset(xoff)
                    amh_resid[i].GetXaxis().SetLabelFont(font)
                    amh_resid[i].GetXaxis().SetLabelSize(size)
                    amh_resid[i].GetXaxis().SetLabelOffset(loff)

                    amh_resid[i].SetYTitle('data / MC')
                    amh_resid[i].GetYaxis().SetTitleFont(font)
                    amh_resid[i].GetYaxis().SetTitleSize(size)
                    amh_resid[i].GetYaxis().SetTitleOffset(yoff)
                    amh_resid[i].GetYaxis().SetLabelFont(font)
                    amh_resid[i].GetYaxis().SetLabelSize(size)
                    amh_resid[i].GetYaxis().SetLabelOffset(loff)

                    # '5' secondary and '05' primary
                    # means 5 divisions will be shown
                    amh_resid[i].GetYaxis().SetNdivisions(505)

                    ### make a line at lrm variable
                    zero = TLine(popts['xrange'][0], lrm, popts['xrange'][1], lrm)
                    amh_line[i] = zero
                    amh_line[i].SetLineColor(kRed)
                    amh_line[i].SetLineWidth(1)

                    amh_resid[i].Draw()
                    amh_line[i].Draw('same')
                
                amh_canv[i].Update()
                amh_canv[i].Print(plotdir+'/a_alphaMultiHit_c'+str(i+1)+'.png')

                
    #-----------------------------------------------------------------
    ### print out the fit results
    if fitting:
        print('\n\n!!!!!  FIT RESULTS  !!!!!\n')
        for key in resultskeys:
            if int(key)+1 in justthese:
                for line in fitresults[key]:
                    print('{0}'.format(line))
    #-----------------------------------------------------------------

    
    # delete the extra crap
    #-------------------------------

    # pre fit
    if not show_prefit:
        try: del pfcanv
        except: pass
        
    # post fit
    if not show_postfit:
        try: del mcanv
        except: pass

    # combined 2x2 plots
    #try: del combPlots
    #except: pass

    if not show_zoom:
        try: del zsh_canv
        except: pass
        try: del zmh_canv
        except: pass

    if not show_alpha:
        try: del ash_canv
        except: pass
        try: del amh_canv
        except: pass

    # individual plots
    if not show_indi:
        try: del canvs
        except: pass

    if fitting:
        print('Time to complete the fit = {0} sec \n'.format(fitTime))
    
    if not batch: input('\n[Enter] to quit\n')
    else: print('Running in batch mode - quitting now! \n')

    return


if __name__ == "__main__":
    main(sys.argv[1:])

