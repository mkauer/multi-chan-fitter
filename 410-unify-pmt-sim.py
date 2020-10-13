#!/usr/bin/env python

############################################################
# Matt Kauer - mkauer@icecube.wisc.edu
#-----------------------------------------------------------
# Make PMT sim truely global and
# and not a collection of pmt pairs

# see CHANGELOG
# version: 2020-08-18
############################################################

import os,sys,re
import shutil
import socket
import copy
import math
import numpy as np
import datetime

# ROOT6 FIX
import ctypes
from ctypes import *

script = os.path.basename(__file__)
V = 'v'+(script.split('-')[0])
print 'INFO: running script -->', script, '('+V+')'

import ROOT
from ROOT import *
#ROOT.gErrorIgnoreLevel = kWarning
ROOT.gErrorIgnoreLevel = kError
vroot = ROOT.__version__
print 'INFO: using root -->', vroot

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs410 import *


### batch job?
#print 'HOSTNAME:', socket.gethostname()
batch = onCup()
#batch = 1

### get the total number of possible crystals (including lsveto)
numx = numX()

### print debug info
debug = 0

# draw options
dopt='hist same'

### ==========  GENERAL INPUTS  ======================================
### note to add to saved plot names?
note = 0
#note = ''

#mcfile = 'backgrounds_404.txt' # pretty solid baseline file
#mcfile = 'backgrounds_405.txt' # changed things to govindas numbers
#mcfile = 'backgrounds_406.txt' # going back to 404 and modifying
#mcfile = 'backgrounds_407.txt' # new expo file format
#mcfile = 'backgrounds_408.txt' # add back other Te cosmogenics
#mcfile = 'backgrounds_409.txt' # can only use "side" components
#mcfile = 'backgrounds_410.txt' # try fixing pmt activities
#mcfile = 'backgrounds_411.txt' # play with teflon surface
#mcfile = 'backgrounds_412.txt' # surface tweaks for root-6

#mcfile = 'backgrounds_413.txt' # new build410 MC

#mcfile = 'backgrounds_414.txt' # try removing lsveto mc
#mcfile = 'backgrounds_414_tweaks.txt' # tweak some activities
#mcfile = 'backgrounds_414_tweaks2.txt' # tweak some activities
#mcfile = 'backgrounds_414_tweaks2_C1.txt' # try C1
#mcfile = 'backgrounds_414_tweaks2_C5.txt' # try C5
#mcfile = 'backgrounds_414_tweaks2_C8.txt' # try C8
#mcfile = 'backgrounds_414_tweaks3.txt' # fit externals with C1

#mcfile = 'backgrounds_414_tweaks2_update.txt' # manually tweaked things
#mcfile = 'backgrounds_414_tweaks3_update.txt' # fit externals with C1
mcfile = 'backgrounds_414_I128.txt' # try fit with I128


#mcfile = 'backgrounds_415.txt' # try removing plastic mc
#mcfile = 'backgrounds_415_tweaks.txt' # tweak some activities
#mcfile = 'backgrounds_415_tweaks2.txt' # tweak some activities


#mcfile = 'testing-expo-surf.txt' # 2020-05-12
#mcfile = 'testing-all-bkgs-for-plotting.txt'
#mcfile = 'testing-expo-pb210.txt' # 2020-03-24
#mcfile = 'testing-data-400.txt' # 2020-03-26
#mcfile = 'testing-alpha-cuts.txt' # 2020-04-13
#mcfile = 'testing-unified-pmt.txt' # 2020-08-07
#mcfile = 'testing-copper-shield.txt' # 2020-08-08
#mcfile = 'testing-lsveto-resolution.txt' # 2020-08-11
#mcfile = 'testing-steel.txt' # 2020-08-13
#mcfile = 'testing-u238-th232.txt' # 2020-08-13
#mcfile = 'testing-na24.txt' # 2020-10-12


print 'INFO: using backgrounds config file -->', mcfile


### ==========  OPTIMIZATION OPTIONS  ================================
### test MC energy shift in bins (12bins = 1keV)
### currently just shifting teflon
### might need different shifting for different xstals
# for xstals 1   2   3   4   5   6   7   8   9
binShift =  [4,  4,  6,  8, -1,  8,  6, -1, -1]
#binShift =  [-1, 0, 11, 11, -1, 11, 11, -1, -1]
shiftWhat = ['teflonbulk','teflonsurf']
#shiftWhat = ['teflon']

### set the fitter step size
### can help to go smaller than default 1e-2
stepSize = 0
stepSize = 1e-4

### test MC smoothing?
### smoothing window in +/- number of bins
smoothing = 0
#smoothing = 5
smoothWhat = ['cushield','steel']


### ==========  FITTING OPTIONS  =====================================
### select channels to fit
#fitchans = 'SM'
fitchans = 'S'  # testing low E fit

### Let's finally try different fit ranges!!
fitranges = [{} for x in range(numx)]
for i in range(numx):

    ### defaults for crystals
    # --------------------------------------------------------------
    fitranges[i]['S0'] = [1,  2,   90]  # singl-hit low-energy
    fitranges[i]['S1'] = [6, 80, 4000]  # singl-hit high-energy
    fitranges[i]['M0'] = [1,  2,   90]  # multi-hit low-energy
    fitranges[i]['M1'] = [6, 80, 4000]  # multi-hit high-energy
    
    # testing low E fit
    fitranges[i]['S0'] = [1, 4, 40]
    fitranges[i]['S1'] = [1, 0,  0]
    fitranges[i]['M0'] = [1, 2, 40]
    #fitranges[i]['M0'] = [1, 20, 40]
    fitranges[i]['M1'] = [1, 0,  0]
    
    # testing high E fit
    #fitranges[i]['S0'] = [1, 0, 0]
    #fitranges[i]['S1'] = [4, 100, 2800]
    #fitranges[i]['M0'] = [1, 40, 100]
    #fitranges[i]['M1'] = [4, 100, 2800]

# special case for C1, C5, C8
#for i in [0, 4, 7]:
#    fitranges[i]['S0'] = [1, 10, 60]
#    fitranges[i]['S1'] = [6, 0, 0]
#    fitranges[i]['M0'] = [1, 10, 90]
#    fitranges[i]['M1'] = [6, 80, 4000]
        
    
### defaults for lsveto
# --------------------------------------------------------------
fitranges[8]['S0'] = [0,0,0]
fitranges[8]['S1'] = [0,0,0]
fitranges[8]['M0'] = [0,0,0]
fitranges[8]['M1'] = [1, 100, 4000] # ext gamma testing


### set bounds separately for some crystals
# --------------------------------------------------------------
# c1 low-energy is messed up below 9 keV
#fitranges[0]['S0'] = [1, 9, 100]



### ==========  EXTRA MC OPTIONS  ====================================
### which MC to fit globally (to all crystals simultaneously)?
globalmc = ['pmt','plastic','lsveto','cushield','innersteel','steel','gamma']

### include bkgs from 'other' pmts and internals?
others    = 1

### plot components in groups? [0,1]
ingroups  = 1

### show the total? [0,1]
showTotal = 1

### show the legends? [0,1]
showlegs  = 1

### plot the total in gray==0, red==1? [0,1]
redtotal  = 1

### combine 'others' into the makePlots() plots?
combine   = 1

### force the reuse of all joined rootfiles in mcfile? [0,1,2]
### very nice for debugging
### [0] default - use whatever is specified in the backgrounds file
### [1] forces reusing of all data/bkgs/sigs
### [2] forces NOT reusing any data/bkgs/sigs
reuse = 0


### ==========  OTHER FITTING OPTIONS  ===============================
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


### ==========  PLOTTING OPTIONS  ====================================
### individual plots for all crystals? [0,1]
indi = 1

### select channels to plot
pltchans = 'SM'

### plotting ranges
loer = [0, 100]
hier = [0, 4000]
#hiYlo = 1e-3 # for 3000 X range
hiYlo = 1e-5 # for 4000 X range
eran = [loer, hier]

### special energy range for lsveto
lsHiE = 4000

### special range for zoomed in plot
zmaxE = 40

### zoomed in residual as dat-mc [0] or dat/mc [1]
zdr = 0

### rebin the final plots [1,inf]
loEplotRebin = 3
hiEplotRebin = 4

### use linear residual scale? [0,1]
linres = 1

### set y scale on resid plot
lrs = [0, 2]

### main plots in linear scale [0,1]
liny = 0


### ==========  CAN EFFECT FIT RESULTS  ==============================
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
    if onCup(): here = '/home/mkauer/mc-fitting'
    else: here = '/home/mkauer/COSINE/CUP/mc-fitting'
    """
    ### for saving the plots...
    plotdir = here+'/plots/c'
    for x in xstals:
        plotdir += str(x)
    if not os.path.exists(plotdir): 
        os.makedirs(plotdir)
    """
    #mcfile = os.path.join(here, mcfile)
    if not os.path.exists(os.path.join(here, mcfile)):
        print '\nERROR: could not find backgrounds file -->', os.path.join(here, mcfile)
        sys.exit()
    
    
    ### where everything gets loaded into dictionary
    #-----------------------------------------------------------------
    #-----------------------------------------------------------------
    allchans = uniqString(fitchans+pltchans)
    data, bkgs, sigs = build410(os.path.join(here, mcfile), others, reuse, allchans, xstals)
    
    ### make of list of crystal runtimes
    runtimes = [0 for x in range(numx)]
    usedkeys = []
    for i in range(numx):
        for key in data:
            #print key
            if 'x'+str(i+1) in key and 'x'+str(i+1) not in usedkeys:
                runtimes[i] = data[key]['runtime']
                usedkeys.append('x'+str(i+1))
                print 'INFO: c'+str(i+1)+' runtime =', round(runtimes[i],1), '(seconds)'
    #print runtimes
    
    datkeys = sortDataKeys92(data)
    if datsumw2:
        for key in datkeys:
            data[key]['hist'].Sumw2()
    
    ### scale into dru units
    if dru:
        data = scaleData70(data, 1)
        bkgs = scaleBkgs410(bkgs)
        sigs = scaleBkgs410(sigs)
    else:
        data = scaleData70(data, 0)
        ### FIX ME - need to fix scaleBkgs to use runtimes[]
        bkgs = scaleBkgs410(bkgs, runtimes)
        sigs = scaleBkgs410(sigs, runtimes)
    
    
    
    #=======================================================================
    ### testing MC energy scaling
    #=======================================================================

    if binShift != 0:
        #print 'energy shifting...'
        sigs = ScaleEnergy410(sigs, binShift, shiftWhat)
        bkgs = ScaleEnergy410(bkgs, binShift, shiftWhat)
    #sys.exit()
    
    #=======================================================================
    #=======================================================================
    
    
    
    ### make plots before combining?
    #makePlots93(bkgs, combine, others)
    #makePlots93(sigs, combine, others)
    #sys.exit()
    
    ### combine after scaling?
    ### FIX ME - combine but don't delete others?
    sigs = combineOthers100(sigs, globalmc)
    bkgs = combineOthers100(bkgs, globalmc)
    
    ### now sort and remove empty histos
    bkgs, bakkeys = sortSimKeys92(bkgs)
    sigs, sigkeys = sortSimKeys92(sigs)
    
    ### make plots after combining?
    #makePlots93(bkgs, combine, others)
    #makePlots93(sigs, combine, others)
    #sys.exit()
    
    ### do histogram smoothing?
    if smoothing:
        bkgs = smooth(bkgs, smoothWhat, smoothing)
        sigs = smooth(sigs, smoothWhat, smoothing)
        print 'INFO: done smoothing histograms'
    
    #-----------------------------------------------------------------
    #-----------------------------------------------------------------
    
    ### only plot crystals that have data
    justthese = []
    for i in range(1, numx+1):
        for key in datkeys:
            if 'x'+str(i) in key and i not in justthese:
                justthese.append(i)
    justthese.sort()
    if len(justthese) == 0:
        print 'ERROR: nothing to plot, quitting...'
        sys.exit()
    print 'INFO: plotting crystals -->', justthese
    
    ### create dir for saving the plots...
    plotdir = here+'/plots/c'
    for x in justthese:
        plotdir += str(x)
    if not os.path.exists(plotdir): 
        os.makedirs(plotdir)

    ### copy the backgrounds file before fitting so it isn't overwritten later - 2020-04-02
    shutil.copyfile(os.path.join(here, mcfile), os.path.join(plotdir, mcfile))

    ### copy the actual script as well just to have it on hand - 2020-04-02
    shutil.copyfile(os.path.join(here, script), os.path.join(plotdir, script))

    
    ### assume all data is using same runs and hist params
    try:    runtag = data[datkeys[0]]['info']['tag']
    except: runtag = 'none'
    try:    params = globalParams(data)
    except: params = globalParams(bkgs)
    
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

    
    ### create unique colors for fit sim
    Nc = len(uniqAll)
    print 'INFO: Total number of unique bkgs and sigs =', Nc
    # ROOT6 FIX
    colors, cis = rainbowSix(uniqAll)
    uniqColor = cis
    
    ### colors for the groups
    gis = {
        'internal':  kBlue,
        'cosmo':     kMagenta+1,
        'surface':   kCyan+1,
        'copper':    kYellow+1,
        'pmts':      kGreen+1,
        'plastic':   kOrange,
        'lsveto':    kOrange+1,
        'steel':     kYellow,
        'ext-gamma': kOrange+4,
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
        
        print 'INFO: global fits to -->', globstr
    
        print 'INFO: building histograms for fit...'

        fitting = 1
        
        fitsigs = {}
        fitbkgs = {}

        ### have separate global signals and backgrounds
        ### 2018-11-15
        fitglobsigs = {}
        fglobsigkeys = []
        fitglobbkgs = {}
        fglobbkgkeys = []

        fitdata = TH1F('globData', 'globData', 1, 0, 1)
        
        ### I think I need to build the fitdata histogram first
        ### store start/stop bins for each crystal
        xstalbins = [[0,0] for x in range(numx)]
        for i in range(numx):
            startbin = fitdata.GetXaxis().GetNbins()
            for C in fitchans:
                for dkey in datkeys:
                    # low Energy
                    if 'x'+str(i+1) in dkey and '-c'+C in dkey and '-e0' in dkey:
                        fitdata = extendHist(fitdata,
                                    fitPrep(data, dkey, fitranges[i][C+'0']))
                    # high Energy
                    if 'x'+str(i+1) in dkey and '-c'+C in dkey and '-e1' in dkey:
                        fitdata = extendHist(fitdata,
                                    fitPrep(data, dkey, fitranges[i][C+'1']))
            stopbin = fitdata.GetXaxis().GetNbins()
            ### this start/stop binning stuff might still be a little wonky
            xstalbins[i] = [startbin-1, stopbin-1]
            if debug: print 'DEBUG: bin range for C'+str(i+1)+' =', xstalbins[i]
        
        ### the number of bins in fitdata should be used...
        fitbins = fitdata.GetNbinsX()

        ftotal = TH1F('globTotal', 'globTotal', fitbins, 0, fitbins)
        ftotal.SetLineColor(kGray+1)
        ftotal.SetMarkerColor(kGray+1)
        ftotal.SetLineWidth(1)

        fresid = TH1F('globResid', 'globResid', fitbins, 0, fitbins)
        fresid.SetLineColor(kBlack)
        fresid.SetMarkerColor(kBlack)
        fresid.SetLineWidth(1)
        
        ### now build the fitsigs and fitbkgs histograms
        for i in range(numx):
            for nc, C in enumerate(fitchans):
                
                # signals
                #---------------------------------------------------------------
                for skey in sigkeys:
                    
                    ### init histograms for global signals (pmt, lsveto, steel, etc.)
                    for gmckey in globalmc:
                        if gmckey in skey:
                            bits = skey.split('-')
                            fgkey = bits[1]+'-'+bits[2]
                            if fgkey not in fglobsigkeys:
                                fglobsigkeys.append(fgkey)
                                fitglobsigs[fgkey] = {}
                                fglob = TH1F(fgkey, fgkey, fitbins, 0, fitbins)
                                fitglobsigs[fgkey]['hist'] = fglob
                    
                    # pad out the blank extra crystals
                    if i+1 in justthese and 'x'+str(i+1) not in skey:
                        if '-c'+C in skey and '-e0' in skey:
                            fskey = skey.split('-c'+C+'-e0')[0]
                            fitsigs = addHistKey(fitsigs, fskey)
                            fitsigs[fskey]['hist'] = extendHist(fitsigs[fskey]['hist'],
                                                fitPrep(sigs, skey, fitranges[i][C+'0'], 1))
                        if '-c'+C in skey and '-e1' in skey:
                            fskey = skey.split('-c'+C+'-e1')[0]
                            fitsigs = addHistKey(fitsigs, fskey)
                            fitsigs[fskey]['hist'] = extendHist(fitsigs[fskey]['hist'],
                                                fitPrep(sigs, skey, fitranges[i][C+'1'], 1))
                    
                    # low Energy
                    if 'x'+str(i+1) in skey and '-c'+C in skey and '-e0' in skey:
                        fskey = skey.split('-c'+C+'-e0')[0]
                        fitsigs = addHistKey(fitsigs, fskey)
                        fitsigs[fskey]['hist'] = extendHist(fitsigs[fskey]['hist'],
                                                fitPrep(sigs, skey, fitranges[i][C+'0']))
                    
                    # high Energy
                    if 'x'+str(i+1) in skey and '-c'+C in skey and '-e1' in skey:
                        fskey = skey.split('-c'+C+'-e1')[0]
                        fitsigs = addHistKey(fitsigs, fskey)
                        fitsigs[fskey]['hist'] = extendHist(fitsigs[fskey]['hist'],
                                                fitPrep(sigs, skey, fitranges[i][C+'1']))

                        
                # backgrounds
                #---------------------------------------------------------------
                for bkey in bakkeys:
                    
                    ### init histograms for global backgrounds (pmt, lsveto, steel, etc.)
                    for gmckey in globalmc:
                        if gmckey in bkey:
                            bits = bkey.split('-')
                            fgbkey = bits[1]+'-'+bits[2]
                            if fgbkey not in fglobbkgkeys:
                                fglobbkgkeys.append(fgbkey)
                                fitglobbkgs[fgbkey] = {}
                                fglob = TH1F(fgbkey, fgbkey, fitbins, 0, fitbins)
                                fitglobbkgs[fgbkey]['hist'] = fglob
                    
                    # pad out the blank extra crystals
                    if i+1 in justthese and 'x'+str(i+1) not in bkey:
                        if '-c'+C in bkey and '-e0' in bkey:
                            fbkey = bkey.split('-c'+C+'-e0')[0]
                            fitbkgs = addHistKey(fitbkgs, fbkey)
                            fitbkgs[fbkey]['hist'] = extendHist(fitbkgs[fbkey]['hist'],
                                                    fitPrep(bkgs, bkey, fitranges[i][C+'0'], 1))
                        if '-c'+C in bkey and '-e1' in bkey:
                            fbkey = bkey.split('-c'+C+'-e1')[0]
                            fitbkgs = addHistKey(fitbkgs, fbkey)
                            fitbkgs[fbkey]['hist'] = extendHist(fitbkgs[fbkey]['hist'],
                                                    fitPrep(bkgs, bkey, fitranges[i][C+'1'], 1))
                    
                    # low Energy
                    if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e0' in bkey:
                        fbkey = bkey.split('-c'+C+'-e0')[0]
                        fitbkgs = addHistKey(fitbkgs, fbkey)
                        fitbkgs[fbkey]['hist'] = extendHist(fitbkgs[fbkey]['hist'],
                                                fitPrep(bkgs, bkey, fitranges[i][C+'0']))
                    # high Energy
                    if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e1' in bkey:
                        fbkey = bkey.split('-c'+C+'-e1')[0]
                        fitbkgs = addHistKey(fitbkgs, fbkey)
                        fitbkgs[fbkey]['hist'] = extendHist(fitbkgs[fbkey]['hist'],
                                                fitPrep(bkgs, bkey, fitranges[i][C+'1']))
            

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
            if debug: print 'DEBUG: zero events so deleting sig key', key
            del fitsigs[key]
        fsigkeys.sort()

        # (2) sort the global sig keys
        fglobsigkeys.sort()

        # (3) build global sigs
        for fskey in fsigkeys:
            for fgskey in fglobsigkeys:
                if fgskey in fskey:
                    fitglobsigs[fgskey]['hist'].Add(fitsigs[fskey]['hist'])
        
        # (4) remove global sigs from fit sigs
        L = len(fsigkeys)-1
        for k, fskey in enumerate(reversed(fsigkeys)):
            for gmckey in globalmc:
                ### make gmckey unique (steel vs innersteel) 2018-10-23
                gmckey = '-'+gmckey+'-'
                if gmckey in fskey:
                    if debug: print 'DEBUG: remove global key from fit sigs', fsigkeys[L-k]
                    del fsigkeys[L-k]

        # (5) now delete the empty sig histos
        delete = []
        fglobsigkeys = []
        for fgskey in fitglobsigs:
            if fitglobsigs[fgskey]['hist'].Integral() > 0:
                fglobsigkeys.append(fgskey)
            else: delete.append(fgskey)
        for key in delete:
            if debug: print 'DEBUG: zero events so deleting global sigs key', key
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
            if debug: print 'DEBUG: zero events so deleting bkg key', key
            del fitbkgs[key]
        fbkgkeys.sort()

        # (2) sort the global bkg keys
        fglobbkgkeys.sort()
        
        # (3) build global bkgs
        for fbkey in fbkgkeys:
            for fgbkey in fglobbkgkeys:
                if fgbkey in fbkey:
                    fitglobbkgs[fgbkey]['hist'].Add(fitbkgs[fbkey]['hist'])
        
        # (4) remove global bkgs from fit bkgs
        L = len(fbkgkeys)-1
        for k, fbkey in enumerate(reversed(fbkgkeys)):
            for gmckey in globalmc:
                ### make gmckey unique (steel vs innersteel) 2018-10-23
                gmckey = '-'+gmckey+'-'
                if gmckey in fbkey:
                    if debug: print 'DEBUG: remove global key from fit bkgs', fbkgkeys[L-k]
                    del fbkgkeys[L-k]

        # (5) now delete the empty bkg histos
        delete = []
        fglobbkgkeys = []
        for fgbkey in fitglobbkgs:
            if fitglobbkgs[fgbkey]['hist'].Integral() > 0:
                fglobbkgkeys.append(fgbkey)
            else: delete.append(fgbkey)
        for key in delete:
            if debug: print 'DEBUG: zero events so deleting global bkgs key', key
            del fitglobbkgs[key]
        fglobbkgkeys.sort()

        ### --------------------------------------------------------------------
        ###         SUBTRACT BACKGROUNDS FROM DATA
        ### --------------------------------------------------------------------

        # (1) subtract normal backgrounds
        for key in fbkgkeys:
            fitdata.Add(fitbkgs[key]['hist'], -1)

        # (2) subtract global backgrounds
        for key in fglobbkgkeys:
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
        
        ### data integral to normalize signal to
        dat_int = fitdata.Integral()
                
        ### set up the fitting object for TFractionFitter
        for i in range(numx):
            fitresults[str(i)] = []
            for fskey in fsigkeys:
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
                        print '\nERROR: No events for --> ',fskey
                        print   '       Remove it from the fit? \n'
                        sys.exit()
                    
                    for C in allchans:
                        for E in range(2):
                            E=str(E)
                            
                            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            
                            newkeys = []
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

                            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            
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
                                        print 'ERROR: do not know what to do with useBounds =',useBounds
                                        sys.exit()
                    
                    bounds.append(sigs[newkeys[0]]['info']['newfbnd'])
                    #---------------------------------------------------------------------
                    
                    ### set errors to zero?
                    fitsigs[fskey]['hist'] = zeroBinError(fitsigs[fskey]['hist'])

            
            fitresults[str(i)].append('Crystal-'+str(i+1)+' fit results')
            fitresults[str(i)].append('runtime = '+str(round(runtimes[i]/60./60./24., 2))+' days')
            if note: fitresults[str(i)].append('note = '+note)
            fitresults[str(i)].append('script = '+script)
            fitresults[str(i)].append('root = '+vroot)
            fitresults[str(i)].append('channels fit = '+fitchans)
            fitresults[str(i)].append('globals = '+str(globalmc))
            fitresults[str(i)].append('use others = '+str(others))
            #fitresults[str(i)].append('hist extend = '+str(extend))
            fitresults[str(i)].append('norm to dru = '+str(dru))
            fitresults[str(i)].append('energy shift in bins = '+str(binShift[i]))
            fitresults[str(i)].append('shifting these = '+str(shiftWhat))
            fitresults[str(i)].append('fitter step size = '+str(stepSize))
            for key in ['S0', 'S1', 'M0', 'M1']:
                fitresults[str(i)].append(key+' fit range = '\
                                          +str(fitranges[i][key][1])\
                                          +' - '\
                                          +str(fitranges[i][key][2])\
                                          +' keV  -->  rebin = '\
                                          +str(fitranges[i][key][0]))
        
        
        ### do the same to the global lsveto
        boundskeys = []
        for fgkey in fglobsigkeys:
            
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
                print '\nERROR: No events for --> ',fgkey
                print   '       Remove it from the fit!\n'
                sys.exit()
            
            for i in range(numx):
                X = str(i+1)
                for C in allchans:
                    for E in range(2):
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
                                    continue
                        
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
                                        print 'ERROR: do not know what to do with useBounds =',useBounds
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


        totalNumFits = len(fsigkeys)+len(fglobsigkeys)
        print 'INFO: total number of hists being fit =', totalNumFits,'\n\n'
        sigObj = TObjArray(totalNumFits)
        
        for fskey in fsigkeys:
            sigObj.append(fitsigs[fskey]['hist'])

        for fgkey in fglobsigkeys:
            sigObj.append(fitglobsigs[fgkey]['hist'])
        
        fit = TFractionFitter(fitdata, sigObj)

        fitter = fit.GetFitter()
        
        ### set fit bounds!!!
        for l in range(len(bounds)):
            #fit.Constrain(l+1, bounds[l][0], bounds[l][1])
            # ROOT6 FIX
            fit.Constrain(l, bounds[l][0], bounds[l][1])
            if stepSize != 0:
                fitter.Config().ParSettings(l).SetStepSize(stepSize)
            
        ### set the fit range
        #fit.SetRangeX(0, fmax*numx)
        fit.SetRangeX(0, fitbins)

        ### if you want to get initial fit values from fitter...
        #fitter = fit.GetFitter()
        #fit_values = fitter.Config().ParamsValues()
        #fit_values = np.asarray(fit_values)
        #print fit_values

        
        #=======================================================================
        #        MACHEN SIE DAS FIT!!!
        #=======================================================================
        fitStartTime = datetime.datetime.now()
        status = fit.Fit()
        #status = 0  # force the plotting for testing
        status = int(status)  # because 6.14.00 returns <ROOT.TFitResultPtr object>
        fitStopTime = datetime.datetime.now()
        fitTime = int((fitStopTime-fitStartTime).total_seconds())
        #=======================================================================

        if status != 0:
            print '\n\n*******************  FIT FAILURE  *******************\n\n'
            print 'Time to fail = '+str(fitTime)+' sec \n'
            sys.exit()
        
        print '\n\n*******************  SUCCESSFUL FIT  *******************\n\n'
        
        
        chi2 = fit.GetChisquare()
        ndf  = fit.GetNDF()
        pval = fit.GetProb()

        #print 'chi2 =', chi2
        #print 'ndf =', ndf
        #print 'pval =', pval
        
        ### get non-zero ndf
        NDF=0
        for n in range(fitbins):
            if fitdata.GetBinContent(n) > 0:
                NDF+=1
        
        #print 'DEBUG: fit chi2/ndf =',chi2,ndf
        #print 'DEBUG: new ndf =',NDF
        ndf=NDF
        
        fitchi2ndf = (chi2/ndf)
        
        count = 0
        
        for i in range(numx):
            
            fitresults[str(i)].append('total number of hists being fit = '+str(totalNumFits))
            #fitresults[str(i)].append('returned fit status = '+str(status))
            fitresults[str(i)].append('chi2/ndf = %.3g/%s = %.3g'%(chi2, ndf, chi2/ndf))
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
                    
                    for C in allchans:
                        for E in range(2):
                            E = str(E)
                            
                            newkey = fskey+'-c'+C+'-e'+E
                            try:
                                ### save the raw scaling factor from the fit
                                sigs[newkey]['hist'].Scale(fscale)
                                
                                ### save converted scaling factor
                                sigs[newkey]['fitscale'] = sigs[newkey]['fitscale'] * fscale
                            except:
                                continue

                            try:
                                ### set error as a percent of the scaling factor
                                sigs[newkey]['fiterror'] = ferror/fscale
                            except:
                                sigs[newkey]['fiterror'] = 1.
                            
                                
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
            
            #newkey=0
            for i in range(numx):
                X = str(i+1)
                for C in allchans:
                    for E in range(2):
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

        
        ### scale the signals to mBq/kg
        if dru: sigs = scaleSigs410(sigkeys, sigs)
        else:   sigs = scaleSigs410(sigkeys, sigs, runtimes)
        
        ### print the fit activities
        for i in range(numx):
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    finit=1
                    for C in allchans:
                        for E in range(2):
                            if finit:
                                E = str(E)
                                newkey = fskey+'-c'+C+'-e'+E
                                try: test = sigs[newkey]['hist']
                                except: continue
                                
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
                                        '%42s = %.2e +/- %.2e mBq  (%.2e, %.2e) %s'
                                        %('fit '+fskey, fitacti, fiterro, lobnd, hibnd, limit))
                                        # show error as %
                                        # oh this messes up my pretty plotting of fit activities
                                        #'%42s = %.2e mBq (%.1f%%) %s'
                                        #%('fit '+fskey, fitacti, 100*fiterro/fitacti, limit))
                                
                                else:
                                    fitresults[str(i)].append(
                                        # show the actual error
                                        '%42s = %.2e +/- %.2e mBq  (%.2e, %.2e)'
                                        %('fit '+fskey, fitacti, fiterro, lobnd, hibnd))
                                        # show error as %
                                        # oh this messes up my pretty plotting of fit activities
                                        #'%42s = %.2e mBq (%.1f%%)'
                                        #%('fit '+fskey, fitacti, 100*fiterro/fitacti))
                                
                                finit = 0

        ### do the same for the globals
        for i in range(numx):
            for fgkey in fglobsigkeys:
                X=str(i+1)
                finit=1
                for C in allchans:
                    for E in range(2):
                        E=str(E)
                        if finit:
                            
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
                                        # show the actual error
                                        '%42s = %.2e +/- %.2e mBq  (%.2e, %.2e) %s'
                                        %('fit '+'x'+X+'-'+fgkey, fitacti, fiterro, lobnd, hibnd, limit))
                                        # show error as %
                                        # oh this messes up my pretty plotting of fit activities
                                        #'%42s = %.2e mBq (%.1f%%) %s'
                                        #%('fit '+'x'+X+'-'+fgkey, fitacti, 100*fiterro/fitacti, limit))
                                
                                else:
                                    fitresults[str(i)].append(
                                        # show the actual error
                                        '%42s = %.2e +/- %.2e mBq  (%.2e, %.2e)'
                                        %('fit '+'x'+X+'-'+fgkey, fitacti, fiterro, lobnd, hibnd))
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
        save += '_globals-'+globstr
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
        if note: save += '_'+note
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
            resultsfile = os.path.join(plotdir, save+'_fit-results.txt')
            outfile = open(resultsfile, 'w')
            for key in resultskeys:
                if int(key)+1 in justthese:
                    for line in fitresults[key]:
                        outfile.write(line+'\n')
            outfile.close()
            
            ### save histograms to a rootfile
            rootoutfile = TFile(plotdir+"/histograms.root", "RECREATE")
            for key in sigkeys:
                sigs[key]['hist'].Write(key)
            for key in bakkeys:
                bkgs[key]['hist'].Write(key)
            for key in data:
                data[key]['hist'].Write(key)
            rootoutfile.Write()
            rootoutfile.Close()
        #-------------------------------------------------------------
        
        
        #=================================================================
        #=================================================================
        #
        #   Plot the multi-chan fit results
        #
        #=================================================================
        #=================================================================

        fcanv = TCanvas('fcanv', 'fcanv', 0, 0, 1400, 900)
        #fcanv.Divide(4,2)
        fcanv.Divide(5,2)

        ftoppad = [0 for x in range(numx)]
        fbotpad = [0 for x in range(numx)]

        sepfitdata = [0 for x in range(numx)]
        sepfitsigs = [0 for x in range(numx)]
        sepfitglobsigs = [0 for x in range(numx)]

        ### have to init the histos
        ### there must be an easier way...
        sepfittotal = [0 for x in range(numx)]
        sepfitresid = [0 for x in range(numx)]
        for x in range(numx):
            sepfittotal[x] = TH1F('sepFitTotal', 'sepFitTotal', fitbins, 0, fitbins)
            sepfittotal[x].SetLineColor(kGray+1)
            sepfittotal[x].SetMarkerColor(kGray+1)
            sepfittotal[x].SetLineWidth(1)
            
            sepfitresid[x] = TH1F('sepFitResid', 'sepFitResid', fitbins, 0, fitbins)
            sepfitresid[x].SetLineColor(kBlack)
            sepfitresid[x].SetMarkerColor(kBlack)
            sepfitresid[x].SetLineWidth(1)
            
        flegs  = [0 for x in range(numx)]
        flegs2 = [0 for x in range(numx)]
        fzeros = [0 for x in range(numx)]

        
        for i in range(numx):
            fcanv.cd(i+1)
            
            gStyle.SetPadTopMargin    (0.07)
            gStyle.SetPadBottomMargin (0.11)
            gStyle.SetPadLeftMargin   (0.12)
            gStyle.SetPadRightMargin  (0.02)
            
            font = 63
            size = 13
            yoff = 4.2
            xoff = 8
            
            
            fraction = 0.3
            pad1 = TPad('pad1','pad1',0,fraction,1,1)
            ftoppad[i] = pad1
            pad2 = TPad('pad2','pad2',0,0,1,fraction)
            fbotpad[i] = pad2

            ftoppad[i].SetBottomMargin(0.01)
            ftoppad[i].SetBorderMode(0)
            ftoppad[i].SetLogy()
            fbotpad[i].SetTopMargin(0.05)
            fbotpad[i].SetBottomMargin(0.3)
            fbotpad[i].SetBorderMode(0)
            fbotpad[i].SetLogy(1)
            ftoppad[i].Draw()
            fbotpad[i].Draw()
            ftoppad[i].cd()

            #sepfitdata[i] = deepcopy(fitdata)
            sepfitdata[i] = TH1F('sepFitData', 'sepFitData', fitbins, 0, fitbins)
            if i+1 in justthese:
                sepfitdata[i] = deepcopy(fitdata)
            newFitTitle = str('Crystal-'+str(i+1)+'   '+'Fit-chans-'+fitchans)
            sepfitdata[i].SetTitle(newFitTitle)
            sepfitdata[i].SetLineColor(kBlack)
            sepfitdata[i].SetMarkerColor(kBlack)
            sepfitdata[i].SetLineWidth(1)

            if dru: sepfitdata[i].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
            else: sepfitdata[i].GetYaxis().SetTitle('arb. counts')

            sepfitdata[i].GetYaxis().SetTitleFont(font)
            sepfitdata[i].GetYaxis().SetTitleSize(size)
            sepfitdata[i].GetYaxis().SetTitleOffset(yoff)
            sepfitdata[i].GetYaxis().SetLabelFont(font)
            sepfitdata[i].GetYaxis().SetLabelSize(size)
            sepfitdata[i].GetYaxis().SetLabelOffset(0.01)
            #sepfitdata[i].GetXaxis().SetTitle('Energy (keV)')
            #sepfitdata[i].GetXaxis().SetLabelFont(font)
            #sepfitdata[i].GetXaxis().SetLabelSize(size)

            if dru:
                sepfitdata[i].SetAxisRange(2e-3, 2e3, 'y')
                if i == 8:
                    sepfitdata[i].SetAxisRange(2e-5, 2e1, 'y')

            sepfitdata[i].SetAxisRange(xstalbins[i][0], xstalbins[i][1], 'x')
            sepfitdata[i].Draw()
            
            
            sepfitsigs[i] = deepcopy(fitsigs)
            sepfitglobsigs[i] = deepcopy(fitglobsigs)
            sepfittotal[i] = deepcopy(ftotal)
            Nfsigs = 0
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    Nfsigs += 1
                    
                    ### find the unique name for color and set color
                    cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
                    sepfitsigs[i][fskey]['hist'].SetMarkerColor(uniqColor[cname])
                    sepfitsigs[i][fskey]['hist'].SetLineColor(uniqColor[cname])
                                        
                    ### draw the sigs
                    if i+1 in justthese:
                        sepfitsigs[i][fskey]['hist'].Draw(dopt)

                    ### add MC to total MC hist
                    #print '-3-', ftotal.GetNbinsX(), sepfitsigs[i][fskey]['hist'].GetNbinsX(), fskey
                    sepfittotal[i].Add(sepfitsigs[i][fskey]['hist'])
                
            for fgkey in fglobsigkeys:
                Nfsigs += 1
                
                ### find the unique name for color and set color
                cname = fgkey
                sepfitglobsigs[i][fgkey]['hist'].SetMarkerColor(uniqColor[cname])
                sepfitglobsigs[i][fgkey]['hist'].SetLineColor(uniqColor[cname])
                
                ### draw the sigs
                if i+1 in justthese:
                    sepfitglobsigs[i][fgkey]['hist'].Draw(dopt)

                ### add MC to total MC hist
                #print '-4-', sepfittotal[i].GetNbinsX(), sepfitglobsigs[i][fgkey]['hist'].GetNbinsX(), fgkey
                sepfittotal[i].Add(sepfitglobsigs[i][fgkey]['hist'])
                
                
            ### select the right range to plot
            ### ===============================================
            if i+1 in justthese:
                sepfittotal[i].SetAxisRange(xstalbins[i][0], xstalbins[i][1], 'x')

            
            ### build legend after getting the numbers of signals
            Nfs = Nfsigs+2
            flnc = 1
            if Nfs >  6: flnc = 2
            if Nfs > 12: flnc = 3
            #if Nfs > 18: flnc = 4
            fxlegstop  = 0.94
            fxlegstart = fxlegstop-(0.2*flnc)
            #fylegstop  = 0.89
            fylegstop  = 0.91
            fylegstart = fylegstop-(0.04*6)
            
            fleg = TLegend(fxlegstart, fylegstart, fxlegstop, fylegstop)
            flegs[i] = (fleg)
            flegs[i].SetNColumns(flnc)
            flegs[i].SetFillColor(0)
            flegs[i].SetBorderSize(0)
            legopt = 'LPE'
            
            flegs[i].AddEntry(sepfitdata[i], 'data - bkgs', legopt)
            
            ### add legend entries in order
            for name in uniqAll:
                for fskey in fsigkeys:
                    if name in fskey and 'x'+str(i+1) in fskey:
                        flegs[i].AddEntry(sepfitsigs[i][fskey]['hist'], fskey, legopt)
                # and for globals
                for fgkey in fglobsigkeys:
                    if name in fgkey:
                        flegs[i].AddEntry(sepfitglobsigs[i][fgkey]['hist'], fgkey, legopt)
            
            
            ### draw to total and legend
            if i+1 in justthese:
                sepfittotal[i].Draw(dopt)
                flegs[i].AddEntry(sepfittotal[i],
                                  'Fit Total (chi2/ndf = '+str(round(fitchi2ndf,2))+')',
                                  legopt)
                flegs[i].Draw('same')
            
            
            ### plot the fit residuals
            #-------------------------------------------------------------
            fbotpad[i].cd()
            leg = TLegend(0.72, 0.78, 0.94, 0.94)
            flegs2[i] = (leg)
            flegs2[i].SetFillColor(0)
            flegs2[i].SetBorderSize(0)
            legopt = 'LPE'

            if i+1 in justthese:
                sepfitresid[i].Divide(sepfitdata[i], sepfittotal[i])
                        
            sepfitresid[i].SetTitle('')
            #sepfitresid[i].SetXTitle('Energy (keVee)')
            sepfitresid[i].SetXTitle('fit bins')
            sepfitresid[i].GetXaxis().SetTitleFont(font)
            sepfitresid[i].GetXaxis().SetTitleSize(size)
            sepfitresid[i].GetXaxis().SetTitleOffset(xoff)
            sepfitresid[i].GetXaxis().SetLabelFont(font)
            sepfitresid[i].GetXaxis().SetLabelSize(size)
            sepfitresid[i].GetXaxis().SetLabelOffset(0.03)
            #sepfitresid[i].SetYTitle('counts / keV')
            sepfitresid[i].SetYTitle('data / MC')
            #sepfitresid[i].SetYTitle('MC / data')
            sepfitresid[i].GetYaxis().SetTitleFont(font)
            sepfitresid[i].GetYaxis().SetTitleSize(size)
            sepfitresid[i].GetYaxis().SetTitleOffset(yoff)
            sepfitresid[i].GetYaxis().SetLabelFont(font)
            sepfitresid[i].GetYaxis().SetLabelSize(size)
            sepfitresid[i].GetYaxis().SetLabelOffset(0.01)
            # '5' secondary and '05' primary
            sepfitresid[i].GetYaxis().SetNdivisions(505)
            
            sepfitresid[i].SetAxisRange(0.1, 10,'y')
            sepfitresid[i].SetAxisRange(xstalbins[i][0], xstalbins[i][1], 'x')

            if linres:
                fbotpad[i].SetLogy(0)
                sepfitresid[i].SetAxisRange(lrs[0], lrs[1], 'y')

            sepfitresid[i].Draw()
            
            #zero = TLine(xstalbins[i][0], 1, xstalbins[i][1], 1)
            zero = TLine(xstalbins[i][0], 1, xstalbins[i][1]+1, 1)
            fzeros[i] = (zero)
            fzeros[i].SetLineColor(kRed)
            fzeros[i].SetLineWidth(1)
            fzeros[i].Draw()
            
            flegs2[i].AddEntry(sepfitresid[i], 'data / MC', legopt)
            #flegs2[i].Draw()
            #-------------------------------------------------------------

            fcanv.Update()
            
            if indi and i+1 in justthese:
                #fisave = 'xstal-'+str(i+1)
                try:
                    ftpad=ftoppad[i].Clone()
                    fbpad=fbotpad[i].Clone()
                    sepFitPlot = TCanvas('ican-fit-'+str(fitchans)+str(i),
                                         'ican-fit-'+str(fitchans)+str(i),
                                         0, 0, 1400, 900)
                    
                    ftpad.Draw()
                    fbpad.Draw()
                    sepFitPlot.Update()
                    fisave  = ''
                    fisave += 'x'+str(i+1)
                    fisave += '_fit'
                    fisave += '_globals-'+globstr
                    #fisave += '_loEfRS'+str(loEfitRebinScale)
                    #fisave += '_hiEfRS'+str(hiEfitRebinScale)
                    fisave += '_dru'+str(dru)
                    fisave += '_cs'+str(fitchans)
                    #fisave += '_ext'+str(extend)
                    fisave += '_oth'+str(others)
                    if note: fisave += '_'+str(note)
                    fisave += '_'+str(V)
                
                    sepFitPlot.Print(plotdir+'/'+fisave+'.png')
                except:
                    print 'WARNING: could not make plot for xstal-'+str(i+1)
                    continue
            
        fcanv.Update()
        
        msave  = ''
        msave += str(runtag)
        msave += '_globals-'+globstr
        #msave += '_loEfit-'+str(int(fLoE[0]))+'-'+str(int(fLoE[1]))
        #msave += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        #msave += '_loEfitRebin-'+str(loEfitRebin)
        #msave += '_hiEfitRebin-'+str(hiEfitRebin)
        #msave += '_loEfitRebinScale'+str(loEfitRebinScale)
        #msave += '_hiEfitRebinScale'+str(hiEfitRebinScale)
        #msave += '_useBounds'+str(useBounds)
        #msave += '_mcsumw2'+str(mcsumw2)
        #msave += '_datsumw2'+str(datsumw2)
        #msave += '_dru'+str(dru)
        #msave += '_loEplotRebin-'+str(loEplotRebin)
        #msave += '_hiEplotRebin-'+str(hiEplotRebin)
        #msave += '_reuse'+str(reuse)
        #msave += '_chans'+str(fitchans)
        #msave += '_extend'+str(extend)
        #msave += '_others'+str(others)
        if note: msave += '_'+note
        msave += '_'+V
        
        #fcanv.Print(plotdir+'/'+save+'.png')
        fcanv.Print(plotdir+'/'+msave+'.png')
        
        
        #=============================================================
        #=============================================================
        #
        #   plot the MEGA combined histograms for fun!
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
        mtoppad.SetLogy()
        
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
        if dru:
            fitdata.SetAxisRange(2e-3, 2e3, 'y')
            #if i == 8:
            #    fitdata.SetAxisRange(2e-5, 2e1, 'y')
        
        fitdata.SetAxisRange(0, fitbins, 'x')
        
        
        Nfsigs=0
        #for i in range(numx):
        for fskey in fsigkeys:
            #if 'x'+str(i+1) in fskey:

            Nfsigs+=1

            # find the unique name for color and set color
            cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
            fitsigs[fskey]['hist'].SetMarkerColor(uniqColor[cname])
            fitsigs[fskey]['hist'].SetLineColor(uniqColor[cname])

            ### draw the sigs
            fitsigs[fskey]['hist'].Draw(dopt)

            ### add MC to total MC hist
            ftotal.Add(fitsigs[fskey]['hist'])

        for fgkey in fglobsigkeys:
            #if i==0:
            Nfsigs+=1

            # find the unique name for color and set color
            #cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
            cname = fgkey
            fitglobsigs[fgkey]['hist'].SetMarkerColor(uniqColor[cname])
            fitglobsigs[fgkey]['hist'].SetLineColor(uniqColor[cname])

            ### draw the sigs
            fitglobsigs[fgkey]['hist'].Draw(dopt)

            # FIX-ME
            # why if i==0???
            #if i==0:
                ### add MC to total MC hist
            ftotal.Add(fitglobsigs[fgkey]['hist'])

                
        ### select the right range
        ftotal.SetAxisRange(0, fitbins, 'x')
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
        for name in uniqAll:
            for fskey in fsigkeys:
                if name in fskey:
                    mleg.AddEntry(fitsigs[fskey]['hist'], fskey, legopt)
            # and for globals
            for fgkey in fglobsigkeys:
                if name in fgkey:
                    mleg.AddEntry(fitglobsigs[fgkey]['hist'], fgkey, legopt)

        #ftotal.Draw('same')

        #if i in wasFit:
        # returned from fit
        mleg.AddEntry(ftotal, 'Fit Total (chi2/ndf = '+str(round(fitchi2ndf,2))+')', legopt)
        # calc by Chi2TestX()
        #flegs[i].AddEntry(ftotal, 'Fit Total (chi2/ndf = '+str(round(fitchi2ndfv2,2))+')', legopt)
        #print 'INFO: Fit total MC from Chi2TestX chi2/ndf = '+str(round(fitchi2ndfv2,2))

        mleg.Draw('same')


        ### plot the fit residuals
        #-------------------------------------------------------------
        mbotpad.cd()
        mleg2 = TLegend(0.72, 0.78, 0.94, 0.94)
        #flegs2.append(leg)
        mleg2.SetFillColor(0)
        mleg2.SetBorderSize(0)
        legopt = 'LPE'
        
        fresid.Divide(fitdata, ftotal)
                
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
        fresid.Draw()

        ### set my reference line to '1'
        #zero = TLine(fmin, 1, fmax, 1)
        #fzeros.append(zero)
        #fzeros[i].SetLineColor(kRed)
        #fzeros[i].SetLineWidth(1)

        mzero = TLine(0, 1, fitbins, 1)
        #fzeros.SetAxisRange(i*fmax, (i+1)*fmax, 'x')
        #fzeros.append(zero)
        mzero.SetLineColor(kRed)
        mzero.SetLineWidth(1)
        mzero.Draw()
        
        mleg2.AddEntry(fresid,'data / MC',legopt)
        #flegs2[i].Draw()
        #-------------------------------------------------------------
        
        mcanv.Update()

        
        #mcanv.Print(plotdir+'/'+msave+'.png')
        mcanv.Print(plotdir+'/'+save+'.png')
        #=============================================================
        #=============================================================

        #raw_input('[Enter] to quit \n')
        
    ### end of fitting bit if you have signals

    
    ### create the updated backgrounds file
    updateBkgsFile300(justthese, os.path.join(plotdir, mcfile), resultsfile, plotdir)
    
    
    # plot the lo and hi energy histograms for all channels
    #=================================================================
    
    # number of energy ranges (lo, hi)
    numE = 2
    # number of channels (single and/or multi-hit)
    numC = len(pltchans)
    
    canvs  = [[0 for x in range(numE)] for x in range(numC)]
    
    ### for separate plots
    sepPlots = [[[0 for x in range(numx)] for x in range(numE)] for x in range(numC)]
    sepTopPad = [[[0 for x in range(numx)] for x in range(numE)] for x in range(numC)]
    sepBotPad = [[[0 for x in range(numx)] for x in range(numE)] for x in range(numC)]
    
    ### seperate memory space for the pads is key!!!!
    toppad = [[0 for x in range(numE)] for x in range(numC)]
    botpad = [[0 for x in range(numE)] for x in range(numC)]
    
    legs   = [[0 for x in range(numE)] for x in range(numC)]
    legs2  = [[0 for x in range(numE)] for x in range(numC)]
    zeros  = [[0 for x in range(numE)] for x in range(numC)]
    
    total  = [[0 for x in range(numE)] for x in range(numC)]
    resid  = [[0 for x in range(numE)] for x in range(numC)]
    
    gbkgs  = [[[{} for x in range(numx)] for x in range(numE)] for x in range(numC)]
    gsigs  = [[[{} for x in range(numx)] for x in range(numE)] for x in range(numC)]
        
    plotRebin = 1
    for C, chan in enumerate(pltchans): 
    
        for E in range(numE):

            if E: plotRebin = hiEplotRebin
            else: plotRebin = loEplotRebin
            
            # have the plotting be seperated out from the 8 crystal loop
            canvs[C][E] = TCanvas('canv'+chan+str(E),
                                  'canv'+chan+str(E),
                                  0, 0, 1400, 900)

            #canvs[C][E].Divide(4,2)
            canvs[C][E].Divide(5,2)
            
            gStyle.SetPadTopMargin    (0.07)
            gStyle.SetPadBottomMargin (0.11)
            gStyle.SetPadRightMargin  (0.02)
            #gStyle.SetPadLeftMargin   (0.12)
            #font = 63
            #size = 13
            #yoff = 4.2
            #xoff = 8
            
            gStyle.SetPadLeftMargin(0.10)
            font = 43
            size = 16
            yoff = 2.5
            xoff = 7
            
            toppad[C][E] = []
            botpad[C][E] = []

            legs[C][E]   = []
            legs2[C][E]  = []
            zeros[C][E]  = []

            total[C][E]  = makeTotal100(chan, E, params[E])
            resid[C][E]  = makeResid100(chan, E, params[E])
            
            for i in range(numx):
                
                canvs[C][E].cd(i+1)
                
                fraction = 0.3
                pad1 = TPad('pad1'+chan+str(E),'pad1'+chan+str(E),0,fraction,1,1)
                toppad[C][E].append(pad1)
                pad2 = TPad('pad2'+chan+str(E),'pad2'+chan+str(E),0,0,1,fraction)
                botpad[C][E].append(pad2)
                
                toppad[C][E][i].SetBottomMargin(0.01)
                toppad[C][E][i].SetBorderMode(0)
                
                toppad[C][E][i].SetLogy(1)
                
                botpad[C][E][i].SetTopMargin(0.05)
                botpad[C][E][i].SetBottomMargin(0.3)
                botpad[C][E][i].SetBorderMode(0)
                
                botpad[C][E][i].SetLogy(1)
                #if linres: botpad[C][E][i].SetLogy(0)
                
                toppad[C][E][i].Draw()
                botpad[C][E][i].Draw()
                toppad[C][E][i].cd()

                if ingroups:
                    #leg = TLegend(0.35, 0.75, 0.94, 0.89)
                    leg = TLegend(0.56, 0.76, 0.98, 0.92)
                    lnc = 3
                else:
                    leg = TLegend(xlegstart, ylegstart, xlegstop, ylegstop)
                legs[C][E].append(leg)
                legs[C][E][i].SetFillColor(0)
                legs[C][E][i].SetBorderSize(0)
                legs[C][E][i].SetNColumns(lnc)
                legopt = 'LPE'

                if redtotal:
                    total[C][E][i].SetMarkerColor(kRed)
                    total[C][E][i].SetLineColor(kRed)
                    #total[C][E][i].SetLineWidth(4)
                    #total[C][E][i].SetMarkerSize(4)
                    
                total[C][E][i].Rebin(plotRebin)
                #total[C][E][i].Scale(1./float(plotRebin))
                #total[C][E][i].Sumw2()
                if dru and i!=8:
                    #total[C][E][i].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
                    if E: total[C][E][i].SetAxisRange(2e-3, 2e1, 'y')
                    else: total[C][E][i].SetAxisRange(2e-3, 3e2, 'y')
                    #if i==8: total[C][E][i].SetAxisRange(2e-4, 2e0, 'y')
                #else: total[C][E][i].GetYaxis().SetTitle('arb. counts')

                
                resid[C][E][i].Rebin(plotRebin)
                #resid[C][E][i].Scale(1./float(plotRebin))
                #resid[C][E][i].Sumw2()
                
                dkey = 0
                datasets = 0
                for key in datkeys:
                    if 'x'+str(i+1) in key and '-c'+chan in key and '-e'+str(E) in key:

                        dkey = key
                        datasets += 1
                        
                        if dru:
                            data[dkey]['hist'].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
                        else:
                            data[dkey]['hist'].GetYaxis().SetTitle('arb. counts')

                        newTitle = str('C'+str(i+1)+'   '
                                    +chanNames(chan)+'   '
                                    +energyNames(E))
                        """
                        newTitle = str('crystal-'+str(i+1)+'   '
                                    +energyNames(E)+'   '
                                    +chanNames(chan)+'   '
                                    +detailNames(i)+'   '
                                    +str(cmass(i))+'kg')
                        """
                        
                        data[dkey]['hist'].SetTitle(newTitle)
                        data[dkey]['hist'].SetTitleFont(font)
                        data[dkey]['hist'].SetTitleSize(size)
                        
                        data[dkey]['hist'].Rebin(plotRebin)
                        data[dkey]['hist'].Scale(1./float(plotRebin))
                        """
                        total[C][E][i].Rebin(plotRebin)
                        total[C][E][i].Scale(1./float(plotRebin))
                        #total[C][E][i].Sumw2()
                        
                        resid[C][E][i].Rebin(plotRebin)
                        resid[C][E][i].Scale(1./float(plotRebin))
                        #resid[C][E][i].Sumw2()
                        """
                        #data[dkey]['hist'].SetMarkerStyle(8)
                        #data[dkey]['hist'].SetMarkerSize(.6)

                        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        #  Multiple data sets?
                        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        """
                        if datasets == 1:
                            data[dkey]['hist'].SetMarkerColor(kBlack)
                            data[dkey]['hist'].SetLineColor(kBlack)
                        if datasets == 2:
                            data[dkey]['hist'].SetMarkerColor(kBlue)
                            data[dkey]['hist'].SetLineColor(kBlue)
                        """
                        data[dkey]['hist'].SetMarkerColor(dcs[datasets-1])
                        data[dkey]['hist'].SetLineColor(dcs[datasets-1])
                        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        
                        data[dkey]['hist'].SetLineWidth(1)
                        
                        #if redtotal:
                            #data[dkey]['hist'].SetMarkerStyle(7)
                            #data[dkey]['hist'].SetMarkerSize(1)
                            
                        data[dkey]['hist'].GetYaxis().SetTitleFont(font)
                        data[dkey]['hist'].GetYaxis().SetTitleSize(size)
                        data[dkey]['hist'].GetYaxis().SetTitleOffset(yoff)
                        data[dkey]['hist'].GetYaxis().SetLabelFont(font)
                        data[dkey]['hist'].GetYaxis().SetLabelSize(size)
                        #data[dkey]['hist'].GetYaxis().SetLabelOffset(0.01)
                        #data[dkey]['hist'].GetXaxis().SetTitle('Energy (keV)')
                        #data[dkey]['hist'].GetXaxis().SetLabelFont(font)
                        #data[dkey]['hist'].GetXaxis().SetLabelSize(size)

                        if E: data[dkey]['hist'].SetAxisRange(hier[0], hier[1], 'x')
                        else: data[dkey]['hist'].SetAxisRange(loer[0], loer[1], 'x')

                        if dru and i!=8:
                            #data[dkey]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            #if E: data[dkey]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            if E: data[dkey]['hist'].SetAxisRange(hiYlo, 2e1, 'y')
                            else: data[dkey]['hist'].SetAxisRange(2e-3, 3e2, 'y')

                        ### LSveto plotting
                        if E and i==8:
                            data[dkey]['hist'].SetAxisRange(0, lsHiE, 'x')
                            data[dkey]['hist'].SetAxisRange(1e-5, 0.5, 'y')
                        
                        data[dkey]['hist'].Draw('same')
                        days = round(data[dkey]['runtime']/86400., 2)
                        
                        legs[C][E][i].AddEntry(data[dkey]['hist'],
                                               data[dkey]['info']['tag']+'  '+data[dkey]['info']['build'],
                                               legopt)
                        
                        #if ingroups: legs[C][E][i].AddEntry(data[dkey]['hist'], 'Data', legopt)
                        #else: legs[C][E][i].AddEntry(data[dkey]['hist'], dkey+' ('+str(days)+' days)', legopt)

                        
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
                        bkgs[key]['hist'].Rebin(plotRebin)
                        bkgs[key]['hist'].Scale(1./float(plotRebin))
                        #bkgs[key]['hist'].Sumw2()
                        
                        if dru:
                            if E: bkgs[key]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            else: bkgs[key]['hist'].SetAxisRange(2e-3, 3e2, 'y')
                        
                        if ingroups:
                            if bkgs[key]['info']['group'] == 'none':
                                try:
                                    gbkgs[C][E][i]['none'][key] = bkgs[key]['hist']
                                except:
                                    gbkgs[C][E][i]['none'] = {}
                                    gbkgs[C][E][i]['none'][key] = bkgs[key]['hist']
                            else:
                                try:
                                    gbkgs[C][E][i][bkgs[key]['info']['group']].Add(bkgs[key]['hist'])
                                except:
                                    gbkgs[C][E][i][bkgs[key]['info']['group']] = bkgs[key]['hist']
                        else:
                            bkgs[key]['hist'].Draw(dopt)

                        ### add MC to total MC hist
                        total[C][E][i].Add(bkgs[key]['hist'])
                        tcount += 1
                        
                        ### create the legend entry for MC
                        #legs[C][E][i].AddEntry(bkgs[key]['hist'], key, legopt)

                skey = 0
                for key in sigkeys:
                    if 'x'+str(i+1) in key and '-c'+chan in key and '-e'+str(E) in key:

                        skey = key
                        
                        #if i in wasFit:
                        # find the unique name for color and set color
                        cname = key.split('-')[1]+'-'+key.split('-')[2]
                        sigs[key]['hist'].SetMarkerColor(uniqColor[cname])
                        sigs[key]['hist'].SetLineColor(uniqColor[cname])

                        #if E and plotRebin:
                        sigs[key]['hist'].Rebin(plotRebin)
                        sigs[key]['hist'].Scale(1./float(plotRebin))
                        #sigs[key]['hist'].Sumw2()

                        ### set range
                        #sigs[key]['hist'].SetAxisRange(1,1000,'y')

                        if ingroups:
                            if sigs[key]['info']['group'] == 'none':
                                try:
                                    gbkgs[C][E][i]['none'][key] = sigs[key]['hist']
                                except:
                                    gbkgs[C][E][i]['none'] = {}
                                    gbkgs[C][E][i]['none'][key] = sigs[key]['hist']
                            else:
                                try:
                                    gbkgs[C][E][i][sigs[key]['info']['group']].Add(sigs[key]['hist'])
                                except:
                                    gbkgs[C][E][i][sigs[key]['info']['group']] = sigs[key]['hist']
                        else:
                            sigs[key]['hist'].Draw(dopt)

                        ### add MC to total MC hist
                        total[C][E][i].Add(sigs[key]['hist'])
                        tcount += 1

                        ### create the legend entry for MC
                        #legs[C][E][i].AddEntry(sigs[key]['hist'], key, legopt)

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
                            legs[C][E][i].AddEntry(gbkgs[C][E][i]['none'][key], key, legopt)
                else:
                    # add legend entries in order
                    #print uniqAll
                    for name in uniqAll:
                        for bkey in bakkeys:
                            #print bkey
                            if bkey == 'x'+str(i+1)+'-'+name+'-c'+chan+'-e'+str(E) \
                               or bkey == 'x'+str(i+1)+'-'+name+'-f'+str(i+1)+'-c'+chan+'-e'+str(E):
                                activ = '(%.2e) '%(bkgs[bkey]['info']['acti'])
                                legs[C][E][i].AddEntry(bkgs[bkey]['hist'], activ+bkey, legopt)
                        for skey in sigkeys:
                            if skey == 'x'+str(i+1)+'-'+name+'-c'+chan+'-e'+str(E) \
                               or skey == 'x'+str(i+1)+'-'+name+'-f'+str(i+1)+'-c'+chan+'-e'+str(E):
                                activ = '(%.2e) '%(sigs[skey]['info']['fitacti'])
                                legs[C][E][i].AddEntry(sigs[skey]['hist'], activ+skey, legopt)
                
                
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
                                total[C][E][i].GetBinError(n+1)/(float(plotRebin)/math.sqrt(2.)))
                
                
                ### chi2 test
                ### get the chi2 of the total mc compared to data
                #=============================================================================
                #-----------------------------------------------------------------------------
                #chi2  = ROOT.Double(-1)
                #ndf   = ROOT.Long(1)
                #igood = ROOT.Long(0)
                # ROOT6 FIX
                chi2  = ctypes.c_double(0.0)
                ndf   = ctypes.c_int(0)
                igood = ctypes.c_int(0)

                #print 'THIS ONE???'
                #print dkey, data[dkey]['hist'].GetNbinsX()
                #print total[C][E][i].GetNbinsX()
                
                if dkey and tcount:
                    if data[dkey]['hist'].GetEntries() > 0:
                        pval  = data[dkey]['hist'].Chi2TestX(total[C][E][i], chi2, ndf, igood, chiopt)
                        #print 'INFO:',dkey,'pval =',pval,'chi2 =',chi2,'ndf =',ndf,'igood =',igood

                        # ROOT6 FIX
                        chi2  = float(chi2.value)
                        ndf   = int(ndf.value)
                        igood = int(igood.value)
                        
                        #print 'INFO:',dkey,'total MC chi2/ndf =',round(chi2/data[dkey]['druScale']/ndf,2)
                        print 'INFO:',dkey,'total MC chi2/ndf =',round(chi2/ndf,2)
                #-----------------------------------------------------------------------------
                #=============================================================================



                #-----------------------------------------------------------------------------
                # still draw total even if nothing so hist box is created
                if showTotal:
                    total[C][E][i].Draw(dopt)
                    if tcount:
                        total[C][E][i].SetLineWidth(1)
                        if ingroups:
                            legs[C][E][i].AddEntry(total[C][E][i], 'Total', legopt)
                        else:
                            total[C][E][i].SetLineColor(kGray+1)
                            legs[C][E][i].AddEntry(total[C][E][i], 'Total', legopt)
                            #legs[C][E][i].AddEntry(total[C][E][i],
                            #    'Total MC (chi2/ndf = '+str(round(chi2/ndf,2))+')', legopt)
                else:
                    # draw an empty total hist to preserve plot box layout?
                    if not tcount and not dkey: total[C][E][i].Draw()
                
                
                ### show the legends?
                if showlegs and (dkey or bkey):
                    legs[C][E][i].Draw('same')


                ### try to get the residuals in!
                #---------------------------------------------------------
                botpad[C][E][i].cd()
                leg = TLegend(0.72, 0.78, 0.94, 0.94)
                legs2[C][E].append(leg)
                legs2[C][E][i].SetFillColor(0)
                legs2[C][E][i].SetBorderSize(0)
                legopt = 'LPE'

                if tcount and dkey:
                    resid[C][E][i].Divide(data[dkey]['hist'], total[C][E][i])
                
                resid[C][E][i].SetTitle('')
                resid[C][E][i].SetXTitle('Energy (keVee)')
                resid[C][E][i].GetXaxis().SetTitleFont(font)
                resid[C][E][i].GetXaxis().SetTitleSize(size)
                resid[C][E][i].GetXaxis().SetLabelFont(font)
                resid[C][E][i].GetXaxis().SetLabelSize(size)
                resid[C][E][i].GetXaxis().SetLabelOffset(0.03)
                resid[C][E][i].GetXaxis().SetTitleOffset(xoff)
                #resid[C][E][i].SetYTitle('counts / keV')
                resid[C][E][i].SetYTitle('data / MC')
                #resid[C][E][i].SetYTitle('MC / data')
                resid[C][E][i].GetYaxis().SetTitleFont(font)
                resid[C][E][i].GetYaxis().SetTitleSize(size)
                resid[C][E][i].GetYaxis().SetTitleOffset(yoff)
                resid[C][E][i].GetYaxis().SetLabelFont(font)
                resid[C][E][i].GetYaxis().SetLabelSize(size)
                resid[C][E][i].GetYaxis().SetLabelOffset(0.01)
                # '5' secondary and '05' primary
                resid[C][E][i].GetYaxis().SetNdivisions(505)

                if E: resid[C][E][i].SetAxisRange(hier[0], hier[1], 'x')
                else: resid[C][E][i].SetAxisRange(loer[0], loer[1], 'x')

                ### LSveto plotting
                if E and i==8:
                    resid[C][E][i].SetAxisRange(0, lsHiE, 'x')
                
                resid[C][E][i].SetAxisRange(0.1,10,'y')
                
                
                ###---------------------------------------------
                if linres:
                    botpad[C][E][i].SetLogy(0)
                    resid[C][E][i].SetAxisRange(lrs[0], lrs[1], 'y')
                ###---------------------------------------------

                
                ###==================================================================
                ### set errors on resid R where R=Data/Total
                ### sigR = R*sqrt((sigD/D**2) + (sigT/T)**2)
                ### but sigT=0 (for now) so...
                ### sigR = R*(sigD/D)
                ###------------------------------------------------------------------
                if dkey:
                    for n in range(resid[C][E][i].GetNbinsX()+1):
                        R  = resid[C][E][i].GetBinContent(n)
                        D  = data[dkey]['hist'].GetBinContent(n)
                        sD = data[dkey]['hist'].GetBinError(n)
                        
                        try:
                            resid[C][E][i].SetBinError(n, R * (sD/D))
                        except:
                            resid[C][E][i].SetBinError(n,0)
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
                ###------------------------------------------------------------------
                ###==================================================================
                
                
                resid[C][E][i].Draw()

                # set my line to '1'
                zero = TLine(eran[E][0], 1, eran[E][1], 1)
                if E and i==8:
                    zero = TLine(0, 1, lsHiE, 1)
                zeros[C][E].append(zero)
                zeros[C][E][i].SetLineColor(kRed)
                zeros[C][E][i].SetLineWidth(1)
                zeros[C][E][i].Draw()
                
                #legs2[C][E][i].AddEntry(resid[C][E][i],'data / MC',legopt)
                #legs2[C][E][i].Draw()
                #---------------------------------------------------------
            

            save = ''
            #if local: save += 'local'
            #else: save += 'on-cup'
            save += str(runtag)
            save += '_globals-'+globstr
            save += '_chan'+chan
            save += '_E'+str(E)
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
            #save += '_chans'+fitchans
            #save += '_chan'+chan
            #save += '_extend'+str(extend)
            #save += '_others'+str(others)
            if note: save += '_'+note
            save += '_'+V
            
            canvs[C][E].Update()
            canvs[C][E].Print(plotdir+'/'+save+'.png')
            
            
            ### Save separate crystal histos?
            #-------------------------------------------------------------
            if indi:
                for i in range(numx):
                    if i+1 in justthese:
                        #isave = 'xstal-'+str(i+1)
                        try:
                            sepTopPad[C][E][i] = toppad[C][E][i].Clone()
                            sepBotPad[C][E][i] = botpad[C][E][i].Clone()
                            sepPlots[C][E][i] = TCanvas('ican-'+str(chan)+str(E)+str(i+1),
                                                        'ican-'+str(chan)+str(E)+str(i+1),
                                                        0, 0, 1200, 700)
                            sepTopPad[C][E][i].Draw()
                            sepBotPad[C][E][i].Draw()
                            sepPlots[C][E][i].Update()
                            isave  = ''
                            isave += 'x'+str(i+1)
                            isave += '-cs'+str(fitchans)
                            isave += '-c'+str(chan)
                            isave += '-e'+str(E)
                            #isave += '-ext'+str(extend)
                            isave += '-oth'+str(others)
                            if note: isave += '_'+str(note)
                            isave += '_'+str(V)

                            sepPlots[C][E][i].Print(plotdir+'/'+isave+'.png')

                        except:
                            print 'WARNING: could not make plot for xstal-'+str(i+1)
                            continue

    
    ### Save 4 plots to one canvas
    #-------------------------------------------------------------
    if indi:
        combPlots = [0 for x in range(numx)]
        tpad = [[0 for x in range(numx)] for x in range(4)]
        bpad = [[0 for x in range(numx)] for x in range(4)]
        
        for i in range(numx):
            if i+1 in justthese:
                try:
                    combPlots[i] = TCanvas('ccan-'+str(i+1),'ccan-'+str(i+1),0,0,1400,900)
                    combPlots[i].Divide(2,2)
                    p=0
                    for C, chan in enumerate(pltchans):
                        for E in range(numE):
                            tpad[p][i] = toppad[C][E][i].Clone()
                            bpad[p][i] = botpad[C][E][i].Clone()
                            combPlots[i].cd(p+1)
                            tpad[p][i].Draw()
                            bpad[p][i].Draw()
                            combPlots[i].Update()
                            p += 1
                    csave  = ''
                    csave += 'x'+str(i+1)
                    csave += '_combined'
                    csave += '_globals-'+globstr
                    #csave += '_ext'+str(extend)
                    csave += '_oth'+str(others)
                    if note: csave += '_'+str(note)
                    csave += '_'+str(V)
                    
                    combPlots[i].Print(plotdir+'/'+csave+'.png')
                    
                except:
                    print 'WARNING: could not make plot for xstal-'+str(i+1)
                    continue


    ### Plot the 0-40 keV region of single hit
    #-------------------------------------------------------------       
    if indi:
        zscanv = [0 for x in range(numx)]
        zshist = [0 for x in range(numx)]
        zsdata = [0 for x in range(numx)]
        zstotal = [0 for x in range(numx)]
        #zsresid = [0 for x in range(numx)]
        zsresid = makeResid100('S', 0, params[0])
        zsline = [0 for x in range(numx)]
        zstop = [0 for x in range(numx)]
        zsbot = [0 for x in range(numx)]
        
        font = 43
        size = 24
        loff = 0.005
        xoff = 3.6
        yoff = 1.1
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:

                zscanv[i] = TCanvas('zscanv-'+str(i+1),'zscanv-'+str(i+1),
                                0, 0, 1200, 700)

                fraction = 0.3
                pad1 = TPad('pad1_'+str(i),'pad1_'+str(i),0,fraction,1,1)
                #zstop[i].append(pad1)
                zstop[i] = (pad1)
                pad2 = TPad('pad2_'+str(i),'pad2_'+str(i),0,0,1,fraction)
                #zsbot[i].append(pad2)
                zsbot[i] = (pad2)
                
                zstop[i].SetBottomMargin(0.01)
                zstop[i].SetBorderMode(0)
                
                #zsbot[i].SetLogy(1)
                zstop[i].SetLogy(0)
                
                zsbot[i].SetTopMargin(0.05)
                zsbot[i].SetBottomMargin(0.3)
                zsbot[i].SetBorderMode(0)
                
                #zsbot[i].SetLogy(1)
                zsbot[i].SetLogy(0)
                
                zstop[i].Draw()
                zsbot[i].Draw()

                zstop[i].cd()

                #gStyle.SetPadLeftMargin(0.09)

                ### just checking residual plot names
                #prim_list = botpad[0][0][i].GetListOfPrimitives()
                #for bla in prim_list:
                #    print bla.GetName(), bla

                prim_list = toppad[0][0][i].GetListOfPrimitives()
                for bla in prim_list:
                    #print bla.GetName(), bla

                    if 'data' in bla.GetName():
                        zsdata[i] = (toppad[0][0][i].GetPrimitive(bla.GetName()))
                    if 'total' in bla.GetName():
                        zstotal[i] = (toppad[0][0][i].GetPrimitive(bla.GetName()))
                    
                    try:
                        zshist[i] = (toppad[0][0][i].GetPrimitive(bla.GetName()))
                        
                        zshist[i].SetXTitle('energy (keV)')
                        zshist[i].GetXaxis().SetTitleFont(font)
                        zshist[i].GetXaxis().SetTitleSize(size)
                        zshist[i].GetXaxis().SetTitleOffset(xoff)
                        zshist[i].GetXaxis().SetLabelFont(font)
                        zshist[i].GetXaxis().SetLabelSize(size)
                        zshist[i].GetXaxis().SetLabelOffset(loff)
                                                
                        zshist[i].SetYTitle('counts (dru)')
                        zshist[i].GetYaxis().SetTitleFont(font)
                        zshist[i].GetYaxis().SetTitleSize(size)
                        zshist[i].GetYaxis().SetTitleOffset(yoff)
                        zshist[i].GetYaxis().SetLabelFont(font)
                        zshist[i].GetYaxis().SetLabelSize(size)
                        zshist[i].GetYaxis().SetLabelOffset(loff)
                        
                        zshist[i].SetAxisRange(0, zmaxE, 'x')

                        if   i+1 == 1: zshist[i].SetAxisRange(0,  8.0, 'y')
                        elif i+1 == 2: zshist[i].SetAxisRange(0,  4.0, 'y')
                        elif i+1 == 3: zshist[i].SetAxisRange(0,  4.0, 'y')
                        elif i+1 == 4: zshist[i].SetAxisRange(0,  4.0, 'y')
                        elif i+1 == 5: zshist[i].SetAxisRange(0, 10.0, 'y')
                        elif i+1 == 6: zshist[i].SetAxisRange(0,  3.5, 'y')
                        elif i+1 == 7: zshist[i].SetAxisRange(0,  3.5, 'y')
                        elif i+1 == 8: zshist[i].SetAxisRange(0, 10.0, 'y')
                        else:          zshist[i].SetAxisRange(0,  5.0, 'y')
                                                
                        if 'data' in bla.GetName():
                            zshist[i].Draw('same')
                        else:
                            zshist[i].Draw(dopt)
                            
                    except: continue
                
                try:
                    #leg = TLegend(toppad[0][0][i].GetPrimitive(bla.GetName()))
                    #leg.Draw("same")
                    legs[0][0][i].Draw('same')
                except: continue

                zsbot[i].cd()
                
                #databins = zsdata[i].GetNbinsX()
                #totalbins = zstotal[i].GetNbinsX()
                #residbins = zsresid[i].GetNbinsX()
                #print databins, totalbins, residbins
                #zsresid[i].Rebin(residbins/databins)
                zsresid[i].Rebin(loEplotRebin)
                
                #zsresid[i] = zsdata[i]/zstotal[i]
                #zsresid[i] = zsdata[i] - zstotal[i]
                
                ### in some cases I may not be plotting the total
                ### so do a 'try'
                try:
                    if zdr:
                        zsresid[i].Divide(zsdata[i], zstotal[i])
                    else:
                        zsresid[i] = zsdata[i] - zstotal[i]
                except:
                    #zsresid[i] = zsdata[i]
                    pass
                
                ### tweak the errors on the residual
                for n in range(zsresid[i].GetNbinsX()+1):
                    R  = zsresid[i].GetBinContent(n)
                    D  = zsdata[i].GetBinContent(n)
                    sD = zsdata[i].GetBinError(n)

                    if zdr: zsresid[i].SetBinError(n, R*(sD/D))
                    else:   zsresid[i].SetBinError(n, sD)
                    
                    #try:
                    #    zsresid[i].SetBinError(n, R*(sD/D))
                    #except:
                    #    zsresid[i].SetBinError(n, 0)
                

                ### now just formatting stuff
                zsresid[i].SetTitle('')

                zsresid[i].SetAxisRange(0, zmaxE, 'x')
                if zdr: zsresid[i].SetAxisRange(0.8, 1.2, 'y')
                else:   zsresid[i].SetAxisRange(-0.4, 0.4, 'y')
                
                zsresid[i].SetXTitle('energy (keV)')
                zsresid[i].GetXaxis().SetTitleFont(font)
                zsresid[i].GetXaxis().SetTitleSize(size)
                zsresid[i].GetXaxis().SetTitleOffset(xoff)
                zsresid[i].GetXaxis().SetLabelFont(font)
                zsresid[i].GetXaxis().SetLabelSize(size)
                zsresid[i].GetXaxis().SetLabelOffset(loff)
                
                if zdr: zsresid[i].SetYTitle('data/mc')
                else:   zsresid[i].SetYTitle('data-mc (dru)')
                zsresid[i].GetYaxis().SetTitleFont(font)
                zsresid[i].GetYaxis().SetTitleSize(size)
                zsresid[i].GetYaxis().SetTitleOffset(yoff)
                zsresid[i].GetYaxis().SetLabelFont(font)
                zsresid[i].GetYaxis().SetLabelSize(size)
                zsresid[i].GetYaxis().SetLabelOffset(loff)
                
                # '5' secondary and '05' primary
                # means 5 divisions will be shown
                zsresid[i].GetYaxis().SetNdivisions(505)

                
                ### make a line at '1' or '0'
                if zdr: zero = TLine(0, 1, zmaxE, 1)
                else:   zero = TLine(0, 0, zmaxE, 0)
                zsline[i] = zero
                zsline[i].SetLineColor(kRed)
                zsline[i].SetLineWidth(1)

                zsresid[i].Draw()
                zsline[i].Draw('same')
                
                zscanv[i].Update()
                zscanv[i].Print(plotdir+'/a_zoomSingleHit_c'+str(i+1)+'.png')

    
    ### Plot the 0-40 keV region of multi hit
    #-------------------------------------------------------------       
    if indi:
        zmcanv = [0 for x in range(numx)]
        zmhist = [0 for x in range(numx)]
        zmdata = [0 for x in range(numx)]
        zmtotal = [0 for x in range(numx)]
        zmresid = makeResid100('M', 0, params[0])
        zmline = [0 for x in range(numx)]
        zmtop = [0 for x in range(numx)]
        zmbot = [0 for x in range(numx)]
        
        font = 43
        size = 24
        loff = 0.005
        xoff = 3.6
        yoff = 1.1
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:
                
                zmcanv[i] = TCanvas('zmcanv-'+str(i+1),'zmcanv-'+str(i+1),
                                0, 0, 1200, 700)

                fraction = 0.3
                pad1 = TPad('pad1_'+str(i),'pad1_'+str(i),0,fraction,1,1)
                #zmtop[i].append(pad1)
                zmtop[i] = (pad1)
                pad2 = TPad('pad2_'+str(i),'pad2_'+str(i),0,0,1,fraction)
                #zmbot[i].append(pad2)
                zmbot[i] = (pad2)
                
                zmtop[i].SetBottomMargin(0.01)
                zmtop[i].SetBorderMode(0)
                
                #zmbot[i].SetLogy(1)
                zmtop[i].SetLogy(0)
                
                zmbot[i].SetTopMargin(0.05)
                zmbot[i].SetBottomMargin(0.3)
                zmbot[i].SetBorderMode(0)
                
                #zmbot[i].SetLogy(1)
                zmbot[i].SetLogy(0)
                
                zmtop[i].Draw()
                zmbot[i].Draw()

                zmtop[i].cd()

                
                prim_list = toppad[1][0][i].GetListOfPrimitives()
                for bla in prim_list:
                    #print bla.GetName(), bla

                    if 'data' in bla.GetName():
                        zmdata[i] = (toppad[1][0][i].GetPrimitive(bla.GetName()))
                    if 'total' in bla.GetName():
                        zmtotal[i] = (toppad[1][0][i].GetPrimitive(bla.GetName()))
                    
                    
                    try:
                        zmhist[i] = (toppad[1][0][i].GetPrimitive(bla.GetName()))

                        zmhist[i].SetXTitle('energy (keV)')
                        zmhist[i].GetXaxis().SetTitleFont(font)
                        zmhist[i].GetXaxis().SetTitleSize(size)
                        zmhist[i].GetXaxis().SetTitleOffset(xoff)
                        zmhist[i].GetXaxis().SetLabelFont(font)
                        zmhist[i].GetXaxis().SetLabelSize(size)
                        zmhist[i].GetXaxis().SetLabelOffset(loff)
                                                
                        zmhist[i].SetYTitle('counts (dru)')
                        zmhist[i].GetYaxis().SetTitleFont(font)
                        zmhist[i].GetYaxis().SetTitleSize(size)
                        zmhist[i].GetYaxis().SetTitleOffset(yoff)
                        zmhist[i].GetYaxis().SetLabelFont(font)
                        zmhist[i].GetYaxis().SetLabelSize(size)
                        zmhist[i].GetYaxis().SetLabelOffset(loff)
                        
                        zmhist[i].SetAxisRange(0, zmaxE, 'x')

                        if   i+1 == 1: zmhist[i].SetAxisRange(0, 1.6, 'y')
                        elif i+1 == 2: zmhist[i].SetAxisRange(0, 4.5, 'y')
                        elif i+1 == 3: zmhist[i].SetAxisRange(0, 2.5, 'y')
                        elif i+1 == 4: zmhist[i].SetAxisRange(0, 2.0, 'y')
                        elif i+1 == 5: zmhist[i].SetAxisRange(0, 2.0, 'y')
                        elif i+1 == 6: zmhist[i].SetAxisRange(0, 1.2, 'y')
                        elif i+1 == 7: zmhist[i].SetAxisRange(0, 1.2, 'y')
                        elif i+1 == 8: zmhist[i].SetAxisRange(0, 0.8, 'y')
                        else:          zmhist[i].SetAxisRange(0, 5.0, 'y')

                        if 'data' in bla.GetName():
                            zmhist[i].Draw('same')
                        else:
                            zmhist[i].Draw(dopt)
                        
                    except: continue
                
                try:
                    #leg = TLegend(toppad[0][0][i].GetPrimitive(bla.GetName()))
                    #leg.Draw("same")
                    legs[1][0][i].Draw('same')
                except: continue

                zmbot[i].cd()
                
                #databins = zmdata[i].GetNbinsX()
                #totalbins = zmtotal[i].GetNbinsX()
                #residbins = zmresid[i].GetNbinsX()
                #print databins, totalbins, residbins
                #zmresid[i].Rebin(residbins/databins)
                zmresid[i].Rebin(loEplotRebin)
                
                #zmresid[i] = zmdata[i]/zmtotal[i]
                #zmresid[i] = zmdata[i] - zmtotal[i]
                
                ### in some cases I may not be plotting the total
                ### so do a 'try'
                try:
                    if zdr:
                        zmresid[i].Divide(zmdata[i], zmtotal[i])
                    else:
                        zmresid[i] = zmdata[i] - zmtotal[i]
                except:
                    #zmresid[i] = zmdata[i]
                    pass
                
                ### tweak the errors on the residual
                for n in range(zmresid[i].GetNbinsX()+1):
                    R  = zmresid[i].GetBinContent(n)
                    D  = zmdata[i].GetBinContent(n)
                    sD = zmdata[i].GetBinError(n)

                    if zdr: zmresid[i].SetBinError(n, R*(sD/D))
                    else:   zmresid[i].SetBinError(n, sD)
                    
                    #try:
                    #    zmresid[i].SetBinError(n, R*(sD/D))
                    #except:
                    #    zmresid[i].SetBinError(n, 0)
                

                ### now just formatting stuff
                zmresid[i].SetTitle('')

                zmresid[i].SetAxisRange(0, zmaxE, 'x')
                if zdr: zmresid[i].SetAxisRange(0.8, 1.2, 'y')
                else:   zmresid[i].SetAxisRange(-0.4, 0.4, 'y')
                
                zmresid[i].SetXTitle('energy (keV)')
                zmresid[i].GetXaxis().SetTitleFont(font)
                zmresid[i].GetXaxis().SetTitleSize(size)
                zmresid[i].GetXaxis().SetTitleOffset(xoff)
                zmresid[i].GetXaxis().SetLabelFont(font)
                zmresid[i].GetXaxis().SetLabelSize(size)
                zmresid[i].GetXaxis().SetLabelOffset(loff)
                
                if zdr: zmresid[i].SetYTitle('data/mc')
                else:   zmresid[i].SetYTitle('data-mc (dru)')
                zmresid[i].GetYaxis().SetTitleFont(font)
                zmresid[i].GetYaxis().SetTitleSize(size)
                zmresid[i].GetYaxis().SetTitleOffset(yoff)
                zmresid[i].GetYaxis().SetLabelFont(font)
                zmresid[i].GetYaxis().SetLabelSize(size)
                zmresid[i].GetYaxis().SetLabelOffset(loff)
                
                # '5' secondary and '05' primary
                # means 5 divisions will be shown
                zmresid[i].GetYaxis().SetNdivisions(505)

                
                ### make a line at '1' or '0'
                if zdr: zero = TLine(0, 1, zmaxE, 1)
                else:   zero = TLine(0, 0, zmaxE, 0)
                zmline[i] = zero
                zmline[i].SetLineColor(kRed)
                zmline[i].SetLineWidth(1)

                zmresid[i].Draw()
                zmline[i].Draw('same')
                
                zmcanv[i].Update()
                zmcanv[i].Print(plotdir+'/a_zoomMultiHit_c'+str(i+1)+'.png')

                
    #-----------------------------------------------------------------
    ### print out the fit results
    if fitting:
        print '\n\n'
        print '!!!!!  FIT RESULTS  !!!!!\n'
        for key in resultskeys:
            if int(key)+1 in justthese:
                for line in fitresults[key]:
                    print line
    #-----------------------------------------------------------------

    
    # delete the extra crap
    #-------------------------------
    
    try: del fcanv
    except: pass
    
    try: del mcanv
    except: pass
    
    try: del sepFitPlot
    except: pass
    
    #try: del combPlots
    #except: pass

    #try: del zscanv
    #except: pass
    
    try: del zmcanv
    except: pass
    
    try: del sepPlots
    except: pass
    
    try: del canvs
    except: pass
    
    if fitting:
        print 'Time to complete the fit = '+str(fitTime)+' sec \n'
    
    if not batch: raw_input('[Enter] to quit \n')
    else: print "Running in batch mode - quitting now! \n"

    return


if __name__ == "__main__":
    main(sys.argv[1:])

