#!/usr/bin/env python

############################################################
# Matt Kauer - mkauer@icecube.wisc.edu
#-----------------------------------------------------------
# Revamp the MC format and hist generation
# 
# version: 2021-01-12
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

# ROOT6 FIX
import ctypes
from ctypes import *

script = os.path.basename(__file__)
V = 'v1'
print('INFO: running script --> {0}'.format(script))

import ROOT
from ROOT import *
#ROOT.gErrorIgnoreLevel = kWarning
ROOT.gErrorIgnoreLevel = kError
try: vroot = ROOT.__version__
except: vroot = '5'
print('INFO: using root --> {0}'.format(vroot))

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from funcs500 import *

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


mcfile = 'backgrounds.txt'



print('INFO: using backgrounds config file --> {0}'.format(mcfile))


### ==========  OPTIMIZATION OPTIONS  ================================
### test MC energy shift in bins (12bins = 1keV)
### currently just shifting teflon
### might need different shifting for different xstals
# for xstals 1   2   3   4   5   6   7   8   9
binShift =  [4,  4,  6,  8, -1,  8,  6, -1, -1]
#binShift =  [-1, 0, 11, 11, -1, 11, 11, -1, -1]
#shiftWhat = ['teflonbulk','teflonsurf']
shiftWhat = ['teflon']

### set the fitter step size
### can help to go smaller than default 1e-2
stepSize = 0
#stepSize = 1e-3
stepSize = 1e-4

### test MC smoothing?
### smoothing window in +/- number of bins
smoothing = 0
#smoothing = 5
smoothWhat = ['cushield', 'steel']


### ==========  FITTING OPTIONS  =====================================
### select channels to fit
fitchans = ['S', 'M']
#fitchans = ['S']  # testing low E fit

### Let's finally try different fit ranges!!
fitranges = [{} for x in range(numx)]
for i in range(numx):

    ### all channels default except alphas
    # --------------------------------------------------------------
    fitranges[i]['S0'] = [1,   2,   90]  # single-hit low-energy
    fitranges[i]['S1'] = [6,  80, 4000]  # single-hit high-energy
    fitranges[i]['S2'] = [1,   0,    0]  # single-hit alphas
    fitranges[i]['M0'] = [1,   2,   90]  # multi-hit low-energy
    fitranges[i]['M1'] = [6,  80, 4000]  # multi-hit high-energy
    fitranges[i]['M2'] = [1,   0,    0]  # multi-hit alphas
    
    ### for just low E fit
    #fitranges[i]['S0'] = [1, 4, 40]
    #fitranges[i]['S1'] = [1, 0,  0]
    #fitranges[i]['M0'] = [1, 2, 40]
    #fitranges[i]['M0'] = [1, 20, 40]
    #fitranges[i]['M1'] = [1, 0,  0]
    
    ### for just high E fit
    fitranges[i]['S0'] = [1, 0, 0]
    fitranges[i]['S1'] = [2, 80, 3500]
    fitranges[i]['M0'] = [1, 40, 100]
    #fitranges[i]['M0'] = [1, 0, 0]
    fitranges[i]['M1'] = [2, 80, 3200]

    ### for just alpha fit
    #fitranges[i]['S0'] = [1, 0, 0]
    #fitranges[i]['S1'] = [1, 0, 0]
    #fitranges[i]['M0'] = [1, 0, 0]
    #fitranges[i]['M1'] = [1, 0, 0]
    #fitranges[i]['S2'] = [1, 2000, 3200]

    ### for just lsveto fit
    #fitranges[i]['S0'] = [1, 0, 0]
    #fitranges[i]['S1'] = [1, 0, 0]
    #fitranges[i]['M0'] = [1, 0, 0]
    #fitranges[i]['M1'] = [1, 0, 0]
    
### special case for C1, C5, C8
#for i in [0, 4, 7]:
#    fitranges[i]['S0'] = [1, 10, 60]
#    fitranges[i]['S1'] = [6, 0, 0]
#    fitranges[i]['M0'] = [1, 10, 90]
#    fitranges[i]['M1'] = [6, 80, 4000]


### defaults for lsveto
# --------------------------------------------------------------
fitranges[8]['S0'] = [1, 0, 0]
fitranges[8]['S1'] = [4, 400, 4000]
fitranges[8]['M0'] = [1, 0, 0]
fitranges[8]['M1'] = [4, 200, 4000]

# not lsveto
#fitranges[8]['S0'] = [1, 0, 0]
#fitranges[8]['S1'] = [1, 0, 0]
#fitranges[8]['M0'] = [1, 0, 0]
#fitranges[8]['M1'] = [1, 0, 0]



### ==========  EXTRA MC OPTIONS  ====================================
### which MC to fit globally (to all crystals simultaneously)?
globalmc = [
    'pmt',
    'plastic',
    'lsveto',
    'cushield',
    'innersteel',
    'steel',
    'gamma',
]

### include bkgs from 'other' internals? [0,1]
others = 1

### plot components in groups? [0,1]
ingroups = 1

### show the total? [0,1]
showTotal = 1

### plot the total in gray==0, red==1? [0,1]
redtotal = 1

### plot alpha channels? [0,1]
showAlphas = 1

### only show q1 quenched peaks? [0,1]
justQ1 = 0

### show the legends? [0,1]
showlegs = 1

### combine 'others' into the makePlots() plots? [0,1]
combine = 1

### print out chi2 of data to total? [0,1]
printChi2 = 0

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
pltchans = ['S','M']

### crystal plotting ranges
loer = [0, 100]
hier = [0, 4000]
aler = [0, 6000]
#hiYlo = 1e-3 # for 3000 X range
hiYlo = 1e-5 # for 4000 X range
eran = [loer, hier, aler]

### lsveto plotting ranges
lsHiE = 8000
lsLoSY = 1e-6
lsHiSY = 10
lsLoMY = 1e-8
lsHiMY = 1

### special range for zoomed in plot
zmaxE = 40

### zoomed in residual as dat-mc [0] or dat/mc [1]
zdr = 0

### rebin the final plots [1,inf]
loEplotRebin = 3
hiEplotRebin = 4
alphaPlotRebin = 4

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
    allchans = list(set(fitchans+pltchans))
    print('INFO: loading/building all the histograms...')
    data, bkgs, sigs = build500(os.path.join(HERE, mcfile), others, reuse, allchans, xstals)

    ### only plot/fit crystals that have data
    justthese = []
    for i in range(1, numx+1):
        for key in data:
            if 'x'+str(i) in key and i not in justthese:
                justthese.append(i)
    justthese.sort()
    if len(justthese) == 0:
        print('ERROR: nothing to plot, quitting...')
        sys.exit()
    print('INFO: plotting/fitting crystals --> {0}'.format(justthese))
    
    ### create dir for saving the plots...
    plotdir = HERE+'/plots/c'
    for x in justthese:
        plotdir += str(x)
    if not os.path.exists(plotdir): 
        os.makedirs(plotdir)

    ### copy the backgrounds file before fitting so it isn't overwritten later
    shutil.copyfile(os.path.join(HERE, mcfile), os.path.join(plotdir, mcfile))

    ### copy the actual script as well just to have it on hand
    shutil.copyfile(os.path.join(HERE, script), os.path.join(plotdir, script))

    ### save all raw histograms to a rootfile so it can be re-used
    rootoutfile = TFile(plotdir+"/histograms-raw.root", "RECREATE")
    for hist_dict in [data, bkgs, sigs]:
        writeHists500(hist_dict)
    rootoutfile.Write()
    rootoutfile.Close()

    
    ### make of list runtimes per crystal, channel, energy
    runtimes = getRuntimes(data)

    ### sort keys for convenience
    datkeys = sortDataKeys92(data)
    if datsumw2:
        for key in datkeys:
            data[key]['hist'].Sumw2()
    
    ### scale into dru or not
    data = scaleData411(data, dru)
    if dru:
        bkgs = scaleBkgs411(bkgs)
        sigs = scaleBkgs411(sigs)
    else:
        bkgs = scaleBkgs411(bkgs, runtimes)
        sigs = scaleBkgs411(sigs, runtimes)
    
    ### testing MC energy scaling
    if binShift != 0:
        #print 'energy shifting...'
        sigs = ScaleEnergy410(sigs, binShift, shiftWhat)
        bkgs = ScaleEnergy410(bkgs, binShift, shiftWhat)
    
    ### make plots before combining?
    #makePlots93(bkgs, combine, others)
    #makePlots93(sigs, combine, others)
    #sys.exit()
    
    ### combine after scaling
    ### FIXME - combine but don't delete others?
    sigs = combineOthers500(sigs, globalmc)
    bkgs = combineOthers500(bkgs, globalmc)
    
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
        print('INFO: done smoothing histograms')

        
    #-----------------------------------------------------------------
    #-----------------------------------------------------------------

    
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
    print('INFO: Total number of unique bkgs and sigs = {0}'.format(Nc))
    # ROOT6 FIX
    colors, cis = rainbowSix(uniqAll)
    uniqColor = cis
    
    ### colors for the groups
    gis = {
        'internal':  kBlue,
        'cosmo':     kMagenta+1,
        'surface':   kCyan+1,
        'cucase':    kYellow+1,
        'xpmts':     kGreen,
        'pmts':      kGreen+1,
        'plastic':   kOrange,
        'lsveto':    kOrange+1,
        'cushield':  kOrange+2,
        'steel':     kYellow,
        'gamma':     kOrange+4,
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
        
        print('INFO: global fits to --> {0}'.format(globstr))
    
        print('INFO: building histograms for fit...')

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
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
            startbin = fitdata.GetXaxis().GetNbins()
            #print 'startbin', startbin
            if startbin == 1:
                startbin = 0
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
                    # alpha Energy
                    if 'x'+str(i+1) in dkey and '-c'+C in dkey and '-e2' in dkey:
                        fitdata = extendHist(fitdata,
                                        fitPrep(data, dkey, fitranges[i][C+'2']))
            
            stopbin = fitdata.GetXaxis().GetNbins()
            ### this start/stop binning stuff might still be a little wonky
            xstalbins[i] = [startbin, stopbin-1]
            if debug: print('DEBUG: bin range for C{0} = {1}'.format(i+1, xstalbins[i]))
            #print fitdata.GetXaxis().GetNbins()
            
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
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
            for nc, C in enumerate(fitchans):
                
                # signals
                #---------------------------------------------------------------
                for skey in sigkeys:
                    
                    ### init histograms for global signals (pmt, lsveto, steel, etc.)
                    bits = skey.split('-')
                    if bits[1] in globalmc:
                    #for gmckey in globalmc:
                    #    if gmckey in skey:
                    #        bits = skey.split('-')
                        fgkey = bits[1]+'-'+bits[2]
                        if fgkey not in fglobsigkeys:
                            fglobsigkeys.append(fgkey)
                            fitglobsigs[fgkey] = {}
                            fglob = TH1F(fgkey, fgkey, fitbins, 0, fitbins)
                            fitglobsigs[fgkey]['hist'] = fglob
                    
                    # pad out the empty histograms
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
                        if '-c'+C in skey and '-e2' in skey:
                            fskey = skey.split('-c'+C+'-e2')[0]
                            fitsigs = addHistKey(fitsigs, fskey)
                            fitsigs[fskey]['hist'] = extendHist(fitsigs[fskey]['hist'],
                                                fitPrep(sigs, skey, fitranges[i][C+'2'], 1))
                    
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

                    # alpha Energy
                    if 'x'+str(i+1) in skey and '-c'+C in skey and '-e2' in skey:
                        bits = skey.split('-c'+C+'-e2')
                        fskey = bits[0]
                        fitsigs = addHistKey(fitsigs, fskey)
                        fitsigs[fskey]['hist'] = extendHist(fitsigs[fskey]['hist'],
                                                    fitPrep(sigs, skey, fitranges[i][C+'2']))

                        #fitsigs[fskey]['hist'].Draw()
                        #input('enter to continue')
                        
                # backgrounds
                #---------------------------------------------------------------
                for bkey in bakkeys:
                    
                    ### init histograms for global backgrounds (pmt, lsveto, steel, etc.)
                    bits = bkey.split('-')
                    if bits[1] in globalmc:
                    #for gmckey in globalmc:
                    #    if gmckey in bkey:
                    #        bits = bkey.split('-')
                        fgbkey = bits[1]+'-'+bits[2]
                        if fgbkey not in fglobbkgkeys:
                            fglobbkgkeys.append(fgbkey)
                            fitglobbkgs[fgbkey] = {}
                            fglob = TH1F(fgbkey, fgbkey, fitbins, 0, fitbins)
                            fitglobbkgs[fgbkey]['hist'] = fglob
                    
                    # pad out the empty histograms
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
                        if '-c'+C in bkey and '-e2' in bkey:
                            fbkey = bkey.split('-c'+C+'-e2')[0]
                            fitbkgs = addHistKey(fitbkgs, fbkey)
                            fitbkgs[fbkey]['hist'] = extendHist(fitbkgs[fbkey]['hist'],
                                                    fitPrep(bkgs, bkey, fitranges[i][C+'2'], 1))
                        
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
                    # alpha Energy
                    if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e2' in bkey:
                        bits = bkey.split('-c'+C+'-e2')
                        fbkey = bits[0]
                        fitbkgs = addHistKey(fitbkgs, fbkey)
                        fitbkgs[fbkey]['hist'] = extendHist(fitbkgs[fbkey]['hist'],
                                                    fitPrep(bkgs, bkey, fitranges[i][C+'2']))
            

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
                if fgskey in fskey:
                    fitglobsigs[fgskey]['hist'].Add(fitsigs[fskey]['hist'])
        
        # (4) remove global sigs from fit sigs
        L = len(fsigkeys)-1
        for k, fskey in enumerate(reversed(fsigkeys)):
            for gmckey in globalmc:
                ### make gmckey unique (steel vs innersteel) 2018-10-23
                gmckey = '-'+gmckey+'-'
                if gmckey in fskey:
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
                if fgbkey in fbkey:
                    fitglobbkgs[fgbkey]['hist'].Add(fitbkgs[fbkey]['hist'])
        
        # (4) remove global bkgs from fit bkgs
        L = len(fbkgkeys)-1
        for k, fbkey in enumerate(reversed(fbkgkeys)):
            for gmckey in globalmc:
                ### make gmckey unique (steel vs innersteel) 2018-10-23
                gmckey = '-'+gmckey+'-'
                if gmckey in fbkey:
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
            if debug: print('DEBUG: adding hist --> {0}'.format(key))
            fitdata.Add(fitbkgs[key]['hist'], -1)

        # (2) subtract global backgrounds
        for key in fglobbkgkeys:
            if debug: print('DEBUG: adding hist --> {0}'.format(key))
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
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
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
                                    print tmpkey
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
                                #print newkey
                                sigs[newkey]['hist'].Scale(dat_int/mc_int)
                                sigs[newkey]['fitscale'] = sigs[newkey]['scale'] * dat_int/mc_int
                                #print sigs[newkey]['fitscale']
                                
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
                                        print('ERROR: do not know what to do with useBounds = {0}'.format(useBounds))
                                        sys.exit()
                    #print '!!!', sigs[newkeys[0]]['info']['newfbnd']
                    bounds.append(sigs[newkeys[0]]['info']['newfbnd'])
                    #---------------------------------------------------------------------
                    
                    ### set errors to zero?
                    fitsigs[fskey]['hist'] = zeroBinError(fitsigs[fskey]['hist'])

            
            fitresults[str(i)].append('Crystal-'+str(i+1)+' fit results')
            fitresults[str(i)].append('runtime = '+str(round(runtimes[i]['S'][0]/60./60./24., 2))+' days')
            if note: fitresults[str(i)].append('note = '+note)
            fitresults[str(i)].append('script = '+script)
            fitresults[str(i)].append('root = '+vroot)
            fitresults[str(i)].append('selected xstals = '+str(justthese))
            fitresults[str(i)].append('channels fit = '+str(fitchans))
            fitresults[str(i)].append('globals = '+str(globalmc))
            fitresults[str(i)].append('use others = '+str(others))
            #fitresults[str(i)].append('hist extend = '+str(extend))
            fitresults[str(i)].append('norm to dru = '+str(dru))
            fitresults[str(i)].append('energy shift in bins = '+str(binShift[i]))
            fitresults[str(i)].append('shifting these = '+str(shiftWhat))
            fitresults[str(i)].append('fitter step size = '+str(stepSize))
            for chan in ['S0', 'S1', 'S2', 'M0', 'M1', 'M2']:
                fitresults[str(i)].append('{0} fit --> rebin = [{1}] --> range = [{2} - {3}] keV'
                                          .format(chan,
                                                  fitranges[i][chan][0],
                                                  fitranges[i][chan][1],
                                                  fitranges[i][chan][2]))
        
        
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
                                        print('ERROR: do not know what to do with useBounds = {0}'.format(useBounds))
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
        print('INFO: total number of hists being fit = {0}'.format(totalNumFits))
        sigObj = TObjArray(totalNumFits)
        
        for fskey in fsigkeys:
            #print 'loading signal to fit -->', fskey
            #print fitsigs[fskey]['hist'].GetNbinsX()
            sigObj.append(fitsigs[fskey]['hist'])

        for fgkey in fglobsigkeys:
            #print 'loading global to fit -->', fskey
            sigObj.append(fitglobsigs[fgkey]['hist'])

        #print 'data', fitdata.GetNbinsX()

        # init the fitter class
        fit = TFractionFitter(fitdata, sigObj)
        fitter = fit.GetFitter()
        #sys.exit()
        
        ### set fit bounds!!!
        for l in range(len(bounds)):
            #fit.Constrain(l+1, bounds[l][0], bounds[l][1])
            # ROOT6 FIX
            fit.Constrain(l, bounds[l][0], bounds[l][1])
            if stepSize != 0:
                fitter.Config().ParSettings(l).SetStepSize(stepSize)
            
        ### set the fit range
        #fit.SetRangeX(0, fmax*numx)
        fit.SetRangeX(0, fitbins-1)

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
            print('\n\n*******************  FIT FAILURE  *******************\n\n')
            print('Time to fail = {} sec \n'.format(fitTime))
            sys.exit()
        
        print('\n\n*******************  SUCCESSFUL FIT  *******************\n\n')
        
        
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
        
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
            
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
            #for i in range(numx):
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
                            if finit:
                                E = str(E)
                                newkey = fskey+'-c'+C+'-e'+E
                                try: test = sigs[newkey]['hist']
                                except: continue

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
            resultsfile = os.path.join(plotdir, save+'.txt')
            outfile = open(resultsfile, 'w')
            for key in resultskeys:
                if int(key)+1 in justthese:
                    for line in fitresults[key]:
                        outfile.write(line+'\n')
            outfile.close()

            ### save scaled histograms to rootfile?
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
        updateBkgsFile300(justthese, os.path.join(plotdir, mcfile), resultsfile, plotdir)

        
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
        for j in range(numx):
            sepfittotal[j] = TH1F('sepFitTotal-'+str(j+1), 'sepFitTotal-'+str(j+1), fitbins, 0, fitbins)
            if redtotal:
                sepfittotal[j].SetLineColor(kRed)
                sepfittotal[j].SetMarkerColor(kRed)
            #else:
            #    sepfittotal[j].SetLineColor(kGray+1)
            #    sepfittotal[j].SetMarkerColor(kGray+1)
            sepfittotal[j].SetLineWidth(1)
            
            sepfitresid[j] = TH1F('sepFitResid-'+str(j+1), 'sepFitResid-'+str(j+1), fitbins, 0, fitbins)
            sepfitresid[j].SetLineColor(kBlack)
            sepfitresid[j].SetMarkerColor(kBlack)
            sepfitresid[j].SetLineWidth(1)
            
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
            newFitTitle = str('Crystal-'+str(i+1)+'   '+'Fit-chans-'+str(fitchans))
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
            """
            if dru:
                sepfitdata[i].SetAxisRange(2e-3, 2e3, 'y')
                if i == 8:
                    sepfitdata[i].SetAxisRange(2e-5, 2e1, 'y')
            """
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

                    ### set errors to 0 for better viewing
                    sepfitsigs[i][fskey]['hist'] = zeroBinError(sepfitsigs[i][fskey]['hist'])
                    
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

                ### set errors to 0 for better viewing
                sepfitglobsigs[i][fgkey]['hist'] = zeroBinError(sepfitglobsigs[i][fgkey]['hist'])
                
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

                if redtotal:
                    sepfittotal[i].SetLineColor(kRed)
                    sepfittotal[i].SetMarkerColor(kRed)

                ### set errors to 0 for better viewing
                sepfittotal[i] = zeroBinError(sepfittotal[i])

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
                    #fisave += '_globals-'+globstr
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
        msave += '_all-xstals'
        msave += '_the-indi-fits'
        #msave += '_globals-'+globstr
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
        """
        if dru:
            fitdata.SetAxisRange(2e-3, 2e3, 'y')
            #if i == 8:
            #    fitdata.SetAxisRange(2e-5, 2e1, 'y')
        """
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

            # FIXME
            # why if i==0???
            #if i==0:
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

        cfsave  = ''
        cfsave += str(runtag)
        cfsave += '_all-xstals'
        cfsave += '_the-comb-fit'
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
        if note: cfsave += '_'+note
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
    if showAlphas:
        numE = 3
    else:
        numE = 2
        
    # number of channels (single and/or multi-hit)
    #numC = len(pltchans)
    numC = 2
    
    canvs  = [[0 for e in range(numE)] for c in range(numC)]
    
    ### for separate plots
    sepPlots = [[[0 for i in range(numx)] for e in range(numE)] for c in range(numC)]
    sepTopPad = [[[0 for i in range(numx)] for e in range(numE)] for c in range(numC)]
    sepBotPad = [[[0 for i in range(numx)] for e in range(numE)] for c in range(numC)]
    
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
        
    plotRebin = 1
    for C, chan in enumerate(pltchans): 
    
        for E in range(numE):
            
            if E == 2: plotRebin = alphaPlotRebin
            elif E == 1: plotRebin = hiEplotRebin
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
            gStyle.SetPadLeftMargin   (0.10)
            
            font = 43
            size = 16
            yoff = 2.5
            xoff = 7
            """
            toppad[C][E] = []
            botpad[C][E] = []

            legs[C][E]   = []
            legs2[C][E]  = []
            zeros[C][E]  = []
            """
            total[C][E] = makeTotal100(chan, E, params[E])
            resid[C][E] = makeResid100(chan, E, params[E])
            
            for i in range(numx):
            #for i in [j-1 for j in justthese]:
                
                canvs[C][E].cd(i+1)
                
                fraction = 0.3
                pad1 = TPad('pad1'+chan+str(E),'pad1'+chan+str(E),0,fraction,1,1)
                toppad[C][E][i] = (pad1)
                pad2 = TPad('pad2'+chan+str(E),'pad2'+chan+str(E),0,0,1,fraction)
                botpad[C][E][i] = (pad2)
                
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
                    
                legs[C][E][i] = (leg)
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
                #if dru and i!=8:
                    #total[C][E][i].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
                    #if E: total[C][E][i].SetAxisRange(2e-3, 2e1, 'y')
                    #else: total[C][E][i].SetAxisRange(2e-3, 3e2, 'y')
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

                        if E == 2: data[dkey]['hist'].SetAxisRange(aler[0], aler[1], 'x')
                        elif E == 1: data[dkey]['hist'].SetAxisRange(hier[0], hier[1], 'x')
                        else: data[dkey]['hist'].SetAxisRange(loer[0], loer[1], 'x')

                        if dru and i!=8:
                            #data[dkey]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            #if E: data[dkey]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            if E: data[dkey]['hist'].SetAxisRange(hiYlo, 2e1, 'y')
                            else: data[dkey]['hist'].SetAxisRange(2e-3, 3e2, 'y')

                        ### LSveto plotting
                        if E and i==8:
                            data[dkey]['hist'].SetAxisRange(0, lsHiE, 'x')
                            if dru:
                                if chan == 'S':
                                    data[dkey]['hist'].SetAxisRange(lsLoSY, lsHiSY, 'y')
                                if chan == 'M':
                                    data[dkey]['hist'].SetAxisRange(lsLoMY, lsHiMY, 'y')

                        
                        data[dkey]['hist'].Draw('same')
                        days = round(data[dkey]['runtime']/86400., 2)
                        dataLegName = data[dkey]['info']['tag']\
                            +'  '+data[dkey]['info']['build']
                        legs[C][E][i].AddEntry(data[dkey]['hist'], dataLegName, legopt)
                        
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
                        """
                        if dru:
                            if E: bkgs[key]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            else: bkgs[key]['hist'].SetAxisRange(2e-3, 3e2, 'y')
                        """
                        if ingroups:
                            if bkgs[key]['info']['group'] == 'none':
                                try:
                                    gbkgs[C][E][i]['none'][key] = deepcopy(bkgs[key]['hist'])
                                except:
                                    gbkgs[C][E][i]['none'] = {}
                                    gbkgs[C][E][i]['none'][key] = deepcopy(bkgs[key]['hist'])
                            else:
                                try:
                                    gbkgs[C][E][i][bkgs[key]['info']['group']].Add(bkgs[key]['hist'])
                                except:
                                    gbkgs[C][E][i][bkgs[key]['info']['group']] = deepcopy(bkgs[key]['hist'])
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
                                    gbkgs[C][E][i]['none'][key] = deepcopy(sigs[key]['hist'])
                                except:
                                    gbkgs[C][E][i]['none'] = {}
                                    gbkgs[C][E][i]['none'][key] = deepcopy(sigs[key]['hist'])
                            else:
                                try:
                                    gbkgs[C][E][i][sigs[key]['info']['group']].Add(sigs[key]['hist'])
                                except:
                                    gbkgs[C][E][i][sigs[key]['info']['group']] = deepcopy(sigs[key]['hist'])
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
                                #legs[C][E][i].AddEntry(bkgs[bkey]['hist'], activ+bkey, legopt)
                                # make legend names smaller
                                bits = bkey.split('-')
                                legname = bits[1]+' '+bits[2]
                                legs[C][E][i].AddEntry(bkgs[bkey]['hist'], legname, legopt)
                        for skey in sigkeys:
                            if skey == 'x'+str(i+1)+'-'+name+'-c'+chan+'-e'+str(E) \
                               or skey == 'x'+str(i+1)+'-'+name+'-f'+str(i+1)+'-c'+chan+'-e'+str(E):
                                activ = '(%.2e) '%(sigs[skey]['info']['fitacti'])
                                #legs[C][E][i].AddEntry(sigs[skey]['hist'], activ+skey, legopt)
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
                                total[C][E][i].GetBinError(n+1)/(float(plotRebin)/math.sqrt(2.)))
                
                
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
                            #total[C][E][i].SetLineColor(kGray+1)
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
                leg2 = TLegend(0.72, 0.78, 0.94, 0.94)
                legs2[C][E][i] = (leg2)
                legs2[C][E][i].SetFillColor(0)
                legs2[C][E][i].SetBorderSize(0)
                legopt = 'LPE'

                if tcount and dkey:
                    total[C][E][i] = cleanSmallBins(total[C][E][i])
                    try:
                        resid[C][E][i].Divide(data[dkey]['hist'], total[C][E][i])
                    except:
                        print 'WARNING: could not divide', dkey, 'by total'
                        pass
                    
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

                if E == 2:   resid[C][E][i].SetAxisRange(aler[0], aler[1], 'x')
                elif E == 1: resid[C][E][i].SetAxisRange(hier[0], hier[1], 'x')
                else:        resid[C][E][i].SetAxisRange(loer[0], loer[1], 'x')

                ### LSveto plotting
                if E and i==8:
                    resid[C][E][i].SetAxisRange(0, lsHiE, 'x')
                
                                
                ###---------------------------------------------
                if linres:
                    botpad[C][E][i].SetLogy(0)
                    resid[C][E][i].SetAxisRange(lrs[0], lrs[1], 'y')
                else:
                    resid[C][E][i].SetAxisRange(0.1,10,'y')
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
                ###------------------------------------------------------------------
                ###==================================================================
                
                
                resid[C][E][i].Draw()

                # set my line to '1'
                zero = TLine(eran[E][0], 1, eran[E][1], 1)
                if E and i==8:
                    zero = TLine(0, 1, lsHiE, 1)
                zeros[C][E][i] = (zero)
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
            save += '_all-xstals'
            #save += '_globals-'+globstr
            save += '_c'+chan
            save += '_e'+str(E)
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
                            #isave += '-cs'+str(fitchans)
                            isave += '-c'+str(chan)
                            isave += '-e'+str(E)
                            #isave += '-ext'+str(extend)
                            isave += '-oth'+str(others)
                            if note: isave += '_'+str(note)
                            isave += '_'+str(V)

                            sepPlots[C][E][i].Print(plotdir+'/'+isave+'.png')

                        except:
                            print('WARNING: could not make plot for xstal-{0}'.format(i+1))
                            continue

    
    ### Save 4 (or 6) plots to one canvas
    #-------------------------------------------------------------
    if indi:
        combPlots = [0 for x in range(numx)]
        #tpad = [[0 for x in range(numx)] for x in range(4)]
        #bpad = [[0 for x in range(numx)] for x in range(4)]
        # extend for alphas
        tpad = [[0 for x in range(numx)] for x in range(6)]
        bpad = [[0 for x in range(numx)] for x in range(6)]
        
        #for i in range(numx):
        for i in [j-1 for j in justthese]:
            #if i+1 in justthese:
                #try:
            combPlots[i] = TCanvas('ccan-'+str(i+1),'ccan-'+str(i+1),0,0,1400,900)
            """
            if showAlphas:
                combPlots[i].Divide(3,2)
                Es = 3
            else:
                combPlots[i].Divide(2,2)
                Es = 2
            """
            combPlots[i].Divide(2,2)
            Es = 2
            
            p=0
            for C, chan in enumerate(pltchans):
                #for E in range(numE):
                # don't plot alphas
                for E in range(Es):
                    tpad[p][i] = toppad[C][E][i].Clone()
                    bpad[p][i] = botpad[C][E][i].Clone()
                    combPlots[i].cd(p+1)
                    tpad[p][i].Draw()
                    bpad[p][i].Draw()
                    combPlots[i].Update()
                    p += 1
            csave  = ''
            csave += 'x'+str(i+1)
            csave += '_all-channels'
            #csave += '_globals-'+globstr
            #csave += '_ext'+str(extend)
            csave += '_oth'+str(others)
            if note: csave += '_'+str(note)
            csave += '_'+str(V)

            combPlots[i].Print(plotdir+'/'+csave+'.png')
                    
                #except:
                #    print('WARNING: could not make plot for xstal-{0}'.format(i+1))
                #    continue


    ### Plot the 0-40 keV region of single hit
    #-------------------------------------------------------------       
    if indi:
        zsh_canv = [0 for x in range(numx)]
        zsh_hist = [0 for x in range(numx)]
        zsh_data = [0 for x in range(numx)]
        zsh_total = [0 for x in range(numx)]
        #zsh_resid = [0 for x in range(numx)]
        zsh_resid = makeResid100('S', 0, params[0])
        zsh_line = [0 for x in range(numx)]
        zsh_top = [0 for x in range(numx)]
        zsh_bot = [0 for x in range(numx)]
        
        font = 43
        size = 24
        loff = 0.005
        xoff = 3.6
        yoff = 1.1
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:

                zsh_canv[i] = TCanvas('zsh_canv-'+str(i+1),'zsh_canv-'+str(i+1),
                                0, 0, 1200, 700)

                fraction = 0.3
                pad1 = TPad('pad1_'+str(i),'pad1_'+str(i),0,fraction,1,1)
                #zsh_top[i].append(pad1)
                zsh_top[i] = (pad1)
                pad2 = TPad('pad2_'+str(i),'pad2_'+str(i),0,0,1,fraction)
                #zsh_bot[i].append(pad2)
                zsh_bot[i] = (pad2)
                
                zsh_top[i].SetBottomMargin(0.01)
                zsh_top[i].SetBorderMode(0)
                
                #zsh_bot[i].SetLogy(1)
                zsh_top[i].SetLogy(0)
                
                zsh_bot[i].SetTopMargin(0.05)
                zsh_bot[i].SetBottomMargin(0.3)
                zsh_bot[i].SetBorderMode(0)
                
                #zsh_bot[i].SetLogy(1)
                zsh_bot[i].SetLogy(0)
                
                zsh_top[i].Draw()
                zsh_bot[i].Draw()

                zsh_top[i].cd()

                #gStyle.SetPadLeftMargin(0.09)

                ### just checking residual plot names
                #prim_list = botpad[0][0][i].GetListOfPrimitives()
                #for bla in prim_list:
                #    print bla.GetName(), bla

                prim_list = toppad[0][0][i].GetListOfPrimitives()
                for bla in prim_list:
                    #print bla.GetName(), bla

                    if 'data' in bla.GetName():
                        zsh_data[i] = (toppad[0][0][i].GetPrimitive(bla.GetName()))
                    if 'total' in bla.GetName():
                        zsh_total[i] = (toppad[0][0][i].GetPrimitive(bla.GetName()))
                    
                    try:
                        zsh_hist[i] = (toppad[0][0][i].GetPrimitive(bla.GetName()))
                        
                        zsh_hist[i].SetXTitle('energy (keV)')
                        zsh_hist[i].GetXaxis().SetTitleFont(font)
                        zsh_hist[i].GetXaxis().SetTitleSize(size)
                        zsh_hist[i].GetXaxis().SetTitleOffset(xoff)
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

                        if   i+1 == 1: zsh_hist[i].SetAxisRange(0,  4.0, 'y')
                        elif i+1 == 2: zsh_hist[i].SetAxisRange(0,  4.0, 'y')
                        elif i+1 == 3: zsh_hist[i].SetAxisRange(0,  4.0, 'y')
                        elif i+1 == 4: zsh_hist[i].SetAxisRange(0,  4.0, 'y')
                        elif i+1 == 5: zsh_hist[i].SetAxisRange(0, 10.0, 'y')
                        elif i+1 == 6: zsh_hist[i].SetAxisRange(0,  3.5, 'y')
                        elif i+1 == 7: zsh_hist[i].SetAxisRange(0,  3.5, 'y')
                        elif i+1 == 8: zsh_hist[i].SetAxisRange(0, 10.0, 'y')
                        else:          zsh_hist[i].SetAxisRange(0, 10.0, 'y')
                                                
                        if 'data' in bla.GetName():
                            zsh_hist[i].Draw('same')
                        else:
                            zsh_hist[i].Draw(dopt)
                            
                    except: continue
                
                try:
                    #leg = TLegend(toppad[0][0][i].GetPrimitive(bla.GetName()))
                    #leg.Draw("same")
                    legs[0][0][i].Draw('same')
                except: continue

                zsh_bot[i].cd()
                
                #databins = zsh_data[i].GetNbinsX()
                #totalbins = zsh_total[i].GetNbinsX()
                #residbins = zsh_resid[i].GetNbinsX()
                #print databins, totalbins, residbins
                #zsh_resid[i].Rebin(residbins/databins)
                zsh_resid[i].Rebin(loEplotRebin)
                
                #zsh_resid[i] = zsh_data[i]/zsh_total[i]
                #zsh_resid[i] = zsh_data[i] - zsh_total[i]
                
                ### in some cases I may not be plotting the total
                ### so do a 'try'
                if zdr:
                    try:
                        zsh_resid[i].Divide(zsh_data[i], zsh_total[i])
                    except:
                        print 'WARNING: could not divide data by total'
                        pass
                else:
                    try:
                        zsh_resid[i] = zsh_data[i] - zsh_total[i]
                    except:
                        print 'WARNING: could not subtract data by total'
                        pass
                                
                ### tweak the errors on the residual
                for n in range(zsh_resid[i].GetNbinsX()+1):
                    R  = zsh_resid[i].GetBinContent(n)
                    D  = zsh_data[i].GetBinContent(n)
                    sD = zsh_data[i].GetBinError(n)

                    if zdr: zsh_resid[i].SetBinError(n, R*(sD/D))
                    else:   zsh_resid[i].SetBinError(n, sD)
                    
                    #try:
                    #    zsh_resid[i].SetBinError(n, R*(sD/D))
                    #except:
                    #    zsh_resid[i].SetBinError(n, 0)
                

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

                
                ### make a line at '1' or '0'
                if zdr: zero = TLine(0, 1, zmaxE, 1)
                else:   zero = TLine(0, 0, zmaxE, 0)
                zsh_line[i] = zero
                zsh_line[i].SetLineColor(kRed)
                zsh_line[i].SetLineWidth(1)

                zsh_resid[i].Draw()
                zsh_line[i].Draw('same')
                
                zsh_canv[i].Update()
                zsh_canv[i].Print(plotdir+'/a_singleHitZoom_c'+str(i+1)+'.png')

    
    ### Plot the 0-40 keV region of multi hit
    #-------------------------------------------------------------       
    if indi:
        zmh_canv = [0 for x in range(numx)]
        zmh_hist = [0 for x in range(numx)]
        zmh_data = [0 for x in range(numx)]
        zmh_total = [0 for x in range(numx)]
        zmh_resid = makeResid100('M', 0, params[0])
        zmh_line = [0 for x in range(numx)]
        zmh_top = [0 for x in range(numx)]
        zmh_bot = [0 for x in range(numx)]
        
        font = 43
        size = 24
        loff = 0.005
        xoff = 3.6
        yoff = 1.1
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:
                
                zmh_canv[i] = TCanvas('zmh_canv-'+str(i+1),'zmh_canv-'+str(i+1),
                                0, 0, 1200, 700)

                fraction = 0.3
                pad1 = TPad('pad1_'+str(i),'pad1_'+str(i),0,fraction,1,1)
                #zmh_top[i].append(pad1)
                zmh_top[i] = (pad1)
                pad2 = TPad('pad2_'+str(i),'pad2_'+str(i),0,0,1,fraction)
                #zmh_bot[i].append(pad2)
                zmh_bot[i] = (pad2)
                
                zmh_top[i].SetBottomMargin(0.01)
                zmh_top[i].SetBorderMode(0)
                
                #zmh_bot[i].SetLogy(1)
                zmh_top[i].SetLogy(0)
                
                zmh_bot[i].SetTopMargin(0.05)
                zmh_bot[i].SetBottomMargin(0.3)
                zmh_bot[i].SetBorderMode(0)
                
                #zmh_bot[i].SetLogy(1)
                zmh_bot[i].SetLogy(0)
                
                zmh_top[i].Draw()
                zmh_bot[i].Draw()

                zmh_top[i].cd()

                
                prim_list = toppad[1][0][i].GetListOfPrimitives()
                for bla in prim_list:
                    #print bla.GetName(), bla

                    if 'data' in bla.GetName():
                        zmh_data[i] = (toppad[1][0][i].GetPrimitive(bla.GetName()))
                    if 'total' in bla.GetName():
                        zmh_total[i] = (toppad[1][0][i].GetPrimitive(bla.GetName()))
                    
                    
                    try:
                        zmh_hist[i] = (toppad[1][0][i].GetPrimitive(bla.GetName()))

                        zmh_hist[i].SetXTitle('energy (keV)')
                        zmh_hist[i].GetXaxis().SetTitleFont(font)
                        zmh_hist[i].GetXaxis().SetTitleSize(size)
                        zmh_hist[i].GetXaxis().SetTitleOffset(xoff)
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

                        if   i+1 == 1: zmh_hist[i].SetAxisRange(0, 4.5, 'y')
                        elif i+1 == 2: zmh_hist[i].SetAxisRange(0, 4.5, 'y')
                        elif i+1 == 3: zmh_hist[i].SetAxisRange(0, 2.5, 'y')
                        elif i+1 == 4: zmh_hist[i].SetAxisRange(0, 2.0, 'y')
                        elif i+1 == 5: zmh_hist[i].SetAxisRange(0, 5.0, 'y')
                        elif i+1 == 6: zmh_hist[i].SetAxisRange(0, 1.2, 'y')
                        elif i+1 == 7: zmh_hist[i].SetAxisRange(0, 1.2, 'y')
                        elif i+1 == 8: zmh_hist[i].SetAxisRange(0, 5.0, 'y')
                        else:          zmh_hist[i].SetAxisRange(0, 5.0, 'y')

                        if 'data' in bla.GetName():
                            zmh_hist[i].Draw('same')
                        else:
                            zmh_hist[i].Draw(dopt)
                        
                    except: continue
                
                try:
                    #leg = TLegend(toppad[0][0][i].GetPrimitive(bla.GetName()))
                    #leg.Draw("same")
                    legs[1][0][i].Draw('same')
                except: continue

                zmh_bot[i].cd()
                
                #databins = zmh_data[i].GetNbinsX()
                #totalbins = zmh_total[i].GetNbinsX()
                #residbins = zmh_resid[i].GetNbinsX()
                #print databins, totalbins, residbins
                #zmh_resid[i].Rebin(residbins/databins)
                zmh_resid[i].Rebin(loEplotRebin)
                
                #zmh_resid[i] = zmh_data[i]/zmh_total[i]
                #zmh_resid[i] = zmh_data[i] - zmh_total[i]
                
                ### in some cases I may not be plotting the total
                ### so do a 'try'
                if zdr:
                    try:
                        zmh_resid[i].Divide(zmh_data[i], zmh_total[i])
                    except:
                        print 'WARNING: could not divide divide data by total'
                        pass
                else:
                    try:
                        zmh_resid[i] = zmh_data[i] - zmh_total[i]
                    except:
                        print 'WARNING: could not subtract data by total'
                        pass
                
                ### tweak the errors on the residual
                for n in range(zmh_resid[i].GetNbinsX()+1):
                    R  = zmh_resid[i].GetBinContent(n)
                    D  = zmh_data[i].GetBinContent(n)
                    sD = zmh_data[i].GetBinError(n)

                    if zdr: zmh_resid[i].SetBinError(n, R*(sD/D))
                    else:   zmh_resid[i].SetBinError(n, sD)
                    
                    #try:
                    #    zmh_resid[i].SetBinError(n, R*(sD/D))
                    #except:
                    #    zmh_resid[i].SetBinError(n, 0)
                

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

                
                ### make a line at '1' or '0'
                if zdr: zero = TLine(0, 1, zmaxE, 1)
                else:   zero = TLine(0, 0, zmaxE, 0)
                zmh_line[i] = zero
                zmh_line[i].SetLineColor(kRed)
                zmh_line[i].SetLineWidth(1)

                zmh_resid[i].Draw()
                zmh_line[i].Draw('same')
                
                zmh_canv[i].Update()
                zmh_canv[i].Print(plotdir+'/a_multiHitZoom_c'+str(i+1)+'.png')

    
    ### Plot the single hit alpha channel
    #-------------------------------------------------------------       
    if indi and showAlphas:
        ash_canv = [0 for x in range(numx)]
        ash_hist = [{} for x in range(numx)]
        ash_data = [0 for x in range(numx)]
        ash_total = makeTotal100('S', 2, params[2])
        ash_resid = makeResid100('S', 2, params[2])
        ash_line = [0 for x in range(numx)]
        ash_top = [0 for x in range(numx)]
        ash_bot = [0 for x in range(numx)]
        ash_leg = [0 for x in range(numx)]
        
        font = 43
        size = 24
        loff = 0.005
        xoff = 3.6
        yoff = 1.1
        
        alphaColors = [
            kBlue    +1,
            kCyan    +1,
            kGreen   +1,
            kOrange  +1,
            kMagenta +1,
            kBlue    -7,
            kCyan    -7,
            kGreen   -7,
            kOrange  -7,
            kMagenta -7,
        ]
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:
            #if i+1 in justthese:
                
                ash_canv[i] = TCanvas('ash_canv-'+str(i+1),'ash_canv-'+str(i+1),
                                0, 0, 1200, 700)
                
                fraction = 0.3
                pad1 = TPad('pad1_'+str(i),'pad1_'+str(i),0,fraction,1,1)
                ash_top[i] = (pad1)
                ash_top[i].SetBottomMargin(0.01)
                ash_top[i].SetBorderMode(0)
                ash_top[i].SetLogy(1)
                #ash_top[i].SetLogy(0)
                ash_top[i].Draw()
                
                pad2 = TPad('pad2_'+str(i),'pad2_'+str(i),0,0,1,fraction)
                ash_bot[i] = (pad2)
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
                        ash_data[i].GetXaxis().SetTitleOffset(xoff)
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
                        
                        ash_data[i].SetAxisRange(aler[0], aler[1], 'x')
                        ash_data[i].SetAxisRange(1e-5, 10, 'y')

                        ash_data[i].Draw('same')
                        ldatname = data[key]['info']['tag']+'  '+data[key]['info']['build']

                # join hists of same type together -fx -qx
                for key in bakkeys:
                    if 'x'+str(i+1) in key and '-cS' in key and '-e2' in key:
                        if justQ1 and '-q1-' not in key: continue
                        if bkgs[key]['hist'].GetEntries() > 0:
                            histkey = re.split('-f(?:[0-9])-', key)[0]
                            if histkey in ash_hist[i]:
                                ash_hist[i][histkey].Add(bkgs[key]['hist'])
                            else:
                                ash_hist[i][histkey] = (bkgs[key]['hist'])
                
                for key in sigkeys:
                    if 'x'+str(i+1) in key and '-cS' in key and '-e2' in key:
                        if justQ1 and '-q1-' not in key: continue
                        if sigs[key]['hist'].GetEntries() > 0:
                            histkey = re.split('-f(?:[0-9])-', key)[0]
                            if histkey in ash_hist[i]:
                                ash_hist[i][histkey].Add(sigs[key]['hist'])
                            else:
                                ash_hist[i][histkey] = (sigs[key]['hist'])
                
                N_hists = len(ash_hist[i])
                # split in half because of q1 and q2 quenching
                #N_hists = N_hists/2
                
                leg = TLegend(0.75, (1-((N_hists+2)*0.08)), 0.98, 0.90)
                ash_leg[i] = (leg)
                ash_leg[i].SetFillColor(0)
                ash_leg[i].SetBorderSize(0)
                legopt = 'LPE'
                
                ash_leg[i].AddEntry(ash_data[i], ldatname, legopt)

                histkeys = sorted([key for key in ash_hist[i]])
                for N, histkey in enumerate(histkeys):
                    ash_hist[i][histkey].SetMarkerColor(alphaColors[N])
                    ash_hist[i][histkey].SetLineColor(alphaColors[N])
                    ash_hist[i][histkey].Draw(dopt)
                    bits = histkey.split('-')
                    legname = bits[1]+' '+bits[2]
                    ash_leg[i].AddEntry(ash_hist[i][histkey], legname, legopt)
                    ash_total[i].Add(ash_hist[i][histkey])


                if showTotal:
                    if redtotal:
                        ash_total[i].SetMarkerColor(kRed)
                        ash_total[i].SetLineColor(kRed)
                    ash_total[i].Draw(dopt)
                    ash_leg[i].AddEntry(ash_total[i], "Total", legopt)
                    
                ash_leg[i].Draw('same')
                
                
                # plot residuals
                ash_bot[i].cd()
                ash_resid[i].Rebin(alphaPlotRebin)
                
                # remove very small bins or the Divide returns 'inf'
                ash_total[i] = cleanSmallBins(ash_total[i])
                ash_resid[i].Divide(ash_data[i], ash_total[i])
                
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

                ash_resid[i].SetAxisRange(aler[0], aler[1], 'x')
                ash_resid[i].SetAxisRange(0.0, 2, 'y')
                                
                ash_resid[i].SetXTitle('energy (keV)')
                ash_resid[i].GetXaxis().SetTitleFont(font)
                ash_resid[i].GetXaxis().SetTitleSize(size)
                ash_resid[i].GetXaxis().SetTitleOffset(xoff)
                ash_resid[i].GetXaxis().SetLabelFont(font)
                ash_resid[i].GetXaxis().SetLabelSize(size)
                ash_resid[i].GetXaxis().SetLabelOffset(loff)
                
                ash_resid[i].SetYTitle('data/mc')
                ash_resid[i].GetYaxis().SetTitleFont(font)
                ash_resid[i].GetYaxis().SetTitleSize(size)
                ash_resid[i].GetYaxis().SetTitleOffset(yoff)
                ash_resid[i].GetYaxis().SetLabelFont(font)
                ash_resid[i].GetYaxis().SetLabelSize(size)
                ash_resid[i].GetYaxis().SetLabelOffset(loff)
                
                # '5' secondary and '05' primary
                # means 5 divisions will be shown
                ash_resid[i].GetYaxis().SetNdivisions(505)

                
                ### make a line at '1'
                zero = TLine(aler[0], 1, aler[1], 1)
                ash_line[i] = zero
                ash_line[i].SetLineColor(kRed)
                ash_line[i].SetLineWidth(1)

                ash_resid[i].Draw()
                ash_line[i].Draw('same')
                
                ash_canv[i].Update()
                ash_canv[i].Print(plotdir+'/a_alphaSingleHit_c'+str(i+1)+'.png')

                
    ### Plot the multi hit alpha channel
    #-------------------------------------------------------------       
    if indi and showAlphas:
        amh_canv = [0 for x in range(numx)]
        amh_hist = [{} for x in range(numx)]
        amh_data = [0 for x in range(numx)]
        amh_total = makeTotal100('M', 2, params[2])
        amh_resid = makeResid100('M', 2, params[2])
        amh_line = [0 for x in range(numx)]
        amh_top = [0 for x in range(numx)]
        amh_bot = [0 for x in range(numx)]
        amh_leg = [0 for x in range(numx)]
        
        font = 43
        size = 24
        loff = 0.005
        xoff = 3.6
        yoff = 1.1
        
        alphaColors = [
            kBlue    +1,
            kCyan    +1,
            kGreen   +1,
            kOrange  +1,
            kMagenta +1,
            kBlue    -7,
            kCyan    -7,
            kGreen   -7,
            kOrange  -7,
            kMagenta -7,
        ]
        
        for i in range(numx):
            if i+1 in justthese and i+1 != 9:
            #if i+1 in justthese:
                
                amh_canv[i] = TCanvas('amh_canv-'+str(i+1),'amh_canv-'+str(i+1),
                                0, 0, 1200, 700)
                
                fraction = 0.3
                pad1 = TPad('pad1_'+str(i),'pad1_'+str(i),0,fraction,1,1)
                amh_top[i] = (pad1)
                amh_top[i].SetBottomMargin(0.01)
                amh_top[i].SetBorderMode(0)
                amh_top[i].SetLogy(1)
                #amh_top[i].SetLogy(0)
                amh_top[i].Draw()
                
                pad2 = TPad('pad2_'+str(i),'pad2_'+str(i),0,0,1,fraction)
                amh_bot[i] = (pad2)
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
                        amh_data[i].GetXaxis().SetTitleOffset(xoff)
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
                        
                        amh_data[i].SetAxisRange(aler[0], aler[1], 'x')
                        amh_data[i].SetAxisRange(1e-7, 1e-2, 'y')

                        amh_data[i].Draw('same')
                        ldatname = data[key]['info']['tag']+'  '+data[key]['info']['build']

                # join hists of same type together -fx -qx
                for key in bakkeys:
                    if 'x'+str(i+1) in key and '-cM' in key and '-e2' in key:
                        if justQ1 and '-q1-' not in key: continue
                        if bkgs[key]['hist'].GetEntries() > 0:
                            histkey = re.split('-f(?:[0-9])-', key)[0]
                            if histkey in amh_hist[i]:
                                amh_hist[i][histkey].Add(bkgs[key]['hist'])
                            else:
                                amh_hist[i][histkey] = (bkgs[key]['hist'])
                
                for key in sigkeys:
                    if 'x'+str(i+1) in key and '-cM' in key and '-e2' in key:
                        if justQ1 and '-q1-' not in key: continue
                        if sigs[key]['hist'].GetEntries() > 0:
                            histkey = re.split('-f(?:[0-9])-', key)[0]
                            if histkey in amh_hist[i]:
                                amh_hist[i][histkey].Add(sigs[key]['hist'])
                            else:
                                amh_hist[i][histkey] = (sigs[key]['hist'])
                
                N_hists = len(amh_hist[i])
                # split in half because of q1 and q2 quenching
                #N_hists = N_hists/2
                
                leg = TLegend(0.75, (1-((N_hists+2)*0.08)), 0.98, 0.90)
                amh_leg[i] = (leg)
                amh_leg[i].SetFillColor(0)
                amh_leg[i].SetBorderSize(0)
                legopt = 'LPE'
                
                amh_leg[i].AddEntry(amh_data[i], ldatname, legopt)

                histkeys = sorted([key for key in amh_hist[i]])
                for N, histkey in enumerate(histkeys):
                    amh_hist[i][histkey].SetMarkerColor(alphaColors[N])
                    amh_hist[i][histkey].SetLineColor(alphaColors[N])
                    amh_hist[i][histkey].Draw(dopt)
                    bits = histkey.split('-')
                    legname = bits[1]+' '+bits[2]
                    amh_leg[i].AddEntry(amh_hist[i][histkey], legname, legopt)
                    amh_total[i].Add(amh_hist[i][histkey])


                if showTotal:
                    if redtotal:
                        amh_total[i].SetMarkerColor(kRed)
                        amh_total[i].SetLineColor(kRed)
                    amh_total[i].Draw(dopt)
                    amh_leg[i].AddEntry(amh_total[i], "Total", legopt)
                    
                amh_leg[i].Draw('same')
                
                
                # plot residuals
                amh_bot[i].cd()
                amh_resid[i].Rebin(alphaPlotRebin)
                
                # remove very small bins or the Divide returns 'inf'
                amh_total[i] = cleanSmallBins(amh_total[i])
                amh_resid[i].Divide(amh_data[i], amh_total[i])
                
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

                amh_resid[i].SetAxisRange(aler[0], aler[1], 'x')
                amh_resid[i].SetAxisRange(0.0, 2, 'y')
                                
                amh_resid[i].SetXTitle('energy (keV)')
                amh_resid[i].GetXaxis().SetTitleFont(font)
                amh_resid[i].GetXaxis().SetTitleSize(size)
                amh_resid[i].GetXaxis().SetTitleOffset(xoff)
                amh_resid[i].GetXaxis().SetLabelFont(font)
                amh_resid[i].GetXaxis().SetLabelSize(size)
                amh_resid[i].GetXaxis().SetLabelOffset(loff)
                
                amh_resid[i].SetYTitle('data/mc')
                amh_resid[i].GetYaxis().SetTitleFont(font)
                amh_resid[i].GetYaxis().SetTitleSize(size)
                amh_resid[i].GetYaxis().SetTitleOffset(yoff)
                amh_resid[i].GetYaxis().SetLabelFont(font)
                amh_resid[i].GetYaxis().SetLabelSize(size)
                amh_resid[i].GetYaxis().SetLabelOffset(loff)
                
                # '5' secondary and '05' primary
                # means 5 divisions will be shown
                amh_resid[i].GetYaxis().SetNdivisions(505)

                
                ### make a line at '1'
                zero = TLine(aler[0], 1, aler[1], 1)
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
    
    try: del fcanv
    except: pass
    
    try: del mcanv
    except: pass
    
    try: del sepFitPlot
    except: pass
    
    #try: del combPlots
    #except: pass

    try: del zsh_canv
    except: pass
    
    try: del zmh_canv
    except: pass

    try: del ash_canv
    except: pass
    
    try: del amh_canv
    except: pass
    
    try: del sepPlots
    except: pass
    
    try: del canvs
    except: pass
    
    if fitting:
        print('Time to complete the fit = {0} sec \n'.format(fitTime))
    
    if not batch: input('[Enter] to quit \n')
    else: print('Running in batch mode - quitting now! \n')

    return


if __name__ == "__main__":
    main(sys.argv[1:])

