#!/usr/bin/env python

######################################################################
# Matt Kauer - mkauer@physics.wisc.edu
######################################################################
# 201-global-bkgs.py

V = 'v201'

# Add global MC that are backgrounds and not signals - deep bug fix
# 
# version: 2018-12-14
# 
# see CHANGELOG for changes
######################################################################

import os,sys,re
import shutil
import socket
import copy
import math
import numpy as np

import ROOT
from ROOT import *
ROOT.gROOT.Reset()
#ROOT.gErrorIgnoreLevel = kWarning
ROOT.gErrorIgnoreLevel = kError

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs101 import *


### get the number of crystals
numx = numX()

### print debug info
debug = 0


### ==========  GENERAL INPUTS  ======================================
### note to add to saved plot names?
note = 0
#note = 'calib-130'

#mcfile = 'backgrounds_200.txt'
#mcfile = 'backgrounds_201.txt'
#mcfile = 'backgrounds_202.txt' # lsveto
#mcfile = 'backgrounds_203.txt' # plastic
mcfile = 'backgrounds_204.txt' # lsveto again

#mcfile = 'backgrounds_202-update.txt'
#mcfile = 'lsveto-testing.txt'

print 'INFO: using backgrounds config file -->', mcfile


### ==========  OPTIMIZATION OPTIONS  ================================
### MC smoothing? Specify a smoothing window in +/- number of bins
#smoothing = 5
smoothing = 0


### ==========  FITTING OPTIONS  =====================================
### select channels to fit
fitchans = 'SM'

### Let's finally try different fit ranges!!
fitranges = [{} for x in range(numx)]
for i in range(numx):
    ### format = [rebin, xmin, xmax]
    fitranges[i]['S0'] = [3,  6,   76]  # single-hit low-energy
    fitranges[i]['S1'] = [4, 70, 2770]  # single-hit high-energy
    fitranges[i]['M0'] = [3,  2,   72]  # multi-hit low-energy
    fitranges[i]['M1'] = [4, 70, 2770]  # multi-hit high-energy
    ### separate range for lsveto
    if i == 8:
        fitranges[i]['S0'] = [0,0,0]
        fitranges[i]['S1'] = [0,0,0]
        fitranges[i]['M0'] = [0,0,0]
        fitranges[i]['M1'] = [4, 300, 2700]

        
### ==========  EXTRA MC OPTIONS  ====================================
### which MC to fit globally (to all crystals simultaneously)?
#globalmc = []
#globalmc = ['plastic', 'lsveto', 'innersteel', 'steel']
globalmc = ['pmt', 'plastic', 'lsveto', 'innersteel', 'steel']
#globalmc = ['copper', 'pmt', 'plastic', 'lsveto', 'innersteel', 'steel']
#globalmc = ['pmt', 'plastic', 'lsveto', 'innersteel', 'steel',
#            'nai-surf-10um', 'teflon-bulk', 'internal', 'copper']

### include bkgs from 'other' pmts and internals?
others  = 1

### plot components in groups? [0,1]
ingroups = 1

### show the total? [0,1]
showTotal = 1

### show the legends? [0,1]
showlegs = 1

### plot the total in red? [0,1]
redtotal = 1

### combine 'others' into the makePlots() plots?
combine = 1

### force the reuse of all joined rootfiles in mcfile? [0,1,2]
### very nice for debugging
### [0] default - use whatever is specified in the backgrounds file
### [1] forces reusing of all data/bkgs/sigs
### [2] forces NOT reusing any data/bkgs/sigs
reuse = 0

### update and save new backgrounds file with fit results
updateMCfile = 1


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
hier = [0, 3000]
eran = [loer, hier]

### special energy range for lsveto
lsHiE = 3500

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
datsumw2 = 1

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
    
    batch = 0
    if onCup(): batch = 1
    
    #-----------------------------------------
    # set as an empty list for default action
    #-----------------------------------------
    xstals = []
    
    if len(argv) > 0:
        batch = 0
        xstals = []
        for c in argv:
            xstals.append(int(c))
    #-----------------------------------------
    #-----------------------------------------
    
    gROOT.SetBatch(batch)
    gStyle.SetPalette (1)
    gStyle.SetOptStat ('')
    gStyle.SetOptFit  (0)
    
    ### for saving the plots...
    if not os.path.exists('./plots'): 
        os.makedirs('./plots')
    
    if not os.path.exists(mcfile):
        print '\nERROR: could not find backgrounds file -->', mcfile
        sys.exit()
    
    
    ### where everything gets loaded into dictionary
    #-----------------------------------------------------------------
    #-----------------------------------------------------------------
    allchans = uniqString(fitchans+pltchans)
    data, bkgs, sigs, runtime = build101(mcfile, others, reuse, allchans, xstals)
    print 'INFO: runtime =', runtime, '(seconds)'

    datkeys = sortDataKeys92(data)
    if datsumw2:
        for key in datkeys:
            data[key]['hist'].Sumw2()
    
    ### scale into dru units
    if dru:
        data = scaleData70(data, 1)
        bkgs = scaleBkgs101(bkgs)
        sigs = scaleBkgs101(sigs)
    else:
        data = scaleData70(data, 0)
        bkgs = scaleBkgs101(bkgs, runtime)
        sigs = scaleBkgs101(sigs, runtime)

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
        bkgs = smooth(bkgs, smoothing)
        sigs = smooth(sigs, smoothing)
        print 'INFO: done smoothing histograms'
    
    #-----------------------------------------------------------------
    #-----------------------------------------------------------------
    
    ### plot all crystals that have data
    justthese = []
    for i in range(1, numx+1):
        for key in datkeys:
            if 'x'+str(i) in key and i not in justthese:
                justthese.append(i)
    justthese.sort()
    print 'INFO: plotting crystals -->', justthese
    
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
    print 'INFO: global fits to -->', globstr
    
    ### Number of colors
    Nc = len(uniqAll)
    #print 'INFO: Total number of unique bkgs and sigs =',Nc
    colors, cis = rainbow(Nc)

    ### Create color dict for unique simulations
    uniqColor = {}
    for i, key in enumerate(uniqAll):
        #print key
        uniqColor[key] = cis[i]

    ### colors for the groups
    gis = {'steel':    kYellow,
           'cosmo':    kMagenta+1,
           'lsveto':   kOrange+1,
           'surface':  kCyan+1,
           'internal': kBlue,
           'pmts':     kGreen+1,
           'copper':   kYellow+1,
           'plastic':  kOrange+2,
           'none':     kRed-1}

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
    
    if len(sigs) > 0:
        
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
            xstalbins[i] = [startbin, stopbin]
        
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
                        fitsigs[fskey]['hist'].Scale(dat_int/mc_int) # scale to data integral
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

                    bounds.append(sigs[newkey]['info']['newfbnd'])
                    #---------------------------------------------------------------------
                    
                    ### set errors to zero?
                    fitsigs[fskey]['hist'] = zeroBinError(fitsigs[fskey]['hist'])

            
            fitresults[str(i)].append('Crystal-'+str(i+1)+' fit results')
            fitresults[str(i)].append('runtime = '+str(round(runtime/60./60./24., 2))+' days')
            if note: fitresults[str(i)].append('note = '+note)
            fitresults[str(i)].append('version = '+V)
            fitresults[str(i)].append('channels fit = '+fitchans)
            fitresults[str(i)].append('global fits = '+globstr)
            fitresults[str(i)].append('other pmts = '+str(others))
            #fitresults[str(i)].append('hist extend = '+str(extend))
            fitresults[str(i)].append('norm to dru = '+str(dru))
            for key in ['S0', 'S1', 'M0', 'M1']:
                fitresults[str(i)].append(key+' fit range = '\
                                          +str(fitranges[i][key][1])\
                                          +' - '\
                                          +str(fitranges[i][key][2])\
                                          +' keV')

        
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
                fitglobsigs[fgkey]['hist'].Scale(dat_int/mc_int) # scale to data integral
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
                                bounds.append(sigs[newkey]['info']['newfbnd'])
            
            
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
        
        ### set fit bounds!!!
        ### l=0 sets all params to the same constrain
        ### set all bounds by default
        #fit.Constrain(0, 0.0, 1.0)
        for l in range(len(bounds)):
            fit.Constrain(l+1, bounds[l][0], bounds[l][1])
        
        ### set the fit range
        #fit.SetRangeX(0, fmax*numx)
        fit.SetRangeX(0, fitbins)
        
        ### do the fit
        status = fit.Fit()
        if status != 0:
            print '\n\n*******************  FIT FAILURE  *******************\n\n'
            sys.exit()
        
        print '\n\n*******************  SUCCESSFUL FIT  *******************\n\n'
        
        chi2 = fit.GetChisquare()
        ndf  = fit.GetNDF()
        pval = fit.GetProb()

        ### get non-zero ndf
        NDF=0
        for n in range(fitbins):
            if fitdata.GetBinContent(n) > 0:
                NDF+=1
        ndf=NDF
        
        count = 0
        for i in range(numx):
            
            fitresults[str(i)].append('total number of hists being fit = '+str(totalNumFits))
            fitresults[str(i)].append('returned fit status = '+str(status))
            fitchi2ndf = (chi2/ndf)
            fitresults[str(i)].append('chi2/ndf = %.3g/%s = %.3g'%(chi2,ndf,chi2/ndf))
            
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    
                    fscale = ROOT.Double(0.0)
                    ferror = ROOT.Double(0.0)
                    #print 'count',count
                    fit.GetResult(count, fscale, ferror)
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

            fscale = ROOT.Double(0.0)
            ferror = ROOT.Double(0.0)
            #print 'count',count
            fit.GetResult(count, fscale, ferror)
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
        if dru: sigs = scaleSigs101(sigkeys, sigs)
        else:   sigs = scaleSigs101(sigkeys, sigs, runtime)
        
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
                                        'fit '+fskey+' = %.2e +/- %.2e mBq  (%.2e, %.2e)  %s'
                                        %(fitacti, fiterro, lobnd, hibnd, limit))
                                else:
                                    fitresults[str(i)].append(
                                        'fit '+fskey+' = %.2e +/- %.2e mBq  (%.2e, %.2e)'
                                        %(fitacti, fiterro, lobnd, hibnd))
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
                                        'fit '+'x'+X+'-'+fgkey+' = %.2e +/- %.2e mBq  (%.2e, %.2e)  %s'
                                        %(fitacti, fiterro, lobnd, hibnd, limit))
                                else:
                                    fitresults[str(i)].append(
                                        'fit '+'x'+X+'-'+fgkey+' = %.2e +/- %.2e mBq  (%.2e, %.2e)'
                                        %(fitacti, fiterro, lobnd, hibnd))
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
        

        ### only write out files if fit is successful
        #-------------------------------------------------------------
        if status == 0:
            
            ### sort the fitresults keys
            resultskeys = []
            for rskey in fitresults:
                resultskeys.append(rskey)
            resultskeys.sort()
                    
            ### write results to file
            resultsfile = './plots/'+save+'_fit-results.txt'
            outfile = open(resultsfile, 'w')
            for key in resultskeys:
                if int(key)+1 in justthese:
                    for line in fitresults[key]:
                        outfile.write(line+'\n')
            outfile.close()
            
            ### create the updated backgrounds file
            #shutil.copyfile(mcfile, './plots/'+mcfile)
            if updateMCfile:
                newbkgs = './plots/'+mcfile[:-4]+'-update.txt'
                updateBkgsFile70(mcfile, resultsfile, newbkgs, BF='BR')
            
            ### save histograms to a rootfile
            rootoutfile = TFile("./plots/histograms.root", "RECREATE")
            for key in sigkeys:
                sigs[key]['hist'].Write(key)
            for key in bakkeys:
                bkgs[key]['hist'].Write(key)
            for key in data:
                data[key]['hist'].Write(key)
            rootoutfile.Write()
            rootoutfile.Close()
            
            ### create the background model table
            #outtable = newbkgs[:-4]+'-table.txt'
            #outputModelTable61(newbkgs, outtable)
        #-------------------------------------------------------------
        
        
        
        # plot the multi-chan fit results
        #=================================================================

        fcanv = TCanvas('fcanv', 'fcanv', 0, 0, 1400, 900)
        #fcanv.Divide(4,2)
        fcanv.Divide(5,2)

        ftoppad = []
        fbotpad = []

        sepfitdata = [0 for x in range(numx)]
        sepfitsigs = [0 for x in range(numx)]
        sepfitglobsigs = [0 for x in range(numx)]
        
        #sepFitPlots = []
        
        flegs = []
        flegs2 = []

        fzeros = []
        
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
            ftoppad.append(pad1)
            pad2 = TPad('pad2','pad2',0,0,1,fraction)
            fbotpad.append(pad2)
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
            
            newFitTitle = str('Crystal-'+str(i+1)+'   '+'Fit-chans-'+fitchans)
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
            #fitdata.GetXaxis().SetTitle('Energy (keV)')
            #fitdata.GetXaxis().SetLabelFont(font)
            #fitdata.GetXaxis().SetLabelSize(size)
            
            fitdata.Draw()
            #flegs[i].AddEntry(fitdata, 'data - bkgs', lopt)
            if dru: fitdata.SetAxisRange(2e-3, 2e3, 'y')
            fitdata.SetAxisRange(xstalbins[i][0], xstalbins[i][1], 'x')
            
            
            Nfsigs=0
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:

                    Nfsigs+=1
                    
                    # find the unique name for color and set color
                    cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
                    fitsigs[fskey]['hist'].SetMarkerColor(uniqColor[cname])
                    fitsigs[fskey]['hist'].SetLineColor(uniqColor[cname])
                    
                    ### draw the sigs
                    fitsigs[fskey]['hist'].Draw('same')

                    ### add MC to total MC hist
                    #print '-3-', ftotal.GetNbinsX(), fitsigs[fskey]['hist'].GetNbinsX(), fskey
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
                fitglobsigs[fgkey]['hist'].Draw('same')

                if i==0:
                    ### add MC to total MC hist
                    #print '-4-', ftotal.GetNbinsX(), fitglobsigs[fgkey]['hist'].GetNbinsX(), fgkey
                    ftotal.Add(fitglobsigs[fgkey]['hist'])
                
                
            ### select the right range
            ftotal.SetAxisRange(xstalbins[i][0], xstalbins[i][1], 'x')
            
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
            flegs.append(fleg)
            flegs[i].SetNColumns(flnc)
            flegs[i].SetFillColor(0)
            flegs[i].SetBorderSize(0)
            lopt = 'LPE'
            
            flegs[i].AddEntry(fitdata, 'data - bkgs', lopt)
            
            ### add legend entries in order
            for name in uniqAll:
                for fskey in fsigkeys:
                    if name in fskey and 'x'+str(i+1) in fskey:
                        flegs[i].AddEntry(fitsigs[fskey]['hist'], fskey, lopt)
                # and for globals
                for fgkey in fglobsigkeys:
                    if name in fgkey:
                        flegs[i].AddEntry(fitglobsigs[fgkey]['hist'], fgkey, lopt)
            
            ### get the chi2 of the total fit mc compared to data
            #-------------------------------------------------------------------
            """
            chi2  = ROOT.Double(0.0)
            ndf   = ROOT.Long(0)
            igood = ROOT.Long(0)
            fitchi2ndfv2 = -1
            #if i in wasFit:
            pval = fitdata.Chi2TestX(ftotal, chi2, ndf, igood, chiopt)
            #fitchi2ndfv2 = chi2/druscale[i]/ndf
            fitchi2ndfv2 = chi2/ndf
            """
            #-------------------------------------------------------------------
            
            ftotal.Draw('same')
            
            #if i in wasFit:
            # returned from fit
            flegs[i].AddEntry(ftotal, 'Fit Total (chi2/ndf = '+str(round(fitchi2ndf,2))+')', lopt)
            # calc by Chi2TestX()
            #flegs[i].AddEntry(ftotal, 'Fit Total (chi2/ndf = '+str(round(fitchi2ndfv2,2))+')', lopt)
            #print 'INFO: Fit total MC from Chi2TestX chi2/ndf = '+str(round(fitchi2ndfv2,2))
            
            flegs[i].Draw('same')
            
            
            ### plot the fit residuals
            #-------------------------------------------------------------
            fbotpad[i].cd()
            leg = TLegend(0.72, 0.78, 0.94, 0.94)
            flegs2.append(leg)
            flegs2[i].SetFillColor(0)
            flegs2[i].SetBorderSize(0)
            lopt = 'LPE'

            fresid.Divide(fitdata, ftotal)
            #fresid.Divide(ftotal, fitdata)

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
            
            fresid.SetAxisRange(0.1,10,'y')

            fresid.SetAxisRange(xstalbins[i][0], xstalbins[i][1], 'x')

            if linres:
                fbotpad[i].SetLogy(0)
                fresid.SetAxisRange(lrs[0], lrs[1], 'y')
            fresid.Draw()
            
            ### set my reference line to '1'
            #zero = TLine(fmin, 1, fmax, 1)
            #fzeros.append(zero)
            #fzeros[i].SetLineColor(kRed)
            #fzeros[i].SetLineWidth(1)
            
            zero = TLine(xstalbins[i][0], 1, xstalbins[i][1], 1)
            #fzeros.SetAxisRange(i*fmax, (i+1)*fmax, 'x')
            fzeros.append(zero)
            fzeros[i].SetLineColor(kRed)
            fzeros[i].SetLineWidth(1)
            fzeros[i].Draw()

            flegs2[i].AddEntry(fresid,'data / MC',lopt)
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
                
                    sepFitPlot.Print(str('./plots/'+fisave+'.png'))
                except:
                    print 'WARNING: could not make plot for xstal-'+str(i+1)
                    continue
            
        fcanv.Update()
        fcanv.Print('./plots/'+save+'.png')

        
        ### plot the MEGA combined histograms for fun!
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
        #flegs[i].AddEntry(fitdata, 'data - bkgs', lopt)
        if dru: fitdata.SetAxisRange(2e-3, 2e3, 'y')
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
            fitsigs[fskey]['hist'].Draw('same')

            ### add MC to total MC hist
            #ftotal.Add(fitsigs[fskey]['hist'])

        for fgkey in fglobsigkeys:
            #if i==0:
            Nfsigs+=1

            # find the unique name for color and set color
            #cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
            cname = fgkey
            fitglobsigs[fgkey]['hist'].SetMarkerColor(uniqColor[cname])
            fitglobsigs[fgkey]['hist'].SetLineColor(uniqColor[cname])

            ### draw the sigs
            fitglobsigs[fgkey]['hist'].Draw('same')

            #if i==0:
                ### add MC to total MC hist
                #ftotal.Add(fitglobsigs[fgkey]['hist'])

                
        ### select the right range
        ftotal.SetAxisRange(0, fitbins, 'x')
        ftotal.Draw('same')
        
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
        lopt = 'LPE'

        mleg.AddEntry(fitdata, 'data - bkgs', lopt)

        ### add legend entries in order
        for name in uniqAll:
            for fskey in fsigkeys:
                if name in fskey:
                    mleg.AddEntry(fitsigs[fskey]['hist'], fskey, lopt)
            # and for globals
            for fgkey in fglobsigkeys:
                if name in fgkey:
                    mleg.AddEntry(fitglobsigs[fgkey]['hist'], fgkey, lopt)

        #ftotal.Draw('same')

        #if i in wasFit:
        # returned from fit
        mleg.AddEntry(ftotal, 'Fit Total (chi2/ndf = '+str(round(fitchi2ndf,2))+')', lopt)
        # calc by Chi2TestX()
        #flegs[i].AddEntry(ftotal, 'Fit Total (chi2/ndf = '+str(round(fitchi2ndfv2,2))+')', lopt)
        #print 'INFO: Fit total MC from Chi2TestX chi2/ndf = '+str(round(fitchi2ndfv2,2))

        mleg.Draw('same')


        ### plot the fit residuals
        #-------------------------------------------------------------
        mbotpad.cd()
        mleg2 = TLegend(0.72, 0.78, 0.94, 0.94)
        #flegs2.append(leg)
        mleg2.SetFillColor(0)
        mleg2.SetBorderSize(0)
        lopt = 'LPE'
        
        #fresid.Divide(fitdata, ftotal)
                
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
        
        mleg2.AddEntry(fresid,'data / MC',lopt)
        #flegs2[i].Draw()
        #-------------------------------------------------------------
        
        mcanv.Update()
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
        
        mcanv.Print('./plots/'+msave+'.png')
        #=============================================================
        #=============================================================

        #raw_input('[Enter] to quit \n')
        
    ### end of fitting bit if you have signals

    
    ### copy the backgrounds file
    shutil.copyfile(mcfile, './plots/'+mcfile)
    
    
    # plot the lo and hi energy histograms for all channels
    #=================================================================
    
    # number of energy ranges (lo, hi)
    numE = 2
    # number of channels (single, multi, lsveto, alpha)
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
            gStyle.SetPadLeftMargin   (0.12)
            gStyle.SetPadRightMargin  (0.02)
            
            font = 63
            size = 13
            yoff = 4.2
            xoff = 8
            
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
                lopt = 'LPE'

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

                # just draw to show plots
                if not showTotal: total[C][E][i].Draw()
                
                resid[C][E][i].Rebin(plotRebin)
                #resid[C][E][i].Scale(1./float(plotRebin))
                #resid[C][E][i].Sumw2()
                
                dkey = 0
                for key in datkeys:
                    if 'x'+str(i+1) in key and '-c'+chan in key and '-e'+str(E) in key:
                        
                        dkey = key
                        
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
                        data[dkey]['hist'].SetMarkerColor(kBlack)
                        data[dkey]['hist'].SetLineColor(kBlack)
                        data[dkey]['hist'].SetLineWidth(1)
                        #if redtotal:
                            #data[dkey]['hist'].SetMarkerStyle(7)
                            #data[dkey]['hist'].SetMarkerSize(1)
                            
                        data[dkey]['hist'].GetYaxis().SetTitleFont(font)
                        data[dkey]['hist'].GetYaxis().SetTitleSize(size)
                        data[dkey]['hist'].GetYaxis().SetTitleOffset(yoff)
                        data[dkey]['hist'].GetYaxis().SetLabelFont(font)
                        data[dkey]['hist'].GetYaxis().SetLabelSize(size)
                        data[dkey]['hist'].GetYaxis().SetLabelOffset(0.01)
                        #data[dkey]['hist'].GetXaxis().SetTitle('Energy (keV)')
                        #data[dkey]['hist'].GetXaxis().SetLabelFont(font)
                        #data[dkey]['hist'].GetXaxis().SetLabelSize(size)

                        if E: data[dkey]['hist'].SetAxisRange(hier[0], hier[1], 'x')
                        else: data[dkey]['hist'].SetAxisRange(loer[0], loer[1], 'x')

                        if dru and i!=8:
                            #data[dkey]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            if E: data[dkey]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            else: data[dkey]['hist'].SetAxisRange(2e-3, 3e2, 'y')

                        ### LSveto plotting
                        if E and i==8:
                            data[dkey]['hist'].SetAxisRange(0, lsHiE, 'x')
                            data[dkey]['hist'].SetAxisRange(1e-5, 0.5, 'y')
                                
                        #popt = 'P E1'
                        #popt = 'HIST'
                        popt = ''
                        data[dkey]['hist'].Draw()
                        #days = round(data[dkey]['runtime']/86400.,2)
                        days = round(runtime/86400., 2)
                        
                        if ingroups:
                            legs[C][E][i].AddEntry(data[dkey]['hist'], 'Data', lopt)
                        else:
                            legs[C][E][i].AddEntry(data[dkey]['hist'], dkey+' ('+str(days)+' days)', lopt)

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
                            bkgs[key]['hist'].Draw('same')

                        ### add MC to total MC hist
                        total[C][E][i].Add(bkgs[key]['hist'])
                        tcount += 1
                        
                        ### create the legend entry for MC
                        #legs[C][E][i].AddEntry(bkgs[key]['hist'], key, lopt)

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
                            sigs[key]['hist'].Draw('same')

                        ### add MC to total MC hist
                        total[C][E][i].Add(sigs[key]['hist'])
                        tcount += 1

                        ### create the legend entry for MC
                        #legs[C][E][i].AddEntry(sigs[key]['hist'], key, lopt)

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
                            gbkgs[C][E][i][group].Draw('same')
                            legs[C][E][i].AddEntry(gbkgs[C][E][i][group], group, lopt)
                    if 'none' in groupSort:
                        for key in gbkgs[C][E][i]['none']:
                            gbkgs[C][E][i]['none'][key].SetMarkerColor(gis['none'])
                            gbkgs[C][E][i]['none'][key].SetLineColor(gis['none'])
                            gbkgs[C][E][i]['none'][key].Draw('same')
                            legs[C][E][i].AddEntry(gbkgs[C][E][i]['none'][key], key, lopt)
                else:
                    # add legend entries in order
                    #print uniqAll
                    for name in uniqAll:
                        for bkey in bakkeys:
                            #print bkey
                            if bkey == 'x'+str(i+1)+'-'+name+'-c'+chan+'-e'+str(E) \
                               or bkey == 'x'+str(i+1)+'-'+name+'-f'+str(i+1)+'-c'+chan+'-e'+str(E):
                                activ = '(%.2e) '%(bkgs[bkey]['info']['acti'])
                                legs[C][E][i].AddEntry(bkgs[bkey]['hist'], activ+bkey, lopt)
                        for skey in sigkeys:
                            if skey == 'x'+str(i+1)+'-'+name+'-c'+chan+'-e'+str(E) \
                               or skey == 'x'+str(i+1)+'-'+name+'-f'+str(i+1)+'-c'+chan+'-e'+str(E):
                                activ = '(%.2e) '%(sigs[skey]['info']['fitacti'])
                                legs[C][E][i].AddEntry(sigs[skey]['hist'], activ+skey, lopt)
                
                
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
                chi2  = ROOT.Double(-1)
                ndf   = ROOT.Long(1)
                igood = ROOT.Long(0)
                #print 'THIS ONE???'
                #print dkey, data[dkey]['hist'].GetNbinsX()
                #print total[C][E][i].GetNbinsX()
                
                if dkey and tcount:
                    if data[dkey]['hist'].GetEntries() > 0:
                        pval  = data[dkey]['hist'].Chi2TestX(total[C][E][i], chi2, ndf, igood, chiopt)
                        #print 'INFO:',dkey,'pval =',pval,'chi2 =',chi2,'ndf =',ndf,'igood =',igood
                        print 'INFO:',dkey,'total MC chi2/ndf =',round(chi2/data[dkey]['druScale']/ndf,2)
                #-----------------------------------------------------------------------------
                #=============================================================================



                #-----------------------------------------------------------------------------
                # still draw total even if nothing so hist box is created
                if showTotal:
                    total[C][E][i].Draw('same')
                    if tcount:
                        total[C][E][i].SetLineWidth(1)
                        if ingroups:
                            legs[C][E][i].AddEntry(total[C][E][i], 'Total', lopt)
                        else:
                            legs[C][E][i].AddEntry(total[C][E][i],
                                'Total MC (chi2/ndf = '+str(round(chi2/ndf,2))+')', lopt)

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
                lopt = 'LPE'

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
                
                #legs2[C][E][i].AddEntry(resid[C][E][i],'data / MC',lopt)
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
            canvs[C][E].Print('./plots/'+save+'.png')
            
            
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
                                                        0, 0, 1400, 900)
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
                            
                            sepPlots[C][E][i].Print(str('./plots/'+isave+'.png'))
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
                    
                    combPlots[i].Print(str('./plots/'+csave+'.png'))
                except:
                    print 'WARNING: could not make plot for xstal-'+str(i+1)
                    continue

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
    
    try: del sepPlots
    except: pass
    
    try: del canvs
    except: pass
    
    
    if not batch:
        raw_input('[Enter] to quit \n')

    return

######################################################################
######################################################################

if __name__ == "__main__":
    main(sys.argv[1:])

