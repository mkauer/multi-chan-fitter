#!/usr/bin/env python

######################################################################
# Matt Kauer - mkauer@physics.wisc.edu
######################################################################
# 101-better-bounds.py

V = 'v101'

# Use better bounds handling for single and multi hit channels as well
#   as different bounds for the LSveto (x9).
# This turned into fixing some bugs with steel and testing out
#   smoothing functions so I didn't actually change anything with
#   with the fit bounds. Will do that in v102.
# 
# version: 2018-11-06
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
ROOT.gErrorIgnoreLevel = kWarning

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
#note = 'default'

#mcfile = 'backgrounds_101.0.txt'
#mcfile = 'backgrounds_101.1.txt'
#mcfile = 'backgrounds_101.1-update.txt'
mcfile = 'backgrounds_200.txt'

print 'INFO: using backgrounds config file -->', mcfile


### ==========  OPTIMIZATION OPTIONS  ================================
### MC smoothing? Specify a smoothing window in +/- number of bins
#smoothing = 5
smoothing = 0


### ==========  FITTING OPTIONS  =====================================
### use pushpa's fitting ranges and binning?
pushpa = 0

### select channels to fit
fitchans = 'SM'

### fitting ranges and rebinnings
if pushpa:
    fLoE = [6,  70]   # Pushpa style
    fHiE = [70, 2000] # Pushpa style
    loEfitRebin = 6   # Pushpa style
    hiEfitRebin = 2   # Pushpa style
else:
    fLoE = [6, 76]
    fHiE = [70, 2770]
    # lsveto testing
    #fLoE = [80, 200]
    #fLoE = [0, 0]
    fHiE = [300, 3300]
    loEfitRebin = 3
    hiEfitRebin = 4


### ==========  EXTRA MC OPTIONS  ====================================
### which MC to fit globally (to all crystals simultaneously)?
globalmc = []
#globalmc = ['lsveto', 'pmt', 'innersteel']
#globalmc = ['lsveto', 'pmt', 'copper']
globalmc = ['lsveto', 'pmt', 'innersteel', 'steel']

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

### decide to stack or extend the channels [0,1]
### [0] for stacking
### [1] for extending
extend = 1


### ==========  PLOTTING OPTIONS  ====================================
### individual plots for all crystals? [0,1]
indi = 1

### select channels to plot
pltchans = 'SM'

### plotting ranges
#loer = [0,   70]   # pushpa style
#hier = [100, 2000] # pushpa style
loer = [0, 100]
hier = [0, 3000]

eran = [loer, hier]

### rebin the final plots [1,inf]
if pushpa:
    loEplotRebin = 6  # Pushpa style
    hiEplotRebin = 2  # Pushpa style
else:
    #loEplotRebin = 6
    #hiEplotRebin = 2
    loEplotRebin = loEfitRebin
    hiEplotRebin = hiEfitRebin

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
mcsumw2  = 0
### set data sumw2()? [0,1]
datsumw2 = 1
### set error on the total? [0,1]
toterr   = 0

### Scale data and mc by 1/fitRebin factor? [0,1]
loEfitRebinScale = 0
hiEfitRebinScale = 0

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
    
    #gStyle.SetPadBottomMargin (0.12)
    #gStyle.SetPadLeftMargin   (0.12)
    #gStyle.SetPadRightMargin  (0.05)
    
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
    #data, bkgs, sigs, runtime = build100(mcfile, others, reuse, allchans, xstals)
    data, bkgs, sigs, runtime = build101(mcfile, others, reuse, allchans, xstals)
    print 'INFO: runtime =', runtime, '(seconds)'

    # 2018-07-01
    # a little debug info on what "others" are being generated for lsveto
    # looks correct...
    """
    keys=[]
    for key in bkgs:
        keys.append(key)
    keys.sort()
    for key in keys:
        print key
    sys.exit()
    """
    
    datkeys = sortDataKeys92(data)
    if datsumw2:
        for key in datkeys:
            data[key]['hist'].Sumw2()


    ### scale into dru units
    if dru:
        data = scaleData70(data, 1)
        bkgs = scaleBkgs100(bkgs)
        sigs = scaleBkgs100(sigs)
    else:
        data = scaleData70(data, 0)
        bkgs = scaleBkgs100(bkgs, runtime)
        sigs = scaleBkgs100(sigs, runtime)

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
    justthese=[]
    for i in range(1, numx+1):
        for key in datkeys:
            if 'x'+str(i) in key and i not in justthese:
                justthese.append(i)
    print 'INFO: plotting crystals -->', justthese
    
    ### assume all data is using same runs and hist params
    try: runtag = data[datkeys[0]]['info']['tag']
    except: runtag = 'none'
    try:
        #allchans  = data[datkeys[0]]['info']['chans']
        params = globalParams(data)
    except:
        #allchans  = bkgs[bakkeys[0]]['info']['chans']
        params = globalParams(bkgs)
    
    #print params
    #sys.exit()
    
    ### find unique names for color scheme?
    ### "internal-K40" for example
    uniqBkgs = []
    uniqSigs = []
    uniqAll = []
    for key in bakkeys:
        uniqBkgs.append(key.split('-')[1]+'-'+key.split('-')[2])
        uniqAll.append(key.split('-')[1]+'-'+key.split('-')[2])
    for key in sigkeys:
        uniqSigs.append(key.split('-')[1]+'-'+key.split('-')[2])
        uniqAll.append(key.split('-')[1]+'-'+key.split('-')[2])

    uniqBkgs = sorted(list(set(uniqBkgs)))
    uniqSigs = sorted(list(set(uniqSigs)))
    uniqAll = sorted(list(set(uniqAll)))
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
    print 'INFO: global fits to -->',globstr
    
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
    

    ### set fit bounds for the fit
    #=================================================================
    # NOTE: These need to be integers!!!
    bpkvLo = params[0][3]
    fLo = [int(fLoE[0]*bpkvLo/loEfitRebin), int(fLoE[1]*bpkvLo/loEfitRebin)]
    fLoBins = fLo[1]-fLo[0]
    try: keVperBinLoE = (fLoE[1]-fLoE[0])/float(fLoBins)
    except: keVperBinLoE = 0
    
    bpkvHi = params[1][3]
    fHi = [int(fHiE[0]*bpkvHi/hiEfitRebin), int(fHiE[1]*bpkvHi/hiEfitRebin)]
    fHiBins = fHi[1]-fHi[0]
    try: keVperBinHiE = (fHiE[1]-fHiE[0])/float(fHiBins)
    except: keVperBinHiE = 0
    
    fbins = int(fLoBins+fHiBins)
    
    fmin = 0
    fmax = fbins
    #nchans = len(fitchans)
    if extend:
        #fmax = fbins*nchans
        fmax = fbins*len(fitchans)
    #=================================================================
    
    
    ### need seperate dicts for rebinned data and MC for plotting to
    ### to work right - they need their own memory space
    #-----------------------------------------------------------------
    ### init dict for rebinned data
    rldata = {}
    rhdata = {}
    for dkey in datkeys:
        #print dkey
        rldata[dkey] = {}
        rhdata[dkey] = {}
        
    ### init dict for rebinned backgrounds
    rlbkgs = {}
    rhbkgs = {}
    for bkey in bakkeys:
        #print bkey
        rlbkgs[bkey] = {}
        rhbkgs[bkey] = {}
        
    ### init dict for rebinned MC/signals 
    rlsigs = {}
    rhsigs = {}
    for skey in sigkeys:
        #print skey
        rlsigs[skey] = {}
        rhsigs[skey] = {}
        
    #-----------------------------------------------------------------
    
    
    ### This part puts the histos together to get ready for the fitting
    ##################################################################
    ### only do the fit if you have signals and data!
    #globstr = ''
    resultsfile = ''
    fitting=0
    if len(sigs) > 0:
        fitting=1
        
        #fitdata = []
        fitsigs = {}
        fitbkgs = {}
        
        fitglob = {}
        fglobkeys = []

        #ftotal = []
        #fresid = []
        #druscale = [1 for x in range(numx)]

        fitdata = TH1F('globData', 'globData', fmax*numx, 0, fmax*numx)

        ftotal = TH1F('globTotal', 'globTotal', fmax*numx, 0, fmax*numx)
        ftotal.SetLineColor(kGray+1)
        ftotal.SetMarkerColor(kGray+1)
        ftotal.SetLineWidth(1)

        fresid = TH1F('globResid', 'globResid', fmax*numx, 0, fmax*numx)
        fresid.SetLineColor(kBlack)
        fresid.SetMarkerColor(kBlack)
        fresid.SetLineWidth(1)

        ### fill fitdata and fitsigs
        for i in range(numx):
            
            # cycle through the channels
            sinit = 1
            binit = 1
            tmpscale1 = 0
            tmpscale2 = 0
            for nc, C in enumerate(fitchans):
                
                # build the histos for fitting
                #-----------------------------------------------------------------------------
                for dkey in datkeys:
                    # lo E
                    if 'x'+str(i+1) in dkey and '-c'+C in dkey and '-e0' in dkey:
                        # skip lsveto in low energy
                        #if i==8:
                        #    print 'skipping',dkey
                        #    continue
                        rldata[dkey]['hist'] = copy.deepcopy(data[dkey]['hist'])
                        rldata[dkey]['hist'].Rebin(loEfitRebin)
                        if loEfitRebinScale:
                            rldata[dkey]['hist'].Scale(1./loEfitRebin)
                        for n in range(fLoBins):
                            if extend:
                                fitdata.SetBinContent(n+1+(fbins*nc)+(i*fmax),
                                                      fitdata.GetBinContent(n+1+(fbins*nc)+(i*fmax))
                                                      + rldata[dkey]['hist'].GetBinContent(fLo[0]+n))
                            else:
                                fitdata.SetBinContent(n+1+(i*fmax),
                                                      fitdata.GetBinContent(n+1+(i*fmax))
                                                      + rldata[dkey]['hist'].GetBinContent(fLo[0]+n))
                    # hi E
                    if 'x'+str(i+1) in dkey and '-c'+C in dkey and '-e1' in dkey:
                        rhdata[dkey]['hist'] = copy.deepcopy(data[dkey]['hist'])
                        rhdata[dkey]['hist'].Rebin(hiEfitRebin)
                        if hiEfitRebinScale:
                            rhdata[dkey]['hist'].Scale(1./hiEfitRebin)
                        r = 0
                        for n in range(fLoBins,fbins):
                            if extend:
                                fitdata.SetBinContent(n+1+(fbins*nc)+(i*fmax),
                                                      fitdata.GetBinContent(n+1+(fbins*nc)+(i*fmax))
                                                      + rhdata[dkey]['hist'].GetBinContent(fHi[0]+r))
                            else:
                                fitdata.SetBinContent(n+1+(i*fmax),
                                                      fitdata.GetBinContent(n+1+(i*fmax))
                                                      + rhdata[dkey]['hist'].GetBinContent(fHi[0]+r))
                            r += 1
                
                
                for skey in sigkeys:

                    ### init lsveto histo to globals
                    for gmckey in globalmc:
                        #if 'lsveto' in skey or 'pmt' in skey:
                        if gmckey in skey:
                            #print skey.split('-')
                            bits = skey.split('-')
                            fgkey = bits[1]+'-'+bits[2]
                            if fgkey not in fglobkeys:
                                fglobkeys.append(fgkey)
                                fitglob[fgkey] = {}
                                fglob = TH1F(fgkey, fgkey, fmax*numx, 0, fmax*numx)
                                fitglob[fgkey]['hist'] = fglob
                            #continue
                    
                    # lo E
                    if 'x'+str(i+1) in skey and '-c'+C in skey and '-e0' in skey:
                        # skip lsveto in low energy
                        #if i==8:
                        #    print 'skipping',skey
                        #    continue
                        fskey = skey.split('-c'+C+'-e0')[0]
                        if sinit:
                            fsig = TH1F(fskey, fskey, fmax*numx, 0, fmax*numx)
                            fitsigs[fskey] = {}
                            fitsigs[fskey]['hist'] = fsig
                        rlsigs[skey]['hist'] = copy.deepcopy(sigs[skey]['hist'])
                        rlsigs[skey]['hist'].Rebin(loEfitRebin)
                        if loEfitRebinScale: rlsigs[skey]['hist'].Scale(1./loEfitRebin)
                        for n in range(fLoBins):
                            try:
                                if extend:
                                    fitsigs[fskey]['hist'].SetBinContent(n+1+(fbins*nc)+(i*fmax),
                                            fitsigs[fskey]['hist'].GetBinContent(n+1+(fbins*nc)+(i*fmax))
                                            + rlsigs[skey]['hist'].GetBinContent(fLo[0]+n))
                                else:
                                    fitsigs[fskey]['hist'].SetBinContent(n+1+(i*fmax),
                                            fitsigs[fskey]['hist'].GetBinContent(n+1+(i*fmax))
                                            + rlsigs[skey]['hist'].GetBinContent(fLo[0]+n))
                            except: pass
                    # hi E
                    if 'x'+str(i+1) in skey and '-c'+C in skey and '-e1' in skey:
                        rhsigs[skey]['hist'] = copy.deepcopy(sigs[skey]['hist'])
                        rhsigs[skey]['hist'].Rebin(hiEfitRebin)
                        if hiEfitRebinScale: rhsigs[skey]['hist'].Scale(1./hiEfitRebin)
                        fskey = skey.split('-c'+C+'-e1')[0]
                        r = 0
                        for n in range(fLoBins,fbins):
                            try:
                                if extend:
                                    fitsigs[fskey]['hist'].SetBinContent(n+1+(fbins*nc)+(i*fmax),
                                            fitsigs[fskey]['hist'].GetBinContent(n+1+(fbins*nc)+(i*fmax))
                                            + rhsigs[skey]['hist'].GetBinContent(fHi[0]+r))
                                else:
                                    fitsigs[fskey]['hist'].SetBinContent(n+1+(i*fmax),
                                            fitsigs[fskey]['hist'].GetBinContent(n+1+(i*fmax))
                                            + rhsigs[skey]['hist'].GetBinContent(fHi[0]+r))
                            except: pass
                            r += 1
                sinit=0
                
                for bkey in bakkeys:
                    # lo E
                    if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e0' in bkey:
                        # skip lsveto in low energy
                        #if i==8:
                        #    print 'skipping',bkey
                        #    continue
                        fbkey = bkey.split('-c'+C+'-e0')[0]
                        if binit:
                            fbak = TH1F(fbkey, fbkey, fmax*numx, 0, fmax*numx)
                            fitbkgs[fbkey] = {}
                            fitbkgs[fbkey]['hist'] = fbak
                        rlbkgs[bkey]['hist'] = copy.deepcopy(bkgs[bkey]['hist'])
                        rlbkgs[bkey]['hist'].Rebin(loEfitRebin)
                        if loEfitRebinScale: rlbkgs[bkey]['hist'].Scale(1./loEfitRebin)
                        for n in range(fLoBins):
                            try:
                                if extend:
                                    fitbkgs[fbkey]['hist'].SetBinContent(n+1+(fbins*nc)+(i*fmax),
                                            fitbkgs[fbkey]['hist'].GetBinContent(n+1+(fbins*nc)+(i*fmax))
                                            + rlbkgs[bkey]['hist'].GetBinContent(fLo[0]+n))
                                else:
                                    fitbkgs[fbkey]['hist'].SetBinContent(n+1+(i*fmax),
                                            fitbkgs[fbkey]['hist'].GetBinContent(n+1+(i*fmax))
                                            + rlbkgs[bkey]['hist'].GetBinContent(fLo[0]+n))
                            except: pass
                    # hi E
                    if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e1' in bkey:
                        rhbkgs[bkey]['hist'] = copy.deepcopy(bkgs[bkey]['hist'])
                        rhbkgs[bkey]['hist'].Rebin(hiEfitRebin)
                        if hiEfitRebinScale: rhbkgs[bkey]['hist'].Scale(1./hiEfitRebin)
                        fbkey = bkey.split('-c'+C+'-e1')[0]
                        r = 0
                        for n in range(fLoBins,fbins):
                            try:
                                if extend:
                                    fitbkgs[fbkey]['hist'].SetBinContent(n+1+(fbins*nc)+(i*fmax),
                                            fitbkgs[fbkey]['hist'].GetBinContent(n+1+(fbins*nc)+(i*fmax))
                                            + rhbkgs[bkey]['hist'].GetBinContent(fHi[0]+r))
                                else:
                                    fitbkgs[fbkey]['hist'].SetBinContent(n+1+(i*fmax),
                                            fitbkgs[fbkey]['hist'].GetBinContent(n+1+(i*fmax))
                                            + rhbkgs[bkey]['hist'].GetBinContent(fHi[0]+r))
                            except: pass
                            r += 1
                binit=0
                

            ### subtract fixed MC from data
            #-----------------------------------------------
            for key in fitbkgs:
                if 'x'+str(i+1) in key:
                    fitdata.Add(fitbkgs[key]['hist'], -1)
            #-----------------------------------------------

        
        ### sort all the fbakkeys
        delete=[]
        fbakkeys=[]
        for fbkey in fitbkgs:
            if fitbkgs[fbkey]['hist'].Integral() > 0:
                fbakkeys.append(fbkey)
            else: delete.append(fbkey)
        for key in delete:
            if debug: print 'DEBUG: zero events so deleting bkg key', key
            del fitbkgs[key]
        fbakkeys.sort()

        ### sort all the fsigkeys
        delete=[]
        fsigkeys=[]
        for fskey in fitsigs:
            if fitsigs[fskey]['hist'].Integral() > 0:
                fsigkeys.append(fskey)
            else: delete.append(fskey)
        for key in delete:
            if debug: print 'DEBUG: zero events so deleting sig key', key
            del fitsigs[key]
        fsigkeys.sort()

        ### sort all the fglobkeys
        """
        delete=[]
        fglobkeys=[]
        for fgkey in fitglob:
            if fitglob[fgkey]['hist'].Integral() > 0:
                fglobkeys.append(fgkey)
            else: delete.append(fgkey)
        for key in delete:
            print 'INFO: deleting fit glob key', key
            del fitglob[key]
        """
        fglobkeys.sort()
        
        
        ### build global sigs
        for fskey in fsigkeys:
            for fgkey in fglobkeys:
                if fgkey in fskey:
                    #print 'INFO: adding',fskey,'to',fgkey
                    fitglob[fgkey]['hist'].Add(fitsigs[fskey]['hist'])
        
        ### remove global sigs from fit sigs
        L = len(fsigkeys)-1
        for k, fskey in enumerate(reversed(fsigkeys)):
            for gmckey in globalmc:
                ### make gmckey unique (steel vs innersteel) 2018-10-23
                gmckey = '-'+gmckey+'-'
                #if debug: print 'DEBUG: is',gmckey,'in',fskey
                #if 'lsveto' in fskey or 'pmt' in fskey:
                if gmckey in fskey:
                    if debug: print 'DEBUG: remove global key from sigs', fsigkeys[L-k]
                    del fsigkeys[L-k]

        ### now delete the empty histos
        delete=[]
        fglobkeys=[]
        for fgkey in fitglob:
            if fitglob[fgkey]['hist'].Integral() > 0:
                fglobkeys.append(fgkey)
            else: delete.append(fgkey)
        for key in delete:
            if debug: print 'DEBUG: zero events so deleting global key', key
            del fitglob[key]
        fglobkeys.sort()

        ### global fit debug infos
        #=============================================================
        """
        tmpcand = TCanvas('tmpcand','tmpcand',1200,600)
        fitdata.Draw()
        fitdata.SetAxisRange(2e-3, 3e2, 'y')
        tmpcand.SetLogy(1)
        tmpcand.Update()
        #raw_input('Enter to continue')
        
        tmpcanb = TCanvas('tmpcanb','tmpcanb',1200,600)
        for bkey in fbakkeys:
            fitbkgs[bkey]['hist'].Draw('same')
            fitbkgs[bkey]['hist'].SetAxisRange(2e-3, 3e2, 'y')
        tmpcanb.SetLogy(1)
        tmpcanb.Update()
        #raw_input('Enter to continue')
        
        tmpcans = TCanvas('tmpcans','tmpcans',1200,600)
        for skey in fsigkeys:
            fitsigs[skey]['hist'].Draw('same')
            fitsigs[skey]['hist'].SetAxisRange(2e-3, 3e2, 'y')
        tmpcans.SetLogy(1)
        tmpcans.Update()
        #raw_input('Enter to continue')
        
        tmpcang = TCanvas('tmpcang','tmpcang',1200,600)
        for gkey in fglobkeys:
            fitglob[gkey]['hist'].Draw('same')
            fitglob[gkey]['hist'].SetAxisRange(2e-3, 3e2, 'y')
        tmpcang.SetLogy(1)
        tmpcang.Update()
        #raw_input('Enter to continue')
        """
        #=============================================================
        
        
        sigObj = []
        fit = []
        fitresults = {}
        fitchi2ndf = []
        wasFit = []
        bounds = []
        
        ### data integral to normalize signal to
        dat_int = fitdata.Integral(0, fmax*numx)

        ### set up the fitting object for TFractionFitter
        for i in range(numx):
            
            fitresults[str(i)] = []
            
            """
            ### Do a per crystal unique of the signals to see how many
            ### signal channels each crystal fit will have
            uniqSig = []
            noEvents = []
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    mc_int = fitsigs[fskey]['hist'].Integral(i*fmax,(i+1)*fmax)
                    if mc_int > 0.0:
                        #print 'Integral of', fskey, '=', mc_int
                        uniqSig.append(fskey.split('-')[1]+'-'+fskey.split('-')[2])
                    else:
                        print '\nWARNING: No events for --> ',fskey
                        print   '         Maybe set it as a background...\n'
                        noEvents.append(fskey)
            
            uniqSig = sorted(list(set(uniqSig)))
            
            ### delete out keys that have noEntries
            if len(noEvents) > 0:
                ### delete the key out of fsigkeys so it does
                ### not get plotted in the fit histos!
                ### go in reverse to preserve index
                L = len(fsigkeys)-1
                for k,fskey in enumerate(reversed(fsigkeys)):
                    if 'x'+str(i+1) in fskey and fskey in noEvents:
                        print 'INFO: deleting fit key',fsigkeys[L-k]
                        del fsigkeys[L-k]

                ### and delete it out of sigs so it does not
                ### get plotted in the lo-E and hi-E histos
                ### go in reverse to preserve index
                L = len(sigkeys)-1
                for k,skey in enumerate(reversed(sigkeys)):
                    for noEvt in noEvents:
                        if 'x'+str(i+1) in skey and noEvt in skey:
                            print 'INFO: deleting signal key',sigkeys[L-k]
                            del sigkeys[L-k]
            
            ### TFractionFitter wants at least 2 MC to converge the fit
            if len(uniqSig) < 2:
                print "\nWARNING: TFractionFitter needs at least 2 MC to converge"
                print   "         Skipping fit to crystal",str(i+1)

                ### delete the key out of fsigkeys so it does
                ### not get plotted in the fit histos!
                ### go in reverse to preserve index
                L = len(fsigkeys)-1
                for k,fskey in enumerate(reversed(fsigkeys)):
                    if 'x'+str(i+1) in fskey:
                        print 'INFO: deleting fit key',fsigkeys[L-k]
                        del fsigkeys[L-k]
                        
                ### and delete it out of sigs so it does not
                ### get plotted in the lo-E and hi-E histos
                ### go in reverse to preserve index
                L = len(sigkeys)-1
                for k,skey in enumerate(reversed(sigkeys)):
                    if 'x'+str(i+1) in skey:
                        print 'INFO: deleting signal key',sigkeys[L-k]
                        del sigkeys[L-k]
                        
                continue
            
            print '\nINFO: Fitting crystal',str(i+1),'with',uniqSig,'\n'
            
            
            #sigObj.append(TObjArray(Nsigs)) # number of MC to fit to
            sigObj.append(TObjArray(len(uniqSig))) # number of MC to fit to
            
            if datsumw2:
                fitdata.Sumw2()
            """
            
            
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:

                    ### 2018-03-05
                    ### remove keys and hists with no events
                    """
                    if fitsigs[fskey]['hist'].Integral(i*fmax,(i+1)*fmax) <= 0:
                        print '\nINFO: No events for --> ',fskey
                        print   '      deleting fitsigs key', fskey
                        del fitsigs[fskey]
                        for k,key in enumerate(fsigkeys):
                            if key == fskey:
                                print   '      deleting fsigkeys',k,key,fskey
                                del fsigkeys[k]
                        continue
                    """
                        
                    mc_int = fitsigs[fskey]['hist'].Integral(i*fmax,(i+1)*fmax) # MC integral

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

                        ### some crystals don't have events from "other" crystals
                        #print '\nINFO: No events for --> ',fskey
                        #print   '      Deleting/ignoring this histogram...\n'
                        #for key in fsigkeys:
                        #del fsigkeys[fskey]
                        #del fitsigs[fskey]
                        
                    for C in allchans:
                        for E in range(2):
                            E=str(E)
                            
                            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            newkeys = []
                            try:
                                tmpkey=fskey+'-c'+C+'-e'+E
                                test = sigs[tmpkey]['hist']
                                newkeys.append(tmpkey)
                            except:
                                for F in range(numx):
                                    F=str(F+1)
                                    try:
                                        tmpkey=fskey+'-f'+F+'-c'+C+'-e'+E
                                        test = sigs[tmpkey]['hist']
                                        newkeys.append(tmpkey)
                                    except:
                                        continue
                            """
                            if 'pmt' not in fskey:
                                tmpkey=fskey+'-c'+C+'-e'+E
                                try:
                                    test = sigs[tmpkey]['hist'].Integral()
                                    #print test,tmpkey
                                    #newkeys.append(tmpkey)
                                except:
                                    print 'INFO: skipping', tmpkey
                                    continue
                                newkeys.append(tmpkey)
                            else:
                                for F in range(numx):
                                    F=str(F+1)
                                    tmpkey=fskey+'-f'+F+'-c'+C+'-e'+E
                                    try:
                                        test = sigs[tmpkey]['hist'].Integral()
                                        #print test,tmpkey
                                        #newkeys.append(tmpkey)
                                    except:
                                        print 'INFO: skipping', tmpkey
                                        continue
                                    newkeys.append(tmpkey)
                            """
                            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            
                            for newkey in newkeys:

                                ### delete if there are no events in the hist
                                #try:
                                #    if sigs[newkey]['hist'].Integral() <= 0:
                                        
                                #try:
                                sigs[newkey]['hist'].Scale(dat_int/mc_int)
                                sigs[newkey]['fitscale'] = sigs[newkey]['scale'] * dat_int/mc_int
                                #except:
                                #    sigs[newkey]['fitscale'] = sigs[newkey]['scale']
                                #    print '\nWARNING: No events for --> ',newkey
                                
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

            ### Print out the fit infos
            """
            ### First make a string list of the globals
            globstr = ''
            if len(globalmc) == 0:
                globstr = 'none'
            else:
                for txt in globalmc:
                    globstr += txt+'-'
                globstr = globstr[:-1]
            """
            
            fitresults[str(i)].append('Crystal-'+str(i+1)+' fit results')
            fitresults[str(i)].append('runtime = '+str(round(runtime/60./60./24., 2))+' days')
            if note: fitresults[str(i)].append('note = '+note)
            fitresults[str(i)].append('version = '+V)
            fitresults[str(i)].append('channels fit = '+fitchans)
            fitresults[str(i)].append('global fits = '+globstr)
            fitresults[str(i)].append('other pmts = '+str(others))
            fitresults[str(i)].append('hist extend = '+str(extend))
            fitresults[str(i)].append('norm to dru = '+str(dru))
            fitresults[str(i)].append('lo-E fit range = '+str(fLoE[0])+' - '+str(fLoE[1])+' keV')
            fitresults[str(i)].append('lo-E fit binning = '+str(keVperBinLoE)+' keV/bin')
            fitresults[str(i)].append('hi-E fit range = '+str(fHiE[0])+' - '+str(fHiE[1])+' keV')
            fitresults[str(i)].append('hi-E fit binning = '+str(keVperBinHiE)+' keV/bin')

        
        ### do the same to the global lsveto
        boundskeys=[]
        for fgkey in fglobkeys:
                
            mc_int = fitglob[fgkey]['hist'].Integral(0, fmax*numx) # total MC integral
            
            ### to weight or not to weight...
            if mcsumw2:
                fitglob[fgkey]['hist'].Sumw2() # set stat weights
            
            ### normalize MC to total data
            ### needed for TFractionFitter to work right
            ### still don't fully understand why
            try:
                fitglob[fgkey]['hist'].Scale(dat_int/mc_int) # scale to data integral
            except:
                print '\nERROR: No events for --> ',fgkey
                print   '       Remove it from the fit!\n'
                sys.exit()
            
            #newkey=0
            #boundskeys=[]
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
            fitglob[fgkey]['hist'] = zeroBinError(fitglob[fgkey]['hist'])
        
        
        ### global fit debug infos
        #=============================================================
        """
        tmpcand2 = TCanvas('tmpcand2','tmpcand2',1200,600)
        fitdata.Draw()
        fitdata.SetAxisRange(2e-3, 3e2, 'y')
        tmpcand2.SetLogy(1)
        tmpcand2.Update()
        #raw_input('Enter to continue')
        
        tmpcanb2 = TCanvas('tmpcanb2','tmpcanb2',1200,600)
        for bkey in fbakkeys:
            fitbkgs[bkey]['hist'].Draw('same')
            fitbkgs[bkey]['hist'].SetAxisRange(2e-3, 3e2, 'y')
        tmpcanb2.SetLogy(1)
        tmpcanb2.Update()
        #raw_input('Enter to continue')
        
        tmpcans2 = TCanvas('tmpcans2','tmpcans2',1200,600)
        for skey in fsigkeys:
            fitsigs[skey]['hist'].Draw('same')
            fitsigs[skey]['hist'].SetAxisRange(2e-3, 3e2, 'y')
        tmpcans2.SetLogy(1)
        tmpcans2.Update()
        #raw_input('Enter to continue')
        
        tmpcang2 = TCanvas('tmpcang2','tmpcang2',1200,600)
        for gkey in fglobkeys:
            fitglob[gkey]['hist'].Draw('same')
            fitglob[gkey]['hist'].SetAxisRange(2e-3, 3e2, 'y')
        tmpcang2.SetLogy(1)
        tmpcang2.Update()
        #raw_input('Enter to continue')
        """
        #=============================================================


        
        ### conflict of interests going on here
        #=========================================================
        if datsumw2:
            fitdata.Sumw2()

        if zeroFitDataError:
            fitdata = zeroBinError(fitdata)
        #=========================================================


        totalNumFits = len(fsigkeys)+len(fglobkeys)
        print 'INFO: total number of hists being fit =', totalNumFits,'\n\n'
        sigObj = TObjArray(totalNumFits)
        
        for fskey in fsigkeys:
            sigObj.append(fitsigs[fskey]['hist'])

        for fgkey in fglobkeys:
            sigObj.append(fitglob[fgkey]['hist'])
        
        fit = TFractionFitter(fitdata, sigObj)
        
        ### set fit bounds!!!
        ### l=0 sets all params to the same constrain
        ### set all bounds by default
        #fit.Constrain(0, 0.0, 1.0)
        for l in range(len(bounds)):
            fit.Constrain(l+1, bounds[l][0], bounds[l][1])
        
        ### set the fit range
        fit.SetRangeX(0, fmax*numx)

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
        for n in range(fmax*numx):
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
        for fgkey in fglobkeys:

            fscale = ROOT.Double(0.0)
            ferror = ROOT.Double(0.0)
            #print 'count',count
            fit.GetResult(count, fscale, ferror)
            count += 1
            
            fitglob[fgkey]['hist'].Scale(fscale)
            
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

        
        ### double check for broken keys
        """
        newkeys=[]
        for skey in sigkeys:
            try:
                test = sigs[skey]
                newkeys.append(skey)
            except:
                pass
        sigkeys = newkeys
        sigkeys.sort()
        """
        
        ### scale the signals to mBq/kg
        if dru: sigs = scaleSigs100(sigkeys, sigs)
        else: sigs = scaleSigs100(sigkeys, sigs, runtime)
        
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
            for fgkey in fglobkeys:
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
        save += '_loEfit-'+str(int(fLoE[0]))+'-'+str(int(fLoE[1]))
        save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        save += '_loEfitRebin-'+str(loEfitRebin)
        save += '_hiEfitRebin-'+str(hiEfitRebin)
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
        """
        gStyle.SetPadTopMargin    (0.06)
        gStyle.SetPadBottomMargin (0.12)
        gStyle.SetPadLeftMargin   (0.12)
        #gStyle.SetPadRightMargin  (0.05)
        gStyle.SetPadRightMargin  (0.02)
        """
        ftoppad=[]
        fbotpad=[]

        sepfitdata = [0 for x in range(numx)]
        sepfitsigs = [0 for x in range(numx)]
        sepfitglob = [0 for x in range(numx)]
        
        #sepFitPlots=[]
        
        flegs=[]
        flegs2=[]

        fzeros=[]
        
        """
        font=63
        size=13
        yoff=4.2
        xoff=8
        """
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
            fitdata.SetAxisRange(i*fmax, (i+1)*fmax, 'x')
            
            
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
                    ftotal.Add(fitsigs[fskey]['hist'])

            for fgkey in fglobkeys:
                #if i==0:
                Nfsigs+=1
                
                # find the unique name for color and set color
                #cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
                cname = fgkey
                fitglob[fgkey]['hist'].SetMarkerColor(uniqColor[cname])
                fitglob[fgkey]['hist'].SetLineColor(uniqColor[cname])
                
                ### draw the sigs
                fitglob[fgkey]['hist'].Draw('same')

                if i==0:
                    ### add MC to total MC hist
                    ftotal.Add(fitglob[fgkey]['hist'])
                
                
            ### select the right range
            ftotal.SetAxisRange(i*fmax, (i+1)*fmax, 'x')
            
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
                for fgkey in fglobkeys:
                    if name in fgkey:
                        flegs[i].AddEntry(fitglob[fgkey]['hist'], fgkey, lopt)
            
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

            fresid.SetAxisRange(i*fmax, (i+1)*fmax, 'x')

            if linres:
                fbotpad[i].SetLogy(0)
                fresid.SetAxisRange(lrs[0], lrs[1], 'y')
            fresid.Draw()
            
            ### set my reference line to '1'
            #zero = TLine(fmin, 1, fmax, 1)
            #fzeros.append(zero)
            #fzeros[i].SetLineColor(kRed)
            #fzeros[i].SetLineWidth(1)
            
            zero = TLine(i*fmax, 1, (i+1)*fmax, 1)
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
                    fisave += '_loEfRS'+str(loEfitRebinScale)
                    fisave += '_hiEfRS'+str(hiEfitRebinScale)
                    fisave += '_dru'+str(dru)
                    fisave += '_cs'+str(fitchans)
                    fisave += '_ext'+str(extend)
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
        fitdata.SetAxisRange(0, fmax*numx, 'x')
        
        
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

        for fgkey in fglobkeys:
            #if i==0:
            Nfsigs+=1

            # find the unique name for color and set color
            #cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
            cname = fgkey
            fitglob[fgkey]['hist'].SetMarkerColor(uniqColor[cname])
            fitglob[fgkey]['hist'].SetLineColor(uniqColor[cname])

            ### draw the sigs
            fitglob[fgkey]['hist'].Draw('same')

            #if i==0:
                ### add MC to total MC hist
                #ftotal.Add(fitglob[fgkey]['hist'])

                
        ### select the right range
        ftotal.SetAxisRange(0, fmax*numx, 'x')
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
            for fgkey in fglobkeys:
                if name in fgkey:
                    mleg.AddEntry(fitglob[fgkey]['hist'], fgkey, lopt)

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
        
        fresid.SetAxisRange(0, fmax*numx, 'x')
        
        if linres:
            mbotpad.SetLogy(0)
            fresid.SetAxisRange(lrs[0], lrs[1], 'y')
        fresid.Draw()

        ### set my reference line to '1'
        #zero = TLine(fmin, 1, fmax, 1)
        #fzeros.append(zero)
        #fzeros[i].SetLineColor(kRed)
        #fzeros[i].SetLineWidth(1)

        mzero = TLine(0, 1, fmax*numx, 1)
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
        msave += '_loEfit-'+str(int(fLoE[0]))+'-'+str(int(fLoE[1]))
        msave += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        msave += '_loEfitRebin-'+str(loEfitRebin)
        msave += '_hiEfitRebin-'+str(hiEfitRebin)
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
    
    #canvs  = [[[] for x in range(numE)] for x in range(numC)]
    canvs  = [[0 for x in range(numE)] for x in range(numC)]
    
    ### for separate plots
    #sepPlots = [[[[] for x in range(numx)] for x in range(numE)] for x in range(numC)]
    sepPlots = [[[0 for x in range(numx)] for x in range(numE)] for x in range(numC)]
    sepTopPad = [[[0 for x in range(numx)] for x in range(numE)] for x in range(numC)]
    sepBotPad = [[[0 for x in range(numx)] for x in range(numE)] for x in range(numC)]
    
    ### seperate memory space for the pads is key!!!!
    #toppad = [[[] for x in range(numE)] for x in range(numC)]
    #botpad = [[[] for x in range(numE)] for x in range(numC)]
    toppad = [[0 for x in range(numE)] for x in range(numC)]
    botpad = [[0 for x in range(numE)] for x in range(numC)]
    
    #legs   = [[[] for x in range(numE)] for x in range(numC)]
    #legs2  = [[[] for x in range(numE)] for x in range(numC)]
    #zeros  = [[[] for x in range(numE)] for x in range(numC)]
    legs   = [[0 for x in range(numE)] for x in range(numC)]
    legs2  = [[0 for x in range(numE)] for x in range(numC)]
    zeros  = [[0 for x in range(numE)] for x in range(numC)]
    
    #total  = [[[] for x in range(numE)] for x in range(numC)]
    #resid  = [[[] for x in range(numE)] for x in range(numC)]
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
                            data[dkey]['hist'].SetAxisRange(0, 3500, 'x')
                            #data[dkey]['hist'].SetAxisRange(7e-6, 2e-1, 'y')
                                
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

                        # don't think I need this anymore
                        #if dru1 or dru2:
                        #    druScale = data['x'+str(i+1)+'-data'+'-e'+str(E)]['druScale']
                        #if dru2:
                        #    sigs[key]['hist'].Scale(druScale)

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
                #total[C][E][i].Sumw2()
                
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
                    resid[C][E][i].SetAxisRange(0, 3500, 'x')
                
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
                    zero = TLine(0, 1, 3500, 1)
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
            save += '_E'+str(E)
            save += '_loEfit-'+str(int(fLoE[0]))+'-'+str(int(fLoE[1]))
            save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
            save += '_loEfitRebin-'+str(loEfitRebin)
            save += '_hiEfitRebin-'+str(hiEfitRebin)
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
            save += '_chan'+chan
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
                            isave += '-ext'+str(extend)
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
                    csave += '_ext'+str(extend)
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

