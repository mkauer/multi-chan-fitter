#!/usr/bin/env python

######################################################################
# Matt Kauer - mkauer@physics.wisc.edu
######################################################################
# 80-global-fit.py

V = 'v80'

# Try to get global fitting to work!
# 
# version: 2017-11-06
#
# see CHANGELOG for changes
######################################################################

import os,sys
import shutil
import socket
import copy
import math
import numpy as np

import ROOT
from ROOT import *
ROOT.gROOT.Reset()
ROOT.gErrorIgnoreLevel = kWarning

#sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
#sys.path.append("/home/mkauer/mc-fitting/")
from funcs80 import *


#xstal = 8
#justthese = [xstal]

### ========== GENERAL INPUTS ==============================
### note to add to saved plot names?
note = 0
#note = 'C'+str(xstal)

### backgrounds file
#mcfile = 'backgrounds800-C'+str(xstal)+'-update.txt'
#mcfile = 'backgrounds800.txt'
mcfile = 'backgrounds801.txt'
#mcfile = 'backgrounds801-C'+str(xstal)+'.txt'


### force the reuse of all joined rootfiles in mcfile? [0,1,2]
### very nice for debugging
### [0] default - use whatever is specified in the backgrounds file
### [1] forces reusing of all data/bkgs/sigs
### [2] forces NOT reusing any data/bkgs/sigs
reuse = 0

### update and save new backgrounds file with fit results
updateMCfile = 0


### ========== FITTING OPTIONS =============================
### select channels to plot and fit
### ['S'] single-hit data channel
### ['M'] multi-hit data channel
### ['SM'] both channels
fitchans = 'SM'

### fitting ranges
### lo and hi energy fit ranges or set [0,0]
fLoE = [2, 100]
#fLoE = [0,0]
#fHiE = [60, 2800]
fHiE = [50, 3100]
#fHiE = [0,0]

### rebin the histos for fitting [1,inf]
loEfitRebin = 4
hiEfitRebin = 10

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


### ========== PLOTTING OPTIONS ============================
### channels to plot
pltchans = 'SM'

### plotting ranges
### lo and hi energy ranges
loer = [0,  100]
#loer = [0,  20]
hier = [0, 3500]
eran = [loer, hier]

### rebin the final plots [1,inf]
loEplotRebin = 6
hiEplotRebin = 10

### show the total? [0,1]
showTotal = 1
### show the legends? [0,1]
showlegs = 1
### plot components in groups? [0,1]
ingroups = 1

### use linear residual scale? [0,1]
linres = 1
### set y scale on resid plot
lrs = [0, 2]

### main plots in linear scale [0,1]
liny = 0

### individual plots for all crystals? [0,1]
indi = 1
### just plot individual for crystals? [1-8]
justthese = [1,2,3,4,5,6,7,8]
#justthese = [1]


### ========== CAN EFFECT FIT RESULTS ======================
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



#batch = 0
#if onCup(): batch = 1
#gROOT.SetBatch(batch)


def myself(argv):

    batch = 0
    if onCup(): batch = 1
    gROOT.SetBatch(batch)
    
    gStyle.SetPalette (1)
    gStyle.SetOptStat ('')
    gStyle.SetOptFit  (0)
    
    gStyle.SetPadBottomMargin (0.12)
    gStyle.SetPadLeftMargin   (0.12)
    gStyle.SetPadRightMargin  (0.05)
    
    ### for saving the plots...
    if not os.path.exists('./plots'): 
        os.makedirs('./plots')
    
    if not os.path.exists(mcfile):
        print '\nERROR: could not find backgrounds file -->', mcfile
        sys.exit()
    
    
    ### where everything gets loaded into dictionary
    #-----------------------------------------------------------------
    allchans = uniqString(fitchans+pltchans)
    data, bkgs, sigs, runtime = build80(mcfile, reuse, allchans)
    #data, bkgs, sigs, runtime = build80(mcfile, reuse, allchans, [5,8])
    datkeys, bakkeys, sigkeys = sortKeys(data, bkgs, sigs)
    if datsumw2:
        for key in datkeys:
            data[key]['hist'].Sumw2()
    #-----------------------------------------------------------------
    
    # assume all data is using same runs and hist params
    try: runtag = data[datkeys[0]]['info']['tag']
    except: runtag = 'none'
    try:
        #allchans  = data[datkeys[0]]['info']['chans']
        params = globalParams(data)
    except:
        #allchans  = bkgs[bakkeys[0]]['info']['chans']
        params = globalParams(bkgs)

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
    print 'INFO: Unique bkgs and sigs =',uniqAll
    
    # scale into dru units
    if dru:
        data = scaleData70(data, 1)
        bkgs = scaleBkgs71(bkgs)
        sigs = scaleBkgs71(sigs)
    else:
        data = scaleData70(data, 0)
        bkgs = scaleBkgs71(bkgs, runtime)
        sigs = scaleBkgs71(sigs, runtime)

    ### Number of colors
    Nc = len(uniqAll)
    print 'INFO: Total number of unique bkgs and sigs =',Nc
    colors, cis = rainbow(Nc)

    ### Create color dict for unique simulations
    uniqColor = {}
    for i, key in enumerate(uniqAll):
        #print key
        uniqColor[key] = cis[i]

    ### colors for groups
    #gcs, gis = rainbow(5)
    gis = [kRed, kOrange, kGreen+1, kBlue, kViolet,
           kRed, kOrange, kGreen+1, kBlue, kViolet]
    
    ### legend length = MC + data + total
    Nlg = Nc+2
    lnc = 1
    if Nlg > 6: lnc = 2
    if Nlg > 12: lnc = 3
    xlegstop  = 0.94
    xlegstart = xlegstop-(0.2*lnc)
    ylegstop  = 0.89
    ylegstart = ylegstop-(0.04*6)
    

    ### set fit bounds for the fit
    #=================================================================
    # NOTE: These need to be integers!!!
    bpkvLo = params[0][3]
    fLo = [int(fLoE[0]*bpkvLo/loEfitRebin), int(fLoE[1]*bpkvLo/loEfitRebin)]
    fLoBins = fLo[1]-fLo[0]
    
    bpkvHi = params[1][3]
    fHi = [int(fHiE[0]*bpkvHi/hiEfitRebin), int(fHiE[1]*bpkvHi/hiEfitRebin)]
    fHiBins = fHi[1]-fHi[0]
    
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
        
    # init dict for rebinned MC/signals 
    rlsigs = {}
    rhsigs = {}
    for skey in sigkeys:
        #print skey
        rlsigs[skey] = {}
        rhsigs[skey] = {}
        
    #-----------------------------------------------------------------
    
    
    ### here's the fitting part!!!
    ##################################################################
    ### only do the fit if you have signals and data!
    resultsfile = ''
    fitting=0
    if len(sigs) > 0:
        fitting=1
        
        fitdata = []
        fitsigs = {}
        fitbkgs = {}
        ftotal = [] # total of fixed MC plus fit MC
        fresid = [] # residual of the data/total
        druscale = [1 for x in range(8)]
        
        ### fill fitdata and fitsigs
        for i in range(8):
            #druscale.append(1)
            
            fdata = TH1F('fdata'+str(i), cnames(i)+' - multi-chan fit', fmax,0,fmax)
            #fdata = TH1F('fdata'+str(i), cnames(i)+' - multi-chan fit', fbins,0,fbins)
            #fdata = TH1I('fdata'+str(i), cnames(i)+' - multi-chan fit', fbins,0,fbins)
            fitdata.append(fdata)
            
            tot = TH1F('ftotal'+str(i), longNames(i), fmax,0,fmax)
            #tot = TH1F('ftotal'+str(i), longNames(i), fbins,0,fbins)
            #tot = TH1I('ftotal'+str(i), longNames(i), fbins,0,fbins)
            tot.SetLineColor(kGray+1)
            tot.SetMarkerColor(kGray+1)
            tot.SetLineWidth(1)
            ftotal.append(tot)
            
            res = TH1F('fresid'+str(i), longNames(i), fmax,0,fmax)
            #res = TH1F('fresid'+str(i), longNames(i), fbins,0,fbins)
            #res = TH1I('fresid'+str(i), longNames(i), fbins,0,fbins)
            res.SetLineColor(kBlack)
            res.SetMarkerColor(kBlack)
            res.SetLineWidth(1)
            fresid.append(res)
            

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
                        tmpscale1 = data[dkey]['druScale']
                        rldata[dkey]['hist'] = copy.deepcopy(data[dkey]['hist'])
                        rldata[dkey]['hist'].Rebin(loEfitRebin)
                        if loEfitRebinScale:
                            rldata[dkey]['hist'].Scale(1./loEfitRebin)
                            tmpscale1 = data[dkey]['druScale'] / float(loEfitRebin)
                        for n in range(fLoBins):
                            if extend:
                                fitdata[i].SetBinContent(n+1+(fbins*nc), fitdata[i].GetBinContent(n+1+(fbins*nc))
                                                         + rldata[dkey]['hist'].GetBinContent(fLo[0]+n))
                            else:
                                fitdata[i].SetBinContent(n+1, fitdata[i].GetBinContent(n+1)
                                                         + rldata[dkey]['hist'].GetBinContent(fLo[0]+n))
                    # hi E
                    if 'x'+str(i+1) in dkey and '-c'+C in dkey and '-e1' in dkey:
                        tmpscale2 = data[dkey]['druScale']
                        rhdata[dkey]['hist'] = copy.deepcopy(data[dkey]['hist'])
                        rhdata[dkey]['hist'].Rebin(hiEfitRebin)
                        if hiEfitRebinScale:
                            rhdata[dkey]['hist'].Scale(1./hiEfitRebin)
                            tmpscale2 = data[dkey]['druScale'] / float(hiEfitRebin)
                        r = 0
                        for n in range(fLoBins,fbins):
                            if extend:
                                fitdata[i].SetBinContent(n+1+(fbins*nc), fitdata[i].GetBinContent(n+1+(fbins*nc))
                                                         + rhdata[dkey]['hist'].GetBinContent(fHi[0]+r))
                            else:
                                fitdata[i].SetBinContent(n+1, fitdata[i].GetBinContent(n+1)
                                                         + rhdata[dkey]['hist'].GetBinContent(fHi[0]+r))
                            r += 1
                
                ### find average dru scale from loE and hiE fit regions
                druscale[i] = ((tmpscale1*fLoBins) + (tmpscale2*(fbins-fLoBins))) / fbins
                
                # DEBUGGING
                #--------------------------------------------------
                """
                if i==6:
                    if nc==0:
                        tmpcan = TCanvas('tmpcan','tmpcan',800,600)
                        fitdata[i].Draw()
                    else:
                        fitdata[i].Draw('same')
                    tmpcan.Update()
                    raw_input('Enter to continue')
                """
                #--------------------------------------------------
                
                for skey in sigkeys:
                    # lo E
                    if 'x'+str(i+1) in skey and '-c'+C in skey and '-e0' in skey:
                        #fskey = skey.split('-e0')[0]
                        #fskey = skey.split('-c')[0]
                        fskey = skey.split('-c'+C+'-e0')[0]
                        if sinit:
                            fsig = TH1F(fskey, fskey, fmax, 0, fmax)
                            fitsigs[fskey] = {}
                            fitsigs[fskey]['hist'] = fsig
                        rlsigs[skey]['hist'] = copy.deepcopy(sigs[skey]['hist'])
                        rlsigs[skey]['hist'].Rebin(loEfitRebin)
                        if loEfitRebinScale: rlsigs[skey]['hist'].Scale(1./loEfitRebin)
                        for n in range(fLoBins):
                            if extend:
                                fitsigs[fskey]['hist'].SetBinContent(n+1+(fbins*nc),
                                                fitsigs[fskey]['hist'].GetBinContent(n+1+(fbins*nc))
                                                + rlsigs[skey]['hist'].GetBinContent(fLo[0]+n))
                            else:
                                fitsigs[fskey]['hist'].SetBinContent(n+1,
                                                fitsigs[fskey]['hist'].GetBinContent(n+1)
                                                + rlsigs[skey]['hist'].GetBinContent(fLo[0]+n))
                    # hi E
                    if 'x'+str(i+1) in skey and '-c'+C in skey and '-e1' in skey:
                        rhsigs[skey]['hist'] = copy.deepcopy(sigs[skey]['hist'])
                        rhsigs[skey]['hist'].Rebin(hiEfitRebin)
                        if hiEfitRebinScale: rhsigs[skey]['hist'].Scale(1./hiEfitRebin)
                        #fskey = skey.split('-e1')[0]
                        #fskey = skey.split('-c')[0]
                        fskey = skey.split('-c'+C+'-e1')[0]
                        r = 0
                        for n in range(fLoBins,fbins):
                            if extend:
                                fitsigs[fskey]['hist'].SetBinContent(n+1+(fbins*nc),
                                                fitsigs[fskey]['hist'].GetBinContent(n+1+(fbins*nc))
                                                + rhsigs[skey]['hist'].GetBinContent(fHi[0]+r))
                            else:
                                fitsigs[fskey]['hist'].SetBinContent(n+1,
                                                fitsigs[fskey]['hist'].GetBinContent(n+1)
                                                + rhsigs[skey]['hist'].GetBinContent(fHi[0]+r))
                            r += 1
                sinit=0
                
                for bkey in bakkeys:
                    # lo E
                    if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e0' in bkey:
                        #fbkey = bkey.split('-e0')[0]
                        #fbkey = bkey.split('-c')[0]
                        fbkey = bkey.split('-c'+C+'-e0')[0]
                        if binit:
                            fbak = TH1F(fbkey, fbkey, fmax, 0, fmax)
                            fitbkgs[fbkey] = {}
                            fitbkgs[fbkey]['hist'] = fbak
                        rlbkgs[bkey]['hist'] = copy.deepcopy(bkgs[bkey]['hist'])
                        rlbkgs[bkey]['hist'].Rebin(loEfitRebin)
                        if loEfitRebinScale: rlbkgs[bkey]['hist'].Scale(1./loEfitRebin)
                        for n in range(fLoBins):
                            if extend:
                                fitbkgs[fbkey]['hist'].SetBinContent(n+1+(fbins*nc),
                                                fitbkgs[fbkey]['hist'].GetBinContent(n+1+(fbins*nc))
                                                + rlbkgs[bkey]['hist'].GetBinContent(fLo[0]+n))
                            else:
                                fitbkgs[fbkey]['hist'].SetBinContent(n+1,
                                                fitbkgs[fbkey]['hist'].GetBinContent(n+1)
                                                + rlbkgs[bkey]['hist'].GetBinContent(fLo[0]+n))
                    # hi E
                    if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e1' in bkey:
                        rhbkgs[bkey]['hist'] = copy.deepcopy(bkgs[bkey]['hist'])
                        rhbkgs[bkey]['hist'].Rebin(hiEfitRebin)
                        if hiEfitRebinScale: rhbkgs[bkey]['hist'].Scale(1./hiEfitRebin)
                        #fbkey = bkey.split('-e1')[0]
                        #fbkey = bkey.split('-c')[0]
                        fbkey = bkey.split('-c'+C+'-e1')[0]
                        r = 0
                        for n in range(fLoBins,fbins):
                            if extend:
                                fitbkgs[fbkey]['hist'].SetBinContent(n+1+(fbins*nc),
                                                fitbkgs[fbkey]['hist'].GetBinContent(n+1+(fbins*nc))
                                                + rhbkgs[bkey]['hist'].GetBinContent(fHi[0]+r))
                            else:
                                fitbkgs[fbkey]['hist'].SetBinContent(n+1,
                                                fitbkgs[fbkey]['hist'].GetBinContent(n+1)
                                                + rhbkgs[bkey]['hist'].GetBinContent(fHi[0]+r))
                            r += 1
                binit=0
                

            ### subtract fixed MC from data
            #-------------------------------------------------------------------
            for key in fitbkgs:
                if 'x'+str(i+1) in key:
                    fitdata[i].Add(fitbkgs[key]['hist'], -1)
            """
            if i==6:
                tmpcan = TCanvas('tmpcan','tmpcan',800,600)
                fitdata[i].Draw()
                tmpcan.Update()
                raw_input('Enter to continue')
            """
            #-------------------------------------------------------------------


        
        ### sort all the fbakkeys
        fbakkeys=[]
        for fbkey in fitbkgs:
            fbakkeys.append(fbkey)
        fbakkeys.sort()

        ### sort all the fsigkeys
        fsigkeys=[]
        for fskey in fitsigs:
            fsigkeys.append(fskey)
        fsigkeys.sort()


        sigObj = []
        fit = []
        #fitresults = []
        fitresults = {}
        fitchi2ndf = []
        wasFit = []
        
        ### set up the fitting object for TFractionFitter
        for i in range(8):

            fitresults[str(i)] = []
            fitchi2ndf.append(-1)
            
            ### Do a per crystal unique of the signals to see how many
            ### signal channels each crystal fit will have
            uniqSig = []
            noEvents = []
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    mc_int = fitsigs[fskey]['hist'].Integral(fmin,fmax)
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
                fitdata[i].Sumw2()
            
            dat_int = fitdata[i].Integral(fmin,fmax) # data integral to normalize to

            bounds = []
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    
                    mc_int = fitsigs[fskey]['hist'].Integral(fmin,fmax) # MC integral

                    ### to weight or not to weight...
                    # what if you weight first?
                    if mcsumw2:
                        fitsigs[fskey]['hist'].Sumw2() # set stat weights

                    ### normalize MC to data
                    ### needed for TFractionFitter to work right
                    try:
                        fitsigs[fskey]['hist'].Scale(dat_int/mc_int) # scale to data integral
                    except:
                        print '\nWARNING: No events for --> ',fskey
                        print   '         Remove it from the fit!\n'
                        sys.exit()
                        
                    for C in allchans:
                        sigs[fskey+'-c'+C+'-e0']['hist'].Scale(dat_int/mc_int) # make sure MC is scaled too
                        sigs[fskey+'-c'+C+'-e1']['hist'].Scale(dat_int/mc_int) # make sure MC is scaled too

                        ### v31 - save the scaling factors so you can convert to mBq/kg later
                        #sigs[fskey+'-e0']['fitscale'] = dat_int/mc_int
                        #sigs[fskey+'-e1']['fitscale'] = dat_int/mc_int

                        ### v42 version of fitscale
                        sigs[fskey+'-c'+C+'-e0']['fitscale'] = sigs[fskey+'-c'+C+'-e0']['scale'] * dat_int/mc_int
                        sigs[fskey+'-c'+C+'-e1']['fitscale'] = sigs[fskey+'-c'+C+'-e1']['scale'] * dat_int/mc_int
                            
                    
                    ### rescale the bounds to the normalized fraction
                    ### and save new values to the sigs info
                    #---------------------------------------------------------------------
                    #useBounds=0
                    for C in allchans:
                        for E in range(2):
                            E=str(E)
                            sigs[fskey+'-c'+C+'-e'+E]['info']['newfbnd'] = [0,0]
                            for k in range(2):
                                renorm = sigs[fskey+'-c'+C+'-e'+E]['scale'] / float(sigs[fskey+'-c'+C+'-e'+E]['fitscale'])
                                
                                if useBounds == 0:
                                    these = [0.00, 1.00]
                                    sigs[fskey+'-c'+C+'-e'+E]['info']['fbnd'][k] = 1./renorm * these[k]
                                    sigs[fskey+'-c'+C+'-e'+E]['info']['newfbnd'] = these
                                elif useBounds == 1:
                                    sigs[fskey+'-c'+C+'-e'+E]['info']['newfbnd'][k] = \
                                        sigs[fskey+'-c'+C+'-e'+E]['info']['fbnd'][k] * renorm
                                elif useBounds == 2:
                                    sigs[fskey+'-c'+C+'-e'+E]['info']['fbnd'] = newBounds
                                    sigs[fskey+'-c'+C+'-e'+E]['info']['newfbnd'][k] = \
                                        newBounds[k] * renorm
                                elif useBounds == 3:
                                    sigs[fskey+'-c'+C+'-e'+E]['info']['fbnd'][k] = 1./renorm * otherBnds[k]
                                    sigs[fskey+'-c'+C+'-e'+E]['info']['newfbnd'] = otherBnds
                                else:
                                    print 'ERROR: do not know what to do with useBounds =',useBounds
                                    sys.exit()
                            
                    bounds.append(sigs[fskey+'-c'+fitchans[0]+'-e0']['info']['newfbnd'])
                    #---------------------------------------------------------------------
                    
                    
                    ### set errors to zero
                    fitsigs[fskey]['hist'] = zeroBinError(fitsigs[fskey]['hist'])
                    
                    ### add the signal object
                    sigObj[-1].Add(fitsigs[fskey]['hist']) # add to the TFractionFitter object

            ### set errors to zero
            fitdata[i] = zeroBinError(fitdata[i])
                
            #fit.append(TFractionFitter(fitdata[i], sigObj[i])) # create the TFF data and MC objects
            fit.append(TFractionFitter(fitdata[i], sigObj[-1])) # create the TFF data and MC objects
            #fit.append(TFractionFitter(fitdata[i], sigObj[-1], "Q")) # create the TFF data and MC objects
            
            ### Print out the fit infos
            fitresults[str(i)].append('Crystal-'+str(i+1)+' fit results')
            if note: fitresults[str(i)].append('note = '+note)
            fitresults[str(i)].append('version = '+V)
            fitresults[str(i)].append('channels fit = '+fitchans)
            fitresults[str(i)].append('hist extend = '+str(extend))
            fitresults[str(i)].append('norm to dru = '+str(dru))
            fitresults[str(i)].append('lo-E fit range = '+str(fLoE[0])+' - '+str(fLoE[1])+' keV')
            fitresults[str(i)].append('hi-E fit range = '+str(fHiE[0])+' - '+str(fHiE[1])+' keV')
            
            ### set fit bounds!!!
            ### l=0 sets all params to the same constrain
            for l in range(len(bounds)):
                fit[-1].Constrain(l+1, bounds[l][0], bounds[l][1])

            ### set the fit range
            fit[-1].SetRangeX(fmin, fmax)

            ### do the fit
            status = fit[-1].Fit()
            ### do MINOS style error analysis?
            #fit[-1].ErrorAnalysis(1)
            wasFit.append(i)

            chi2 = fit[-1].GetChisquare()
            chi2 = chi2 / druscale[i]
            ndf  = fit[-1].GetNDF()
            pval = fit[-1].GetProb()
            
            fitresults[str(i)].append('returned fit status = '+str(status))
            fitchi2ndf[-1] = (chi2/ndf)
            fitresults[str(i)].append('chi2/ndf = %.3g/%s = %.3g'%(chi2,ndf,chi2/ndf))
            
            count = 0
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    fscale = ROOT.Double(0.0)
                    ferror = ROOT.Double(0.0)
                    fit[-1].GetResult(count, fscale, ferror)
                    fitsigs[fskey]['hist'].Scale(fscale)
                    
                    for C in allchans:
                        for E in range(2):
                            E = str(E)

                            ### save the raw scaling factor from the fit
                            sigs[fskey+'-c'+C+'-e'+E]['hist'].Scale(fscale)
                            
                            ### save converted scaling factor
                            sigs[fskey+'-c'+C+'-e'+E]['fitscale'] = sigs[fskey+'-c'+C+'-e'+E]['fitscale'] * fscale
                            
                            ### set error as a percent of the scaling factor
                            try:
                                sigs[fskey+'-c'+C+'-e'+E]['fiterror'] = ferror/fscale
                            except:
                                sigs[fskey+'-c'+C+'-e'+E]['fiterror'] = 1.
                            
                    count += 1

        ### scale the signals to mBq/kg
        if dru: sigs = scaleSigs71(sigkeys, sigs)
        else: sigs = scaleSigs71(sigkeys, sigs, runtime)
        
        ### print the fit activities
        for i in range(8):
            if i not in wasFit:
                continue
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    finit=1
                    for C in allchans:
                        for E in range(2):
                            if finit:
                                E = str(E)
                                ### print out activity and error and bounds
                                fitacti = sigs[fskey+'-c'+C+'-e'+E]['info']['fitacti']
                                fiterro = sigs[fskey+'-c'+C+'-e'+E]['info']['fiterro']
                                lobnd = sigs[fskey+'-c'+C+'-e'+E]['info']['acti'] * \
                                        sigs[fskey+'-c'+C+'-e'+E]['info']['fbnd'][0]
                                hibnd = sigs[fskey+'-c'+C+'-e'+E]['info']['acti'] * \
                                        sigs[fskey+'-c'+C+'-e'+E]['info']['fbnd'][1]
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
            fitresults[str(i)].append('\n')

            
        save = ''
        #if local: save += 'local'
        #else:     save += 'on-cup'
        save += str(runtag)
        save += '_Nchan-fit'
        save += '_loEfit-'+str(int(fLoE[0]))+'-'+str(int(fLoE[1]))
        save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        save += '_loEfitRebin-'+str(loEfitRebin)
        save += '_hiEfitRebin-'+str(hiEfitRebin)
        save += '_loEfitRebinScale'+str(loEfitRebinScale)
        save += '_hiEfitRebinScale'+str(hiEfitRebinScale)
        save += '_useBounds'+str(useBounds)
        save += '_mcsumw2'+str(mcsumw2)
        save += '_datsumw2'+str(datsumw2)
        save += '_dru'+str(dru)
        save += '_loEplotRebin-'+str(loEplotRebin)
        save += '_hiEplotRebin-'+str(hiEplotRebin)
        save += '_reuse'+str(reuse)
        save += '_chans'+str(fitchans)
        save += '_extend'+str(extend)
        if note: save += '_'+note
        save += '_'+V
        

        ### print out the results and update the backgrounds file
        #-------------------------------------------------------------
        ### sort all the fitresults
        resultskeys=[]
        for rskey in fitresults:
            resultskeys.append(rskey)
        resultskeys.sort()

        ### write results to file
        resultsfile = './plots/'+save+'_fit-results.txt'
        outfile = open(resultsfile, 'w')
        for key in resultskeys:
            for line in fitresults[key]:
                outfile.write(line+'\n')
        outfile.close()
        
        ### create the updated backgrounds file
        if updateMCfile:
            shutil.copyfile(mcfile, './plots/'+mcfile)
            newbkgs = './plots/'+mcfile[:-4]+'-update.txt'
            updateBkgsFile70(mcfile, resultsfile, newbkgs, BF='BR')
        
        ### create the background model table
        #outtable = newbkgs[:-4]+'-table.txt'
        #outputModelTable61(newbkgs, outtable)
        #-------------------------------------------------------------
        
        
        
        # plot the multi-chan fit results
        #=================================================================

        fcanv = TCanvas('fcanv', 'fcanv', 0, 0, 1400, 900)
        fcanv.Divide(4,2)

        ftoppad=[]
        fbotpad=[]

        #sepFitPlots=[]
        
        flegs=[]
        flegs2=[]
        fzeros=[]

        font=63
        size=13
        yoff=4.2
        
        for i in range(8):
            fcanv.cd(i+1)
            
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
            
            if dru: fitdata[i].SetAxisRange(2e-3, 2e3, 'y')
            
            newFitTitle = str('Crystal-'+str(i+1)+'   '+'Fit-chans-'+fitchans)
            fitdata[i].SetTitle(newFitTitle)

            fitdata[i].SetLineColor(kBlack)
            fitdata[i].SetMarkerColor(kBlack)
            fitdata[i].SetLineWidth(1)

            if dru: fitdata[i].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
            else: fitdata[i].GetYaxis().SetTitle('arb. counts')

            fitdata[i].GetYaxis().SetTitleFont(font)
            fitdata[i].GetYaxis().SetTitleSize(size)
            fitdata[i].GetYaxis().SetTitleOffset(yoff)
            fitdata[i].GetYaxis().SetLabelFont(font)
            fitdata[i].GetYaxis().SetLabelSize(size)
            fitdata[i].GetYaxis().SetLabelOffset(0.01)
            #fitdata[i].GetXaxis().SetTitle('Energy (keV)')
            #fitdata[i].GetXaxis().SetLabelFont(font)
            #fitdata[i].GetXaxis().SetLabelSize(size)
            fitdata[i].Draw()
            #flegs[i].AddEntry(fitdata[i], 'data - bkgs', lopt)

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
                    ftotal[i].Add(fitsigs[fskey]['hist'])

            
            ### build legend after getting the numbers of signals
            Nfs = Nfsigs+2
            flnc = 1
            if Nfs > 6: flnc = 2
            if Nfs > 12: flnc = 3
            fxlegstop  = 0.94
            fxlegstart = fxlegstop-(0.2*flnc)
            fylegstop  = 0.89
            fylegstart = fylegstop-(0.04*6)
            
            fleg = TLegend(fxlegstart, fylegstart, fxlegstop, fylegstop)
            flegs.append(fleg)
            flegs[i].SetNColumns(flnc)
            flegs[i].SetFillColor(0)
            flegs[i].SetBorderSize(0)
            lopt = 'LPE'
            
            flegs[i].AddEntry(fitdata[i], 'data - bkgs', lopt)
            
            ### add legend entries in order
            for name in uniqAll:
                for fskey in fsigkeys:
                    if name in fskey and 'x'+str(i+1) in fskey:
                        flegs[i].AddEntry(fitsigs[fskey]['hist'], fskey, lopt)
            
            
            ### get the chi2 of the total fit mc compared to data
            #-------------------------------------------------------------------
            chi2  = ROOT.Double(0.0)
            ndf   = ROOT.Long(0)
            igood = ROOT.Long(0)
            fitchi2ndfv2 = -1
            if i in wasFit:
                pval = fitdata[i].Chi2TestX(ftotal[i], chi2, ndf, igood, chiopt)
                fitchi2ndfv2 = chi2/druscale[i]/ndf
            #-------------------------------------------------------------------
            
            ftotal[i].Draw('same')
            
            if i in wasFit:
                # returned from fit
                flegs[i].AddEntry(ftotal[i], 'Fit Total (chi2/ndf = '+str(round(fitchi2ndf[i],2))+')', lopt)
                # calc by Chi2TestX()
                #flegs[i].AddEntry(ftotal[i], 'Fit Total (chi2/ndf = '+str(round(fitchi2ndfv2,2))+')', lopt)
                print 'INFO: Fit total MC from Chi2TestX chi2/ndf = '+str(round(fitchi2ndfv2,2))
                
                flegs[i].Draw('same')
            
            
            ### plot the fit residuals
            #-------------------------------------------------------------
            fbotpad[i].cd()
            leg = TLegend(0.72, 0.78, 0.94, 0.94)
            flegs2.append(leg)
            flegs2[i].SetFillColor(0)
            flegs2[i].SetBorderSize(0)
            lopt = 'LPE'

            fresid[i].Divide(fitdata[i], ftotal[i])
            #fresid[i].Divide(ftotal[i], fitdata[i])

            fresid[i].SetTitle('')
            #fresid[i].SetXTitle('Energy (keVee)')
            fresid[i].SetXTitle('fit bins')
            fresid[i].GetXaxis().SetTitleFont(font)
            fresid[i].GetXaxis().SetTitleSize(size)
            fresid[i].GetXaxis().SetTitleOffset(8)
            fresid[i].GetXaxis().SetLabelFont(font)
            fresid[i].GetXaxis().SetLabelSize(size)
            fresid[i].GetXaxis().SetLabelOffset(0.03)
            #fresid[i].SetYTitle('counts / keV')
            fresid[i].SetYTitle('data / MC')
            #fresid[i].SetYTitle('MC / data')
            fresid[i].GetYaxis().SetTitleFont(font)
            fresid[i].GetYaxis().SetTitleSize(size)
            fresid[i].GetYaxis().SetTitleOffset(yoff)
            fresid[i].GetYaxis().SetLabelFont(font)
            fresid[i].GetYaxis().SetLabelSize(size)
            fresid[i].GetYaxis().SetLabelOffset(0.01)
            # '5' secondary and '05' primary
            fresid[i].GetYaxis().SetNdivisions(505)
            
            fresid[i].SetAxisRange(0.1,10,'y')
            if linres:
                fbotpad[i].SetLogy(0)
                fresid[i].SetAxisRange(lrs[0], lrs[1], 'y')
            fresid[i].Draw()
            
            ### set my reference line to '1'
            zero = TLine(fmin, 1, fmax, 1)
            fzeros.append(zero)
            fzeros[i].SetLineColor(kRed)
            fzeros[i].SetLineWidth(1)
            fzeros[i].Draw()

            flegs2[i].AddEntry(fresid[i],'data / MC',lopt)
            #flegs2[i].Draw()
            #-------------------------------------------------------------

            fcanv.Update()
            if indi and i+1 in justthese:
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
                fisave += '_loEfRS'+str(loEfitRebinScale)
                fisave += '_hiEfRS'+str(hiEfitRebinScale)
                fisave += '_dru'+str(dru)
                fisave += '_cs'+str(fitchans)
                fisave += '_ext'+str(extend)
                if note: fisave += '_'+str(note)
                fisave += '_'+str(V)
                
                try: sepFitPlot.Print(str('./plots/'+fisave+'.png'))
                except: pass
            try: del sepFitPlot
            except: pass
            
        fcanv.Update()
        fcanv.Print('./plots/'+save+'.png')
    
    ### end of fitting bit if you have signals
    
    
    
    # plot the lo and hi energy histograms for all channels
    #=================================================================
    
    # number of energy ranges (lo, hi)
    numE = 2
    # number of channels (all, single, multi, combos)
    numC = len(pltchans)
    
    #canvs  = [[[] for x in range(numE)] for x in range(numC)]
    canvs  = [[0 for x in range(numE)] for x in range(numC)]
    
    ### for separate plots
    #sepPlots = [[[[] for x in range(8)] for x in range(numE)] for x in range(numC)]
    sepPlots = [[[0 for x in range(8)] for x in range(numE)] for x in range(numC)]
        
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
    
    gbkgs  = [[[{} for x in range(8)] for x in range(numE)] for x in range(numC)]
    gsigs  = [[[{} for x in range(8)] for x in range(numE)] for x in range(numC)]
        
    plotRebin = 1
    for C, chan in enumerate(pltchans): 
    
        for E in range(numE):

            if E: plotRebin = hiEplotRebin
            else: plotRebin = loEplotRebin
            
            # have the plotting be seperated out from the 8 crystal loop
            canvs[C][E] = TCanvas('canv'+chan+str(E),
                                  'canv'+chan+str(E),
                                  0, 0, 1400, 900)
            canvs[C][E].Divide(4,2)

            toppad[C][E] = []
            botpad[C][E] = []

            legs[C][E]   = []
            legs2[C][E]  = []
            zeros[C][E]  = []

            total[C][E]  = makeTotal64(chan, E, params[E])
            resid[C][E]  = makeResid64(chan, E, params[E])
            
            font=63
            size=13
            yoff=4.2

            for i in range(8):
                
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
                    leg = TLegend(0.35, 0.75, 0.94, 0.89)
                    lnc = 3
                else:
                    leg = TLegend(xlegstart, ylegstart, xlegstop, ylegstop)
                legs[C][E].append(leg)
                legs[C][E][i].SetFillColor(0)
                legs[C][E][i].SetBorderSize(0)
                legs[C][E][i].SetNColumns(lnc)
                lopt = 'LPE'
                
                total[C][E][i].Rebin(plotRebin)
                #total[C][E][i].Scale(1./float(plotRebin))
                #total[C][E][i].Sumw2()
                if dru:
                    #total[C][E][i].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
                    if E: total[C][E][i].SetAxisRange(2e-3, 2e1, 'y')
                    else: total[C][E][i].SetAxisRange(2e-3, 3e2, 'y')
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

                        if dru:
                            #data[dkey]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            if E: data[dkey]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            else: data[dkey]['hist'].SetAxisRange(2e-3, 3e2, 'y')
                        
                        #popt = 'P E1'
                        #popt = 'HIST'
                        popt = ''
                        data[dkey]['hist'].Draw()
                        #days = round(data[dkey]['runtime']/86400.,2)
                        days = round(runtime/86400.,2)
                        
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
                        
                        if i in wasFit:
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
                    c=0
                    for group in gbkgs[C][E][i]:
                        if group != 'none':
                            gbkgs[C][E][i][group].SetMarkerColor(gis[c])
                            gbkgs[C][E][i][group].SetLineColor(gis[c])
                            gbkgs[C][E][i][group].Draw('same')
                            legs[C][E][i].AddEntry(gbkgs[C][E][i][group], group, lopt)
                            c+=1
                    for group in gbkgs[C][E][i]:
                        if group == 'none':
                            for j, key in enumerate(gbkgs[C][E][i]['none']):
                                gbkgs[C][E][i]['none'][key].SetMarkerColor(gis[c+j])
                                gbkgs[C][E][i]['none'][key].SetLineColor(gis[c+j])
                                gbkgs[C][E][i]['none'][key].Draw('same')
                                legs[C][E][i].AddEntry(gbkgs[C][E][i]['none'][key], key, lopt)
                else:
                    # add legend entries in order
                    for name in uniqAll:
                        for bkey in bakkeys:
                            if bkey == 'x'+str(i+1)+'-'+name+'-c'+chan+'-e'+str(E):
                                activ = '(%.2e) '%(bkgs[bkey]['info']['acti'])
                                legs[C][E][i].AddEntry(bkgs[bkey]['hist'], activ+bkey, lopt)
                        for skey in sigkeys:
                            if skey == 'x'+str(i+1)+'-'+name+'-c'+chan+'-e'+str(E):
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
                resid[C][E][i].GetXaxis().SetTitleOffset(8)
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
            save += '_E'+str(E)
            save += '_loEfit-'+str(int(fLoE[0]))+'-'+str(int(fLoE[1]))
            save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
            save += '_loEfitRebin-'+str(loEfitRebin)
            save += '_hiEfitRebin-'+str(hiEfitRebin)
            save += '_loEfitRebinScale'+str(loEfitRebinScale)
            save += '_hiEfitRebinScale'+str(hiEfitRebinScale)
            save += '_useBounds'+str(useBounds)
            save += '_mcsumw2'+str(mcsumw2)
            save += '_datsumw2'+str(datsumw2)
            save += '_dru'+str(dru)
            save += '_loEplotRebin-'+str(loEplotRebin)
            save += '_hiEplotRebin-'+str(hiEplotRebin)
            save += '_reuse'+str(reuse)
            save += '_chans'+fitchans
            save += '_chan'+chan
            save += '_extend'+str(extend)
            if note: save += '_'+note
            save += '_'+V
            
            canvs[C][E].Update()
            canvs[C][E].Print('./plots/'+save+'.png')
            
            
            ### Save separate crystal histos?
            #-------------------------------------------------------------
            if indi:
                for i in range(8):
                    if i+1 in justthese:
                        tpad=toppad[C][E][i].Clone()

                        #tpad.SetLogy(0)
                        #tpad.SetAxisRange(0,10,'y')
                        
                        bpad=botpad[C][E][i].Clone()
                        sepPlots[C][E][i] = TCanvas('ican-'+str(chan)+str(E)+str(i+1),
                                                    'ican-'+str(chan)+str(E)+str(i+1),
                                                    0, 0, 1400, 900)
                        #sepPlots[C][E][i].cd()
                        tpad.Draw()
                        bpad.Draw()
                        sepPlots[C][E][i].Update()
                        isave  = ''
                        isave += 'x'+str(i+1)
                        isave += '-cs'+str(fitchans)
                        isave += '-c'+str(chan)
                        isave += '-e'+str(E)
                        isave += '-ext'+str(extend)
                        if note: isave += '_'+str(note)
                        isave += '_'+str(V)

                        try: sepPlots[C][E][i].Print(str('./plots/'+isave+'.png'))
                        except: pass

    
    ### Save 4 plots to one canvas
    #-------------------------------------------------------------
    if indi:
        combPlots = [0 for x in range(8)]
        for i in range(8):
            if i+1 in justthese:
                combPlots[i] = TCanvas('ccan-'+str(i+1),'ccan-'+str(i+1),0,0,1400,900)
                combPlots[i].Divide(2,2)
                tpad = [0 for x in range(4)]
                bpad = [0 for x in range(4)]
                p=0
                for C, chan in enumerate(pltchans):
                    for E in range(numE):
                        tpad[p]=toppad[C][E][i].Clone()
                        bpad[p]=botpad[C][E][i].Clone()
                        combPlots[i].cd(p+1)
                        tpad[p].Draw()
                        bpad[p].Draw()
                        combPlots[i].Update()
                        p += 1
                csave  = ''
                csave += 'x'+str(i+1)
                csave += '-combined'
                if note: csave += '_'+str(note)
                csave += '_'+str(V)

                try: combPlots[i].Print(str('./plots/'+csave+'.png'))
                except: pass

    #-----------------------------------------------------------------
    ### print out the fit results
    if fitting:
        print '\n\n'
        print '!!!!!  FIT RESULTS  !!!!!\n'
        for key in resultskeys:
            for line in fitresults[key]:
                print line
    #-----------------------------------------------------------------

    # delete the extra crap
    #-------------------------------
    #try: del combPlots
    #except: pass

    try: del sepPlots
    except: pass

    try: del canvs
    except: pass
        
    
    if not batch:
        raw_input('[Enter] to quit \n')


######################################################################
######################################################################

if __name__ == "__main__":
    myself(sys.argv[1:])

