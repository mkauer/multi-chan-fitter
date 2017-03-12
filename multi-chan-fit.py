#!/usr/bin/env python

######################################################################
# Matt Kauer - mkauer@physics.wisc.edu
######################################################################
# 52-whats-next.py

V = 'v52'

# not sure what's next, probably a lot of bug fixing that will present
#  itself when I start adding more MC
# 
# version: 2017-03-12
#
# note: run 1616 is the first run after calibration-campaign-2
# 
# see CHANGELOG for changes
######################################################################

import os,sys
import socket
import copy
import math
import numpy as np
import ROOT
from ROOT import *

from funcs50 import *


### user inputs
#================================================================

### extra notes to add to the saved plot file names? [0, 'something']
note=0
#note = 'find-data'

### backgrounds file to use?
#mcfile = 'backgrounds41.txt'
mcfile = 'backgrounds50.txt'
#mcfile = 'back-just-data.txt'
#mcfile = 'back-just-c3.txt'

### force reuse of all joined rootfiles in mcfile? [0,1,2]
### nice for debugging
reuse = 0
### [0] default - use whatever is specified in the backgrounds file
### [1] forces reusing of all data/bkgs/sigs
### [2] forces NOT reusing any data/bkgs/sigs
reuse = 1

### force a particular set of hit chan data? [0,1,2,3]
### nice for debugging
mychans = 0
### [ 0 ] default - use whatever is specified in the backgrounds file
### ['A'] force all-hit data selection channel
### ['S'] force single-hit data selection channel
### ['M'] force multi-hi data selection channel
# seems to help doing 'MS' over 'SM' ???
mychans = 'MS'

### plotting ranges
### lo and hi energy ranges
loer = [0, 200]
hier = [0, 3000]

### individual plots for all crystals? [0,1]
indi = 1
### just plot individual for crystals? [1-8]
justthese = [3]

### rebin the hi-E final plots [1,inf]
hiEplotRebin = 10

### rebin the hi-E histo for fitting [1,inf]
hiEfitRebin = 10

### scale to dru?
dru = 1

### This seems to be a must for TFF to weight all MC equally
### pre scale the MC? [0,1]
mcscale = 1

### This doesn't seem to effect the fit results at all
### set MC sumw2()? [0,1]
mcsumw2 = 0
### set data sumw2()? [0,1]
datsumw2 = 0

### I don't understand why this is effecting the fit so much
### But it does help the fit converage and especially with multi-chan
### scale data and mc hi-E hists to 1/hiEfitRebin factor? [0,1]
fitRebinScale = 1
#================================================================


### automated selections...
#==================================================
### figure out if running on laptop or at CUP
### master.cunpa.ibs
local = amLocal()
### 0 = show-plots -- 1 = dont-show-plots
batch = 0
if not local:
    batch = 1
gROOT.SetBatch(batch)
#==================================================


def _myself_(argv):
    
    gROOT.Reset()
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
        print 'ERROR: could not find backgrounds file -->', mcfile
        sys.exit()
    
    data, bkgs, sigs = build50(mcfile, reuse, mychans)
    datkeys, bakkeys, sigkeys = sortKeys2(data, bkgs, sigs)
    
    # assume all data is using same run and channels
    runNum = data[datkeys[0]]['info']['run']
    chans  = data[datkeys[0]]['info']['chans']
    
    ### find unique names for color scheme?
    ### "internal-K40" for example
    uniqAll = []
    for key in bakkeys:
        uniqAll.append(key.split('-')[1]+'-'+key.split('-')[2])
    for key in sigkeys:
        uniqAll.append(key.split('-')[1]+'-'+key.split('-')[2])
    uniqAll = sorted(list(set(uniqAll)))
    print 'INFO: Unique bkgs and sigs =',uniqAll
    
    # scale into dru units
    if dru: data = dataDRU40(data)
    bkgs = scaleBkgs40(bkgs)
    sigs = scaleBkgs40(sigs)
    
    ### Number of colors
    Nc = len(uniqAll)
    print 'INFO: Total number of unique bkgs and sigs =',Nc
    colors, cis = rainbow(Nc)

    ### Create color dict for unique simulations
    uniqColor = {}
    for i,key in enumerate(uniqAll):
        uniqColor[key] = cis[i]
        
    ### legend length = MC + data + total
    Nlg = Nc+2
    
    
    ### set fit bounds for the fit
    #=================================================================
    #=================================================================
    
    #fLo = [10, 100]
    fLo = [20, 180]
    fLoBins = fLo[1]-fLo[0]
    
    #fHiE = [400, 2000]
    fHiE = [300, 2400]
    fHi = [fHiE[0]/hiEfitRebin, fHiE[1]/hiEfitRebin]
    fHiBins = fHi[1]-fHi[0]
    
    fbins = fLoBins+fHiBins
    
    fmin=0
    fmax=fbins
    #=================================================================
    #=================================================================
    

    ### need seperate dicts for rebinned data and MC for plotting to
    ### to work right - they need their own memory space
    #-----------------------------------------------------------------
    ### init dict for rebinned data
    rdata = {}
    for dkey in datkeys:
        #print dkey
        rdata[dkey] = {}
        
    ### init dict for rebinned backgrounds
    rbkgs = {}
    for bkey in bakkeys:
        #print bkey
        rbkgs[bkey] = {}

    # init dict for rebinned MC/signals 
    rsigs = {}
    for skey in sigkeys:
        #print skey
        rsigs[skey] = {}

    #-----------------------------------------------------------------
    
    
    ### here's the fitting part!!!
    ##################################################################
    ### only do the fit if you have signals and data!
    if len(sigs) > 0:
        
        fitdata = []
        fitsigs = {}
        fitbkgs = {}
        ftotal = [] # total of fixed MC plus fit MC
        fresid = [] # residual of the data/total
        
        ### fill fitdata and fitsigs
        for i in range(8):
            
            fdata = TH1F('fdata'+str(i), cnames(i)+' - multi-chan fit', fbins,0,fbins)
            fitdata.append(fdata)
            
            tot = TH1F('ftotal'+str(i), longNames(i), fbins,0,fbins)
            tot.SetLineColor(kGray+1)
            tot.SetMarkerColor(kGray+1)
            tot.SetLineWidth(1)
            ftotal.append(tot)
            
            res = TH1F('fresid'+str(i), longNames(i), fbins,0,fbins)
            res.SetLineColor(kBlack)
            res.SetMarkerColor(kBlack)
            res.SetLineWidth(1)
            fresid.append(res)
            

            # cycle through the channels
            sinit = 1
            binit = 1
            for C in chans:

                # build the histos for fitting
                #-----------------------------------------------------------------------------
                for dkey in datkeys:
                    if 'x'+str(i+1) in dkey and '-c'+C in dkey and '-e0' in dkey:
                        for n in range(fLoBins):
                            fitdata[i].SetBinContent(n+1, fitdata[i].GetBinContent(n+1)
                                                     + data[dkey]['hist'].GetBinContent(fLo[0]+n))
                    if 'x'+str(i+1) in dkey and '-c'+C in dkey and '-e1' in dkey:
                        rdata[dkey]['hist'] = copy.deepcopy(data[dkey]['hist'])
                        rdata[dkey]['hist'].Rebin(hiEfitRebin)
                        if fitRebinScale: rdata[dkey]['hist'].Scale(1./hiEfitRebin)
                        r = 0
                        for n in range(fLoBins,fbins):
                            fitdata[i].SetBinContent(n+1, fitdata[i].GetBinContent(n+1)
                                                     + rdata[dkey]['hist'].GetBinContent(fHi[0]+r))
                            r += 1
                
                for skey in sigkeys:
                    if 'x'+str(i+1) in skey and '-c'+C in skey and '-e0' in skey:
                        #fskey = skey.split('-e0')[0]
                        fskey = skey.split('-c')[0]
                        if sinit:
                            fsig = TH1F(fskey, fskey, fbins, 0, fbins)
                            fitsigs[fskey] = {}
                            fitsigs[fskey]['hist'] = fsig
                        for n in range(fLoBins):
                            fitsigs[fskey]['hist'].SetBinContent(n+1, fitsigs[fskey]['hist'].GetBinContent(n+1)
                                                                 + sigs[skey]['hist'].GetBinContent(fLo[0]+n))
                    if 'x'+str(i+1) in skey and '-c'+C in skey and '-e1' in skey:
                        rsigs[skey]['hist'] = copy.deepcopy(sigs[skey]['hist'])
                        rsigs[skey]['hist'].Rebin(hiEfitRebin)
                        if fitRebinScale: rsigs[skey]['hist'].Scale(1./hiEfitRebin)
                        #fskey = skey.split('-e1')[0]
                        fskey = skey.split('-c')[0]
                        r = 0
                        for n in range(fLoBins,fbins):
                            fitsigs[fskey]['hist'].SetBinContent(n+1, fitsigs[fskey]['hist'].GetBinContent(n+1)
                                                                 + rsigs[skey]['hist'].GetBinContent(fHi[0]+r))
                            r += 1
                sinit=0
                
                for bkey in bakkeys:
                    if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e0' in bkey:
                        #fbkey = bkey.split('-e0')[0]
                        fbkey = bkey.split('-c')[0]
                        if binit:
                            fbak = TH1F(fbkey, fbkey, fbins, 0, fbins)
                            fitbkgs[fbkey] = {}
                            fitbkgs[fbkey]['hist'] = fbak
                        for n in range(fLoBins):
                            fitbkgs[fbkey]['hist'].SetBinContent(n+1, fitbkgs[fbkey]['hist'].GetBinContent(n+1)
                                                                 + bkgs[bkey]['hist'].GetBinContent(fLo[0]+n))
                    if 'x'+str(i+1) in bkey and '-c'+C in bkey and '-e1' in bkey:
                        rbkgs[bkey]['hist'] = copy.deepcopy(bkgs[bkey]['hist'])
                        rbkgs[bkey]['hist'].Rebin(hiEfitRebin)
                        if fitRebinScale: rbkgs[bkey]['hist'].Scale(1./hiEfitRebin)
                        #fbkey = bkey.split('-e1')[0]
                        fbkey = bkey.split('-c')[0]
                        r = 0
                        for n in range(fLoBins,fbins):
                            fitbkgs[fbkey]['hist'].SetBinContent(n+1, fitbkgs[fbkey]['hist'].GetBinContent(n+1)
                                                                 + rbkgs[bkey]['hist'].GetBinContent(fHi[0]+r))
                            r += 1
                binit=0
                
                
                ### subtract fixed MC from data
                #-------------------------------------------------------------------
                for key in fitbkgs:
                    if 'x'+str(i+1) in key:
                        fitdata[i].Add(fitbkgs[key]['hist'], -1)
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
        
        ### set up the fitting object for TFractionFitter
        for i in range(8):

            fitresults[str(i)] = []
            
            ### Do a per crystal unique of the signals to see how many
            ### signal channels each crystal fit will have
            uniqSig = []
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    uniqSig.append(fskey.split('-')[1]+'-'+fskey.split('-')[2])
            uniqSig = sorted(list(set(uniqSig)))
            
            ### only fit a crystal that has 2 or more signals
            ### TFF wants at least 2 MC to converge the fit
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

                    ### to scale or not to scale...
                    if mcscale:
                        fitsigs[fskey]['hist'].Scale(dat_int/mc_int) # scale to data integral

                        for C in chans:
                            sigs[fskey+'-c'+C+'-e0']['hist'].Scale(dat_int/mc_int) # make sure MC is scaled too
                            sigs[fskey+'-c'+C+'-e1']['hist'].Scale(dat_int/mc_int) # make sure MC is scaled too
                            
                            ### v31 - save the scaling factors so you can convert to mBq/kg later
                            #sigs[fskey+'-e0']['fitscale'] = dat_int/mc_int
                            #sigs[fskey+'-e1']['fitscale'] = dat_int/mc_int
                            
                            ### v42 version of fitscale
                            sigs[fskey+'-c'+C+'-e0']['fitscale'] = sigs[fskey+'-c'+C+'-e0']['scale'] * dat_int/mc_int
                            sigs[fskey+'-c'+C+'-e1']['fitscale'] = sigs[fskey+'-c'+C+'-e1']['scale'] * dat_int/mc_int
                            
                    else:
                        ### why doesn't this work??
                        ### maybe it does now??
                        for C in chans:
                            sigs[fskey+'-c'+C+'-e0']['fitscale'] = sigs[fskey+'-c'+C+'-e0']['scale']
                            sigs[fskey+'-c'+C+'-e1']['fitscale'] = sigs[fskey+'-c'+C+'-e1']['scale']


                    #renorm = sigs[fskey+'-c'+C+'-e0']['scale'] / sigs[fskey+'-c'+C+'-e0']['fitscale']
                    renorm = sigs[fskey+'-c'+chans[0]+'-e0']['scale'] / sigs[fskey+'-c'+chans[0]+'-e0']['fitscale']
                    ### use bounds from backgrounds files
                    #bounds.append([renorm * sigs[fskey+'-c'+C+'-e0']['info']['fbnd'][0],
                    #               renorm * sigs[fskey+'-c'+C+'-e0']['info']['fbnd'][1]])
                    bounds.append([renorm * sigs[fskey+'-c'+chans[0]+'-e0']['info']['fbnd'][0],
                                   renorm * sigs[fskey+'-c'+chans[0]+'-e0']['info']['fbnd'][1]])
                    
                    #sigObj[i].Add(fitsigs[fskey]['hist']) # add to the TFractionFitter object
                    sigObj[-1].Add(fitsigs[fskey]['hist']) # add to the TFractionFitter object
            

            #fit.append(TFractionFitter(fitdata[i], sigObj[i])) # create the TFF data and MC objects
            fit.append(TFractionFitter(fitdata[i], sigObj[-1])) # create the TFF data and MC objects
            #fit.append(TFractionFitter(fitdata[i], sigObj[-1], "Q")) # create the TFF data and MC objects
            
            ### for some reason on CUP the fit status gets returned as a pointer
            ### like <ROOT.TFitResultPtr object at 0x37246c0>
            #fitresults.append('NaI-C'+str(i+1)+' fit completed with status '+str(status))
            #fitresults.append('NaI-C'+str(i+1)+' fit results')
            fitresults[str(i)].append('NaI-C'+str(i+1)+' fit results')
            #fitresults.append('MC fit constrained to 0.0001 - 10.0')
            
            ### set fit bounds!!!
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #for l in range(len(uniqSig)):
            ### seems like Constrain starts with param=1
            ### l=0 sets all params to the same constrain
            ### very confusing
            for l in range(len(bounds)):
                ### use numbers from backgrounds file
                #fit[-1].Constrain(l+1, bounds[l][0], bounds[l][1])
                ### or force other less strict contraints
                fit[-1].Constrain(l+1, 0.01, 10.0)
            
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            
            #fit[i].SetRangeX(fmin, fmax) # set global range, should be the same as fmin and fmax?
            fit[-1].SetRangeX(fmin, fmax) # set global range, should be the same as fmin and fmax?


            ### try doing the fit
            #------------------------------
            #status = fit[i].Fit()
            status = fit[-1].Fit()
            #------------------------------


            #chi2 = fit[i].GetChisquare()
            #ndf  = fit[i].GetNDF()
            #pval = fit[i].GetProb()
            chi2 = fit[-1].GetChisquare()
            ndf  = fit[-1].GetNDF()
            pval = fit[-1].GetProb()

            fitchi2ndf.append(chi2/ndf)

            ### for some reason on CUP the fit status gets returned as a pointer
            ### like <ROOT.TFitResultPtr object at 0x37246c0>
            #fitresults.append('chi2/ndf = '+str(round(chi2,2))+'/'+str(ndf)+' = '+str(round(chi2/float(ndf),2)))
            fitresults[str(i)].append('chi2/ndf = '+str(round(chi2,2))+'/'+str(ndf)+' = '+str(round(chi2/float(ndf),2)))

            ### p-val is always 0 - need to look into this
            #fitresults.append('p-value = '+str(pval))
            fitresults[str(i)].append('p-value = '+str(pval))

            count = 0
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    fscale = ROOT.Double(0.0)
                    ferror = ROOT.Double(0.0)
                    fit[-1].GetResult(count, fscale, ferror)
                    #fitresults.append(fskey+' = '+str(round(fscale,4))+' +/- '+str(round(ferror,4)))
                    fitresults[str(i)].append('fit scale  -  '+fskey+' = '+str(round(fscale,4))+' +/- '+str(round(ferror,4)))
                    fitsigs[fskey]['hist'].Scale(fscale)

                    for C in chans:
                        for E in range(2):
                            E = str(E)
                            #fscale = ROOT.Double(0.0)
                            #ferror = ROOT.Double(0.0)
                            #fit[-1].GetResult(count, fscale, ferror)
                            #fitresults.append(fskey+' = '+str(round(fscale,4))+' +/- '+str(round(ferror,4)))
                            #fitsigs[fskey]['hist'].Scale(fscale)
                            
                            sigs[fskey+'-c'+C+'-e'+E]['hist'].Scale(fscale)
                            #sigs[fskey+'-c'+C+'-e1']['hist'].Scale(fscale)
                            
                            ### v31 - save the scaling factors so you can convert to mBq/kg later
                            sigs[fskey+'-c'+C+'-e'+E]['fitscale'] = sigs[fskey+'-c'+C+'-e'+E]['fitscale'] * fscale
                            #sigs[fskey+'-c'+C+'-e1']['fitscale'] = sigs[fskey+'-c'+C+'-e1']['fitscale'] * fscale
                            
                            sigs[fskey+'-c'+C+'-e'+E]['fiterror'] = ferror
                            #sigs[fskey+'-c'+C+'-e1']['fiterror'] = ferror

                    count += 1

            #fitresults.append('\n')
        
        
        ### print the fit results
        #print '\n\n'
        #for line in fitresults:
        #    print line
        
        ### scale the signals to mBq/kg
        sigs = scaleSigs40(sigkeys, sigs)
        
        ### print the fit activities
        for i in range(8):
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:
                    finit=1
                    for C in chans:
                        for E in range(2):
                            if finit:
                                E = str(E)
                                #print fskey, sigs[fskey+'-c'+C+'-e'+E]['info']['acti']
                                fitresults[str(i)].append('activity   -  '+fskey+' '
                                                          +str(sigs[fskey+'-c'+C+'-e'+E]['info']['acti'])+' mBq')
                                finit=0
            #print '\n'
            fitresults[str(i)].append('\n')

            
        #-------------------------------------------------------------
        ### sort all the fitresults
        #resultskeys=[]
        #for rskey in fitresults:
        #    resultskeys.append(rskey)
        #resultskeys.sort()

        ### print out the goods
        #for key in resultskeys:
        #    for line in fitresults[key]:
        #        print line
        #-------------------------------------------------------------
        

        save = ''
        if local:         save += 'local'
        else:             save += 'cup'
        save += '_'+str(runNum)
        save += '_Nchan-fit'
        save += '_loEfit-'+str(int(fLo[0]))+'-'+str(int(fLo[1]))
        save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        save += '_hiEfitRebin-'+str(hiEfitRebin)
        if fitRebinScale: save += '_rebinscale'
        if mcscale:       save += '_mcscale'
        if mcsumw2:       save += '_mcsumw2'
        if datsumw2:      save += '_datsumw2'
        if dru:           save += '_dru'
        if hiEplotRebin:  save += '_hiEplotRebin-'+str(hiEplotRebin)
        save += '_reuse'+str(reuse)
        save += '_chans'+str(chans)
        if note:          save += '_'+note
        save += '_'+V

        
        #-------------------------------------------------------------
        ### sort all the fitresults
        resultskeys=[]
        for rskey in fitresults:
            resultskeys.append(rskey)
        resultskeys.sort()

        ### write results to file
        outfile = open('./plots/'+save+'_fit-results.txt', 'w')
        #for line in fitresults:
        #    outfile.write(str(line)+'\n')
        for key in resultskeys:
            for line in fitresults[key]:
                outfile.write(line+'\n')
        outfile.close()    
        #-------------------------------------------------------------
        
        
        
        # plot the multi-chan fit results
        #=================================================================

        fcanv = TCanvas('fcanv', 'fcanv', 0, 0, 1400, 900)
        fcanv.Divide(4,2)

        ftoppad=[]
        fbotpad=[]

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
            fbotpad[i].SetLogy()
            ftoppad[i].Draw()
            fbotpad[i].Draw()
            ftoppad[i].cd()
            
            ylegstart = 0.88
            ylegend = (ylegstart-(Nlg*0.035))
            space = '  '
            leg = TLegend(0.55, ylegend, 0.94, ylegstart)
            flegs.append(leg)
            flegs[i].SetFillColor(0)
            flegs[i].SetBorderSize(0)
            lopt = 'LPE'

            if dru: fitdata[i].SetAxisRange(2e-2, 3e2, 'y')
            
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
            flegs[i].AddEntry(fitdata[i], space+'data - bkgs', lopt)

            
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:

                    # find the unique name for color and set color
                    cname = fskey.split('-')[1]+'-'+fskey.split('-')[2]
                    fitsigs[fskey]['hist'].SetMarkerColor(uniqColor[cname])
                    fitsigs[fskey]['hist'].SetLineColor(uniqColor[cname])
                    
                    ### draw the sigs
                    fitsigs[fskey]['hist'].Draw('same')

                    ### add MC to total MC hist
                    ftotal[i].Add(fitsigs[fskey]['hist'])

            
            # add legend entries in order
            for name in uniqAll:
                #for fbkey in fbakkeys:
                #    if name in fbkey and 'x'+str(i+1) in fbkey:
                #        flegs[i].AddEntry(fitbkgs[fbkey]['hist'], space+fbkey, lopt)
                for fskey in fsigkeys:
                    if name in fskey and 'x'+str(i+1) in fskey:
                        flegs[i].AddEntry(fitsigs[fskey]['hist'], space+fskey, lopt)
                        #activ = ' ('+str(sigs[name]['acti'])+')  '
                        #flegs[i].AddEntry(fitsigs[fskey]['hist'], activ+fskey, lopt)

            
            ### get the chi2 of the total fit mc compared to data
            #---------------------------------------------------------
            chi2  = ROOT.Double(0.0)
            ndf   = ROOT.Long(0)
            igood = ROOT.Long(0)
            pval  = fitdata[i].Chi2TestX(ftotal[i], chi2, ndf, igood, "WW")
            print 'INFO:','fit = crystal-'+str(i+1),'pval =',pval,'chi2 =',chi2,'ndf =',ndf,'igood =',igood
            fitchi2ndfv2=chi2/ndf
            #---------------------------------------------------------

            
            ftotal[i].Draw('same')
            #flegs[i].AddEntry(ftotal[i], space+'Fit Total', lopt)
            ### chi2/ndf from the fit results
            #flegs[i].AddEntry(ftotal[i], space+'Fit Total (chi2/ndf = '+str(round(fitchi2ndf[i],2))+')', lopt)
            ### chi2/ndf from the Chi2TestX function
            flegs[i].AddEntry(ftotal[i], space+'Fit Total (chi2/ndf = '+str(round(fitchi2ndfv2,2))+')', lopt)
            flegs[i].Draw('same')
            

            ### fit residuals plots
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
            fresid[i].GetYaxis().SetNdivisions(505) # '5' secondary and '05' primary
            
            fresid[i].SetAxisRange(0.1,10,'y')
            fresid[i].Draw()
            
            ### set my line to '1'
            zero = TLine(fmin, 1, fmax, 1)
            fzeros.append(zero)
            fzeros[i].SetLineColor(kRed)
            fzeros[i].SetLineWidth(1)
            fzeros[i].Draw()

            flegs2[i].AddEntry(fresid[i],space+'data / MC',lopt)
            #flegs2[i].Draw()
            #-------------------------------------------------------------
            
            
        fcanv.Update()
        fcanv.Print('./plots/'+save+'.png')
        
    ### end of fitting bit if you have signals
    
    
    # plot the lo and hi energy histograms for all channels
    #=================================================================

    # number of energy ranges (lo, hi)
    numE = 2
    # number of channels (all, single, multi, combos)
    numC = len(chans)
    
    canvs  = [[[] for x in range(numE)] for x in range(numC)]

    ### for separate plots
    #sepPlots = [[] for x in range(numE)]
    #sepPlots = [[[] for x in range(8)] for x in range(numE)]
    sepPlots = [[[[] for x in range(8)] for x in range(numE)] for x in range(numC)]
    
    ### seperate memory space for the pads is key!!!!
    toppad = [[[] for x in range(numE)] for x in range(numC)]
    botpad = [[[] for x in range(numE)] for x in range(numC)]

    legs   = [[[] for x in range(numE)] for x in range(numC)]
    legs2  = [[[] for x in range(numE)] for x in range(numC)]
    zeros  = [[[] for x in range(numE)] for x in range(numC)]

    total  = [[[] for x in range(numE)] for x in range(numC)]
    resid  = [[[] for x in range(numE)] for x in range(numC)]


    for C, chan in enumerate(chans): 
        #print '!!!!!!!!!!!!!!!!!!!',C,chan
    
        for E in range(numE):

            # have the plotting be seperated out from the loop
            canvs[C][E] = TCanvas('canv'+chan+str(E), 'canv'+chan+str(E), 0, 0, 1400, 900)
            canvs[C][E].Divide(4,2)

            toppad[C][E] = []
            botpad[C][E] = []

            legs[C][E]   = []
            legs2[C][E]  = []
            zeros[C][E]  = []

            total[C][E]  = makeTotal50(chan,E)
            resid[C][E]  = makeResid50(chan,E)

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
                toppad[C][E][i].SetLogy()
                botpad[C][E][i].SetTopMargin(0.05)
                botpad[C][E][i].SetBottomMargin(0.3)
                botpad[C][E][i].SetBorderMode(0)
                botpad[C][E][i].SetLogy()
                toppad[C][E][i].Draw()
                botpad[C][E][i].Draw()
                toppad[C][E][i].cd()

                ylegstart = 0.88
                ylegend = (ylegstart-(Nlg*0.035))
                space = '  '
                leg = TLegend(0.55, ylegend, 0.94, ylegstart)
                legs[C][E].append(leg)
                legs[C][E][i].SetFillColor(0)
                legs[C][E][i].SetBorderSize(0)
                lopt = 'LPE'

                for dkey in datkeys:
                    if 'x'+str(i+1) in dkey and '-c'+chan in dkey and '-e'+str(E) in dkey:
                        data[dkey]['hist'].GetYaxis().SetTitle('arb. counts')
                        if dru:
                            data[dkey]['hist'].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')

                #---------------------------------------------------------
                if E and hiEplotRebin:
                    for dkey in datkeys:
                        if 'x'+str(i+1) in dkey and '-c'+chan in dkey and '-e'+str(E) in dkey:
                            data[dkey]['hist'].Rebin(hiEplotRebin)
                            data[dkey]['hist'].Scale(1./float(hiEplotRebin))

                    total[C][E][i].Rebin(hiEplotRebin)
                    total[C][E][i].Scale(1./float(hiEplotRebin))
                    #total[C][E][i].Sumw2()

                    resid[C][E][i].Rebin(hiEplotRebin)
                    resid[C][E][i].Scale(1./float(hiEplotRebin))
                    #resid[C][E][i].Sumw2()
                #---------------------------------------------------------
                for key in datkeys:
                    if 'x'+str(i+1) in key and '-c'+chan in key and '-e'+str(E) in key:

                        ### needed for Chi2TestX()
                        ### and needed later too
                        dkey = key

                        #if dru1 or dru2:
                            #druScale = data['x'+str(i+1)+'-data'+'-e'+str(E)]['druScale']
                            #druScale = data[key]['druScale']

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
                            if E: data[dkey]['hist'].SetAxisRange(2e-3, 2e1, 'y')
                            else: data[dkey]['hist'].SetAxisRange(2e-2, 3e2, 'y')

                        data[dkey]['hist'].Draw()
                        days = round(data[dkey]['runtime']/86400.,2)
                        #legs[C][E][i].AddEntry(data[key]['hist'], 'data', lopt)
                        #legs[C][E][i].AddEntry(data[dkey]['hist'], 'data ('+str(days)+' days)', lopt)
                        legs[C][E][i].AddEntry(data[dkey]['hist'], dkey+' ('+str(days)+' days)', lopt)

                for key in bakkeys:
                    if 'x'+str(i+1) in key and '-c'+chan in key and '-e'+str(E) in key:

                        # find the unique name for color and set color
                        cname = key.split('-')[1]+'-'+key.split('-')[2]
                        bkgs[key]['hist'].SetMarkerColor(uniqColor[cname])
                        bkgs[key]['hist'].SetLineColor(uniqColor[cname])

                        # don't think I need this anymore
                        #if dru1 or dru2:
                        #    druScale = data['x'+str(i+1)+'-data'+'-e'+str(E)]['druScale']
                        #if dru2:
                        #    bkgs[key]['hist'].Scale(druScale)

                        if E and hiEplotRebin:
                            bkgs[key]['hist'].Rebin(hiEplotRebin)
                            bkgs[key]['hist'].Scale(1./float(hiEplotRebin))
                            #bkgs[key]['hist'].Sumw2()

                        ### temp - don't draw MC if you just want to show data
                        bkgs[key]['hist'].Draw('same')

                        ### add MC to total MC hist
                        total[C][E][i].Add(bkgs[key]['hist'])

                        ### create the legend entry for MC
                        #legs[C][E][i].AddEntry(bkgs[key]['hist'], space+key, lopt)

                for key in sigkeys:
                    if 'x'+str(i+1) in key and '-c'+chan in key and '-e'+str(E) in key:

                        # find the unique name for color and set color
                        cname = key.split('-')[1]+'-'+key.split('-')[2]
                        sigs[key]['hist'].SetMarkerColor(uniqColor[cname])
                        sigs[key]['hist'].SetLineColor(uniqColor[cname])

                        # don't think I need this anymore
                        #if dru1 or dru2:
                        #    druScale = data['x'+str(i+1)+'-data'+'-e'+str(E)]['druScale']
                        #if dru2:
                        #    sigs[key]['hist'].Scale(druScale)

                        if E and hiEplotRebin:
                            sigs[key]['hist'].Rebin(hiEplotRebin)
                            sigs[key]['hist'].Scale(1./float(hiEplotRebin))
                            #sigs[key]['hist'].Sumw2()

                        ### set range
                        #sigs[key]['hist'].SetAxisRange(1,1000,'y')

                        ### draw sigs
                        sigs[key]['hist'].Draw('same')

                        ### add MC to total MC hist
                        total[C][E][i].Add(sigs[key]['hist'])

                        ### create the legend entry for MC
                        #legs[C][E][i].AddEntry(sigs[key]['hist'], space+key, lopt)


                # add legend entries in order
                for name in uniqAll:
                    for bkey in bakkeys:
                        if name in bkey and 'x'+str(i+1) in bkey and '-c'+chan in bkey and '-e'+str(E) in bkey:
                            #legs[C][E][i].AddEntry(bkgs[bkey]['hist'], space+bkey, lopt)
                            activ = '('+str(bkgs[bkey]['info']['acti'])+') '
                            legs[C][E][i].AddEntry(bkgs[bkey]['hist'], activ+bkey, lopt)
                    for skey in sigkeys:
                        if name in skey and 'x'+str(i+1) in skey and '-c'+chan in skey and '-e'+str(E) in skey:
                            #legs[C][E][i].AddEntry(sigs[skey]['hist'], space+skey, lopt)
                            activ = '('+str(sigs[skey]['info']['acti'])+') '
                            legs[C][E][i].AddEntry(sigs[skey]['hist'], activ+skey, lopt)


                ### you need to scale the error by the dru scaling and/or the rebinning
                #-----------------------------------------------------------------------------
                #total[C][E][i].Sumw2()

                if not dru:
                    for n in range(total[C][E][i].GetNbinsX()):
                        total[C][E][i].SetBinError(n+1, total[C][E][i].GetBinError(n+1)/(float(hiEplotRebin)/math.sqrt(2.)))

                else:
                    for n in range(total[C][E][i].GetNbinsX()):
                        #total[C][E][i].SetBinError(n+1, total[C][E][i].GetBinError(n+1)*druScale)
                        total[C][E][i].SetBinError(n+1, total[C][E][i].GetBinError(n+1)*data[dkey]['druScale'])

                #if dru2:
                #    for n in range(total[C][E][i].GetNbinsX()):
                #        total[C][E][i].SetBinError(n+1, total[C][E][i].GetBinError(n+1)*druScale)


                ### chi2 test
                ### get the chi2 of the total mc compared to data
                #=============================================================================
                #-----------------------------------------------------------------------------
                chi2  = ROOT.Double(0.0)
                ndf   = ROOT.Long(0)
                igood = ROOT.Long(0)
                pval  = data[dkey]['hist'].Chi2TestX(total[C][E][i], chi2, ndf, igood, "WW")
                print 'INFO:',dkey,'pval =',pval,'chi2 =',chi2,'ndf =',ndf,'igood =',igood
                print 'INFO: Total MC chi2/ndf =',chi2/ndf
                #-----------------------------------------------------------------------------
                #=============================================================================



                #-----------------------------------------------------------------------------
                if len(bkgs) + len(sigs) > 0:        
                    total[C][E][i].Draw('same')
                    #legs[C][E][i].AddEntry(total[C][E][i], 'Total MC', lopt)
                    legs[C][E][i].AddEntry(total[C][E][i], 'Total MC (chi2/ndf = '+str(round(chi2/ndf,2))+')', lopt)

                legs[C][E][i].Draw('same')


                ### try to get the residuals in!
                #---------------------------------------------------------
                botpad[C][E][i].cd()
                leg = TLegend(0.72, 0.78, 0.94, 0.94)
                legs2[C][E].append(leg)
                legs2[C][E][i].SetFillColor(0)
                legs2[C][E][i].SetBorderSize(0)
                lopt = 'LPE'

                #resid[C][E][i].Divide(data['x'+str(i+1)+'-data'+'-e'+str(E)]['hist'], total[C][E][i])
                resid[C][E][i].Divide(data[dkey]['hist'], total[C][E][i])
                #resid[C][E][i].Divide(total[C][E][i], data['x'+str(i+1)+'-data'+'-e'+str(E)]['hist'])

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
                resid[C][E][i].GetYaxis().SetNdivisions(505) # '5' secondary and '05' primary

                if E: resid[C][E][i].SetAxisRange(hier[0], hier[1], 'x')
                else: resid[C][E][i].SetAxisRange(loer[0], loer[1], 'x')

                resid[C][E][i].SetAxisRange(0.1,10,'y')
                resid[C][E][i].Draw()

                par = histparam(E)
                # set my line to '1'
                zero = TLine(par[1], 1, par[2], 1)
                zeros[C][E].append(zero)
                zeros[C][E][i].SetLineColor(kRed)
                zeros[C][E][i].SetLineWidth(1)
                zeros[C][E][i].Draw()

                #legs2[C][E][i].AddEntry(resid[C][E][i],space+'data / MC',lopt)
                #legs2[C][E][i].Draw()
                #---------------------------------------------------------


                ### get the chi2 of the total mc compared to data
                #---------------------------------------------------------
                #chi2  = ROOT.Double(0.0)
                #ndf   = ROOT.Long(0)
                #igood = ROOT.Long(0)
                #pval  = data[dkey]['hist'].Chi2TestX(total[C][E][i], chi2, ndf, igood, "WW")
                #print 'INFO:',dkey,'pval =',pval,'chi2 =',chi2,'ndf =',ndf,'igood =',igood
                #print 'INFO: Total MC chi2/ndf =',chi2/ndf
                #---------------------------------------------------------


            save = ''
            if local:         save += 'local'
            else:             save += 'cup'
            save += '_'+str(runNum)
            if E:             save += '_hiE'
            else:             save += '_loE'
            save += '_loEfit-'+str(int(fLo[0]))+'-'+str(int(fLo[1]))
            save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
            save += '_hiEfitRebin-'+str(hiEfitRebin)
            if fitRebinScale: save += '_rebinscale'
            if mcscale:       save += '_mcscale'
            if mcsumw2:       save += '_mcsumw2'
            if datsumw2:      save += '_datsumw2'
            if dru:           save += '_dru'
            if hiEplotRebin:  save += '_hiEplotRebin-'+str(hiEplotRebin)
            save += '_reuse'+str(reuse)
            save += '_chans'+chans
            save += '_chan'+chan
            if note:          save += '_'+note
            save += '_'+V
            
            canvs[C][E].Update()
            canvs[C][E].Print('./plots/'+save+'.png')


            ### Save separate crystal histos?
            #-------------------------------------------------------------
            if indi:
                for i in range(8):
                    if i+1 in justthese:
                        tpad=toppad[C][E][i].Clone()
                        bpad=botpad[C][E][i].Clone()
                        sepPlots[C][E][i] = TCanvas('ican'+chan+str(E)+str(i), 'ican'+chan+str(E)+str(i), 0, 0, 1400, 900)
                        #sepPlots[C][E][i].cd()
                        tpad.Draw()
                        bpad.Draw()
                        sepPlots[C][E][i].Update()
                        isave  = ''
                        isave += 'x'+str(i+1)
                        isave += '-cs'+chans
                        isave += '-c'+chan
                        isave += '-e'+str(E)
                        if note: isave += '_'+note
                        isave += '_'+V

                        sepPlots[C][E][i].Print('./plots/'+isave+'.png')

    ### but don't show all those plots
    if indi: del sepPlots

    
    #-----------------------------------------------------------------
    ### print out the fit results
    print '\n\n'
    print '!!!!!  FIT RESULTS  !!!!!\n'
    for key in resultskeys:
        for line in fitresults[key]:
            print line
    #-----------------------------------------------------------------
    
    
    if not batch:
        raw_input('[Enter] to quit \n')


######################################################################
######################################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])

