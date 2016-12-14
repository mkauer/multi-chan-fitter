#!/usr/bin/env python

######################################################################
# Matt Kauer - mkauer@physics.wisc.edu
######################################################################
# 21-get-activity-values.py

V = 'v21'

# Scale the MC by real activity values
#
# version: 2016-12-13
# 
# see CHANGELOG for changes
######################################################################

import os,sys
import socket
import numpy as np
from ROOT import *
import ROOT
from funcs2 import *
import copy
import math


### user inputs
#================================================================
### I'm pretty sure you want to convert to dru after the fit
### More stats for the fit will give bettet errors
### ie use dru2 option!!
### dru1 has some bugs still... the fits won't converge
dru1 = 1   ### convert data to dru before fit? [0,1]
dru2 = 0   ### convert data and mc to dru after fit? [0,1]

# dru safty check...
if dru2: dru1 = 0

hiEplotRebin = 10  ### rebin the hi-E final plots [1,inf]

# something weird with the fit rebinning - need to FIX this ASAP!
# use rebin=1 for now
hiEfitRebin = 10   ### rebin the hi-E histo for fitting

# fill the fitdata and fitmc starting with n plus number?
# i think this should be 1 because bin 0 is underflows right?
# yes, just checked, this should be 1
np = 1

reuse = 1   ### use joined rootfile data? [0,1]

mcscale = 1   ### pre scale the MC? [0,1]
mcweight = 1  ### set MC weights? [0,1]

# this doesn't work at all, not sure why...
fitweight = 0   ### set fit weights? [0,1]

# still working on this - I don't understand why
# this is effecting the fit so much
fitRebinScale = 0   ### scale data and mc hi-E hists to 1/hiEfitRebin factor? [0,1]
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
    
    #runNum = 1324
    runNum = 1544
    mcfile = 'backgrounds.txt'
    
    if reuse:
        rootfile = './root-join-read/join2-'+str(runNum)+'-master.root'
        data, bkgs, sigs = readROOT2(rootfile, mcfile)
    
    else:
        data = getData2(runNum)
        bkgs, sigs = buildMC2(mcfile)
    
    datkeys, bakkeys, sigkeys = sortKeys2(data, bkgs, sigs)
    
    if dru1:
        data = dataDRU2(data)
            
    
    # create a color scheme for MC
    # color list and indexes
    # divide by 2 because e0 and e1 of same component
    # divide by 8 because there's 8 crystals
    Nbkgs = len(bkgs)/2/8
    Nsigs = len(sigs)/2/8
    # number of colors
    Nc = Nbkgs + Nsigs
    print 'number of bkgs and sigs =',Nc
    colors, cis = rainbow(Nc)
    
    # legend length = MC + data + total
    Nlg = Nc+2
    
    ### fit the MC to data
    #-----------------------------------------------------------------
    
    fLo = [10, 100]
    fLoBins = fLo[1]-fLo[0]
    
    fHiE = [400, 2000]
    fHi = [fHiE[0]/hiEfitRebin, fHiE[1]/hiEfitRebin]
    fHiBins = fHi[1]-fHi[0]
    print fHi
    print fHiBins
    
    fbins = fLoBins+fHiBins
    
    fmin=0
    fmax=fbins

        
    ### need seperate dicts for rebinned data and MC for plotting to
    ### to work right - they need their own memory space
    #-----------------------------------------------------------------
    ### init dict for rebinned data
    rdata = {}
    for key in datkeys:
        #print key
        rdata[key] = {}
        
    ### init dict for rebinned backgrounds
    rbkgs = {}
    for key in bakkeys:
        #print key
        rbkgs[key] = {}

    # init dict for rebinned MC/signals 
    rsigs = {}
    for key in sigkeys:
        #print key
        rsigs[key] = {}

    #-----------------------------------------------------------------
    
    
    # pre-scale all backgrounds for testing purposes
    # generally needed if not using tons of data
    if not reuse:
        for key in bakkeys:
            bkgs[key]['hist'].Scale(0.001)

    
    ### only do the fit if you have signals!
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


            # fill the low energy part
            #------------------------------------------------
            for key in datkeys:
                if 'x'+str(i+1) in key and '-e0' in key:
                    for n in range(fLoBins):
                        fitdata[i].SetBinContent(n+np, data[key]['hist'].GetBinContent(fLo[0]+n))
                        fitdata[i].SetBinError(n+np, data[key]['hist'].GetBinError(fLo[0]+n))

            for key in sigkeys:
                if 'x'+str(i+1) in key and '-e0' in key:
                    fskey = key.split('-e0')[0]
                    fsig = TH1F(fskey, fskey, fbins, 0, fbins)
                    fitsigs[fskey] = {}
                    fitsigs[fskey]['hist'] = fsig
                    for n in range(fLoBins):
                        fitsigs[fskey]['hist'].SetBinContent(n+np, sigs[key]['hist'].GetBinContent(fLo[0]+n))

            for key in bakkeys:
                if 'x'+str(i+1) in key and '-e0' in key:
                    fbkey = key.split('-e0')[0]
                    fbak = TH1F(fbkey, fbkey, fbins, 0, fbins)
                    fitbkgs[fbkey] = {}
                    fitbkgs[fbkey]['hist'] = fbak
                    for n in range(fLoBins):
                        fitbkgs[fbkey]['hist'].SetBinContent(n+np, bkgs[key]['hist'].GetBinContent(fLo[0]+n))


            # fill the high energy part
            #-------------------------------------------------------------

            for key in datkeys:
                if 'x'+str(i+1) in key and '-e1' in key:
                    rdata[key]['hist'] = copy.deepcopy(data[key]['hist'])
                    rdata[key]['hist'].Rebin(hiEfitRebin)
                    if fitRebinScale: rdata[key]['hist'].Scale(1./hiEfitRebin)

            for key in sigkeys:
                if 'x'+str(i+1) in key and '-e1' in key:
                    rsigs[key]['hist'] = copy.deepcopy(sigs[key]['hist'])
                    rsigs[key]['hist'].Rebin(hiEfitRebin)
                    if fitRebinScale: rsigs[key]['hist'].Scale(1./hiEfitRebin)

            for key in bakkeys:
                if 'x'+str(i+1) in key and '-e1' in key:
                    rbkgs[key]['hist'] = copy.deepcopy(bkgs[key]['hist'])
                    rbkgs[key]['hist'].Rebin(hiEfitRebin)
                    if fitRebinScale: rbkgs[key]['hist'].Scale(1./hiEfitRebin)


            for key in datkeys:
                if 'x'+str(i+1) in key and '-e1' in key:
                    r = 0
                    for n in range(fLoBins,fbins):
                        fitdata[i].SetBinContent(n+np, rdata[key]['hist'].GetBinContent(fHi[0]+r))
                        fitdata[i].SetBinError(n+np, rdata[key]['hist'].GetBinError(fHi[0]+r))
                        r += 1

            for key in sigkeys:
                if 'x'+str(i+1) in key and '-e1' in key:
                    fskey = key.split('-e1')[0]
                    r = 0
                    for n in range(fLoBins,fbins):
                        fitsigs[fskey]['hist'].SetBinContent(n+np, rsigs[key]['hist'].GetBinContent(fHi[0]+r))
                        r += 1

            for key in bakkeys:
                if 'x'+str(i+1) in key and '-e1' in key:
                    fbkey = key.split('-e1')[0]
                    r = 0
                    for n in range(fLoBins,fbins):
                        fitbkgs[fbkey]['hist'].SetBinContent(n+np, rbkgs[key]['hist'].GetBinContent(fHi[0]+r))
                        r += 1


            ### subtract fixed MC from data
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ### this works I think, but then this plots the data-fixedMC
            for key in fitbkgs:
                if 'x'+str(i+1) in key:
                    fitdata[i].Add(fitbkgs[key]['hist'], -1)
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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


        ### make weights histo
        ### testing with TFractionFitter
        ### this doesn't seem to work yet...
        weights = TH1F('weights','weights',fbins,0,fbins)
        for n in range(0,fbins):
            weights.SetBinContent(n, 1)
        weights.Sumw2()
        #-----------------------------------------------------------------


        sigObj = []
        fit = []
        fitresults = []

        ### set up the fitting object for TFractionFitter
        for i in range(8):

            sigObj.append(TObjArray(Nsigs)) # number of MC to fit to

            fitdata[i].Sumw2()
            dat_int = fitdata[i].Integral(fmin,fmax) # data integral to normalize to

            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:

                    mc_int = fitsigs[fskey]['hist'].Integral(fmin,fmax) # MC integral

                    ### to weight or not to weight...
                    # what if you weight first?
                    if mcweight:
                        fitsigs[fskey]['hist'].Sumw2() # set stat weights

                    ### to scale or not to scale...
                    if mcscale:
                        fitsigs[fskey]['hist'].Scale(dat_int/mc_int) # scale to data integral
                        sigs[fskey+'-e0']['hist'].Scale(dat_int/mc_int) # make sure MC is scaled too
                        sigs[fskey+'-e1']['hist'].Scale(dat_int/mc_int) # make sure MC is scaled too


                    sigObj[i].Add(fitsigs[fskey]['hist']) # add to the TFractionFitter object

            fit.append(TFractionFitter(fitdata[i], sigObj[i])) # create the TFF data and MC objects

            ### for some reason on CUP the fit status gets returned as a pointer
            ### like <ROOT.TFitResultPtr object at 0x37246c0>
            #fitresults.append('NaI-C'+str(i+1)+' fit completed with status '+str(status))
            fitresults.append('NaI-C'+str(i+1)+' fit results')
            fitresults.append('MC fit constrained to 0.0001 - 10.0')

            for l in range(Nsigs):
                # set fit constraints on the MC put into sigObj TObject
                fit[i].Constrain(l, 0.0001, 10.0)

                # set fit weights to 1
                # this doesn't work yet, not sure why, on the to-do list...
                if fitweight:
                    fit[i].SetWeight(l, weights)
                    #fit[i].SetWeight(l, 1)

            fit[i].SetRangeX(fmin, fmax) # set global range, should be the same as fmin and fmax?


            ### try doing the fit
            #------------------------------
            status = fit[i].Fit()
            #------------------------------


            chi2 = fit[i].GetChisquare()
            ndf  = fit[i].GetNDF()
            pval = fit[i].GetProb()

            ### for some reason on CUP the fit status gets returned as a pointer
            ### like <ROOT.TFitResultPtr object at 0x37246c0>
            fitresults.append('chi2/ndf = '+str(round(chi2,2))+'/'+str(ndf)+' = '+str(round(chi2/float(ndf),2)))

            ### p-val is always 0 - need to look into this
            #fitresults.append('p-value = '+str(pval))

            count = 0
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:

                    value = ROOT.Double(0.0)
                    error = ROOT.Double(0.0)
                    fit[i].GetResult(count, value, error)
                    fitresults.append(fskey+' = '+str(round(value,4))+' +/- '+str(round(error,4)))
                    fitsigs[fskey]['hist'].Scale(value)
                    sigs[fskey+'-e0']['hist'].Scale(value)
                    sigs[fskey+'-e1']['hist'].Scale(value)

                    count += 1

            fitresults.append('\n')

        print '\n\n'
        for line in fitresults:
            print line

        ### write fit results to text file
        #-----------------------------------------------------------------
        save = ''
        if local:        save += 'local'
        else:            save += 'cup'
        save += '_'+str(runNum)
        save += '_Nchan-fit'
        save += '_loEfit-'+str(int(fLo[0]))+'-'+str(int(fLo[1]))
        save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        save += '_hiEfitRebin-'+str(hiEfitRebin)
        if fitRebinScale:   save += '_rebinscale'
        if mcscale:      save += '_mcscale'
        if mcweight:     save += '_mcweight'
        if fitweight:    save += '_fitweight'
        if dru1:         save += '_dru1'
        if dru2:         save += '_dru2'
        if hiEplotRebin: save += '_hiEplotRebin-'+str(hiEplotRebin)
        if reuse:        save += '_reuse'
        save += '_'+V

        outfile = open(save+'_fit-results.txt', 'w')
        for line in fitresults:
            outfile.write(str(line)+'\n')
        outfile.close()    



        # plot the multi-chan fit results
        #=================================================================

        flegs=[]
        flegs2=[]
        fzeros=[]

        fcanv = TCanvas('fcanv', 'fcanv', 0, 0, 1400, 900)
        fcanv.Divide(4,2)

        ftoppad=[]
        fbotpad=[]

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

            if dru1:
                fitdata[i].SetAxisRange(10**-1, 10**2, 'y')
            else:
                fitdata[i].SetAxisRange(10**1, 10**4, 'y')

            fitdata[i].SetLineColor(kBlack)
            fitdata[i].SetMarkerColor(kBlack)
            fitdata[i].SetLineWidth(1)

            if dru1:
                fitdata[i].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
            else:
                fitdata[i].GetYaxis().SetTitle('arb. counts')

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
            flegs[i].AddEntry(fitdata[i], space+'data', lopt)

            color = 0
            #for loc in locs:
            #    for iso in isos:
            for fbkey in fbakkeys:
                if 'x'+str(i+1) in fbkey:

                    ### set bak colors
                    fitbkgs[fbkey]['hist'].SetMarkerColor(cis[color])
                    fitbkgs[fbkey]['hist'].SetLineColor(cis[color])
                    color += 1

                    #fitsigs[loc][iso]['hist'][i].SetAxisRange(1,1000,'y')

                    ### draw the bkgs
                    fitbkgs[fbkey]['hist'].Draw('same')

                    ### add MC to total MC hist
                    ftotal[i].Add(fitbkgs[fbkey]['hist'])

                    ### create the legend entry for MC
                    flegs[i].AddEntry(fitbkgs[fbkey]['hist'], space+fbkey, lopt)

            color = Nbkgs
            for fskey in fsigkeys:
                if 'x'+str(i+1) in fskey:

                    ### set mc colors
                    fitsigs[fskey]['hist'].SetMarkerColor(cis[color])
                    fitsigs[fskey]['hist'].SetLineColor(cis[color])
                    color += 1

                    #fitsigs[loc][iso]['hist'][i].SetAxisRange(1,1000,'y')

                    ### draw the sigs
                    fitsigs[fskey]['hist'].Draw('same')

                    ### add MC to total MC hist
                    ftotal[i].Add(fitsigs[fskey]['hist'])

                    ### create the legend entry for MC
                    flegs[i].AddEntry(fitsigs[fskey]['hist'], space+fskey, lopt)

            ftotal[i].Draw('same')
            flegs[i].AddEntry(ftotal[i], space+'Total MC', lopt)
            flegs[i].Draw('same')


            ### fit residuals plots
            #-------------------------------------------------------------
            fbotpad[i].cd()
            leg = TLegend(0.72, 0.78, 0.94, 0.94)
            flegs2.append(leg)
            flegs2[i].SetFillColor(0)
            flegs2[i].SetBorderSize(0)
            lopt = 'LPE'

            ### subtract total-mc from data
            #fresid[i].Add(data[E][i], total[E][i], 1, -1)
            ### divide data by total-mc
            fresid[i].Divide(fitdata[i], ftotal[i])

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
        fcanv.Print(save+'.png')

    ### end of fitting bit if you have signals
    
    
    # plot the lo and hi energy histograms
    #=================================================================
    
    canvs  = [None for x in range(2)]
    if dru2:
        data = dataDRU2(data)
    
    total  = [None for x in range(2)]
    resid  = [None for x in range(2)]
    
    legs   = [None for x in range(2)]
    legs2  = [None for x in range(2)]
    zeros  = [None for x in range(2)]


    
    for E in range(2):
        
        legs[E]=[]
        legs2[E]=[]
        zeros[E]=[]

        total[E] = makeTotal(E)
        resid[E] = makeResid(E)
                    
        # have the plotting be seperated out from the loop
        canvs[E] = TCanvas('canv'+str(E), 'canv'+str(E), 0, 0, 1400, 900)
        canvs[E].Divide(4,2)
        
        toppad=[]
        botpad=[]

        font=63
        size=13
        yoff=4.2
        
        for i in range(8):
            canvs[E].cd(i+1)

            fraction = 0.3
            pad1 = TPad('pad1','pad1',0,fraction,1,1)
            toppad.append(pad1)
            pad2 = TPad('pad2','pad2',0,0,1,fraction)
            botpad.append(pad2)
            toppad[i].SetBottomMargin(0.01)
            toppad[i].SetBorderMode(0)
            toppad[i].SetLogy()
            botpad[i].SetTopMargin(0.05)
            botpad[i].SetBottomMargin(0.3)
            botpad[i].SetBorderMode(0)
            botpad[i].SetLogy()
            toppad[i].Draw()
            botpad[i].Draw()
            toppad[i].cd()

            ylegstart = 0.88
            ylegend = (ylegstart-(Nlg*0.035))
            space = '  '
            leg = TLegend(0.55, ylegend, 0.94, ylegstart)
            legs[E].append(leg)
            legs[E][i].SetFillColor(0)
            legs[E][i].SetBorderSize(0)
            lopt = 'LPE'
            
            for key in datkeys:
                if 'x'+str(i+1) in key and '-e'+str(E) in key:
                    data[key]['hist'].GetYaxis().SetTitle('arb. counts')
                    if dru1 or dru2:
                        data[key]['hist'].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
                            
            #---------------------------------------------------------
            if E and hiEplotRebin:
                for key in datkeys:
                    if 'x'+str(i+1) in key and '-e'+str(E) in key:
                        data[key]['hist'].Rebin(hiEplotRebin)
                        data[key]['hist'].Scale(1./float(hiEplotRebin))
                                        
                total[E][i].Rebin(hiEplotRebin)
                total[E][i].Scale(1./float(hiEplotRebin))
                #total[E][i].Sumw2()
                
                resid[E][i].Rebin(hiEplotRebin)
                resid[E][i].Scale(1./float(hiEplotRebin))
                #resid[E][i].Sumw2()
            #---------------------------------------------------------
            for key in datkeys:
                if 'x'+str(i+1) in key and '-e'+str(E) in key:
                    data[key]['hist'].GetYaxis().SetTitleFont(font)
                    data[key]['hist'].GetYaxis().SetTitleSize(size)
                    data[key]['hist'].GetYaxis().SetTitleOffset(yoff)
                    data[key]['hist'].GetYaxis().SetLabelFont(font)
                    data[key]['hist'].GetYaxis().SetLabelSize(size)
                    data[key]['hist'].GetYaxis().SetLabelOffset(0.01)
                    #data[key]['hist'].GetXaxis().SetTitle('Energy (keV)')
                    #data[key]['hist'].GetXaxis().SetLabelFont(font)
                    #data[key]['hist'].GetXaxis().SetLabelSize(size)
                    
                    if dru1 or dru2:
                        if not E:
                            data[key]['hist'].SetAxisRange(10**-1, 10**3, 'y')
                        else:
                            data[key]['hist'].SetAxisRange(10**-2, 10**2, 'y')
                    
                    data[key]['hist'].Draw()
                    legs[E][i].AddEntry(data[key]['hist'], space+'data', lopt)
                    
            color = 0
            for key in bakkeys:
                if 'x'+str(i+1) in key and '-e'+str(E) in key:
                    ### set mc colors
                    bkgs[key]['hist'].SetMarkerColor(cis[color])
                    bkgs[key]['hist'].SetLineColor(cis[color])
                    color += 1
                    
                    if dru1 or dru2:
                        druScale = data['x'+str(i+1)+'-data'+'-e'+str(E)]['druScale']
                    if dru2:
                        bkgs[key]['hist'].Scale(druScale)
                    
                    if E and hiEplotRebin:
                        bkgs[key]['hist'].Rebin(hiEplotRebin)
                        bkgs[key]['hist'].Scale(1./float(hiEplotRebin))
                        #bkgs[key]['hist'].Sumw2()
                    
                    ### temp - don't draw MC if you just want to show data
                    bkgs[key]['hist'].Draw('same')

                    ### add MC to total MC hist
                    total[E][i].Add(bkgs[key]['hist'])
                    
                    ### create the legend entry for MC
                    legs[E][i].AddEntry(bkgs[key]['hist'], space+key, lopt)

            for key in sigkeys:
                if 'x'+str(i+1) in key and '-e'+str(E) in key:
                    ### set mc colors
                    sigs[key]['hist'].SetMarkerColor(cis[color])
                    sigs[key]['hist'].SetLineColor(cis[color])
                    color += 1
                    
                    if dru1 or dru2:
                        druScale = data['x'+str(i+1)+'-data'+'-e'+str(E)]['druScale']
                    if dru2:
                        sigs[key]['hist'].Scale(druScale)
                    
                    if E and hiEplotRebin:
                        sigs[key]['hist'].Rebin(hiEplotRebin)
                        sigs[key]['hist'].Scale(1./float(hiEplotRebin))
                        #sigs[key]['hist'].Sumw2()
                    
                    ### set range
                    #sigs[key]['hist'].SetAxisRange(1,1000,'y')

                    ### draw sigs
                    sigs[key]['hist'].Draw('same')

                    ### add MC to total MC hist
                    total[E][i].Add(sigs[key]['hist'])
                    
                    ### create the legend entry for MC
                    legs[E][i].AddEntry(sigs[key]['hist'], space+key, lopt)

            
            total[E][i].Sumw2()

            
            ### you need to scale the error by the dru scaling and/or the rebinning
            #-----------------------------------------------------------------------------
            if not dru1 and not dru2:
                for n in range(total[E][i].GetNbinsX()):
                    total[E][i].SetBinError(n, total[E][i].GetBinError(n)/(float(hiEplotRebin)/math.sqrt(2.)))
            if dru1:
                for n in range(total[E][i].GetNbinsX()):
                    total[E][i].SetBinError(n, total[E][i].GetBinError(n)*druScale)
                    
            if dru2:
                for n in range(total[E][i].GetNbinsX()):
                    total[E][i].SetBinError(n, total[E][i].GetBinError(n)*druScale)
                    
            #-----------------------------------------------------------------------------
                    
                    
            total[E][i].Draw('same')
            legs[E][i].AddEntry(total[E][i], space+'Total MC', lopt)
            legs[E][i].Draw('same')


            ### try to get the residuals in!
            #---------------------------------------------------------
            botpad[i].cd()
            leg = TLegend(0.72, 0.78, 0.94, 0.94)
            legs2[E].append(leg)
            legs2[E][i].SetFillColor(0)
            legs2[E][i].SetBorderSize(0)
            lopt = 'LPE'

            resid[E][i].Divide(data['x'+str(i+1)+'-data'+'-e'+str(E)]['hist'], total[E][i])

            resid[E][i].SetTitle('')
            resid[E][i].SetXTitle('Energy (keVee)')
            resid[E][i].GetXaxis().SetTitleFont(font)
            resid[E][i].GetXaxis().SetTitleSize(size)
            resid[E][i].GetXaxis().SetLabelFont(font)
            resid[E][i].GetXaxis().SetLabelSize(size)
            resid[E][i].GetXaxis().SetLabelOffset(0.03)
            resid[E][i].GetXaxis().SetTitleOffset(8)
            #resid[E][i].SetYTitle('counts / keV')
            resid[E][i].SetYTitle('data / MC')
            resid[E][i].GetYaxis().SetTitleFont(font)
            resid[E][i].GetYaxis().SetTitleSize(size)
            resid[E][i].GetYaxis().SetTitleOffset(yoff)
            resid[E][i].GetYaxis().SetLabelFont(font)
            resid[E][i].GetYaxis().SetLabelSize(size)
            resid[E][i].GetYaxis().SetLabelOffset(0.01)
            resid[E][i].GetYaxis().SetNdivisions(505) # '5' secondary and '05' primary

            resid[E][i].SetAxisRange(0.1,10,'y')
            resid[E][i].Draw()
            
            par = histparam(E)
            # set my line to '1'
            zero = TLine(par[1], 1, par[2], 1)
            zeros[E].append(zero)
            zeros[E][i].SetLineColor(kRed)
            zeros[E][i].SetLineWidth(1)
            zeros[E][i].Draw()

            #legs2[E][i].AddEntry(resid[E][i],space+'data / MC',lopt)
            #legs2[E][i].Draw()
            #---------------------------------------------------------

        save = ''
        if local:        save += 'local'
        else:            save += 'cup'
        save += '_'+str(runNum)
        if E:            save += '_hiE'
        else:            save += '_loE'
        save += '_loEfit-'+str(int(fLo[0]))+'-'+str(int(fLo[1]))
        save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        save += '_hiEfitRebin-'+str(hiEfitRebin)
        if fitRebinScale:   save += '_rebinscale'
        if mcscale:      save += '_mcscale'
        if mcweight:     save += '_mcweight'
        if fitweight:    save += '_fitweight'
        if dru1:         save += '_dru1'
        if dru2:         save += '_dru2'
        if hiEplotRebin: save += '_hiEplotRebin-'+str(hiEplotRebin)
        if reuse:        save += '_reuse'
        save += '_'+V
        
        canvs[E].Update()
        canvs[E].Print(save+'.png')

    if not batch:
        raw_input('[Enter] to quit \n')


######################################################################
######################################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])

