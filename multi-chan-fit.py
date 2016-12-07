#!/usr/bin/env python

######################################################################
# Matt Kauer - mkauer@physics.wisc.edu
######################################################################
# 14-set-bkg-activity.py

V = 'v14'

# Scale some MC to a fixed activity
#
# version: 2016-12-06
# 
# see CHANGELOG for changes
######################################################################
#
# Where is the MC I'm using?
# /data/MC/KIMS-NaI/user-scratch/sim/processed/K40/set2
# /data/MC/KIMS-NaI/user-scratch/sim/processed/U238/set2
# /data/MC/KIMS-NaI/user-scratch/sim/processed/Th232/set2
#
# Where is raw data I'm using?
# /data/KIMS/COSINE/PHY_RUN
#
# My processed data is currently in
# /home/mkauer/temp
#
######################################################################

import os,sys
import socket
import numpy as np
from ROOT import *
import ROOT
import funcs
from funcs import *
import copy
import math


### user inputs
#================================================================
### I'm pretty sure you want to convert to dru after the fit
### More stats for the fit will give bettet errors
### ie use dru2 option!!
dru1 = 0   ### convert data to dru before fit? [0,1]
dru2 = 1   ### convert data and mc to dru after fit? [0,1]

# dru safty check...
if dru2: dru1 = 0

hiEplotRebin = 10  ### rebin the hi-E final plots [1,inf]

# something weird with the fit rebinning - need to FIX this ASAP!
# use rebin=1 for now
hiEfitRebin = 1   ### rebin the hi-E histo for fitting [1,10]

# fill the fitdata and fitmc starting with n plus number?
# i think this should be 1 because bin 0 is underflows right?
# yes, just checked, this should be 1
np = 1

reuse = 1   ### use joined rootfile data? [0,1]

# this 'energy' variable doesn't matter anymore...
energy = 0   ### [0] = low-energy -- [1] = high-energy

mcscale = 1   ### pre scale the MC? [0,1]
mcweight = 1   ### set MC weights? [0,1]

# this doesn't work at all, not sure why...
fitweight = 0   ### set fit weights? [0,1]

# still working on this - I don't understand why this is
# effecting the fit so much - rebinscale=0 seems to be
# the right thing to do at this time...
rebinscale = 0   ### scale data and mc hi-E hists to 1/hiEfitRebin factor? [0,1]
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
    gStyle.SetPalette(1)
    gStyle.SetOptStat("")
    gStyle.SetOptFit(0)
    
    gStyle.SetPadBottomMargin(0.12)
    gStyle.SetPadLeftMargin  (0.12)
    gStyle.SetPadRightMargin (0.05)
    
    
    if reuse:
        rootfile = "./root-join-read/joined-master.root"
        #if local:
        #    rootfile = "/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/joined-master.root"
            #rootfile = "/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/join-test.root"
        #else:
        #    rootfile = "/home/mkauer/mc-fitting/root-join-read/joined-master.root"
        dataLo, mcLo, locs, isos = readROOT(rootfile, 0)
        dataHi, mcHi, locs, isos = readROOT(rootfile, 1)
        
    else:
        if local: path = "/home/mkauer/COSINE/CUP/mc-fitting/data/phys/"
        else:     path = "/home/mkauer/temp/"
        dataLo = getData([path+"*1324*root*"], 0)
        dataHi = getData([path+"*1324*root*"], 1)
        
        ### define what MC you want and from where
        locs = ['internal','pmt']
        isos = ['K40','U238','Th232']
        mcLo = buildMC(locs, isos, 0)
        mcHi = buildMC(locs, isos, 1)

    if dru1:
        dataLo, scalesLo = dataDRU(dataLo)
        dataHi, scalesHi = dataDRU(dataHi)

        ### this should be the same factor...
        ### checked, it is the same factor!!
        #print scalesLo
        #print scalesHi
    
    
    ### set the order so the color of the histos is consistent
    locs = ['internal','pmt']
    isos = ['K40','Th232','U238']
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    scaletest=0.005
    for i in range(8):
        mcLo['internal']['K40']['act'][i] = scaletest
        mcLo['internal']['K40']['hist'][i].Scale(scaletest)
        mcHi['internal']['K40']['act'][i] = scaletest
        mcHi['internal']['K40']['hist'][i].Scale(scaletest)
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # create a color scheme for MC
    # color list and indexes
    Nmc = int((len(locs)*len(isos)))
    colors, cis = rainbow(Nmc)
    
    # legend length = MC + data + total
    Nlg = Nmc+2
    
    ### fit the MC to data
    #-----------------------------------------------------------------
    
    fLo = [10, 100]
    fLoBins = fLo[1]-fLo[0]
    
    ### rebinning the hi-E part to have close to equal bins as the lo-E
    ### part will help to weight both parts more equally - I think...
    ### hiEfitRebin variable is now at the top
    
    fHiE = [400, 2000]
    fHi = [fHiE[0]/hiEfitRebin, fHiE[1]/hiEfitRebin]
    fHiBins = fHi[1]-fHi[0]
    
    print '!!!',fLo,fHi
    
    fbins = fLoBins+fHiBins
    
    fmin=0
    fmax=fbins

    print '!!!', fLoBins, fHiBins, fbins
    
    mcObj = []
    fit = []
    fitresults = []
    
    #-----------------------------------------------------------------
    ### create the multi-chan histos for fitting

    fitdata = []
    fitmc = {}
    
    ### init new fitmc dict
    for loc in locs:
        fitmc[loc] = {}
        for iso in isos:
            fitmc[loc][iso] = {}
            fitmc[loc][iso]["hist"] = []
            fitmc[loc][iso]["weight"] = []

            
    ### need seperate dicts for rebinned data and MC for plotting to
    ### to work right - they need their own memory space
    #-----------------------------------------------------------------
    ### init dict for rebinned data
    rdata = []

    # init dict for rebinned MC 
    rmc = {}
    for loc in locs:
        rmc[loc] = {}
        for iso in isos:
            rmc[loc][iso] = {}
            rmc[loc][iso]["hist"] = []
            #rmc[loc][iso]["weight"] = []
    #-----------------------------------------------------------------
    
    fmcsum=[] # sum of the fixed activity MC
    ftotal=[] # total of fixed MC plus fit MC
    fresid=[] # residual of the data/total
    ### fill fitdata and fitmc
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
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        mcsum = TH1F('fmcsum'+str(i), longNames(i), fbins,0,fbins)
        mcsum.SetLineColor(kBlack)
        mcsum.SetMarkerColor(kBlack)
        mcsum.SetLineWidth(1)
        fmcsum.append(mcsum)
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        
        # fill the low energy part
        #------------------------------------------------
        for n in range(fLoBins):
        #for n in range(1, fLoBins+1):
            fitdata[i].SetBinContent(n+np, dataLo[i].GetBinContent(fLo[0]+n))
            fitdata[i].SetBinError(n+np, dataLo[i].GetBinError(fLo[0]+n))
            #fitdata.append(fdata)
        for loc in locs:
            for iso in isos:
                
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                if mcLo[loc][iso]['act'][i] > 0:
                    for n in range(fLoBins):
                    #for n in range(1, fLoBins+1):
                        #fmcsum[i].SetBinContent(n+np, mcLo[loc][iso]['hist'][i].GetBinContent(fLo[0]+n)*mcLo[loc][iso]['act'][i])
                        fmcsum[i].SetBinContent(n+np, mcLo[loc][iso]['hist'][i].GetBinContent(fLo[0]+n))
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                else:
                    key = loc+'_'+iso+'_'+str(i)
                    #key = 'fmc'+loc+'_'+iso+'_'+str(i)
                    fmc = TH1F(key, key, fbins, 0, fbins)
                    #fmc = TH1F(key, longNames(i), fbins, 0, fbins)
                    fitmc[loc][iso]["hist"].append(fmc)
                    for n in range(fLoBins):
                    #for n in range(1, fLoBins+1):
                        fitmc[loc][iso]["hist"][i].SetBinContent(n+np, mcLo[loc][iso]['hist'][i].GetBinContent(fLo[0]+n))
                    #print 'filling lo-E with',key
        

        # fill the high energy part
        #-------------------------------------------------------------
        
        rdata.append(copy.deepcopy(dataHi[i]))
        rdata[i].Rebin(hiEfitRebin)
        
        ### I think we need to scale the counts by hiEfitRebin factor
        ### not sure that's true anymore...
        if rebinscale: rdata[i].Scale(1./hiEfitRebin)
        for n in range(fLoBins,fbins):
        #for n in range(fLoBins+1, fbins+1):
            fitdata[i].SetBinContent(n+np, rdata[i].GetBinContent(fHi[0]+n))
            fitdata[i].SetBinError(n+np, rdata[i].GetBinError(fHi[0]+n))
        for loc in locs:
            for iso in isos:
                rmc[loc][iso]['hist'].append(copy.deepcopy(mcHi[loc][iso]['hist'][i]))
                rmc[loc][iso]['hist'][i].Rebin(hiEfitRebin)

                ### I think we need to scale the counts by hiEfitRebin factor
                ### not sure that's true anymore...
                if rebinscale: rmc[loc][iso]['hist'][i].Scale(1./hiEfitRebin)

                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                if mcHi[loc][iso]['act'][i] > 0:
                    #for n in range(fLoBins): # was this
                    for n in range(fLoBins,fbins): # fixed 2016-12-05
                        # doesn't seem to make much difference in the fit...
                    #for n in range(fLoBins+1, fbins+1):
                        #fmcsum[i].SetBinContent(n+np, rmc[loc][iso]['hist'][i].GetBinContent(fHi[0]+n)*mcHi[loc][iso]['act'][i])
                        fmcsum[i].SetBinContent(n+np, rmc[loc][iso]['hist'][i].GetBinContent(fHi[0]+n))
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                else:
                    for n in range(fLoBins,fbins):
                    #for n in range(fLoBins+1, fbins+1):
                        fitmc[loc][iso]["hist"][i].SetBinContent(n+np, rmc[loc][iso]['hist'][i].GetBinContent(fHi[0]+n))
                    #print 'filling hi-E with',key

        ### subtract fixed MC from data
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ### this works I think, but then this plots the data-fixedMC
        fitdata[i].Add(fmcsum[i], -1)
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        
    ### make weights histo
    ### testing with TFractionFitter
    ### this doesn't seem to work yet...
    weights = TH1F('weights','weights',fbins,0,fbins)
    for n in range(0,fbins):
        weights.SetBinContent(n, 1)
    weights.Sumw2()
    #-----------------------------------------------------------------
    
    for i in range(8):
        
        mcObj.append(TObjArray(Nmc)) # number of MC to fit to
        
        fitdata[i].Sumw2()
        dat_int = fitdata[i].Integral(fmin,fmax) # data integral to normalize to

        for loc in locs:
            for iso in isos:
                
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                if mcHi[loc][iso]['act'][i] > 0: continue
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                mc_int = fitmc[loc][iso]["hist"][i].Integral(fmin,fmax) # MC integral
                #print iso,loc,'mc_int =',mc_int
                
                
                ### to weight or not to weight...
                # what if you weight first?
                if mcweight:
                    fitmc[loc][iso]["hist"][i].Sumw2() # set stat weights
                                    
                ### to scale or not to scale...
                if mcscale:
                    fitmc[loc][iso]["hist"][i].Scale(dat_int/mc_int) # scale to data integral
                    mcLo[loc][iso]["hist"][i].Scale(dat_int/mc_int) # make sure MC is scaled too
                    mcHi[loc][iso]["hist"][i].Scale(dat_int/mc_int) # make sure MC is scaled too
                    
                                
                mcObj[i].Add(fitmc[loc][iso]["hist"][i]) # add to the TFractionFitter object
                
        fit.append(TFractionFitter(fitdata[i], mcObj[i])) # create the TFF data and MC objects

        ### for some reason on CUP the fit status gets returned as a pointer
        ### like <ROOT.TFitResultPtr object at 0x37246c0>
        #fitresults.append("NaI-C"+str(i+1)+" fit completed with status "+str(status))
        fitresults.append("NaI-C"+str(i+1)+" fit results")
        fitresults.append("MC fit constrained to 0.0001 - 10.0")
        
        for l in range(Nmc):
            # set fit constraints on the MC put into mcObj TObject
            fit[i].Constrain(l, 0.0001, 10.0)

            # set fit weights to 1
            # this doesn't work yet, not sure why, on the to-do list...
            if fitweight:
                fit[i].SetWeight(l, weights)
                #fit[i].SetWeight(l, 1)
        
        fit[i].SetRangeX(fmin, fmax) # set global range, should be the same as fmin and fmax?
        
        status = fit[i].Fit() # do the fit!
        
        
        #-------------------------------------------------------------
        ### hum, didn't know I could do this...
        ### keep it in my pocket, probably use this later...
        #-------------------------------------------------------------
        #result = TH1F("result", "result", par[0], par[1], par[2])
        #result = fit[i].GetPlot()
        #test = TCanvas('test', 'test', 0, 0, 800, 600)
        #result.Draw()
        #data[i].Draw('same')
        #test.Update()
        #-------------------------------------------------------------
        
        
        chi2 = fit[i].GetChisquare()
        ndf  = fit[i].GetNDF()
        pval = fit[i].GetProb()
        
        ### for some reason on CUP the fit status gets returned as a pointer
        ### like <ROOT.TFitResultPtr object at 0x37246c0>
        #fitresults.append("NaI-C"+str(i+1)+" fit completed with status "+str(status))
        #fitresults.append("NaI-C"+str(i+1)+" fit completed")
        fitresults.append("chi2/ndf = "+str(round(chi2,2))+"/"+str(ndf)+" = "+str(round(chi2/float(ndf),2)))

        ### p-val is always 0 - need to look into this
        #fitresults.append('p-value = '+str(pval))
        
        count = 0
        for loc in locs:
            for iso in isos:
                
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                if mcHi[loc][iso]['act'][i] > 0: continue
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                value = ROOT.Double(0.0)
                error = ROOT.Double(0.0)
                #l = j+k
                #print l
                fit[i].GetResult(count, value, error)
                #print iso,loc,'=',value,' +/- ',error
                fitresults.append(str(loc)+' '+str(iso)+' = '+str(round(value,4))+' +/- '+str(round(error,4)))
                fitmc[loc][iso]["hist"][i].Scale(value)
                mcLo[loc][iso]["hist"][i].Scale(value)
                mcHi[loc][iso]["hist"][i].Scale(value)
                ### need to scale the hiE back by the hiEfitRebinning factor?
                #mcHi[loc][iso]["hist"][i].Scale(1./hiEfitRebin)
                #mcHi[loc][iso]["hist"][i].Scale(hiEfitRebin)
                
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
    save += '_Nchan-fit'
    save += '_loEfit-'+str(int(fLo[0]))+'-'+str(int(fLo[1]))
    save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
    save += '_hiEfitRebin-'+str(hiEfitRebin)
    if rebinscale:   save += '_rebinscale'
    if mcscale:      save += '_mcscale'
    if mcweight:     save += '_mcweight'
    if fitweight:    save += '_fitweight'
    if dru1:         save += '_dru1'
    if dru2:         save += '_dru2'
    if hiEplotRebin: save += '_hiEplotRebin-'+str(hiEplotRebin)
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
        pad1 = TPad("pad1","pad1",0,fraction,1,1)
        ftoppad.append(pad1)
        pad2 = TPad("pad2","pad2",0,0,1,fraction)
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
        ylegend = (ylegstart-(Nlg*0.06))
        space = '  '
        leg = TLegend(0.65, ylegend, 0.94, ylegstart)
        flegs.append(leg)
        flegs[i].SetFillColor(0)
        flegs[i].SetBorderSize(0)
        lopt = 'LPE'

        """
        if local and not reuse:
            fitdata[i].SetAxisRange(1,1000,"y")
        else:
            fitdata[i].SetAxisRange(1,1000000,"y")
        """
        
        fitdata[i].SetLineColor(kBlack)
        fitdata[i].SetMarkerColor(kBlack)
        fitdata[i].SetLineWidth(1)
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
        for loc in locs:
            for iso in isos:
                
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                if mcHi[loc][iso]['act'][i] > 0: continue
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                #+++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                ### set mc colors
                fitmc[loc][iso]["hist"][i].SetMarkerColor(cis[color])
                fitmc[loc][iso]["hist"][i].SetLineColor(cis[color])
                color += 1

                #fitmc[loc][iso]["hist"][i].SetAxisRange(1,1000,"y")

                ### temp - don't draw MC if you just want to show data
                fitmc[loc][iso]["hist"][i].Draw("same")

                ### add MC to total MC hist
                ftotal[i].Add(fitmc[loc][iso]["hist"][i])

                ### create the legend entry for MC
                flegs[i].AddEntry(fitmc[loc][iso]["hist"][i], space+loc+'-'+iso, lopt)

        ftotal[i].Draw("same")
        flegs[i].AddEntry(ftotal[i], space+'Total MC', lopt)
        flegs[i].Draw("same")

        
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
        #fresid[i].SetXTitle("Energy (keVee)")
        fresid[i].SetXTitle("fit bins")
        fresid[i].GetXaxis().SetTitleFont(font)
        fresid[i].GetXaxis().SetTitleSize(size)
        fresid[i].GetXaxis().SetTitleOffset(8)
        fresid[i].GetXaxis().SetLabelFont(font)
        fresid[i].GetXaxis().SetLabelSize(size)
        fresid[i].GetXaxis().SetLabelOffset(0.03)
        #fresid[i].SetYTitle("counts / keV")
        fresid[i].SetYTitle("data / MC")
        fresid[i].GetYaxis().SetTitleFont(font)
        fresid[i].GetYaxis().SetTitleSize(size)
        fresid[i].GetYaxis().SetTitleOffset(yoff)
        fresid[i].GetYaxis().SetLabelFont(font)
        fresid[i].GetYaxis().SetLabelSize(size)
        fresid[i].GetYaxis().SetLabelOffset(0.01)
        fresid[i].GetYaxis().SetNdivisions(505) # '5' secondary and '05' primary

        fresid[i].SetAxisRange(0.1,10,"y")
        #if local and not reuse:
        #    fresid[i].SetAxisRange(-100,100,"y")
        #else:
        #    fresid[i].SetAxisRange(-400,400,"y")

        fresid[i].Draw()

        #par = histparam(E)
        ### set my line to "1"
        zero = TLine(fmin, 1, fmax, 1)
        fzeros.append(zero)
        fzeros[i].SetLineColor(kRed)
        fzeros[i].SetLineWidth(1)
        fzeros[i].Draw()

        flegs2[i].AddEntry(fresid[i],space+"data / MC",lopt)
        #flegs2[i].Draw()
        #-------------------------------------------------------------
    """
    save = ''
    if local:      save += 'local'
    else:          save += 'cup'
    save += '_Nchan-fit'
    save += '_loEfit-'+str(int(fLo[0]))+'-'+str(int(fLo[1]))
    save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
    save += '_hiEfitRebin-'+str(hiEfitRebin)
    if rebinscale: save += '_rebinscale'
    if mcscale:    save += '_mcscale'
    if mcweight:   save += '_mcweight'
    if fitweight:  save += '_fitweight'
    if dru1:       save += '_dru1'
    if dru2:       save += '_dru2'
    if hiEplotRebin:   save += '_hiEplotRebin-'+str(hiEplotRebin)
    """
    fcanv.Update()
    fcanv.Print(save+'.png')

    
    # plot the lo and hi energy histograms
    #=================================================================
    
    canvs  = [None for x in range(2)]
    data   = [None for x in range(2)]
    mc     = [None for x in range(2)]
    scales = [None for x in range(2)]
    
    total  = [None for x in range(2)]
    resid  = [None for x in range(2)]
    
    legs   = [None for x in range(2)]
    legs2  = [None for x in range(2)]
    zeros  = [None for x in range(2)]
    
    for E in range(2):
        
        legs[E]=[]
        legs2[E]=[]
        zeros[E]=[]
        
        if E:
            #data[E] = dataHi
            data[E] = copy.deepcopy(dataHi)
            #mc.append(mcHi)
            #mc[E] = mcHi
            mc[E] = copy.deepcopy(mcHi)
            
        else:
            #data[E] = dataLo
            data[E] = copy.deepcopy(dataLo)
            #mc.append(mcLo)
            #mc[E] = mcLo
            mc[E] = copy.deepcopy(mcLo)

        if dru2:
            data[E], scales[E] = dataDRU(data[E])

        #total = makeTotal(E) # make a set of total MC histos
        #resid = makeResid(E) # make a set of resid histos
        total[E] = makeTotal(E)
        resid[E] = makeResid(E)
                    
        # have the plotting be seperated out from the loop
        canvs[E] = TCanvas('canv'+str(E), 'canv'+str(E), 0, 0, 1400, 900)
        #canvs.append(canv)
        canvs[E].Divide(4,2)
        
        toppad=[]
        botpad=[]
        #legs=[]
        #legs2=[]
        #zeros=[]

        font=63
        size=13

        for i in range(8):
            canvs[E].cd(i+1)
            #canvs[E].cd(i+1).SetLogy()

            fraction = 0.3
            pad1 = TPad("pad1","pad1",0,fraction,1,1)
            toppad.append(pad1)
            pad2 = TPad("pad2","pad2",0,0,1,fraction)
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
            ylegend = (ylegstart-(Nlg*0.06))
            space = '  '
            leg = TLegend(0.65, ylegend, 0.94, ylegstart)
            legs[E].append(leg)
            legs[E][i].SetFillColor(0)
            legs[E][i].SetBorderSize(0)
            lopt = 'LPE'
            
            #if dru1 or dru2:
            #    data[E][i].SetAxisRange(.1, 100, "y")
            
            """
            if local and not reuse:
                data[E][i].SetAxisRange(1,1000,"y")
            else:
                data[E][i].SetAxisRange(1,1000000,"y")
            """
            
            #data[E][i].SetLineColor(kBlack)
            #data[E][i].SetMarkerColor(kBlack)
            #data[E][i].SetLineWidth(1)
            data[E][i].GetYaxis().SetTitle('arb. counts')
            if dru1 or dru2:
                data[E][i].GetYaxis().SetTitle('counts / day / kg / keV  (dru)')
                            
            #---------------------------------------------------------
            if E and hiEplotRebin:
                data[E][i].Rebin(hiEplotRebin)
                data[E][i].Scale(1./float(hiEplotRebin))
                #data[E][i].Sumw2()
                #data[E][i].SetAxisRange(.1, 10, "y")
                
                total[E][i].Rebin(hiEplotRebin)
                total[E][i].Scale(1./float(hiEplotRebin))
                #total[E][i].Sumw2()
                
                resid[E][i].Rebin(hiEplotRebin)
                resid[E][i].Scale(1./float(hiEplotRebin))
                #resid[E][i].Sumw2()
            #---------------------------------------------------------
            
            data[E][i].GetYaxis().SetTitleFont(font)
            data[E][i].GetYaxis().SetTitleSize(size)
            data[E][i].GetYaxis().SetTitleOffset(yoff)
            data[E][i].GetYaxis().SetLabelFont(font)
            data[E][i].GetYaxis().SetLabelSize(size)
            data[E][i].GetYaxis().SetLabelOffset(0.01)
            #data[E][i].GetXaxis().SetTitle('Energy (keV)')
            #data[E][i].GetXaxis().SetLabelFont(font)
            #data[E][i].GetXaxis().SetLabelSize(size)
            data[E][i].Draw()
            legs[E][i].AddEntry(data[E][i], space+'data', lopt)
            
            color = 0
            for loc in locs:
                for iso in isos:
                    ### set mc colors
                    mc[E][loc][iso]["hist"][i].SetMarkerColor(cis[color])
                    mc[E][loc][iso]["hist"][i].SetLineColor(cis[color])
                    color += 1

                    # scale to dru
                    # dru1 is already scaled... (don't scale it again)
                    #if dru1:
                    #    mc[E][loc][iso]["hist"][i].Scale(scalesHi[i])
                    if dru2:
                        mc[E][loc][iso]["hist"][i].Scale(scales[E][i])
                    
                    if E and hiEplotRebin:
                        mc[E][loc][iso]["hist"][i].Rebin(hiEplotRebin)
                        mc[E][loc][iso]["hist"][i].Scale(1./float(hiEplotRebin))
                        #mc[E][loc][iso]["hist"][i].Sumw2()
                    
                    ### set range
                    #mc[E][loc][iso]["hist"][i].SetAxisRange(1,1000,"y")

                    ### temp - don't draw MC if you just want to show data
                    mc[E][loc][iso]["hist"][i].Draw("same")

                    ### add MC to total MC hist
                    total[E][i].Add(mc[E][loc][iso]["hist"][i])
                    
                    ### create the legend entry for MC
                    legs[E][i].AddEntry(mc[E][loc][iso]["hist"][i], space+loc+'-'+iso, lopt)

            
            total[E][i].Sumw2()

            
            ### you need to scale the error by the dru scaling and/or the rebinning
            #-----------------------------------------------------------------------------
            if not dru1 and not dru2:
                for n in range(total[E][i].GetNbinsX()):
                    #total[E][i].SetBinError(n, total[E][i].GetBinError(n)/float(hiEplotRebin))
                    #total[E][i].SetBinError(n, total[E][i].GetBinError(n)/np.sqrt(float(hiEplotRebin)))
                    #total[E][i].SetBinError(n, total[E][i].GetBinError(n)/math.sqrt(float(hiEplotRebin)))
                    total[E][i].SetBinError(n, total[E][i].GetBinError(n)/(float(hiEplotRebin)/math.sqrt(2.)))
            if dru1:
                for n in range(total[E][i].GetNbinsX()):
                    total[E][i].SetBinError(n, total[E][i].GetBinError(n)*scalesHi[i])
            if dru2:
                for n in range(total[E][i].GetNbinsX()):
                    total[E][i].SetBinError(n, total[E][i].GetBinError(n)*scales[E][i])
            #-----------------------------------------------------------------------------
                    
                    
            total[E][i].Draw("same")
            legs[E][i].AddEntry(total[E][i], space+'Total MC', lopt)
            legs[E][i].Draw("same")


            ### try to get the residuals in!
            #---------------------------------------------------------
            botpad[i].cd()
            leg = TLegend(0.72, 0.78, 0.94, 0.94)
            legs2[E].append(leg)
            legs2[E][i].SetFillColor(0)
            legs2[E][i].SetBorderSize(0)
            lopt = 'LPE'

            #resid[E][i].Sumw2()
            #resid[E][i].Add(data[E][i], total[E][i], 1, -1)
            resid[E][i].Divide(data[E][i], total[E][i])

            #---------------------------------------------------------
            #if E and hiEplotRebin:
                #resid[E][i].Rebin(hiEplotRebin)
                #resid[E][i].Scale(1./hiEplotRebin)
                #resid[E][i].Sumw2()
            #---------------------------------------------------------
            
            resid[E][i].SetTitle('')
            resid[E][i].SetXTitle("Energy (keVee)")
            resid[E][i].GetXaxis().SetTitleFont(font)
            resid[E][i].GetXaxis().SetTitleSize(size)
            resid[E][i].GetXaxis().SetLabelFont(font)
            resid[E][i].GetXaxis().SetLabelSize(size)
            resid[E][i].GetXaxis().SetLabelOffset(0.03)
            resid[E][i].GetXaxis().SetTitleOffset(8)
            #resid[E][i].SetYTitle("counts / keV")
            resid[E][i].SetYTitle("data / MC")
            resid[E][i].GetYaxis().SetTitleFont(font)
            resid[E][i].GetYaxis().SetTitleSize(size)
            resid[E][i].GetYaxis().SetTitleOffset(yoff)
            resid[E][i].GetYaxis().SetLabelFont(font)
            resid[E][i].GetYaxis().SetLabelSize(size)
            resid[E][i].GetYaxis().SetLabelOffset(0.01)
            resid[E][i].GetYaxis().SetNdivisions(505) # '5' secondary and '05' primary

            resid[E][i].SetAxisRange(0.1,10,"y")
            #if local and not reuse:
            #    resid[E][i].SetAxisRange(-100,100,"y")
            #else:
            #    resid[E][i].SetAxisRange(-400,400,"y")
            
            resid[E][i].Draw()
            
            par = histparam(E)
            # set my line to "1"
            zero = TLine(par[1], 1, par[2], 1)
            zeros[E].append(zero)
            zeros[E][i].SetLineColor(kRed)
            zeros[E][i].SetLineWidth(1)
            zeros[E][i].Draw()

            #legs2[E][i].AddEntry(resid[E][i],space+"data / MC",lopt)
            #legs2[E][i].Draw()
            #---------------------------------------------------------

        save = ''
        if local:        save += 'local'
        else:            save += 'cup'
        if E:            save += '_hiE'
        else:            save += '_loE'
        save += '_loEfit-'+str(int(fLo[0]))+'-'+str(int(fLo[1]))
        save += '_hiEfit-'+str(int(fHiE[0]))+'-'+str(int(fHiE[1]))
        save += '_hiEfitRebin-'+str(hiEfitRebin)
        if rebinscale:   save += '_rebinscale'
        if mcscale:      save += '_mcscale'
        if mcweight:     save += '_mcweight'
        if fitweight:    save += '_fitweight'
        if dru1:         save += '_dru1'
        if dru2:         save += '_dru2'
        if hiEplotRebin: save += '_hiEplotRebin-'+str(hiEplotRebin)
        save += '_'+V
        
        canvs[E].Update()
        canvs[E].Print(save+'.png')

    if not batch:
        raw_input("[Enter] to quit \n")


######################################################################
######################################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])

