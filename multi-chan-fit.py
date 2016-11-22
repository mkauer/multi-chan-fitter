#!/usr/bin/env python

######################################################################
# 12-scale-hi-energy.py
#
# Try rebinning the the hi-energy part to about the same number of
# bins as the low energy part. This should help with the fit
# weighting issue. If the hi-E part has 2x more bins, the fit is
# going to weight that 2x more than the lo-E part. 
#
# version: 2016-11-15
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ I'm probably missing another factor here but it acutally seems
#   that by NOT scaling the rebinned data/mc by 1/rebin, the fit
#   does a better job - I don't understand this yet...
# + scale the data and mc by 1/rebin factor or the lo-E to hi-E scale
#   will bias the higher stats from the rebinned histos
# + I kinda fixed the rebinning bug. still needs more testing...
#   The fit still freaks out sometimes with some rebin factors
#   Not sure why but I need to investigate that more...
# + fixed another big bug with the combined fit histo gap
# + fixed a major bug with the filling of the fit-histos
# ~ let's make the residual show data/MC not data-MC
# + add more lists for total,resid,legs,zeros so they have their
#   own memory space and aren't lost when modifying the canvases
# ~ resolution function should be divided by 2 (see funcs.py)
#   I generated a new joined rootfile with the MC reso / 2
# + added TFractionFitter weight histo but it doesn't work? 
#   Still need to work on that. Don't understand why that doesn't work.
# + plot both loE and hiE histos that are scaled from the multi-fit
# + combine loE and hiE histos into 1 histo for fitting
# + need copy.copy or copy.deepcopy when reading histos from rootfile
#   and adding to a list or dict. Going to use deepcopy just to be sure.
# ~ now use buildMC() to get your MC
# ~ all histograms should have a unique key
# ~ fixed hostname selection
# + NOTE: The funcs.py import will only work for 10-mc-fit.py and later!!!
# + I guess importing other personal python modules now works
#   with the cup cluster? It didn't work in the past. Weird.
#   Well okay, now I have a seperate funcs.py to import and it seems
#   to be working with the cup cluster as well. 
# ~ big re-write - set mc as mc[loc][iso] - just seems more logical
#   just to be clear - it was mc[iso][loc] - now it's mc[loc][iso]
# ~ shit! fixed a bug with enumerate [iso][loc] j+k thing
#   that doesn't work!!
#   that also effected the mc scaling and plotting... 
# + get p-value from fit - doesn't seem right?
# ~ tweaked the plotting cosmetics a bit to make it look better,
#   good grief, what a bitch...
# + put residuals on the same canvas! awesome!
# ~ add options to scale and/or weight the MC
# ~ change plot names if generated on CUP or not
# ~ low energy fit needs to be >= 5keV? not 0keV!
# ~ that was a lot of if-else statements... good grief...
# + option for CUP testing - just 1 root file for testing
# + option for hi/lo energy fitting
# ~ tweak for low energy fit
# ~ fixed mc+data+tot+resid length issue
# + another canvas showing the residuals
# + show the total mc in gray after the fit
# + get hostname and select data path based off hostname
# + try to fit high energy stuff - based off my dm17 fitting code
# ~ tweak legend size to make a little more consistent,
#   but it's still not quite right, i don't know why..
# ~ plot low energy things - ah, use pmtXX.qc5 not qc_5 - oops
# + added low energy calib and resolution functions
# ~ use assumption that res~srqt(E) with intercept=0
# + working on resolution params, but need at least a linear
#   fit to the resolution points to get anything reasonable
# + working on hiE parameters
# + shit! finally got colors to work. What a pain in the ass.
#   The color index can't be too big > 10,000 and you have to
#   save the list of colors and indexes or they get dumped
#   out of memory and then it doesn't work right.
# + urgh... finally got hists into dictionary
#   hist name must be same as hist var-name (don't know why)
# + start adding in external backgrounds too
# ~ use Pushpa's data calibrations - wrong but a start
# ~ just use this as the default so I don't have multiple
#   versions of scripts around
# ~ for some weird reason importing my "funcs.py" doesn't work
#   on the cup cluster? Very weird. So I've put all the extra
#   functions in this script and now it works on CUP. ???
# ~ moved sim stuff to under CUP
# ~ tweaked spacings and labels on the canvas
# + got legends working for the plots
# ~ tweaked the layout to be 4x2 as real layout
# + seperate funcs.py import for extra functions
# + adding in U238 with K40
# + functionized some things in prep for the future fitting
# + ROOT batch settings [0 or 1]
# + canvas as a python list
# + histos as a python list
# + hacking at multi-plot-test.py
# + hacking at plot-test.py
#
# email me: mkauer@physics.wisc.edu
######################################################################
#
# Where is MC?
# /data/MC/KIMS-NaI/user-scratch/sim/processed/K40/set2
# /data/MC/KIMS-NaI/user-scratch/sim/processed/U238/set2
# /data/MC/KIMS-NaI/user-scratch/sim/processed/Th232/set2
#
# Where is raw data?
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


### user inputs
#================================================================
reuse      = 1   ### use joined rootfile data?
energy     = 0   ### 0 = low-energy  -- 1 = high-energy
mcscale    = 1   ### pre scale the MC?
mcweight   = 1   ### set MC weights?
fitweight  = 0   ### set fit weights?

### still working on this - I don't understand why this is
### effecting the fit so much - rebinscale=0 seems to be
### the right thing to do at this time...
rebinscale = 0   ### scale hi-E hist to 1/rebin factor?
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
        if local:
            rootfile = "/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/joined-master.root"
            #rootfile = "/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/join-test.root"
        else:
            rootfile = "/home/mkauer/mc-fitting/root-join-read/joined-master.root"
        dataLo, mcLo, locs, isos = readROOT(rootfile, 0)
        dataHi, mcHi, locs, isos = readROOT(rootfile, 1)
    else:
        if local: path = "/home/mkauer/COSINE/CUP/mc-fitting/data/phys/"
        else:     path = "/home/mkauer/temp/"
        data = getData([path+"*1324*root*"], energy)
        
        ### define what MC you want and from where
        locs = ['internal','pmt']
        isos = ['K40','U238','Th232']
        mc = buildMC(locs, isos, energy)

    
    ### just for testing and making sure plotting still works
    #if energy:
    #    data = dataHi
    #    mc = mcHi
    #else:
    #    data = dataLo
    #    mc = mcLo
    
    
    # maintain the same fitting order and what not?
    # will need to FIX this eventually...
    # seems the order you load the MC into the fit object effects things greatly
    locs = ['internal','pmt']
    isos = ['K40','U238','Th232']
    #isos = ['Th232','U238','K40']
    
    #locs = ['internal']
    #isos = ['K40','U238','Th232']
    #isos = ['Th232','U238','K40']
    
    
    # create a color scheme for MC
    # color list and indexes
    Nmc = int((len(locs)*len(isos)))
    colors, cis = rainbow(Nmc)
    
    # legend length = MC + data + total
    Nlg = Nmc+2
    
    ### fit the MC to data
    #=================================================================
    #=================================================================
    
    fLo = [10, 100]
    fLoBins = fLo[1]-fLo[0]
    
    ### rebinning the hi-E part to have close to equal bins as the lo-E
    ### part will help to weight both parts more equally - I think...
    rebin = 10
    
    fHi = [400/rebin, 2000/rebin]
    fHiBins = fHi[1]-fHi[0]

    print '!!!',fLo,fHi
    
    fbins = fLoBins+fHiBins
    
    fmin=0
    fmax=fbins

    print '!!!', fLoBins, fHiBins, fbins
    
    mcObj = []
    fit = []
    results = []
    
    ##############################################
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

    ftotal=[]
    fresid=[]
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
        
        ### low energy hists
        for n in range(fLoBins):
            fitdata[i].SetBinContent(n, dataLo[i].GetBinContent(fLo[0]+n))
            fitdata[i].SetBinError(n, dataLo[i].GetBinError(fLo[0]+n))
            #fitdata.append(fdata)
        for loc in locs:
            for iso in isos:
                key = loc+'_'+iso+'_'+str(i)
                fmc = TH1F(key,key,fbins,0,fbins)
                fitmc[loc][iso]["hist"].append(fmc)
                for n in range(fLoBins):
                    #key = loc+'_'+iso+'_'+str(i)
                    #fmc = TH1F(key,key,fbins,0,fbins)
                    #fmc.SetBinContent(n, mcLo[loc][iso]['hist'][i].GetBinContent(fLo[0]+n))
                    #fitmc[loc][iso]["hist"].append(fmc)
                    fitmc[loc][iso]["hist"][i].SetBinContent(n, mcLo[loc][iso]['hist'][i].GetBinContent(fLo[0]+n))
                #print 'filling lo-E with',key
                
        ### high energy hists
        dtemp = copy.deepcopy(dataHi[i])
        dtemp.Rebin(rebin)
        ### I'm pretty sure I need to scale the counts by rebin factor
        if rebinscale: dtemp.Scale(1./rebin)
        for n in range(fLoBins,fbins):
            #dtemp = copy.deepcopy(dataHi[i])
            #dtemp.Rebin(rebin)
            # I'm pretty sure I need to scale the counts by rebin factor
            #dtemp.Scale(1./rebin)
            fitdata[i].SetBinContent(n, dtemp.GetBinContent(fHi[0]+n))
            fitdata[i].SetBinError(n, dtemp.GetBinError(fHi[0]+n))
            #fitdata[i].SetBinContent(n, dataHi[i].GetBinContent(fHi[0]+n))
            #fitdata[i].SetBinError(n, dataHi[i].GetBinError(fHi[0]+n))
        for loc in locs:
            for iso in isos:
                mctemp = copy.deepcopy(mcHi[loc][iso]['hist'][i])
                mctemp.Rebin(rebin)
                ### I'm pretty sure I need to scale the counts by rebin factor
                if rebinscale: mctemp.Scale(1./rebin)
                for n in range(fLoBins,fbins):
                    #mctemp = copy.deepcopy(mcHi[loc][iso]['hist'][i])
                    #mctemp.Rebin(rebin)
                    # I'm pretty sure I need to scale the counts by rebin factor
                    #mctemp.Scale(1./rebin)
                    fitmc[loc][iso]["hist"][i].SetBinContent(n, mctemp.GetBinContent(fHi[0]+n))
                    #fitmc[loc][iso]["hist"][i].SetBinContent(n, mcHi[loc][iso]['hist'][i].GetBinContent(fHi[0]+n))
                #print 'filling hi-E with',key
                
    ### make weights histo
    ### testing with TFracFit
    weights = TH1F('weights','weights',fbins,0,fbins)
    for n in range(0,fbins):
        weights.SetBinContent(n, 1)
    weights.Sumw2()
    ##############################################
    
    for i in range(8):
        
        mcObj.append(TObjArray(Nmc)) # number of MC to fit to
        
        fitdata[i].Sumw2()
        dat_int = fitdata[i].Integral(fmin,fmax) # data integral to normalize to

        for loc in locs:
            for iso in isos:
                mc_int = fitmc[loc][iso]["hist"][i].Integral(fmin,fmax) # MC integral
                #print iso,loc,'mc_int =',mc_int
                
                
                #-----------------------------------------------------------------------
                #=======================================================================
                ########################################################################
                
                ### to scale or not to scale...
                if mcscale:
                    fitmc[loc][iso]["hist"][i].Scale(dat_int/mc_int) # scale to data integral
                    mcLo[loc][iso]["hist"][i].Scale(dat_int/mc_int) # make sure MC is scaled too
                    mcHi[loc][iso]["hist"][i].Scale(dat_int/mc_int) # make sure MC is scaled too
                    
                ### to weight or not to weight...
                # this doesn't work yet, not sure why, on the to-do list...
                if mcweight:
                    fitmc[loc][iso]["hist"][i].Sumw2() # set stat weights
                
                ########################################################################
                #=======================================================================
                #-----------------------------------------------------------------------
                
                
                mcObj[i].Add(fitmc[loc][iso]["hist"][i]) # add to the TFractionFitter object
                
        fit.append(TFractionFitter(fitdata[i], mcObj[i])) # create the TFF data and MC objects
        
        for l in range(Nmc):
            # set bounds on the MC put into mcObj TObject
            fit[i].Constrain(l, 0.0001, 10.0)
            # set fit weights to 1
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
        
        #print "Fit done with status ",status," -  chi2/ndf =",chi2,"/",ndf
        # for some reason on CUP the fit status gets returned as a pointer
        # like <ROOT.TFitResultPtr object at 0x37246c0>
        #results.append("NaI-C"+str(i+1)+" fit completed with status "+str(status))
        results.append("NaI-C"+str(i+1)+" fit completed")
        results.append("chi2/ndf = "+str(round(chi2,2))+"/"+str(ndf)+" = "+str(round(chi2/float(ndf),2)))
        results.append('p-value = '+str(pval))
        
        count = 0
        for loc in locs:
            for iso in isos:
                value = ROOT.Double(0.0)
                error = ROOT.Double(0.0)
                #l = j+k
                #print l
                fit[i].GetResult(count, value, error)
                #print iso,loc,'=',value,' +/- ',error
                results.append(str(loc)+' '+str(iso)+' = '+str(round(value,4))+' +/- '+str(round(error,4)))
                fitmc[loc][iso]["hist"][i].Scale(value)
                mcLo[loc][iso]["hist"][i].Scale(value)
                mcHi[loc][iso]["hist"][i].Scale(value)
                count += 1
        results.append('\n')

    print '\n\n'
    for line in results:
        print line

    #=================================================================
    # plot the multi-chan fit results
    #=================================================================
    
    flegs=[]
    flegs2=[]
    fzeros=[]

    fcanv = TCanvas('fcanv', 'fcanv', 0, 0, 1400, 900)
    fcanv.Divide(4,2)

    toppad=[]
    botpad=[]
    
    font=63
    size=13

    for i in range(8):
        fcanv.cd(i+1)
        
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
        flegs.append(leg)
        flegs[i].SetFillColor(0)
        flegs[i].SetBorderSize(0)
        lopt = 'LPE'
        
        if local and not reuse:
            fitdata[i].SetAxisRange(1,1000,"y")
        else:
            fitdata[i].SetAxisRange(1,1000000,"y")

        fitdata[i].SetLineColor(kBlack)
        fitdata[i].SetMarkerColor(kBlack)
        fitdata[i].SetLineWidth(1)
        fitdata[i].GetYaxis().SetTitle('arb. counts (for now)')
        fitdata[i].GetYaxis().SetTitleFont(font)
        fitdata[i].GetYaxis().SetTitleSize(size)
        fitdata[i].GetYaxis().SetTitleOffset(4.2)
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
                ### set mc colors
                fitmc[loc][iso]["hist"][i].SetMarkerColor(cis[color])
                fitmc[loc][iso]["hist"][i].SetLineColor(cis[color])
                color += 1

                fitmc[loc][iso]["hist"][i].SetAxisRange(1,1000,"y")

                ### temp - don't draw MC if you just want to show data
                fitmc[loc][iso]["hist"][i].Draw("same")

                ### add MC to total MC hist
                ftotal[i].Add(fitmc[loc][iso]["hist"][i])

                ### create the legend entry for MC
                flegs[i].AddEntry(fitmc[loc][iso]["hist"][i], space+loc+'-'+iso, lopt)

        #ftotal[i].Sumw2()
        ftotal[i].Draw("same")
        flegs[i].AddEntry(ftotal[i], space+'Total MC', lopt)
        flegs[i].Draw("same")

        
        ### fit residuals plots
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        botpad[i].cd()
        leg = TLegend(0.72, 0.78, 0.94, 0.94)
        flegs2.append(leg)
        flegs2[i].SetFillColor(0)
        flegs2[i].SetBorderSize(0)
        lopt = 'LPE'

        #fresid[i].Sumw2()
        #fresid[i].Add(data[E][i], total[E][i],1,-1)
        fresid[i].Divide(fitdata[i], ftotal[i])

        fresid[i].SetTitle('')
        #fresid[i].SetXTitle("Energy (keVee)")
        fresid[i].SetXTitle("fit bins")
        fresid[i].GetXaxis().SetTitleFont(font)
        fresid[i].GetXaxis().SetTitleSize(size)
        fresid[i].GetXaxis().SetLabelFont(font)
        fresid[i].GetXaxis().SetLabelSize(size)
        fresid[i].GetXaxis().SetLabelOffset(0.03)
        fresid[i].GetXaxis().SetTitleOffset(8)
        #fresid[i].SetYTitle("counts / keV")
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
        # set my line to "1"
        zero = TLine(fmin, 1, fmax, 1)
        fzeros.append(zero)
        fzeros[i].SetLineColor(kRed)
        fzeros[i].SetLineWidth(1)
        fzeros[i].Draw()

        flegs2[i].AddEntry(fresid[i],space+"data / MC",lopt)
        flegs2[i].Draw()
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    save = ''
    if local:      save += 'local'
    else:          save += 'cup'
    save += '_Nchan-fit'
    save += '_fit-'+str(int(fLo[0]))+'-'+str(int(fLo[1]))
    save += '_fit-'+str(int(fHi[0]))+'-'+str(int(fHi[1]))
    save += '_rebin-'+str(rebin)
    if rebinscale: save += '_rebinscale'
    if mcscale:    save += '_mcscale'
    if mcweight:   save += '_mcweight'
    if fitweight:  save += '_fitweight'

    fcanv.Update()
    fcanv.Print(save+'.png')

    
    #=================================================================
    # plot the lo and hi energy histograms
    #=================================================================
    
    canvs=[None for x in range(2)]
    data=[None for x in range(2)]
    
    #mc=[]
    
    total=[None for x in range(2)]
    resid=[None for x in range(2)]
    
    legs=[None for x in range(2)]
    legs2=[None for x in range(2)]
    zeros=[None for x in range(2)]
    
    for E in range(2):
        
        legs[E]=[]
        legs2[E]=[]
        zeros[E]=[]
        
        if E:
            data[E] = dataHi
            #mc.append(mcHi)
            mc = mcHi
        else:
            data[E] = dataLo
            #mc.append(mcLo)
            mc = mcLo
            
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

            if local and not reuse:
                data[E][i].SetAxisRange(1,1000,"y")
            else:
                data[E][i].SetAxisRange(1,1000000,"y")

            #data[E][i].SetLineColor(kBlack)
            #data[E][i].SetMarkerColor(kBlack)
            #data[E][i].SetLineWidth(1)
            data[E][i].GetYaxis().SetTitle('arb. counts (for now)')
            data[E][i].GetYaxis().SetTitleFont(font)
            data[E][i].GetYaxis().SetTitleSize(size)
            data[E][i].GetYaxis().SetTitleOffset(4.2)
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
                    mc[loc][iso]["hist"][i].SetMarkerColor(cis[color])
                    mc[loc][iso]["hist"][i].SetLineColor(cis[color])
                    color += 1

                    mc[loc][iso]["hist"][i].SetAxisRange(1,1000,"y")

                    ### temp - don't draw MC if you just want to show data
                    mc[loc][iso]["hist"][i].Draw("same")

                    ### add MC to total MC hist
                    total[E][i].Add(mc[loc][iso]["hist"][i])

                    ### create the legend entry for MC
                    legs[E][i].AddEntry(mc[loc][iso]["hist"][i], space+loc+'-'+iso, lopt)

            total[E][i].Sumw2()
            total[E][i].Draw("same")
            legs[E][i].AddEntry(total[E][i], space+'Total MC', lopt)
            legs[E][i].Draw("same")


            ### try to get the residuals in!
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            botpad[i].cd()
            leg = TLegend(0.72, 0.78, 0.94, 0.94)
            legs2[E].append(leg)
            legs2[E][i].SetFillColor(0)
            legs2[E][i].SetBorderSize(0)
            lopt = 'LPE'

            #resid[E][i].Sumw2()
            #resid[E][i].Add(data[E][i], total[E][i],1,-1)
            resid[E][i].Divide(data[E][i], total[E][i])
                        
            resid[E][i].SetTitle('')
            resid[E][i].SetXTitle("Energy (keVee)")
            resid[E][i].GetXaxis().SetTitleFont(font)
            resid[E][i].GetXaxis().SetTitleSize(size)
            resid[E][i].GetXaxis().SetLabelFont(font)
            resid[E][i].GetXaxis().SetLabelSize(size)
            resid[E][i].GetXaxis().SetLabelOffset(0.03)
            resid[E][i].GetXaxis().SetTitleOffset(8)
            #resid[E][i].SetYTitle("counts / keV")
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

            legs2[E][i].AddEntry(resid[E][i],space+"data / MC",lopt)
            legs2[E][i].Draw()
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        save = ''
        if local:      save += 'local'
        else:          save += 'cup'
        if E:          save += '_hiE'
        else:          save += '_loE'
        save += '_fit-'+str(int(fLo[0]))+'-'+str(int(fLo[1]))
        save += '_fit-'+str(int(fHi[0]))+'-'+str(int(fHi[1]))
        save += '_rebin-'+str(rebin)
        if rebinscale: save += '_rebinscale'
        if mcscale:    save += '_mcscale'
        if mcweight:   save += '_mcweight'
        if fitweight:  save += '_fitweight'
        
        canvs[E].Update()
        canvs[E].Print(save+'.png')

    if not batch:
        raw_input("[Enter] to quit \n")


######################################################################
######################################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])

