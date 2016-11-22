#!/usr/bin/env python

######################################################################
# 10-mc-func-fit.py
#
# Fit the MC to data!
# Note: This version and higher will use the funcs.py!
# 
# version: 2016-11-01
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ resolution function should be divided by 2 (see funcs.py)
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
# + get p-value from fit
# ~ tweaked the plotting cosmetics a bit to make it look better,
#   good grief, what a bitch...
# + put residuals on the same canvas! awesome!
# ~ add options to scale and/or weight
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
# ~ plot low energy things - ah, use pmtXX.qc5 no qc_5 - oops
# + added low energy calib and resolution functions
# ~ use assumption that res~srqt(E) with intercept=0
# + working on resolution params, but need at least a linear
#   fit to the resolution points to get anything reasonable
# + working on hiE parameters
# + shit! finally got colors to work. What a pain in the butt.
#   The color index can't be too big > 10,000 and you have to
#   save the list of colors and indexes or they get dumped
#   out of memory and then it doesn't work right.
# + urgh... finally got hists into dictionary
#   hist name must be same as hist var-name (don't know why)
# + start adding in external backgrounds too
# ~ use Pushpa's data calibrations - wrong but a start
# ~ just use this as the default so I don't have multiple
#   versions of scripts around
# ~ for some weird reason importing my "extra.py" doesn't work
#   on the cup cluster? Very weird. So I've put all the extra
#   functions in this script and now it works on CUP. ???
# ~ moved sim stuff to under CUP
# ~ tweaked spacings and labels on the canvas
# + got legends working for the plots
# ~ tweaked the layout to be 4x2 as real layout
# + seperate functions import for extra functions
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


### user inputs
#================================================================
reuse   = 0   ### use joined rootfile data?
energy  = 0   ### 0 = low-energy  -- 1 = high-energy
scale   = 1   ### pre scale the MC?
weight  = 1   ### set MC weights?
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
        rootfile = "/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/joined-master.root"
        #rootfile = "/home/mkauer/COSINE/CUP/mc-fitting/root-join-read/join-test.root"
        data, mc, locs, isos = readROOT(rootfile, energy)
    else:
        if local: path = "/home/mkauer/COSINE/CUP/mc-fitting/data/phys/"
        else:     path = "/home/mkauer/temp/"
        data = getData([path+"*1324*root*"], energy)
        
        ### define what MC you want and from where
        locs = ['internal','pmt']
        isos = ['K40','U238','Th232']
        mc = buildMC(locs, isos, energy)
    
    # maintain the same fitting order and what not?
    # will need to FIX this eventually...
    locs = ['internal','pmt']
    isos = ['K40','U238','Th232']
    #locs = ['pmt','internal']
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
    if energy:
        fmin = 100  # fit min E window
        fmax = 3000 # fit max E window
    else:
        fmin = 15   # fit min E window
        fmax = 100  # fit max E window
        
    fmc = []
    fit = []
    results = []
    for i in range(8):
        fmc.append(TObjArray(Nmc)) # number of MC to fit to
        dat_int = data[i].Integral(fmin,fmax) # data integral to normalize to
        #print i,'dat_int =',dat_int
        for loc in locs:
            for iso in isos:
                mc_int = mc[loc][iso]["hist"][i].Integral(fmin,fmax) # MC integral
                #print iso,loc,'mc_int =',mc_int
                
                
                #-----------------------------------------------------------------------
                #=======================================================================
                ########################################################################
                
                ### to scale or not to scale...
                if scale:
                    mc[loc][iso]["hist"][i].Scale(dat_int/mc_int) # scale to data integral
                
                ### to weight or not to weight...
                if weight:
                    mc[loc][iso]["hist"][i].Sumw2() # set stat weights
                
                ########################################################################
                #=======================================================================
                #-----------------------------------------------------------------------
                
                
                fmc[i].Add(mc[loc][iso]["hist"][i]) # add to the TFractionFitter object
                
        fit.append(TFractionFitter(data[i], fmc[i])) # create the TFF data and MC objects
        
        for l in range(Nmc):
            # set bounds on the MC put into fmc TObject
            fit[i].Constrain(l, 0.0001, 10.0)
        
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
                mc[loc][iso]["hist"][i].Scale(value)
                count += 1
        results.append('\n')

    print '\n\n'
    for line in results:
        print line
    #=================================================================
    #=================================================================
    
    # have the plotting be seperated out from the loop
    canv = TCanvas('canv', 'canv', 0, 0, 1400, 900)
    canv.Divide(4,2)
    toppad=[]
    botpad=[]
    legs=[]
    legs2=[]
    zeros=[]
    
    font=63
    size=13
    
    tot = makeTotal(energy) # make a set of total MC histos
    resid = makeResid(energy) # make a set of resid histos
    
    for i in range(8):
        canv.cd(i+1)
        #canv.cd(i+1).SetLogy()
        
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
        toppad[i].Draw()
        botpad[i].Draw()
        toppad[i].cd()
        
        ylegstart = 0.88
        ylegend = (ylegstart-(Nlg*0.06))
        space = '  '
        leg = TLegend(0.65, ylegend, 0.94, ylegstart)
        legs.append(leg)
        legs[i].SetFillColor(0)
        legs[i].SetBorderSize(0)
        lopt = 'LPE'
        
        if local and not reuse:
            data[i].SetAxisRange(1,1000,"y")
        else:
            data[i].SetAxisRange(1,1000000,"y")
        
        data[i].GetYaxis().SetTitle('arb. counts (for now)')
        data[i].GetYaxis().SetTitleFont(font)
        data[i].GetYaxis().SetTitleSize(size)
        data[i].GetYaxis().SetTitleOffset(4.2)
        data[i].GetYaxis().SetLabelFont(font)
        data[i].GetYaxis().SetLabelSize(size)
        data[i].GetYaxis().SetLabelOffset(0.01)
        #data[i].GetXaxis().SetTitle('Energy (keV)')
        #data[i].GetXaxis().SetLabelFont(font)
        #data[i].GetXaxis().SetLabelSize(size)
        data[i].Draw()
        legs[i].AddEntry(data[i], space+'data', lopt)
        
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
                tot[i].Add(mc[loc][iso]["hist"][i])
                
                ### create the legend entry for MC
                legs[i].AddEntry(mc[loc][iso]["hist"][i], space+loc+'-'+iso, lopt)
        
        tot[i].Sumw2()
        tot[i].Draw("same")
        legs[i].AddEntry(tot[i], space+'Total MC', lopt)
        legs[i].Draw("same")
        
        
        ### try to get the residuals in!
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        botpad[i].cd()
        leg = TLegend(0.72, 0.78, 0.94, 0.94)
        legs2.append(leg)
        legs2[i].SetFillColor(0)
        legs2[i].SetBorderSize(0)
        lopt = 'LPE'
        
        resid[i].Add(data[i],tot[i],1,-1)
        resid[i].SetTitle('')
        resid[i].SetXTitle("Energy (keVee)")
        resid[i].GetXaxis().SetTitleFont(font)
        resid[i].GetXaxis().SetTitleSize(size)
        resid[i].GetXaxis().SetLabelFont(font)
        resid[i].GetXaxis().SetLabelSize(size)
        resid[i].GetXaxis().SetLabelOffset(0.03)
        resid[i].GetXaxis().SetTitleOffset(8)
        #resid[i].SetYTitle("counts / keV")
        resid[i].GetYaxis().SetLabelFont(font)
        resid[i].GetYaxis().SetLabelSize(size)
        resid[i].GetYaxis().SetLabelOffset(0.01)
        resid[i].GetYaxis().SetNdivisions(505) # '5' secondary and '05' primary
        
        if local and not reuse:
            resid[i].SetAxisRange(-100,100,"y")
        else:
            resid[i].SetAxisRange(-400,400,"y")
        
        resid[i].Draw()

        par = histparam(energy)
        zero = TLine(par[1],0,par[2],0)
        zeros.append(zero)
        zeros[i].SetLineColor(kRed)
        zeros[i].SetLineWidth(1)
        zeros[i].Draw()
        
        legs2[i].AddEntry(resid[i],space+"residual",lopt)
        legs2[i].Draw()
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    save = ''
    if local:  save += 'local'
    else:      save += 'cup'
    save += '_fit-'+str(int(fmin))+'-'+str(int(fmax))
    #if energy: save += '_hi-energy'
    #else:      save += '_lo-energy'
    if scale:  save += '_scaled'
    if weight: save += '_weighted'
    canv.Update()
    canv.Print(save+'.png')
    
    if not batch:
        raw_input("[Enter] to quit \n")


######################################################################
######################################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])

