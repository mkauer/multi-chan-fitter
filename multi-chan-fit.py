#!/usr/bin/env python

######################################################################
# Fit the MC to data! 
#
# version: 2016-10-26
#
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
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
# ~ for some weird reason importing my "extra.py" does work
#   on the cup cluster? Very weird. So I've put all the extra
#   functions in this script - now it works on CUP. ???
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


### user inputs
#================================================================
energy  = 0   ### 0 = low-energy  -- 1 = high-energy
testing = 0   ### 0 = not-testing -- 1 = testing = one-rootfile
scale   = 1   ### pre scale the MC?
weight  = 1   ### set MC weights?
#================================================================


### automated selections...
#==================================================
### figure out if running on laptop or at CUP
local = 0
host = str(socket.gethostname())
print 'hostname =',host
if 'local' in host:
    local = 1
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
    #gStyle.SetTitleOffset(1.2,"x")
    #gStyle.SetTitleOffset(1.3,"y")
    #gStyle.SetTitleSize(0.050,"xy")
    #gStyle.SetLabelSize(0.045,"xy")

    
    ### get your data!
    if local:
        if energy:
            if testing: data = getHiEdata(["./data/phys/*1324*root.001"])
            else: data = getHiEdata(["./data/phys/*1324*root*"])
        else:
            if testing: data = getLoEdata(["./data/phys/*1324*root.001"])
            else: data = getLoEdata(["./data/phys/*1324*root*"])
    else:
        if energy:
            if testing: data = getHiEdata(["/home/mkauer/temp/*1324*root.001"])
            else: data = getHiEdata(["/home/mkauer/temp/*1324*root*"])
        else:
            if testing: data = getLoEdata(["/home/mkauer/temp/*1324*root.001"])
            else: data = getLoEdata(["/home/mkauer/temp/*1324*root*"])

    
    ### define what MC you want and from where
    locs = ['internal','pmt']
    isos = ['K40','U238','Th232']
    #locs = ['internal']
    #isos = ['U238','Th232']
    
    
    # create a color scheme for MC
    # color list and indexes
    Nmc = int((len(locs)*len(isos)))
    colors, cis = rainbow(Nmc)

    # legend length = MC + data + total
    Nlg = Nmc+2
    
    mc = {}
    for loc in locs:
        mc[loc] = {}
        for iso in isos:
            mc[loc][iso] = {}
            chain=TChain("MC","")
            if local:
                chain.Add('./sim/'+iso+'/'+'*'+loc+'*root')
            else:
                if testing:
                    chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/'+iso+'/set2/'+'set2_*'+loc+'*-1-*root')
                else:
                    chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/'+iso+'/set2/'+'set2_*'+loc+'*root')
            mc[loc][iso]["chain"] = chain
            mc[loc][iso]["hist"] = []

    
    # just to show I can iterate over this thing
    for key in mc:
        #print key
        for nkey in mc[key]:
            #print key,nkey
            for thing in mc[key][nkey]:
                print key,nkey,thing
                
    
    # get histogram parameters for MC
    if energy: par = hiEhistparam()
    else:      par = loEhistparam()

        
    for i in range(0,8): # is primary crystal of origin
        for loc in locs:
            for iso in isos:
                key = str(loc+'-'+iso+'-C'+str(i+1))
                histo = TH1F('histo', key, par[0], par[1], par[2])
                cut1 = TCut('edep['+str(i)+']*1000. > 0')
                if loc == 'internal':
                    cut2 = TCut('primVolumeName == "'+volumeNames(i)+'"')
                if loc == 'pmt':
                    cut2 = TCut('primVolumeName == "phys_pmt" ')
                
                ### test out resolution smearing
                ###-------------------------------------------------------------------------------
                ### using Box-Muller? method here for "rng" alias
                mc[loc][iso]["chain"].SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')

                ### set the resolution function - right now it's linear with intercept of 0
                if energy:
                    mc[loc][iso]["chain"].SetAlias('sigma', str(hiEres(i))+' / sqrt(edep['+str(i)+']*1000.)')
                else:
                    mc[loc][iso]["chain"].SetAlias('sigma', str(loEres(i))+' / sqrt(edep['+str(i)+']*1000.)')
                                
                ### then draw the shit
                mc[loc][iso]["chain"].Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> histo', cut1+cut2)
                ###-------------------------------------------------------------------------------

                mc[loc][iso]["hist"].append(histo)
                #mc[loc][iso]["hist"][i].Sumw2()
                #print 'this',mc[loc][iso]["hist"][i]
                

    ### fit the MC to data
    #=================================================================
    #=================================================================
    if energy:
        fmin = 100  # fit min E window
        fmax = 3000 # fit max E window
    else:
        #fmin = 5    # fit min E window
        fmin = 20    # fit min E window
        fmax = 100  # fit max E window
        
    fmc = []
    fit = []
    results = []
    for i in range(0,8):
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
    
    tot = makeTotal(par) # make a set of total MC histos
    resid = makeResid(par) # make a set of resid histos
            
    for i in range(0,8):
        canv.cd(i+1)
        canv.cd(i+1).SetLogy()

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

        if local or testing: data[i].SetAxisRange(1,1000,"y")
        else:                data[i].SetAxisRange(1,1000000,"y")
        
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
        
        if local or testing: resid[i].SetAxisRange(-100,100,"y")
        else:                resid[i].SetAxisRange(-400,400,"y")
            
        resid[i].Draw()
        
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






def getData(rootfiles=['data/dryRun/*.root']):
    """
    Build histograms of the data from all crystals
    """
    
    chain = TChain("ntp","")
    for rfile in rootfiles:
        print rfile
        chain.Add(rfile)
    data = []
    par = histparam()
    for i in range(0,8):
        dat = TH1F("dat", longNames(i), par[0], par[1], par[2])
        chain.Draw('crystal'+str(i+1)+'.qc5 *'+str(calib(i))+' >> dat')
        dat.Sumw2()
        dat.SetLineColor(kBlack)
        dat.SetMarkerColor(kBlack)
        dat.SetLineWidth(1)
        data.append(dat)
    return data


def getHiEdata(rootfiles=['data/dryRun/*.root']):
    """
    Build histograms of the data from all crystals
    """
    ### use pmtXX.rqcD1_5 for high energy
    chain = TChain("ntp","")
    for rfile in rootfiles:
        chain.Add(rfile)
    data = []
    par = hiEhistparam()
    #par = histparam()
    for i in range(0,8):
        dat = TH1F("dat", longNames(i), par[0], par[1], par[2])
        chain.Draw('(pmt'+str(i+1)+'1.rqcD1_5 *'+str(hiEcalib(i,0))
                  +' + pmt'+str(i+1)+'2.rqcD1_5 *'+str(hiEcalib(i,1))+')/2.'
                  +' >> dat')
        dat.Sumw2()
        dat.SetLineColor(kBlack)
        dat.SetMarkerColor(kBlack)
        dat.SetLineWidth(1)
        data.append(dat)
    return data


def getLoEdata(rootfiles=['data/dryRun/*.root']):
    """
    Build histograms of the data from all crystals
    """
    ### use pmtXX.qc5 for low energy
    chain = TChain("ntp","")
    for rfile in rootfiles:
        chain.Add(rfile)
    data = []
    par = loEhistparam()
    #par = histparam()
    for i in range(0,8):
        dat = TH1F("dat", longNames(i), par[0], par[1], par[2])
        chain.Draw('(pmt'+str(i+1)+'1.qc5 *'+str(loEcalib(i,0))
                  +' + pmt'+str(i+1)+'2.qc5 *'+str(loEcalib(i,1))+')/2.'
                  +' >> dat')
        dat.Sumw2()
        dat.SetLineColor(kBlack)
        dat.SetMarkerColor(kBlack)
        dat.SetLineWidth(1)
        data.append(dat)
    return data


def makeTotal(par):
    total = []
    #par = hiEhistparam()
    for i in range(0,8):
        tot = TH1F("tot", longNames(i), par[0], par[1], par[2])
        tot.SetLineColor(kGray+1)
        tot.SetMarkerColor(kGray+1)
        tot.SetLineWidth(1)
        total.append(tot)
    return total


def makeResid(par):
    resid = []
    #par = hiEhistparam()
    for i in range(0,8):
        res = TH1F("res", longNames(i), par[0], par[1], par[2])
        res.SetLineColor(kBlack)
        res.SetMarkerColor(kBlack)
        res.SetLineWidth(1)
        resid.append(res)
    return resid


def calib(i):
    """
    Return the calibrations for the crystals
    """
    # Pushba 2016-09-27 slides - WRONG!!!
    calibrations = [0.0001093478,
                    0.0001056760,
                    0.0001148900,
                    0.0001127000,
                    0.0002670000,
                    0.0001175340,
                    0.0001119000,
                    0.0003201700]
    return calibrations[int(i)]


def hiEcalib(i,j):
    """
    Return the high energy calibrations >100 keV for each PMT
    """
    # from Estella
    # (i) is the crystal and (j) is the pmt
    calibrations = [
        [0.0061480859, 0.0156416027],
	[0.0081384425, 0.0087327399],
	[0.0069518946, 0.0086127910],
	[0.0084301584, 0.0092437824],
	[0.0126542583, 0.0228144999],
	[0.0115685020, 0.0107546633],
	[0.0068681361, 0.0093653723],
	[0.0452269698, 0.0375218755]
        ]
    return calibrations[int(i)][int(j)]


def loEcalib(i,j):
    """
    Return the low energy calibrations <100 keV for each PMT
    """
    # from Estella
    # (i) is the crystal and (j) is the pmt
    calibrations = [
        [0.0002204598, 0.0002192131],
	[0.0002127388, 0.0002207598],
	[0.0002123268, 0.0002505941],
	[0.0002233761, 0.0002295871],
	[0.0006665573, 0.0005003395],
	[0.0002315343, 0.0002327635],
	[0.0002294907, 0.0002272016],
	[0.0007645054, 0.0008078809]
        ]
    return calibrations[int(i)][int(j)]


def hiEres(i):
    """
    Return the high energy resolution for combined spectra
    """
    # from Estella
    # (i) is the crystal number
    # res = sig/E
    # slope = res * sqrt(E)
    # res(E) = slope / sqrt(E) + 0
    resolution = [1.3976069981,
	          1.3796169490,
	          1.2337377995,
	          1.1778559545,
	          2.4846296071,
	          1.1813408142,
	          1.3249767662,
	          2.7941471926]
    return resolution[int(i)]


def loEres(i):
    """
    Return the high energy resolution for combined spectra
    """
    # from Estella
    # (i) is the crystal number
    # res = sig/E
    # slope = res * sqrt(E)
    # res(E) = slope / sqrt(E) + 0
    resolution = [0.8174646681,
	          0.7927112730,
	          0.7274639316,
	          0.6471664710,
	          1.2783412699,
	          0.7201966764,
	          0.7240077873,
	          1.5955441464]
    return resolution[int(i)]


def histparam():
    """
    Histogram parameters for number of bins, min, and max in keV
    """
    hmin = 0
    hmax = 1800
    bins = (hmax-hmin)/4
    return [bins, hmin, hmax]


def loEhistparam():
    """
    Histogram parameters for the low energy fitting
    """
    hmin = 0
    hmax = 100
    bins = (hmax-hmin)
    return [bins, hmin, hmax]


def hiEhistparam():
    """
    Histogram parameters for the high energy fitting
    """
    hmin = 0
    hmax = 3000
    bins = (hmax-hmin)
    return [bins, hmin, hmax]


def names(i):
    """
    The names of the crystals
    """
    
    crystals = ['NaI-01',
                'NaI-02',
                'NaI-03',
                'NaI-04',
                'NaI-05',
                'NaI-06',
                'NaI-07',
                'NaI-08']
    return crystals[int(i)]


def longNames(i):
    """
    Full name and specs for the crystal
    """
    crystals = ['C1  NaI-001  Sample-B  8.3kg',
                'C2  NaI-002  Sample-C  9.2kg',
                'C3  NaI-007  WimpScint-2  9.2kg',
                'C4  AS-3  WimpScint-2  18.5kg',
                'C5  AS-1  Sample-C  18.5kg',
                'C6  NaI-011  WimpScint-3  12.5kg',
                'C7  NaI-012  WimpScint-3  12.5kg',
                'C8  AS-2  Sample-C  18.5kg']
    
    return crystals[int(i)]


def volumeNames(i):
    """
    For the simulation primary volume names
    """
    volumeNames = ['NaIDet01Crystal',
                   'NaIDet02Crystal',
                   'NaIDet03Crystal',
                   'NaIDet04Crystal',
                   'NaIDet05Crystal',
                   'NaIDet06Crystal',
                   'NaIDet07Crystal',
                   'NaIDet08Crystal']
    return volumeNames[int(i)]

    
def rainbow(N):
    from random import randint
    colors=[]
    cis=[]
    # ci cannot be too large > 10,000!
    ci = randint(1000, 5000)
    for h in range(N):
        H = float(h)/float(N)
        ci += 1
        if H <= 1/5. :
            R=1.
            G=1.*5*H
            B=0.
        elif H > 1/5. and H <= 2/5. :
            R=1.-(1*5*(H-1/5.))
            G=1.
            B=0.
        elif H > 2/5. and H <= 3/5. :
            R=0.
            G=1.
            B=1.*5*(H-2/5.)
        elif H > 3/5. and H <= 4/5. :
            R=0.
            G=1.-(1*5*(H-3/5.))
            B=1.
        elif H > 4/5. and H <= 1. :
            R=1.*5*(H-4/5.)
            G=0.
            B=1.
        elif H > 1. :
            R=1.
            G=1.
            B=1.
        
        color = TColor(ci, R, G, B)
        # must keep the color and the ci in memory
        colors.append(color)
        cis.append(ci)
        
    return colors, cis


######################################################################
######################################################################


if __name__ == "__main__":
    _myself_(sys.argv[1:])

