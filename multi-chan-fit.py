#!/usr/bin/env python

######################################################################
# Fit the MC to data! 
#
# version: 2016-10-25
#
# Change Log
#---------------------------------------------------------------------
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
# mkauer@physics.wisc.edu
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
######################################################################

import os,sys
import socket
import numpy as np
from ROOT import *
import ROOT


### figure out if running on laptop or at CUP
local = 0
host = str(socket.gethostname())
if 'local' in host:
    local = 1

### 0 == show-plots  1 == dont-show-plots
batch = 0
if not local:
    batch = 1
gROOT.SetBatch(batch)


def _myself_(argv):

    gROOT.Reset()
    gStyle.SetPalette(1)
    gStyle.SetOptStat("")
    gStyle.SetOptFit(0)
    
    gStyle.SetPadBottomMargin(0.12)
    gStyle.SetPadLeftMargin  (0.12)
    gStyle.SetPadRightMargin (0.05)
    gStyle.SetTitleOffset(1.2,"x")
    gStyle.SetTitleOffset(1.3,"y")
    gStyle.SetTitleSize(0.050,"xy")
    gStyle.SetLabelSize(0.045,"xy")


    

    ### get your data!
    if local:
        data = getHiEdata(["./data/phys/*1324*root*"])
        #data = getLoEdata(["./data/phys/*1324*root*"])
    else:
        ### ALL files
        data = getHiEdata(["/home/mkauer/temp/*1324*root*"])
        ### Just one file
        #data = getHiEdata(["/home/mkauer/temp/*1324*root.001"])
        #data = getLoEdata(["/home/mkauer/temp/*1324*root*"])

    ### define what MC you want and from where
    isos = ['K40','U238','Th232']
    #locs = ['internal','pmt']
    #isos = ['U238','Th232']
    locs = ['internal']
    
    # legend length = MC + data + total-mc
    Nlg = int((len(isos)*len(locs))+2)
    #print 'legend length =',ln

    # create a color scheme for MC
    # color list and indexes
    Nmc = int((len(isos)*len(locs)))
    colors, cis = rainbow(Nmc)
    
    mc = {}
    for iso in isos:
        mc[iso] = {}
        for loc in locs:
            mc[iso][loc] = {}
            chain=TChain("MC","")
            if local:
                chain.Add('./sim/'+iso+'/'+'*'+loc+'*root')
            else:
                ### ALL files
                chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/'+iso+'/set2/'+'set2_*'+loc+'*root')
                ### Just one file
                #chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/'+iso+'/set2/'+'set2_*'+loc+'*-1-*root')
            mc[iso][loc]["chain"] = chain
            mc[iso][loc]["hist"] = []

    
    # just to show I can iterate over this thing
    for key in mc:
        #print key
        for nkey in mc[key]:
            #print key,nkey
            for thing in mc[key][nkey]:
                print key,nkey,thing
                
    
    # get histogram parameters for MC
    par = hiEhistparam()
    
    for i in range(0,8): # is primary crystal of origin
        for iso in isos:
            for loc in locs:
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
                mc[iso][loc]["chain"].SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')

                ### set the resolution function - right now it's linear with intercept of 0
                #mc[iso][loc]["chain"].SetAlias('sigma', str(loEres(i))+' / sqrt(edep['+str(i)+']*1000.)')
                mc[iso][loc]["chain"].SetAlias('sigma', str(hiEres(i))+' / sqrt(edep['+str(i)+']*1000.)')
                                
                ### then draw the shit
                mc[iso][loc]["chain"].Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> histo', cut1+cut2)
                ###-------------------------------------------------------------------------------

                mc[iso][loc]["hist"].append(histo)
                #mc[iso][loc]["hist"][i].Sumw2()
                #print 'this',mc[iso][loc]["hist"][i]
                

    ### fit the MC to data
    #=================================================================
    #=================================================================
    fmin = 100 # fit min E window
    fmax = 3000 # fit max E window

    fmc = []
    fit = []
    results = []
    for i in range(0,8):
        fmc.append(TObjArray(Nmc)) # number of MC to fit to
        dat_int = data[i].Integral(fmin,fmax) # data integral to normalize to
        #print i,'dat_int =',dat_int
        for iso in isos:
            for loc in locs:
                mc_int = mc[iso][loc]["hist"][i].Integral(fmin,fmax) # MC integral
                #print iso,loc,'mc_int =',mc_int

                ### to scale or not to scale...
                #---------------------------------------------------------------------
                mc[iso][loc]["hist"][i].Scale(dat_int/mc_int) # scale to data integral
                mc[iso][loc]["hist"][i].Sumw2() # reset weights
                #---------------------------------------------------------------------
                
                fmc[i].Add(mc[iso][loc]["hist"][i]) # add to the TFractionFitter object
                
        fit.append(TFractionFitter(data[i], fmc[i])) # create the TFF data and MC objects
        
        for l in range(Nmc):
            # set bounds on the MC put into fmc TObject
            fit[i].Constrain(l, 0.01, 10.0)
        
        fit[i].SetRangeX(fmin, fmax) # set global range, should be the same as fmin and fmax?
        
        status = fit[i].Fit() # do the fit!
        chi2 = fit[i].GetChisquare()
        ndf = fit[i].GetNDF()
        #print "Fit done with status ",status," -  chi2/ndf =",chi2,"/",ndf
        results.append("NaI-C"+str(i+1)+" fit completed with status "+str(status))
        results.append("fit chi2/ndf = "+str(round(chi2,2))+"/"+str(ndf)+" = "+str(round(chi2/float(ndf),2)))
        
        for j, iso in enumerate(isos):
            for k, loc in enumerate(locs):
                value = ROOT.Double(0.0)
                error = ROOT.Double(0.0)
                l = j+k
                #print l
                fit[i].GetResult(l, value, error)
                #print iso,loc,'=',value,' +/- ',error
                results.append(str(iso)+' '+str(loc)+' = '+str(round(value,4))+' +/- '+str(round(error,4)))
                mc[iso][loc]["hist"][i].Scale(value)
        results.append('\n')
    
    print '\n\n'
    for line in results:
        print line
    #=================================================================
    #=================================================================
        
    # have the plotting be seperated out from the loop
    canv = TCanvas('canv', 'canv', 0, 0, 1600, 900)
    canv.Divide(4,2)
    legs=[]
    tot = makeHiEtotal() # make a set of total MC histos
    
    for i in range(0,8):
        canv.cd(i+1)
        canv.cd(i+1).SetLogy()
        #yleg = (1.-(6*ln*0.01))
        ylegstart = 0.88
        ylegend = (ylegstart-(Nlg*0.04))
        space = '  '
        leg = TLegend(0.60, ylegend, 0.94, ylegstart)
        legs.append(leg)
        legs[i].SetFillColor(0)
        legs[i].SetBorderSize(0)
        lopt = 'LPE'
        
        data[i].GetYaxis().SetTitle('arb. counts')
        data[i].GetXaxis().SetTitle('Energy (keV)')
        data[i].Draw()
        legs[i].AddEntry(data[i], space+'data', lopt)

        for j, iso in enumerate(isos):
            for k, loc in enumerate(locs):
                ### set mc colors
                mc[iso][loc]["hist"][i].SetMarkerColor(cis[j+k])
                mc[iso][loc]["hist"][i].SetLineColor(cis[j+k])
                
                ### temp - don't draw MC if you just want to show data
                mc[iso][loc]["hist"][i].Draw("same")

                ### add MC to total MC hist
                tot[i].Add(mc[iso][loc]["hist"][i])

                ### create the legend entry for MC
                legs[i].AddEntry(mc[iso][loc]["hist"][i], space+loc+'-'+iso, lopt)

        tot[i].Sumw2()
        tot[i].Draw("same")
        legs[i].AddEntry(tot[i], space+'Total MC', lopt)
        legs[i].Draw("same")
        
        # put data back on top
        #data[i].Draw('same')
            
    canv.Update()
    canv.Print('testing-hi-energy-fit.png')
    
    ### draw the residuals
    #-----------------------------------------------------------------
    canv2 = TCanvas('canv2', 'canv2', 0, 0, 1600, 900)
    canv2.Divide(4,2)
    legs2=[]
    resid = makeHiEresid() # make a set of resid histos
    zeros=[]
    
    for i in range(0,8):
        canv2.cd(i+1)
        ylegstart = 0.88
        ylegend = (ylegstart-(0.05))
        leg = TLegend(0.55, ylegend, 0.94, ylegstart)
        legs2.append(leg)
        legs2[i].SetFillColor(0)
        legs2[i].SetBorderSize(0)
        lopt = 'LPE'
        
        resid[i].Add(data[i],tot[i],1,-1)
        resid[i].SetXTitle("Energy keV")
        resid[i].SetYTitle("counts / keV")
        resid[i].SetAxisRange(-100,100,"y")
        resid[i].Draw()
        
        zero = TLine(par[1],0,par[2],0)
        zeros.append(zero)
        zeros[i].SetLineColor(kRed)
        zeros[i].SetLineWidth(1)
        zeros[i].Draw()
        
        legs2[i].AddEntry(resid[i],"data-MC residual",lopt)
        legs2[i].Draw()

    canv2.Update()
    canv2.Print('testing-hi-energy-fit-resid.png')
    #-----------------------------------------------------------------
    
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


def makeHiEtotal():
    total = []
    par = hiEhistparam()
    for i in range(0,8):
        tot = TH1F("tot", longNames(i), par[0], par[1], par[2])
        tot.SetLineColor(kGray+1)
        tot.SetMarkerColor(kGray+1)
        tot.SetLineWidth(1)
        total.append(tot)
    return total


def makeHiEresid():
    resid = []
    par = hiEhistparam()
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
    bins = (hmax-hmin)*2
    return [bins, hmin, hmax]


def hiEhistparam():
    """
    Histogram parameters for the high energy fitting
    """
    
    hmin = 0
    hmax = 3000
    # fit energy needs to be normalized to bins???
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

if __name__ == "__main__":
    _myself_(sys.argv[1:])
