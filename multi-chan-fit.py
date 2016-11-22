#!/usr/bin/env python

############################################################
# Try to get everything into a dictionary or list
#
# version: 2016-10-24
#
# Change Log
#-----------------------------------------------------------
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
############################################################

import os,sys
import numpy as np
from ROOT import *
import ROOT


noShow = 0  ### 0 == show-plots  1 == dont-show-plots
gROOT.SetBatch(noShow)

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


    #data = getData(["data/dryRun/*root"])
    #data = getHiEdata(["data/phys/*root*"])
    #data = getData(["../temp/*1324*root*"])
    data = getHiEdata(["data/phys/*1324*root*"])
    
    isos = ['K40','U238','Th232']
    locs = ['internal','pmt']
    #isos = ['K40']
    #locs = ['internal']
    
    # legend length = MC + data
    ln = int((len(isos)*len(locs))+1)
    # color list and indexes
    colors, cis = rainbow(ln-1)
    
    mc = {}
    for iso in isos:
        mc[iso] = {}
        for loc in locs:
            mc[iso][loc] = {}
            chain=TChain("MC","")
            chain.Add('sim/'+iso+'/*'+loc+'*root')
            mc[iso][loc]["chain"] = chain
            mc[iso][loc]["hist"] = []

    
    # just to show I can iterate over this thing
    for key in mc:
        #print key
        for nkey in mc[key]:
            #print key,nkey
            for thing in mc[key][nkey]:
                print key,nkey,thing
                
    
    # get standard histogram parameters
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
                ### ch1.SetAlias("rng","sin(2.*pi*rndm)*sqrt(-2.*log(rndm))")
                mc[iso][loc]["chain"].SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                ### the above seems to work fine...
                
                ### this goes by my res formula of res = p0 * e^(-E+p1/p2) - E/p3 + p4
                ### /home/mkauer/sandbox/calibration/plots-resolution/
                ### ch1.SetAlias("sigma","0.435866*exp(-(EnergyDep+0)/0.0291485)-EnergyDep/205.408+0.0530099")
                #mc[iso][loc]["chain"].SetAlias('sigma', '0.435866 * exp(-(edep['+str(i)+'] + 0.) / 0.0291485) - edep['+str(i)+'] / 205.408 + 0.0530099')
                ### the above works! just not right for this case!

                ### and this should be linear right?
                ### ah, but still need an intercept, this assumes the intercept is 0
                ### crap, add some offset based off dmice17
                ### just guessing here.... actually doesn't look too bad... :)
                ### ah, but energy resolution is going up with energy?? That's not right...
                ### acoroding to estella's numbers, this should be 8e-5 or 0.00008
                ### CRAP - NO! I need a at lest a line, two points, one point won't do it...
                ### ch1.SetAlias("sigma","0.435866*exp(-(EnergyDep+0)/0.0291485)-EnergyDep/205.408+0.0530099")
                #mc[iso][loc]["chain"].SetAlias('sigma', '(-0.0001*edep['+str(i)+']*1000.) + 0.05')
                #mc[iso][loc]["chain"].SetAlias('sigma', str(hiEres(i, 'edep['+str(i)+']*1000.')))
                mc[iso][loc]["chain"].SetAlias('sigma', str(hiEres(i))+' / sqrt(edep['+str(i)+']*1000.)')
                
                ### maybe E needs to be defined?
                #mc[iso][loc]["chain"].SetAlias('E', 'edep['+str(i)+']*1000.')
                #mc[iso][loc]["chain"].SetAlias('sigma', str(hiEres(i, 'E')))
                
                
                ### then draw the shit
                ### ch1.Draw("EnergyDep*1000. + sigma*EnergyDep*1000.*rng >> h_sim1","EnergyDep > 0")
                mc[iso][loc]["chain"].Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> histo', cut1+cut2)
                ###-------------------------------------------------------------------------------

                ### skip this guy below for the E res smearing testing above ^^
                #mc[iso][loc]["chain"].Draw('edep['+str(i)+']*1000. >> histo', cut1+cut2)
                
                mc[iso][loc]["hist"].append(histo)
                mc[iso][loc]["hist"][i].Sumw2()

    # have the plotting be seperated out from the loop
    canv = TCanvas('canv', 'canv', 0, 0, 1600, 900)
    canv.Divide(4,2)
    legs=[]

    for i in range(0,8):
        canv.cd(i+1)
        canv.cd(i+1).SetLogy()
        yleg = (1.-(6*ln*0.01))
        space = '  '
        leg = TLegend(0.60, yleg, 0.94, 0.88)
        legs.append(leg)
        
        data[i].GetYaxis().SetTitle('arb. counts')
        data[i].GetXaxis().SetTitle('Energy (keV)')
        data[i].Draw()

        legs[i].SetFillColor(0)
        legs[i].SetBorderSize(0)
        lopt = 'LPE'
        legs[i].AddEntry(data[i], space+'data', lopt)

        count = 0
        for j, loc in enumerate(locs):
            for k, iso in enumerate(isos):
                mc[iso][loc]["hist"][i].SetMarkerColor(cis[count])
                mc[iso][loc]["hist"][i].SetLineColor(cis[count])
                
                ### temp don't draw MC if you just want to show data
                mc[iso][loc]["hist"][i].Draw("same")

                legs[i].AddEntry(mc[iso][loc]["hist"][i], space+loc+'-'+iso, lopt)
                count += 1
        legs[i].Draw("same")
        
        # put data back on top
        data[i].Draw('same')
        canv.Update()
    
    canv.Update()
    canv.Print('testing-hi-energy.png')
    
    if not noShow:
        if raw_input("[Enter] to continue, [q] to quit \n") == 'q':
            sys.exit()


############################################################


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

    # use pmtXX.rqcD1_5 for high energy
    
    chain = TChain("ntp","")
    for rfile in rootfiles:
        chain.Add(rfile)
    data = []
    par = hiEhistparam()
    #par = histparam()
    for i in range(0,8):
        dat = TH1F("dat", longNames(i), par[0], par[1], par[2])
        #chain.Draw('crystal'+str(i+1)+'.qc5 *'+str(calib(i))+' >> dat')
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

    # use pmtXX.qc_5 for low energy
    
    chain = TChain("ntp","")
    for rfile in rootfiles:
        chain.Add(rfile)
    data = []
    par = loEhistparam()
    #par = histparam()
    for i in range(0,8):
        dat = TH1F("dat", longNames(i), par[0], par[1], par[2])
        #chain.Draw('crystal'+str(i+1)+'.qc5 *'+str(calib(i))+' >> dat')
        chain.Draw('(pmt'+str(i+1)+'1.qc_5 *'+str(loEcalib(i,0))
                  +' + pmt'+str(i+1)+'2.qc_5 *'+str(loEcalib(i,1))+')/2.'
                  +' >> dat')
        dat.Sumw2()
        dat.SetLineColor(kBlack)
        dat.SetMarkerColor(kBlack)
        dat.SetLineWidth(1)
        data.append(dat)
    return data


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
    bins = (hmax-hmin)/4
    return [bins, hmin, hmax]


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


############################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])
