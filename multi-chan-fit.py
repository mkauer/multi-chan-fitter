#!/usr/bin/env python

############################################################
#
# version: 2016-09-14
#
# Change Log
#-----------------------------------------------------------
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
    gStyle.SetPadLeftMargin  (0.11)
    gStyle.SetPadRightMargin (0.05)
    gStyle.SetTitleOffset(1.2,"x")
    gStyle.SetTitleOffset(1.0,"y")
    gStyle.SetTitleSize(0.045,"xy")
    gStyle.SetLabelSize(0.040,"xy")
    
    ch1=TChain("MC","")
    #ch1.Add("sim/data/K40/set2_internalK40-1-nEv200000.root")
    ch1.Add("sim/data/K40/set2_internalK40-2-nEv200000.root")
    #ch1.Add("sim/data/K40/*.root")
    
    volnames = ['NaIDet01Crystal', 'NaIDet02Crystal', 'NaIDet03Crystal', 'NaIDet04Crystal',
                'NaIDet05Crystal', 'NaIDet06Crystal', 'NaIDet07Crystal', 'NaIDet08Crystal']
    
    par = histparam()
    canvs = []
    data = getData(["data/dryRun/*.root"])
    
    #for i in range(0,8): # is primary crystal of origin
    for i in range(0,1): # is primary crystal of origin
        
        sims = []
        for j in range(0,8): # is surounding crystals and self
            sim = TH1F("sim", names(j), par[0], par[1], par[2])
            cut1 = TCut('event_MC.edepResol['+str(j)+']*1000. > 0')
            cut2 = TCut('event_MC.primVolumeName == "'+volnames[i]+'"')
            ch1.Draw('event_MC.edepResol['+str(j)+']*1000. >> sim',cut1+cut2)
            sim.Sumw2()
            sims.append(sim)
        
        canv = TCanvas(names(i), names(i), 0, 0, 800, 800)
        canvs.append(canv)
        canvs[i].Divide(3,3)
        for k,sim in enumerate(sims):
            canvs[i].cd(k+1)
            canvs[i].cd(k+1).SetLogy()
            if k == i:
                sim.SetMarkerColor(kRed+1)
                sim.SetLineColor(kRed+1)
            else:
                sim.SetMarkerColor(kBlue+1)
                sim.SetLineColor(kBlue+1)
            sim.SetLineWidth(1)
            #sim.Draw("same")
            #sim.Draw()
            #dat.Draw("same")
            data[k].Draw()
            sim.Draw("same")
            
        canvs[i].Update()
        canvs[i].Print('testing-'+names(i)+'.png')
                
    # pause here?
    if not noShow:
        raw_input("[Enter] to continue...\n")



def getData(rootfiles=['data/dryRun/*.root']):
    chain = TChain("ntp","")
    for rfile in rootfiles:
        chain.Add(rfile)
    data = []
    par = histparam()
    for i in range(0,8):
        dat = TH1F("dat", names(i), par[0], par[1], par[2])
        chain.Draw('crystal'+str(i+1)+'.qc *'+str(calib(i))+' >> dat')
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
    calibrations = [0.000324,
                    0.000324,
                    0.000324,
                    0.000324,
                    0.000324,
                    0.000324,
                    0.000324,
                    0.000324]
    return calibrations[int(i)]


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


def histparam():
    """
    Histogram parameters for number of bins, min, and max in keV
    """
    hmin = 0
    hmax = 1800
    bins = (hmax-hmin)/2
    return [bins, hmin, hmax]


if __name__ == "__main__":
    _myself_(sys.argv[1:])
