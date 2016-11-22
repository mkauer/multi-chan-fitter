#!/usr/bin/env python

import os,sys
import numpy as np
from ROOT import *
import ROOT

def _myself_(argv):
    
    gROOT.Reset()
    gStyle.SetPalette(1)
    gStyle.SetOptStat("")
    gStyle.SetOptFit(0)
    
    gStyle.SetPadBottomMargin(0.12)
    gStyle.SetPadLeftMargin  (0.11)
    #gStyle.SetPadTopMargin   (0.02)
    gStyle.SetPadRightMargin (0.05)
    
    gStyle.SetTitleOffset(1.2,"x")
    gStyle.SetTitleOffset(1.0,"y")
    gStyle.SetTitleSize(0.045,"xy")
    gStyle.SetLabelSize(0.040,"xy")
    
    ch1=TChain("MC","")

    ch1.Add("sim/data/K40/set2_internalK40-1-nEv200000.root")
    
    hmin = 0
    hmax = 1800
    R = 2
    bins = (hmax-hmin)/R

    volnames = ['NaIDet01Crystal','NaIDet02Crystal','NaIDet03Crystal','NaIDet04Crystal',
                'NaIDet05Crystal','NaIDet06Crystal','NaIDet07Crystal','NaIDet08Crystal']
    names = ['NaI-01','NaI-02','NaI-03','NaI-04',
             'NaI-05','NaI-06','NaI-07','NaI-08']
    sims=[]
    for i in range(0,8):
        sim = TH1F("sim", names[i], bins, hmin, hmax)
        ch1.Draw('edepResol['+str(i)+']*1000. >> sim','edepResol['+str(i)+']*1000. > 0 && primVolumeName == "'+volnames[0]+'"')
        sims.append(sim)
    
    c1 = TCanvas("c1","c1",0,0,800,800)
    c1.Divide(3,3)
    for i,sim in enumerate(sims):
        i=i+1
        c1.cd(i)
        c1.cd(i).SetLogy()
        
        sim.SetMarkerColor(kRed+1)
        sim.SetLineColor(kRed+1)
        sim.SetLineWidth(1)
        #sim.Draw("same")
        sim.Draw()
        
    c1.Update()
    #c1.Print("pytest-fitting.png")
    

    
    # see plots hack
    raw_input("[Enter] to quit\n")

    
if __name__ == "__main__":
    _myself_(sys.argv[1:])
