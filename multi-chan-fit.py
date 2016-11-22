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

    sims=[]
    for i in range(0,1):
        #del.sim
        h_sim1 = TH1F('h_sim1','h_sim1',bins,hmin,hmax)
        #sim = TH1F('sim','sim',bins,hmin,hmax)
        #sim.append(h_sim1)
        ch1.Draw('edepResol[0]*1000. >> h_sim1','edepResol[0]*1000. > 0 && primVolumeName == "NaIDet01Crystal"')
        #ch1.Draw('edepResol[0]*1000. >> sim','edepResol[0]*1000. > 0 && primVolumeName == "NaIDet01Crystal"')
        #sims.append(sim)
    
    c1 = TCanvas("c1","c1",0,0,800,800)
    c1.Divide(4,4)
    c1.cd(1)
    c1_1.SetLogy()

    h_sim1.SetMarkerColor(kRed+1)
    h_sim1.SetLineColor(kRed+1)
    h_sim1.SetLineWidth(1)
    h_sim1.Draw("same")
    
    #sim.SetMarkerColor(kRed+1)
    #sim.SetLineColor(kRed+1)
    #sim.SetLineWidth(1)
    #sim.Draw("same")
    
    c1.Update()
    #c1.Print("pytest-fitting.png")
    

    
    # see plots hack
    raw_input("[Enter] to quit\n")

    
if __name__ == "__main__":
    _myself_(sys.argv[1:])
