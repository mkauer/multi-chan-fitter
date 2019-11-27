#!/usr/bin/env python

import os,sys,re
import shutil
import socket
from copy import deepcopy
import math
import numpy as np
import datetime

import ROOT
from ROOT import *
ROOT.gROOT.Reset()
ROOT.gErrorIgnoreLevel = kError


# write new rootfile?
WRITE = 0


def main():
    
    base  = './root-join-read/'
    build = 'build311'

    allhistos = []
    
    for x in ['1','2','3','4','6','7']:
    #for x in ['1']:
        
        X = '-c'+x+'-'
        
        rootfile01 = base+build+'/'+build+X+'nai-surf-0.01umPb210_GRND.root'
        rootfile1  = base+build+'/'+build+X+'nai-surf-1umPb210_GRND.root'
        rootfile10 = base+build+'/'+build+X+'nai-surf-10umPb210_GRND.root'
        rfile01 = TFile(rootfile01, "READ")
        rfile1  = TFile(rootfile1,  "READ")
        rfile10 = TFile(rootfile10, "READ")

        if WRITE: newrootfile = base+build+'/'+build+X+'nai-surf-expoPb210_GRND.root'
        if WRITE: newrfile = TFile(newrootfile, 'RECREATE')
        allhistos = []
        
        for e in ['0','1']:
            for c in ['S','M']:
                for f in ['1','2','3','4','5','6','7','8']:
                    
                    # define the hist keys
                    key01 = 'x'+x+'-naisurf0.01um-Pb210_GRND-f'+f+'-c'+c+'-e'+e
                    key1  = 'x'+x+'-naisurf1um-Pb210_GRND-f'+f+'-c'+c+'-e'+e
                    key10 = 'x'+x+'-naisurf10um-Pb210_GRND-f'+f+'-c'+c+'-e'+e

                    #print key1

                    # copy the histograms
                    hist01 = deepcopy(TH1F(rfile01.Get(key01)))
                    hist1  = deepcopy(TH1F(rfile1.Get(key1)))
                    hist10 = deepcopy(TH1F(rfile10.Get(key10)))

                    # get the number of generated events
                    gen01 = TH1F(rfile01.Get(key01+'_generated')).GetEntries()
                    gen1  = TH1F(rfile1.Get(key1+'_generated')).GetEntries()
                    gen10 = TH1F(rfile10.Get(key10+'_generated')).GetEntries()

                    #print gen01, gen1, gen10
                    #print gen01/gen1
                    #print gen01/gen10
                    
                    # normalize the counts to per 1um and scale by expo weight
                    # except: pass - for divide by zero errors of empty histos
                    expo01 = 0.28
                    expo1  = 0.16
                    expo10 = 0.005
                    #try: hist01.Scale((gen1/gen01) * expo01)
                    try: hist01.Scale((gen1/gen01) * expo01 * 0.01)
                    except: pass
                    #try: hist1.Scale((gen1/gen1) * expo1)
                    try: hist1.Scale((gen1/gen1) * expo1 * 1.0)
                    except: pass
                    #try: hist10.Scale((gen1/gen10) * expo10)
                    try: hist10.Scale((gen1/gen10) * expo10 * 10.0)
                    except: pass
                    
                    newkey = 'x'+x+'-naisurfexpo-Pb210_GRND-f'+f+'-c'+c+'-e'+e
                    nbins = int(hist01.GetNbinsX())
                    xmin = int(hist01.GetBinLowEdge(1))
                    xmax = int(hist01.GetBinLowEdge(nbins+1))

                    print newkey
                    print nbins, xmin, xmax

                    hist = TH1F(newkey, newkey, nbins, xmin, xmax)
                    hist.Add(hist01)
                    hist.Add(hist1)
                    hist.Add(hist10)

                    #canv1 = TCanvas('canv1', 'canv1', 0, 0, 900, 600)
                    #canv1.SetLogy(1)
                    #hist.Draw()
                    #hist01.Draw()
                    #hist1.Draw()
                    #hist10.Draw()

                    hist_gen = TH1F(newkey+'_generated','generated',1,0,1)
                    hist_gen.SetBinContent(1, gen01+gen1+gen10)
                    hist_gen.SetEntries(gen01+gen1+gen10)

                    #canv2 = TCanvas('canv2', 'canv2', 0, 0, 900, 600)
                    #hist_gen.Draw()

                    allhistos.append(hist)
                    allhistos.append(hist_gen)
                    
                    #raw_input('[Enter] to quit \n')

        if WRITE: newrfile.Write()
        if WRITE: newrfile.Close()
            
    #raw_input('[Enter] to quit \n')
    return


if __name__ == "__main__":
    main()

