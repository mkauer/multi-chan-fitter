#!/usr/bin/env python

############################################################
#
# version: 2016-10-05
#
# Change Log
#-----------------------------------------------------------
# + start adding in external backgrounds too
# ~ use Pushpa's data calibrations - wrong
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


noShow = 1  ### 0 == show-plots  1 == dont-show-plots
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
    data = getData(["data/phys/*root*"])
    #data = getData(["../phy_run/*root*"])
    #data = getData(["../temp/*root*"])
    
    ch1=TChain("MC","")
    ch1.Add("sim/K40/*internal*root")
    #ch1.Add("/data/MC/KIMS-NaI/user-scratch/sim/processed/K40/set2/set2_internalK40-1-nEv200000.root")
    
    ch2=TChain("MC","")
    ch2.Add("sim/U238/*internal*root")
    #ch2.Add("/data/MC/KIMS-NaI/user-scratch/sim/processed/U238/set2/set2_internalU238-1-nEv100000.root")

    ch3=TChain("MC","")
    ch3.Add("sim/Th232/*internal*root")
    #ch3.Add("/data/MC/KIMS-NaI/user-scratch/sim/processed/Th232/set2/set2_internalTh232-1-nEv100000.root")
    
    ch4=TChain("MC","")
    ch4.Add("sim/K40/*pmt*root")
    #ch4.Add("/data/MC/KIMS-NaI/user-scratch/sim/processed/K40/set2/set2_pmtK40-1-nEv200000.root")
    
    par = histparam()
    canvs = []
    #legs = []

    
    sims4 = []
    for i in range(0,8): # is primary crystal of origin
    #for i in range(0,1): # is primary crystal of origin
        
        sims1 = []
        sims2 = []
        sims3 = []
        #sims4 = []
        for j in range(0,8): # is surounding crystals and self
            sim1 = TH1F("sim1", names(j), par[0], par[1], par[2])
            sim2 = TH1F("sim2", names(j), par[0], par[1], par[2])
            sim3 = TH1F("sim3", names(j), par[0], par[1], par[2])
            #cut1 = TCut('event_MC.edepResol['+str(j)+']*1000. > 0')
            #cut2 = TCut('event_MC.primVolumeName == "'+volumeNames(i)+'"')
            #ch1.Draw('event_MC.edepResol['+str(j)+']*1000. >> sim1',cut1+cut2)
            #ch2.Draw('event_MC.edepResol['+str(j)+']*1000. >> sim2',cut1+cut2)
            cut1 = TCut('MC.edepResol['+str(j)+']*1000. > 0')
            cut2 = TCut('MC.primVolumeName == "'+volumeNames(i)+'"')
            ch1.Draw('MC.edepResol['+str(j)+']*1000. >> sim1',cut1+cut2)
            ch2.Draw('MC.edepResol['+str(j)+']*1000. >> sim2',cut1+cut2)
            ch3.Draw('MC.edepResol['+str(j)+']*1000. >> sim3',cut1+cut2)
            sim1.Sumw2()
            sims1.append(sim1)
            sim2.Sumw2()
            sims2.append(sim2)
            sim3.Sumw2()
            sims3.append(sim3)

            # only need to generate external backgrounds one time
            if i == 0:
                sim4 = TH1F("sim4", names(j), par[0], par[1], par[2])
                cut1 = TCut('MC.edepResol['+str(j)+']*1000. > 0')
                cut2 = TCut('MC.primVolumeName == "phys_pmt" ')
                ch4.Draw('MC.edepResol['+str(j)+']*1000. >> sim4',cut1+cut2)
                sim4.Sumw2()
                sims4.append(sim4)

        
        canv = TCanvas(names(i), names(i), 0, 0, 1600, 900)
        canvs.append(canv)
        canvs[i].Divide(4,2)
        legs = []
        for k in range(0,8):
            
            canvs[i].cd(k+1)
            canvs[i].cd(k+1).SetLogy()
            
            leg = TLegend(0.60, 0.65, 0.92, 0.88)
            legs.append(leg)
            
            if k == i:
                sims1[k].SetMarkerColor(kRed+1)
                sims1[k].SetLineColor(kRed+1)
                sims2[k].SetMarkerColor(kBlue+1)
                sims2[k].SetLineColor(kBlue+1)
                sims3[k].SetMarkerColor(kGreen+1)
                sims3[k].SetLineColor(kGreen+1)
                
            else:
                sims1[k].SetMarkerColor(kRed+2)
                sims1[k].SetLineColor(kRed+2)
                sims2[k].SetMarkerColor(kBlue+2)
                sims2[k].SetLineColor(kBlue+2)
                sims3[k].SetMarkerColor(kGreen+2)
                sims3[k].SetLineColor(kGreen+2)

            sims4[k].SetMarkerColor(kYellow+2)
            sims4[k].SetLineColor(kYellow+2)
            
            sims1[k].SetLineWidth(1)
            sims2[k].SetLineWidth(1)
            sims3[k].SetLineWidth(1)
            sims4[k].SetLineWidth(1)
            
            data[k].GetYaxis().SetTitle('arb. counts')
            data[k].GetXaxis().SetTitle('Energy (keV)')
            
            data[k].Draw()
            sims1[k].Draw("same")
            sims2[k].Draw("same")
            sims3[k].Draw("same")
            sims4[k].Draw("same")
            
            legs[k].SetFillColor(0)
            legs[k].SetBorderSize(0)
            lopt='LPE'
            legs[k].AddEntry(data[k],  '  data',     lopt)
            legs[k].AddEntry(sims1[k], '  xtal K40',   lopt)
            legs[k].AddEntry(sims2[k], '  xtal U238',  lopt)
            legs[k].AddEntry(sims3[k], '  xtal Th232', lopt)
            legs[k].AddEntry(sims4[k], '  pmt K40', lopt)
            legs[k].Draw("same")


        canvs[i].Update()
        canvs[i].Print('testing-'+names(i)+'.png')

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


def calib(i):
    """
    Return the calibrations for the crystals
    """
    # Pushba 2016-09-27 slides
    calibrations = [0.0001093478,
                    0.000105676,
                    0.00011489,
                    0.0001127,
                    0.000267,
                    0.000117534,
                    0.0001119,
                    0.00032017]
    #return calibrations[int(i)]
    return 0.0004

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
    bins = (hmax-hmin)/2
    return [bins, hmin, hmax]


############################################################

if __name__ == "__main__":
    _myself_(sys.argv[1:])
