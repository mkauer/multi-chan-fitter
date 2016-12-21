#!/usr/bin/env python
######################################################################
# Import funcs2.py but now adding support for broken chains
#
# Works with v30 and later versions
#
# version: 2016-12-20
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + v3 for broken chain support
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
# /data/COSINE/PHYS
# But just use the processed data...
# 
# Where is the processed data?
# Right now it's in:
# /data/COSINE/NTP/phys/V00-00-02
# Currently I'm using run 1544
# 
######################################################################

import os,sys
import numpy as np
from ROOT import *
import ROOT
from funcs2 import *


def buildMC3(fileName='backgrounds.txt', version=1):
    """
    Build the mega MC dictionary!
    """
    
    local = amLocal()
    
    bkgs = {}
    sigs = {}

    with open(fileName) as infile:
        for line in infile:
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue
            bits = line.split()
            bits = [x.strip() for x in bits]
            bsdr = str(bits[0])
            # skip the data lines
            if 'D' in bsdr: continue
            # save everything else
            xstl = int(bits[1])
            loca = str(bits[2])
            isot = str(bits[3])
            acti = float(bits[4])
            erro = float(bits[5])
            fbnd = [float(bits[6]), float(bits[7])]
            
            chain = TChain("MC","")
            nfiles=0
            if local:
                nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/sim/'+isot+'/'+'*'+loca+'*root')
            else:
                #nfiles = chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/'+isot+'/set2/'+'set2_*'+loca+'*-1-*root')
                nfiles = chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/'+isot+'/set2/'+'set2_*'+loca+'*root')
            
            if nfiles > 0:
                print nfiles,'files found for',loca,isot
                for E in range(2):
                    
                    i = xstl-1
                    
                    ### see if i can use unique names for all the histograms
                    key  = 'x'+str(xstl)
                    key += '-'+loca
                    key += '-'+isot
                    key += '-e'+str(E)
                    
                                      
                    par = histparam(E)
                    #histo = TH1F(key, key, par[0], par[1], par[2])
                    #histo = TH1F(key, longNames(i), par[0], par[1], par[2])
                    histo = TH1F(key, longNames(i), par[0], par[1], par[2])
                    
                    cut1 = TCut('edep['+str(i)+']*1000. > 0')
                    if loca == 'internal':
                        cut2 = TCut('primVolumeName == "'+volumeNames(i)+'"')
                    if loca == 'pmt':
                        cut2 = TCut('primVolumeName == "phys_pmt" ')
                    
                    key2 = key+'_generated'
                    temp = TH1F(key2, 'generated',1,0,1)
                    chain.Draw('primVolumeName >> '+key2, cut2)
                    generated = temp.GetEntries()
                    
                    ### test out resolution smearing
                    ###-------------------------------------------------------------------------------
                    ### using Box-Muller? method here for "rng" alias
                    #mc[loc][iso]["chain"].SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                    chain.SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                    
                    ### set the resolution function
                    if version == 1:
                        ### from Estella
                        ### assume reso = p0/sqrt(energy)
                        p0 = resol(i,E)
                        chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.)')
                    if version == 2:
                        ### from Pushpa
                        ### assume reso = p0/sqrt(energy) + p1
                        p0, p1 = resol2(i,E)
                        chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.) + '+str(p1))
                    
                    ### then draw the chain
                    #mc[loc][iso]["chain"].Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> histo', cut1+cut2)
                    #mc[loc][iso]["chain"].Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> '+key, cut1+cut2)
                    chain.Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> '+key, cut1+cut2)
                    detected = histo.GetEntries()
                    ###-------------------------------------------------------------------------------
                    
                    histo.SetLineColor(kBlack)
                    histo.SetMarkerColor(kBlack)
                    histo.SetLineWidth(1)
                    
                    print 'Eff =',detected/generated
                    
                    if 'B' in bsdr:
                        bkgs[key] = {}
                        bkgs[key]['hist'] = histo
                        bkgs[key]['generated'] = temp
                        bkgs[key]['acti'] = acti
                        bkgs[key]['erro'] = erro
                                            
                    if 'S' in bsdr:
                        sigs[key] = {}
                        sigs[key]['hist'] = histo
                        sigs[key]['generated'] = temp
                        sigs[key]['fbnd'] = fbnd
                        
                    
            else:
                print 'Warning:', loc, iso, 'not found...'
    
    infile.close()
    
    return bkgs, sigs

