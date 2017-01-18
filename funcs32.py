#!/usr/bin/env python
######################################################################
# Import funcs3.py but now adding support for "newGeometry" simulation
#
# Works with v32 and later versions
#
# version: 2017-01-17
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + readROOT32() to read in rootfiles dynamically 
# ~ run 1546 data is under /data/COSINE/NTP/phys/V00-02-00
# + funcs32.py for newGeometry simulation and run 1546
# 
# email me: mkauer@physics.wisc.edu
######################################################################
# 
# Where is MC?
# /data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/
# 
# Where is the processed data?
# Right now I'm using:
# /data/COSINE/NTP/phys/V00-02-00
# and run ntp_I001546*
# 
######################################################################

import os,sys
import numpy as np
from ROOT import *
import ROOT

from funcs3 import *


def buildMC32(fileName='backgrounds3.txt', calib=2):
    """
    Build the mega MC dictionary!
    calib=1 is Estella's calib and resol
    calib=2 is Pushpa's calib and resol
    """
    
    local = amLocal()
    
    bkgs = {}
    sigs = {}
    
    for line in readFile(fileName):
        bits = line.split()
        #bits = [x.strip() for x in bits]
        bsdr = str(bits[0])
        
        ### skip the data lines
        if 'D' in bsdr: continue

        ### crystal number
        xstl = int(bits[1])
        ### background location
        loca = str(bits[2])
        ### top level isotope file name
        isof = str(bits[3])
        ### isotope chain break start
        chst = str(bits[4])
        ### isotope chain break stop
        chsp = str(bits[5])
        ### activity
        acti = float(bits[6])
        ### error on activity
        erro = float(bits[7])
        ### fit bounds
        fbnd = [float(bits[8]), float(bits[9])]

        chain = TChain("MC","")
        nfiles=0
        if local:
            nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/sim/newGeometry/'+isof+'/set2/'+'*'+loca+'*root')
            
        else:
            nfiles = chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'+isof+'/set2/'+'*'+loca+'*root')

        if nfiles > 0:
            print nfiles,'files found for',loca,isof
            for E in range(2):

                i = xstl-1

                ### see if i can use unique names for all the histograms
                key  = 'x'+str(xstl)
                key += '-'+loca
                if chst == chsp: key += '-'+chst
                else: key += '-'+chst+'_'+chsp
                key += '-e'+str(E)

                print key

                par = histparam(E)
                histo = TH1F(key, longNames(i), par[0], par[1], par[2])

                cut1 = TCut('edep['+str(i)+']*1000. > 0')

                if loca == 'internal':
                    cut2 = TCut('primVolumeName == "'+volumeNames(i)+'"')
                if loca == 'pmt':
                    cut2 = TCut('primVolumeName == "phys_pmt"')

                cut3 = 0
                if chst == 'U238' and chsp == 'Rn222':
                    cut3 = TCut('(groupNo >= 11) && (groupNo <= 14)')
                if chst == 'Pb210' and chsp == 'Pb210':
                    cut3 = TCut('groupNo == 15')
                
                key2 = key+'_generated'
                temp = TH1F(key2, 'generated',1,0,1)

                ### v31 - maybe need cut3 for generated events
                ### to get the eff/normalization right?
                #chain.Draw('primVolumeName >> '+key2, cut2)
                if cut3:
                    chain.Draw('primVolumeName >> '+key2, cut2+cut3)
                else:
                    chain.Draw('primVolumeName >> '+key2, cut2)
                    
                generated = temp.GetEntries()

                ### test out resolution smearing
                ###-------------------------------------------------------------------------------
                ### using Box-Muller? method here for "rng" alias
                chain.SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')

                ### set the resolution function
                if calib == 1:
                    ### from Estella
                    ### assume reso = p0/sqrt(energy)
                    p0 = resol(i,E)
                    chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.)')
                if calib == 2:
                    ### from Pushpa
                    ### assume reso = p0/sqrt(energy) + p1
                    p0, p1 = resol2(i,E)
                    chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.) + '+str(p1))

                if cut3:
                    chain.Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> '+key, cut1+cut2+cut3)
                else:
                    chain.Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> '+key, cut1+cut2)

                detected = histo.GetEntries()
                ###-------------------------------------------------------------------------------
                
                histo.SetLineColor(kBlack)
                histo.SetMarkerColor(kBlack)
                histo.SetLineWidth(1)
                
                #print 'Eff =',detected/generated
                
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

        else: print 'WARNING:', loca, isof, 'not found...'
        
    return bkgs, sigs


def getData32(runNum=1546, build='V00-02-00'):
    """
    Build histograms of the data from all crystals
    """
    
    local = amLocal()
    
    data = {}
    
    chain = TChain("ntp","")
    nfiles=0
    if local:
        nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/data/phys/'+build+'/ntp_I*'+str(runNum)+'*root*')
    else:
        nfiles = chain.Add('/data/COSINE/NTP/phys/'+build+'/ntp_I*'+str(runNum)+'*root*')
        
    if nfiles > 0:
        print nfiles,'data files found'
        for E in range(2):

            ### use crystalX.energyD for high energy
            ### use crystalX.qc5 for low energy
            if E:
                edep = 'energyD'
            else:
                edep = 'qc5'
            
            par = histparam(E)
            for i in range(8):
                key  = 'x'+str(i+1)
                key += '-data'
                key += '-e'+str(E)
                
                data[key] = {}
                
                #histo = TH1F(key, key, par[0], par[1], par[2])
                histo = TH1F(key, longNames(i), par[0], par[1], par[2])
                #chain.Draw('(pmt'+str(i+1)+'1.'+str(edep)+'*'+str(calib(i,0,E))
                #           +' + pmt'+str(i+1)+'2.'+str(edep)+'*'+str(calib(i,1,E))+')/2.'
                #           +' >> '+str(key))
                chain.Draw('crystal'+str(i+1)+'.'+str(edep)+'*'+str(calib22(i,E))
                           +' >> '+str(key))
                
                histo.Sumw2()
                histo.SetLineColor(kBlack)
                histo.SetMarkerColor(kBlack)
                histo.SetLineWidth(1)
                
                data[key]['hist'] = histo
                
                key2 = key+'_subruns'
                temp = TH1F(key2,'subruns',1,0,1)
                #temp.SetBinContent(0, nfiles)
                for i in range(nfiles):
                    temp.Fill(0.5)
                data[key]['subruns'] = temp
                
    else:
        print 'ERROR: no data files found... quitting...'
        sys.exit()
    
    return data


def readROOT32(fileName='backgrounds3.txt'):
    
    #if not os.path.exists(rootfile):
    #    print 'rootfile [',rootfile,'] file not found'
    #    sys.exit()
    
    if not os.path.exists(fileName):
        print 'backgrounds [',fileName,'] file not found'
        sys.exit()

    #----------------------------------------
    # still need to deepcopy the histos!!!
    from copy import deepcopy
    #----------------------------------------
    
    #rfile = TFile(rootfile, "READ")
    # switch to different dir in rootfile
    #rfile.cd('')
    # to print out all the keys
    #rfile.ls()
    
    #rnames = [key.GetName() for key in gDirectory.GetListOfKeys()]
    #rnames.sort()
    #print rnames
    
    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(fileName):
        bits = line.split()
        #bits = [x.strip() for x in bits]
        
        rootfile = str(bits[-1])
        if not os.path.exists(rootfile):
            print 'rootfile [',rootfile,'] file not found'
            sys.exit()
        rfile = TFile(rootfile, "READ")

        
        bsdr = str(bits[0])

        if 'D' in bsdr:
            xstl = int(bits[1])
            loca = str(bits[2])

            for e in range(2):
                key = 'x'+str(xstl)+'-'+loca+'-e'+str(e)
                data[key] = {}
                data[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                data[key]['subruns'] = deepcopy(TH1F(rfile.Get(key+'_subruns')))

        if 'B' in bsdr:
            xstl = int(bits[1])
            loca = str(bits[2])
            isof = str(bits[3])
            chst = str(bits[4])
            chsp = str(bits[5])
            acti = float(bits[6])
            erro = float(bits[7])

            for E in range(2):
                key  = 'x'+str(xstl)
                key += '-'+loca
                if chst == chsp: key += '-'+chst
                else: key += '-'+chst+'_'+chsp
                key += '-e'+str(E)

                bkgs[key] = {}
                bkgs[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                bkgs[key]['generated'] = deepcopy(TH1F(rfile.Get(key+'_generated')))
                bkgs[key]['acti'] = acti
                bkgs[key]['erro'] = erro


        if 'S' in bsdr:
            xstl = int(bits[1])
            loca = str(bits[2])
            isof = str(bits[3])
            chst = str(bits[4])
            chsp = str(bits[5])
            fbnd = [float(bits[8]), float(bits[9])]

            for E in range(2):
                key  = 'x'+str(xstl)
                key += '-'+loca
                if chst == chsp: key += '-'+chst
                else: key += '-'+chst+'_'+chsp
                key += '-e'+str(E)

                sigs[key] = {}
                sigs[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                sigs[key]['generated'] = deepcopy(TH1F(rfile.Get(key+'_generated')))
                sigs[key]['fbnd'] = fbnd

        rfile.Close()
    
    return data, bkgs, sigs

