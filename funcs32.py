#!/usr/bin/env python
######################################################################
# Import funcs3.py but now adding support for "newGeometry" simulation
#
# Works with v32 and later versions
#
# version: 2017-01-24
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ cleaned up scaleSigs32() and scaleBkgs32() 
# + add scaleSigs32() and doesn't need data input anymore
# + add scaleBkgs32() and doesn't need data input anymore
# + TCut on primParticleName! This was critial to get generated
#   event count correct.
# ~ throw a warning in readROOT32() if the histo isn't found
# + add more selection criteria for the extra bkgs in buildMC32()
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

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs3 import *


def buildMC32(fileName='backgrounds32.txt', calib=2):
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
                elif loca == 'pmt':
                    cut2 = TCut('primVolumeName == "phys_pmt"')
                elif loca == 'lsveto':
                    cut2 = TCut('primVolumeName == "lsveto"')
                elif loca == 'lsvetoair':
                    cut2 = TCut('(primVolumeName == "DetPMTCover") || (primVolumeName == "DetPMTEnvelope")')
                elif loca == 'airshield':
                    cut2 = TCut('primVolumeName == "LSVetoAirRoom"')
                elif loca == 'steel':
                    # not sure this works the way I think it should?
                    cut2 = TCut('primVolumeName != ""')
                else:
                    print "WARNING: No selection criteria for  --> ", loca
                    continue
                
                ### this is needed to get the generated event numbers right!
                cut3 = TCut('primParticleName == "'+chst+'"')
                
                ### this is needed? to cut out crap events
                ### (event_info.Type) is new
                ### (evt_Type) is old
                #cut4 = TCut('(event_info.Type > 10) || (evt_Type > 10)')
                cut4 = TCut('event_info.Type > 10')
                
                cut100 = 0
                if chst == 'U238' and chsp == 'Rn222':
                    cut100 = TCut('(groupNo >= 11) && (groupNo <= 14)')
                if chst == 'Pb210' and chsp == 'Pb210':
                    cut100 = TCut('groupNo == 15')
                
                key2 = key+'_generated'
                temp = TH1F(key2, 'generated',1,0,1)
                
                ### v31 - maybe need cut100 for generated events
                ### to get the eff/normalization right?
                #chain.Draw('primVolumeName >> '+key2, cut2)
                if cut100:
                    #chain.Draw('primVolumeName >> '+key2, cut2+cut3+cut100)
                    chain.Draw('primVolumeName >> '+key2, cut2+cut3)
                else:
                    chain.Draw('primVolumeName >> '+key2, cut2+cut3)
                
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

                if cut100:
                    chain.Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> '+key, cut1+cut2+cut4+cut100)
                else:
                    chain.Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> '+key, cut1+cut2+cut4)

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
        print 'ERROR: No data files found... quitting...'
        sys.exit()
    
    return data


def readROOT32(fileName='backgrounds32.txt'):
    
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

                try:
                    data[key] = {}
                    data[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                    data[key]['subruns'] = deepcopy(TH1F(rfile.Get(key+'_subruns')))
                except:
                    del data[key]
                    print "WARNING: could not find hist -->",key
                    continue
                
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

                try:
                    bkgs[key] = {}
                    bkgs[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                    bkgs[key]['generated'] = deepcopy(TH1F(rfile.Get(key+'_generated')))
                    bkgs[key]['acti'] = acti
                    bkgs[key]['erro'] = erro
                except:
                    del bkgs[key]
                    print "WARNING: could not find histo -->",key
                    continue
                

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
                
                try:
                    sigs[key] = {}
                    sigs[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                    sigs[key]['generated'] = deepcopy(TH1F(rfile.Get(key+'_generated')))
                    sigs[key]['fbnd'] = fbnd
                except:
                    del sigs[key]
                    print "WARNING: could not find histo -->",key
                    continue
        
        rfile.Close()
    
    return data, bkgs, sigs


def scaleBkgs32(bkgs):
    
    for name in bkgs:
        x = name.split('-')[0]
        loca = name.split('-')[1]
        e = name.split('-')[-1]
        
        #druscale = data[x+'-data-'+e]['druScale']
        #runtime  = data[x+'-data-'+e]['runtime']
        
        kev  = 1.     # keV/bin
        day  = 86400. # in seconds
        
        xkgs = cmass(int(x[-1])-1)
        pmts = 16.
        lskg = 1800.
        
        generated = float(bkgs[name]['generated'].GetEntries())
        if generated < 1:
            print "WARNING: 0 events generated for -->", name
            continue
        
        if loca == 'internal':
            #scale = bkgs[name]['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[name]['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[name]['hist'].Scale(scale)
        elif loca == 'pmt':
            #scale = bkgs[name]['acti'] * (pmts) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[name]['acti'] * (pmts) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[name]['hist'].Scale(scale)
        elif loca == 'lsveto':
            #scale = bkgs[name]['acti'] * (lskg) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[name]['acti'] * (lskg) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[name]['hist'].Scale(scale)
        elif loca == 'lsvetoair':
            #scale = bkgs[name]['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[name]['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[name]['hist'].Scale(scale)
        elif loca == 'airshield':
            #scale = bkgs[name]['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[name]['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[name]['hist'].Scale(scale)
        elif loca == 'steel':
            #scale = bkgs[name]['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[name]['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[name]['hist'].Scale(scale)
        else:
            print "WARNING: No background scaling for  --> ", loca
            continue


    return bkgs


def scaleSigs32(sigkeys, sigs):

    for name in sigkeys:
        x = name.split('-')[0]
        loca = name.split('-')[1]
        e = name.split('-')[-1]
        
        #druscale = data[x+'-data-'+e]['druScale']
        #runtime = data[x+'-data-'+e]['runtime']

        kev  = 1.     # keV/bin
        day  = 86400. # in seconds
        
        xkgs = cmass(int(x[-1])-1)
        pmts = 16.
        lskg = 1800.
        
        generated = float(sigs[name]['generated'].GetEntries())
        if generated < 1:
            print "WARNING: 0 events generated for -->", name
            continue
        
        # verbose?
        V = 1
        E = int(e[-1])
        
        if loca == 'internal':
            #fitActivity = sigs[name]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[name]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[name]['acti'] = fitActivity
            if V and E: print '!!!!', name, sigs[name]['acti'],'mBq/kg'
        if loca == 'pmt':
            #fitActivity = sigs[name]['fitscale'] * (1./pmts) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[name]['fitscale'] * (1./pmts) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[name]['acti'] = fitActivity
            if V and E: print '!!!!', name, sigs[name]['acti'],'mBq/pmt'
        if loca == 'lsveto':
            #fitActivity = sigs[name]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[name]['fitscale'] * (1./lskg) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[name]['acti'] = fitActivity
            if V and E: print '!!!!', name, sigs[name]['acti'],'mBq/pmt'
        if loca == 'lsvetoair':
            #fitActivity = sigs[name]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[name]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[name]['acti'] = fitActivity
            if V and E: print '!!!!', name, sigs[name]['acti'],'mBq/pmt'
        if loca == 'airshield':
            #fitActivity = sigs[name]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[name]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[name]['acti'] = fitActivity
            if V and E: print '!!!!', name, sigs[name]['acti'],'mBq/pmt'
        if loca == 'steel':
            #fitActivity = sigs[name]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[name]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[name]['acti'] = fitActivity
            if V and E: print '!!!!', name, sigs[name]['acti'],'mBq/pmt'


    return sigs

