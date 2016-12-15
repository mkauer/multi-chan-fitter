#!/usr/bin/env python
######################################################################
# Import old funcs.py but now add more for reading in the
# universal backgrounds file.
# 
# version: 2016-12-15
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ tweaked buildMC2() to select either resol() or resol2()
# + added resol2() and Pushpa's resols
# + added getData22() to use new root branches for calib22()
# + added calib22() and Pushpa's calibs
# + added key data[key]['runtime'] to dataDRU2()
# + added scaleBkgs() to scale into real mBq/kg units
# ~ use backgrounds2.txt for v21
# ~ don't use bkgs scaling in readROOT2() for v21
# ~ run 1544 has 2 hour subruns
# ~ use data run 1544
# ~ fixed dataDRU2() for the new data key format
# + added sortKeys2()
# + added dataDRU2()
# ~ build MC now returns bkgs and sigs
# ~ still need deepcopy to get hists from rootfile
# ~ change the hist naming format for data and MC
# + adding new functions to build data and mc from the
#   backgrounds file
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
from funcs import *


def buildMC2(fileName='backgrounds.txt', version=1):
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
                    if version == 2:
                        ### assume reso = p0/sqrt(energy) + p1
                        p0, p1 = resol2(i,E)
                        chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.) + '+str(p1))
                    else:
                        ### assume reso = p0/sqrt(energy)
                        p0 = resol(i,E)
                        chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.)')

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


def getData2(runNum=1544):
    """
    Build histograms of the data from all crystals
    """
    
    local = amLocal()
    
    data = {}
    
    chain = TChain("ntp","")
    nfiles=0
    if local:
        nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/data/phys/ntp_I*'+str(runNum)+'*root*')
    else:
        #nfiles = chain.Add('/home/mkauer/temp/*'+str(runNum)+'*root*')
        ### testing...
        #nfiles = chain.Add('/data/COSINE/NTP/phys/V00-00-02/ntp_I001544.root.001')
        nfiles = chain.Add('/data/COSINE/NTP/phys/V00-00-02/ntp_I*'+str(runNum)+'*root*')
        
    if nfiles > 0:
        print nfiles,'data files found'
        for E in range(2):

            ### use pmtXX.rqcD1_5 for high energy
            ### use pmtXX.qc5 for low energy
            if E:
                edep = 'rqcD1_5'
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
                chain.Draw('(pmt'+str(i+1)+'1.'+str(edep)+'*'+str(calib(i,0,E))
                           +' + pmt'+str(i+1)+'2.'+str(edep)+'*'+str(calib(i,1,E))+')/2.'
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
        print 'Warning: no data files found...'
    
    return data


def readROOT2(rootfile, fileName='backgrounds.txt', scale=0):
    
    if not os.path.exists(rootfile):
        print 'rootfile [',rootfile,'] file not found'
        sys.exit()
    
    if not os.path.exists(fileName):
        print 'backgrounds [',fileName,'] file not found'
        sys.exit()

    #----------------------------------------
    # still need to deepcopy the histos!!!
    from copy import deepcopy
    #----------------------------------------
    
    rfile = TFile(rootfile, "READ")
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
    with open(fileName) as infile:
        for line in infile:
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue
            #print line
            bits = line.split()
            bits = [x.strip() for x in bits]
            bsdr = str(bits[0])
            
            if 'D' in bsdr:
                xstl = int(bits[1])
                loca = str(bits[2])
                
                for e in range(2):
                    name = 'x'+str(xstl)+'-'+loca+'-e'+str(e)
                    data[name] = {}
                    data[name]['hist'] = deepcopy(TH1F(rfile.Get(name)))
                    data[name]['subruns'] = deepcopy(TH1F(rfile.Get(name+'_subruns')))
            
            if 'B' in bsdr:
                xstl = int(bits[1])
                loca = str(bits[2])
                isot = str(bits[3])
                acti = float(bits[4])
                erro = float(bits[5])
                
                for e in range(2):
                    name = 'x'+str(xstl)+'-'+loca+'-'+isot+'-e'+str(e)
                    bkgs[name] = {}
                    bkgs[name]['hist'] = deepcopy(TH1F(rfile.Get(name)))
                    bkgs[name]['generated'] = deepcopy(TH1F(rfile.Get(name+'_generated')))
                    bkgs[name]['acti'] = acti
                    bkgs[name]['erro'] = erro

                    ### scale the bkgs right away?
                    # I think scaling should be a seperate function
                    # Keep this in for v20 comparison
                    if scale:
                        bkgs[name]['hist'].Scale(acti)
                    
                    
            if 'S' in bsdr:
                xstl = int(bits[1])
                loca = str(bits[2])
                isot = str(bits[3])
                fbnd = [float(bits[6]), float(bits[7])]
                
                for e in range(2):
                    name = 'x'+str(xstl)+'-'+loca+'-'+isot+'-e'+str(e)
                    sigs[name] = {}
                    sigs[name]['hist'] = deepcopy(TH1F(rfile.Get(name)))
                    sigs[name]['generated'] = deepcopy(TH1F(rfile.Get(name+'_generated')))
                    sigs[name]['fbnd'] = fbnd
    
    rfile.Close()
    infile.close()
    
    return data, bkgs, sigs


def dataDRU2(data):
    """
    Scale data into DRU
    """
        
    # this will have to get fixed eventually if we want finer binning in the lo-E region
    # or more coarse binning in the hi-E
    keVperBin = float(1)
    
    # run 1324 has 1 hour subruns
    # run 1544 has 2 hour subruns
    subrunTime = float(2) # in hours
    
    #scales = []
    for key in data:
    #for i in range(8):

        i = int(key.split('-')[0].split('x')[-1]) - 1
        nfiles = data[key]['subruns'].GetEntries()
        days = float((subrunTime*nfiles)/24.)
        scale = float(1./(days*cmass(i)*keVperBin))

        data[key]['hist'].Scale(scale)

        # runtime is in sec
        data[key]['runtime'] = nfiles*subrunTime*60*60

        data[key]['druScale'] = scale
        
    return data


def scaleBkgs(bkgs, data):
    
    ### only works for internal backgrounds at this point
    
    for name in bkgs:
        x = name.split('-')[0]
        loca = name.split('-')[1]
        e = name.split('-')[-1]
        druscale = data[x+'-data-'+e]['druScale']
        print name,'druscale = ',druscale
        runtime = data[x+'-data-'+e]['runtime']
        #print 'runtime = ',runtime

        if loca == 'internal':
            mbqkg = bkgs[name]['acti']
            kgs = cmass(int(x[-1])-1)
            mbq = mbqkg*kgs
            bq = mbq / 1000.0
            print name,'bq = ',bq
            
            generated = float(bkgs[name]['generated'].GetEntries())
            #print 'generated = ',generated
            
            detected = float(bkgs[name]['hist'].GetEntries())
            #print 'detected = ',detected
            
            eff = detected / generated
            #print 'eff = ',eff
            
            # number of events needed from MC
            events = bq * runtime * eff
            #print 'events = ',events
            
            ### my first thought
            #scale = events / generated
            ### i think this should be normalized to the number of detected events
            ### because the eff scaling is already included...
            #scale = events / detected
            ### which then becomes
            #scale = bq * runtime * (1./generated)
            
            ### not making sence, just try some things
            #scale = (events / generated) * druscale
            
            ### 2016-12-15
            ### this looks about right compared Estella's plots
            ### slightly lower because of resol smearing
            scale = (events / detected) * druscale
            ### same as doing
            #scale = bq * runtime * (1./generated) * druscale
            print name,'scale = ',scale
            
            bkgs[name]['hist'].Scale(scale)
            
        if loca == 'pmt':
            mbqpmt = bkgs[name]['acti']
            # not sure what to do for pmts
            pmts = 8*2
            mbq = mbqpmt*pmts
            bq = mbq / 1000.0
            print name,'bq = ',bq
            
            generated = float(bkgs[name]['generated'].GetEntries())
            #print 'generated = ',generated
            
            detected = float(bkgs[name]['hist'].GetEntries())
            #print 'detected = ',detected
            
            eff = detected / generated
            #print 'eff = ',eff
            
            # number of events needed from MC
            events = bq * runtime * eff
            #print 'events = ',events
            
            ### my first thought
            #scale = events / generated
            
            ### i think this should be normalized to the number of detected events
            ### because the eff scaling is already included...
            #scale = events / detected
            ### which then becomes
            #scale = bq * runtime * (1./generated)
            
            ### not making sence, just try some things
            #scale = (events / generated) * druscale
            scale = (events / detected) * druscale
            print name,'scale = ',scale
            
            bkgs[name]['hist'].Scale(scale)
            
        #sys.exit()
        
    return bkgs


def sortKeys2(data, bkgs, sigs):
    datkeys=[]
    bakkeys=[]
    sigkeys=[]

    for key in data:
        datkeys.append(key)
    datkeys.sort()

    for key in bkgs:
        bakkeys.append(key)
    bakkeys.sort()

    for key in sigs:
        sigkeys.append(key)
    sigkeys.sort()

    return datkeys, bakkeys, sigkeys


def calib22(i, E=0):
    """
    Return the crystal calibrations
    """
    # from Pushpa
    # (i) is the crystal 
    hiEcalib = [
        1.000,
	0.992,
	1.000,
	1.000,
	0.966,
	1.014,
	0.984,
	1.000
        ]
    loEcalib = [
        1.06000e-4,
	1.05676e-4,
	1.15900e-4,
	1.12700e-4,
	2.86907e-4,
	1.17634e-4,
	1.11900e-4,
	3.84263e-4
        ]
    if E:
        return hiEcalib[int(i)]
    else:
        return loEcalib[int(i)]


def getData22(runNum=1544):
    """
    Build histograms of the data from all crystals
    """
    
    local = amLocal()
    
    data = {}
    
    chain = TChain("ntp","")
    nfiles=0
    if local:
        nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/data/phys/ntp_I*'+str(runNum)+'*root*')
    else:
        #nfiles = chain.Add('/home/mkauer/temp/*'+str(runNum)+'*root*')
        ### testing...
        #nfiles = chain.Add('/data/COSINE/NTP/phys/V00-00-02/ntp_I001544.root.001')
        nfiles = chain.Add('/data/COSINE/NTP/phys/V00-00-02/ntp_I*'+str(runNum)+'*root*')
        
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
        print 'Warning: no data files found...'
    
    return data


def resol2(i, E=0):
    """
    Return the crystal resolutions
    """
    # from Pushpa
    # (i) is the crystal
    # res = p[0]/sqrt(x) + p[1]
    hiEresol = [
        [0.7127, 0.004879],
	[0.7127, 0.004879],
	[0.7127, 0.004879],
	[0.7127, 0.004879],
	[0.7127, 0.004879],
	[0.7127, 0.004879],
	[0.7127, 0.004879],
	[0.7127, 0.004879]
    ]
    
    # C5 = C8 == C4
    # C6 = C7 == C1
    loEresol = [
        [0.2349,  0.018600],
	[0.2699,  0.014920],
	[0.2459,  0.015020],
	[0.3797, -0.003919],
	[0.3797, -0.003919],
	[0.2349,  0.018600],
	[0.2349,  0.018600],
	[0.3797, -0.003919]
    ]
    
    if E:
        return hiEresol[int(i)]
    else:
        return loEresol[int(i)]

