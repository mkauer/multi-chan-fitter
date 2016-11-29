#!/usr/bin/env python
######################################################################
# Bunch of functions to use
# Works with "10-mc-fit.py" and later versions, local and CUP
# Trying to make sure this stays backwards compitable
#
# version: 2016-11-21
#
# Change Log (key == [+] added, [-] removed, [~] changed)
#=====================================================================

#---------------------------------------------------------------------
# should hists be filled and saved as key,key or key,longName?
#---------------------------------------------------------------------

# + added return of scalings to dataDRU()
# + added dataDRU() to scale data into dru units
# + added cmass() to get crystal mass in kg
# + added function cnames() to give C1,C2,etc.
# ~ combined resolution should be divided by 2!
#   see the toy-MC ipython notebook
# ~ changed the histogram naming format
# ~ changed a lot of functions to use energy parameter
# + need copy or deepcopy for reading histos from rootfile
# + make unique keys for ALL the histograms
# ~ makeTotal and makeResid now use energy param to get hist values
# + buildMC() to put together all the MC
# ~ changed hiE and loE data hist names - prep for the root joiner
# + huh, importing personal python moduales now works I guess?
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
# /data/KIMS/COSINE/PHY_RUN
#
# My processed data is currently in
# /home/mkauer/temp
#
######################################################################

import os,sys
import numpy as np
from ROOT import *
import ROOT


def amLocal():
    """
    Check to see what machine you're on so you can get right paths
    """
    import socket
    ### master.cunpa.ibs
    local = 1
    host = str(socket.gethostname())
    #print 'hostname =',host
    if 'cunpa' in host:
        local = 0
    return local


def buildMC(locs=['internal'], isos=['K40'], energy=0):
    """
    Build the mega MC dictionary!
    """
    
    local = amLocal()
    
    locs.sort()
    isos.sort()
    
    mc = {}
    for loc in locs:
        mc[loc] = {}
        for iso in isos:
            mc[loc][iso] = {}
            chain = TChain("MC","")
            if local:
                nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/sim/'+iso+'/'+'*'+loc+'*root')
            else:
                #nfiles = chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/'+iso+'/set2/'+'set2_*'+loc+'*-1-*root')
                nfiles = chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/'+iso+'/set2/'+'set2_*'+loc+'*root')

            if nfiles > 0:
                print nfiles,'files found for',loc,iso
                mc[loc][iso]["chain"] = chain
                mc[loc][iso]["hist"] = []
            else:
                print 'Warning:', loc, iso, 'not found...'
                mc[loc][iso]["chain"] = None
                mc[loc][iso]["hist"] = None
    
    # get histogram parameters for MC
    par = histparam(energy)
    
    for i in range(8): # is primary crystal of origin
        for loc in locs:
            for iso in isos:
                
                if not mc[loc][iso]["chain"]: continue
                
                ### see if i can use unique names for all the histograms
                key  = 'e'+str(energy)
                key += '_'+loc
                key += '_'+iso
                key += '_x'+str(i)
                #histo = TH1F(key, key, par[0], par[1], par[2])
                histo = TH1F(key, longNames(i), par[0], par[1], par[2])
                
                cut1 = TCut('edep['+str(i)+']*1000. > 0')
                if loc == 'internal':
                    cut2 = TCut('primVolumeName == "'+volumeNames(i)+'"')
                if loc == 'pmt':
                    cut2 = TCut('primVolumeName == "phys_pmt" ')
                
                ### test out resolution smearing
                ###-------------------------------------------------------------------------------
                ### using Box-Muller? method here for "rng" alias
                mc[loc][iso]["chain"].SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                
                ### set the resolution function - right now it's linear with intercept of 0
                mc[loc][iso]["chain"].SetAlias('sigma', str(resol(i,energy))+' / sqrt(edep['+str(i)+']*1000.)')
                
                ### then draw the shit
                #mc[loc][iso]["chain"].Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> histo', cut1+cut2)
                mc[loc][iso]["chain"].Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> '+key, cut1+cut2)
                ###-------------------------------------------------------------------------------

                histo.SetLineColor(kBlack)
                histo.SetMarkerColor(kBlack)
                histo.SetLineWidth(1)
                mc[loc][iso]["hist"].append(histo)
                
    return mc


def getData(rootfiles=['./data/phys/*root*'], energy=0):
    """
    Build histograms of the data from all crystals
    """
    ### use pmtXX.rqcD1_5 for high energy
    ### use pmtXX.qc5 for low energy
    if energy:
        edep = 'rqcD1_5'
    else:
        edep = 'qc5'
    chain = TChain("ntp","")
    nfiles=0
    for rfile in rootfiles:
        nfiles += chain.Add(rfile)
    print nfiles,'data files found'
    data = []
    par = histparam(energy)
    for i in range(8):
        key  = 'e'+str(energy)
        key += '_data'
        key += '_x'+str(i)
        #histo = TH1F(key, key, par[0], par[1], par[2])
        histo = TH1F(key, longNames(i), par[0], par[1], par[2])
        chain.Draw('(pmt'+str(i+1)+'1.'+str(edep)+'*'+str(calib(i,0,energy))
                   +' + pmt'+str(i+1)+'2.'+str(edep)+'*'+str(calib(i,1,energy))+')/2.'
                   +' >> '+str(key))
        histo.Sumw2()
        histo.SetLineColor(kBlack)
        histo.SetMarkerColor(kBlack)
        histo.SetLineWidth(1)
        data.append(histo)
    return data


def dataDRU(data, nfiles=78):
    """
    Scale data into DRU
    """
        
    nfiles = float(nfiles)
    # FIX - will eventually need to count the runs or add the run times 
    
    # this will have to get fixed eventually if we want finer binning in the lo-E region
    # or more coarse binning in the hi-E
    keVperBin = float(1)
    
    # each subrun-file is 2 hours long
    # estella says the 1324 runs are 1 hour subruns...
    subrun = float(1) # in hours
    days = float((subrun*nfiles)/24.)

    scales = []
    for i in range(8):
        scale = float(1./(days*cmass(i)*keVperBin))
        data[i].Scale(scale)
        scales.append(scale)
        
    return data, scales


def makeTotal(energy):
    total = []
    par = histparam(energy)
    for i in range(8):
        key = 'e'+str(energy)
        key += '_total_x'+str(i)
        tot = TH1F(key, longNames(i), par[0], par[1], par[2])
        tot.SetLineColor(kGray+1)
        tot.SetMarkerColor(kGray+1)
        tot.SetLineWidth(1)
        total.append(tot)
    return total


def makeResid(energy):
    resid = []
    par = histparam(energy)
    for i in range(8):
        key = 'e'+str(energy)
        key += '_resid_x'+str(i)
        res = TH1F(key, longNames(i), par[0], par[1], par[2])
        res.SetLineColor(kBlack)
        res.SetMarkerColor(kBlack)
        res.SetLineWidth(1)
        resid.append(res)
    return resid


def calib(i, j, energy=0):
    """
    Return the calibrations for the crystals
    """
    # from Estella
    # (i) is the crystal and (j) is the pmt
    hiEcalib = [
        [0.0061480859, 0.0156416027],
	[0.0081384425, 0.0087327399],
	[0.0069518946, 0.0086127910],
	[0.0084301584, 0.0092437824],
	[0.0126542583, 0.0228144999],
	[0.0115685020, 0.0107546633],
	[0.0068681361, 0.0093653723],
	[0.0452269698, 0.0375218755]
        ]
    loEcalib = [
        [0.0002204598, 0.0002192131],
	[0.0002127388, 0.0002207598],
	[0.0002123268, 0.0002505941],
	[0.0002233761, 0.0002295871],
	[0.0006665573, 0.0005003395],
	[0.0002315343, 0.0002327635],
	[0.0002294907, 0.0002272016],
	[0.0007645054, 0.0008078809]
        ]
    if energy:
        return hiEcalib[int(i)][int(j)]
    else:
        return loEcalib[int(i)][int(j)]


def resol(i, energy=0):
    """
    Return the combined PMT crystal resolution
    """
    
    ##################################################################
    ### 2016-11-01
    ### from toy-MC, combined resolution should also be divided by 2
    ### tested and mc compared to data looks much better!!
    ##################################################################
    
    # from Estella
    # (i) is the crystal number
    # res = sig/E
    # slope = res * sqrt(E)
    # res(E) = slope / sqrt(E) + 0
    hiEresol = [
        1.3976069981,
	1.3796169490,
	1.2337377995,
	1.1778559545,
	2.4846296071,
	1.1813408142,
	1.3249767662,
	2.7941471926
    ]
    loEresol = [
        0.8174646681,
	0.7927112730,
	0.7274639316,
	0.6471664710,
	1.2783412699,
	0.7201966764,
	0.7240077873,
	1.5955441464
    ]
    # 2016-11-01 - now divide by 2
    if energy:
        return hiEresol[int(i)] / 2.
    else:
        return loEresol[int(i)] / 2.


def histparam(energy=0):
    """
    Histogram parameters for number of bins, min, and max in keV
    """
    if energy:
        hmin = 0
        hmax = 3000
        bins = (hmax-hmin)
    else:
        hmin = 0
        hmax = 100
        bins = (hmax-hmin)
    return [bins, hmin, hmax]


def names(i):
    """
    The names of the crystals
    """
    crystals = [
        'NaI-01',
        'NaI-02',
        'NaI-03',
        'NaI-04',
        'NaI-05',
        'NaI-06',
        'NaI-07',
        'NaI-08'
    ]
    return crystals[int(i)]


def cnames(i):
    """
    The C names of the crystals
    """
    crystals = [
        'C1',
        'C2',
        'C3',
        'C4',
        'C5',
        'C6',
        'C7',
        'C8'
    ]
    return crystals[int(i)]


def longNames(i):
    """
    Full name and specs for the crystal
    """
    crystals = [
        'C1  NaI-001  Sample-B  8.3kg',
        'C2  NaI-002  Sample-C  9.2kg',
        'C3  NaI-007  WimpScint-2  9.2kg',
        'C4  AS-3  WimpScint-2  18.5kg',
        'C5  AS-1  Sample-C  18.5kg',
        'C6  NaI-011  WimpScint-3  12.5kg',
        'C7  NaI-012  WimpScint-3  12.5kg',
        'C8  AS-2  Sample-C  18.5kg'
    ]
    return crystals[int(i)]


def cmass(i):
    """
    crystal masses in kg
    """
    mass = [
        8.3,
        9.2,
        9.2,
        18.5,
        18.5,
        12.5,
        12.5,
        18.5
        ]
    return mass[int(i)]


def volumeNames(i):
    """
    For the simulation primary volume names
    """
    volumeNames = [
        'NaIDet01Crystal',
        'NaIDet02Crystal',
        'NaIDet03Crystal',
        'NaIDet04Crystal',
        'NaIDet05Crystal',
        'NaIDet06Crystal',
        'NaIDet07Crystal',
        'NaIDet08Crystal'
    ]
    return volumeNames[int(i)]


def rainbow(N):
    """
    Generate an array of colors for MC plotting
    """
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


def readROOT(rootfile, energy=0):
    
    if not os.path.exists(rootfile):
        print 'rootfile [',rootfile,'] not found'
        sys.exit()
    
    # have to copy or deepcopy the histos!!!
    # going to deepcopy just to make sure
    import copy
    
    rfile = TFile(rootfile, "READ")
    # switch to different dir in rootfile
    #rfile.cd('')
    # to print out all the keys
    #rfile.ls()

    names = [key.GetName() for key in gDirectory.GetListOfKeys()]
    
    mc = {}
    data = [0 for i in range(8)]
    for name in names:
        bits = name.split('_')
        
        if int(bits[0][-1]) != energy: continue
        
        if len(bits) < 4:
            xtal = int(bits[2][-1])
            histo = TH1F(rfile.Get(name))
            data[xtal] = copy.deepcopy(histo)
            #data[xtal] = copy.copy(histo)
            #data[xtal] = histo.Clone()
            #data[xtal] = data[xtal].Clone()
            continue
        
        loc = str(bits[1])
        iso = str(bits[2])
        xtal = int(bits[3][-1])
        
        try: mc[loc]
        except: mc[loc] = {}
        
        try: mc[loc][iso]
        except: mc[loc][iso] = {}
        
        try: mc[loc][iso]['hist']
        except: mc[loc][iso]['hist'] = [0 for i in range(8)]
                
        histo = TH1F(rfile.Get(name))
        mc[loc][iso]['hist'][xtal] = copy.deepcopy(histo)
        #mc[loc][iso]['hist'][xtal] = copy.copy(histo)
        #mc[loc][iso]['hist'][xtal] = histo.Clone()
        #mc[loc][iso]['hist'][xtal] = mc[loc][iso]['hist'][xtal].Clone()
        
    rfile.Close()
    
    # create the locs and isos lists
    # not sure i need this but might as well
    locs=[]
    isos=[]
    for loc in mc:
        if loc not in locs:
            locs.append(loc)
        for iso in mc[loc]:
            if iso not in isos:
                isos.append(iso)
    locs.sort()
    isos.sort()
   
    return data, mc, locs, isos

