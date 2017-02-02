#!/usr/bin/env python
######################################################################
# funcs40.py
# 
# Get single-hit data into the processing stream
# 
# Works with v40 and later versions
# 
# version: 2017-02-01
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + move to V00.02.03 processing of the 1546 data
# + add force reuse options in funcs40.py
# + reuse rootfiles works now in funcs40.py
# + forgot to add "from copy import deepcopy" in funcs40.py - oops
# + many tweaks to get "reuse" working
# + add scaleBkgs40() in funcs40.py
# + add scaleSigs40() in funcs40.py
# + add info{} to data,bkgs,sigs in funcs40.py
# + add dataDRU40() in funcs40.py
# + add build() in funcs40.py
# + add buildMC40() in funcs40.py
# + add buildData40() in funcs40.py
# + add getInfo() in funcs40.py
# + funcs40.py clean slate
# 
# email me: mkauer@physics.wisc.edu
######################################################################
# 
# Where is MC?
# /data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/
# 
# Where is the processed data?
# Right now I'm using:
# /data/COSINE/NTP/phys/V00-02-03
# and run ntp_I001546*
# 
######################################################################

import os,sys
from copy import deepcopy
import numpy as np
from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs32 import *


def getInfo40(line, freuse=0):
    
    info={}

    info['line'] = line
    
    bits = line.split()
    #bits = [x.strip() for x in bits]
        
    # data type
    info['type'] = str(bits[0])

    #-----------------------------------------------------------------
    # reuse a joined rootfile
    if 'R' in info['type']: info['reuse'] = 1
    else: info['reuse'] = 0
    # force reuse of everything - nice for debugging
    if freuse == 1: info['reuse'] = 1
    # force not reuse of everything - nice for debugging
    if freuse == 2: info['reuse'] = 0
    #-----------------------------------------------------------------
    
    # channel
    info['chan'] = str(bits[1])
    
    # crystal
    info['xstl'] = int(bits[2])
    
    # reuse rootfile name
    if 'root' in str(bits[-1]):
        info['rootfile'] = str(bits[-1])
    else:
        info['rootfile'] = 0


    if 'D' in info['type']:
        
        ### run number
        info['run'] = str(bits[4])

        ### processing version
        info['build'] = str(bits[5])

        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-data'
        key += '-c'+info['chan']
        info['key'] = key
        
    else:
        
        ### background location
        info['loca'] = str(bits[3])

        ### top level isotope file name
        info['isof'] = str(bits[4])

        ### isotope chain break start
        info['chst'] = str(bits[5])

        ### isotope chain break stop
        info['chsp'] = str(bits[6])

        ### activity
        info['acti'] = float(bits[7])

        ### error on activity
        info['erro'] = float(bits[8])

        ### fit bounds
        info['fbnd'] = [float(bits[9]), float(bits[10])]

        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-'+info['loca']
        if info['chst'] == info['chsp']: key += '-'+info['chst']
        else: key += '-'+info['chst']+'_'+info['chsp']
        key += '-c'+info['chan']
        info['key'] = key
        
    return info


def buildData40(info, data):
    
    if info['reuse']:
        
        if info['rootfile']:
            rootfile = info['rootfile']
            #print 'INFO: using rootfile -->',rootfile
        else:
            print 'WARNING: no rootfile specified in backgrounds file'
            return data
        
        if not os.path.exists(rootfile):
            print 'WARNING: rootfile not found -->', rootfile
            return data
        
        rfile = TFile(rootfile, "READ")

        for e in range(2):
            key = info['key']+'-e'+str(e)
            try:
                data[key] = {}
                data[key]['info'] = info
                data[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
            except:
                print "WARNING: could not find hist -->",key

            try:
                data[key]['subruns_hist'] = deepcopy(TH1F(rfile.Get(key+'_subruns')))
                data[key]['subruns'] = data[key]['subruns_hist'].GetBinContent(1)
            except:
                print "WARNING: could not find subruns_hist -->",key+'_subruns'
                
            try:
                data[key]['runtime_hist'] = deepcopy(TH1F(rfile.Get(key+'_runtime')))
                data[key]['runtime'] = data[key]['runtime_hist'].GetBinContent(1)
            except:
                print "WARNING: could not find runtime_hist -->",key+'_runtime'
            
            #except:
            #    del data[key]
            #    print "WARNING: could not find hist -->",key
            #    return data
        
        return data
    
    else:
        local = amLocal()

        ### will need to do this dynamically eventually
        ### run 1324 has 1 hour subruns
        ### run 1544 has 2 hour subruns
        ### run 1546 has 2 hour subruns
        subrunTime = float(2.) # in hours
        
        chain = TChain("ntp","")
        nfiles=0
        if local:
            nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/data/phys/'+info['build']+'/ntp_I*'+info['run']+'*root*')
        else:
            nfiles = chain.Add('/data/COSINE/NTP/phys/'+info['build']+'/ntp_I*'+info['run']+'*root*')

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
                #for i in range(8):
                i = info['xstl'] - 1

                #key  = 'x'+info['xstl']
                #key += '-data'
                #key += '-c'+info['chan']
                #key += '-e'+str(E)
                key = info['key'] + '-e'+str(E)
                
                data[key] = {}
                data[key]['info'] = info

                histo = TH1F(key, longNames(i), par[0], par[1], par[2])
                
                ### reminder - [A]All-hits, [S]Single-hits, [M]Multi-hits
                #-------------------------------------------------------------------------
                if info['chan'] == 'A':
                    chain.Draw('crystal'+str(i+1)+'.'+str(edep)+'*'+str(calib22(i,E))+' >> '+str(key))
                #elif info['chan'] == 'S':
                #    TCut = 
                else:
                    print 'ERROR: I do not know what to do with channel -->',info['chan']
                    print 'Available channels are [A]All-hits, [S]Single-hits, [M]Multi-hits'
                    sys.exit()
                
                histo.Sumw2()
                histo.SetLineColor(kBlack)
                histo.SetMarkerColor(kBlack)
                histo.SetLineWidth(1)

                data[key]['hist'] = histo

                subruns = nfiles
                key2 = key+'_subruns'
                temp2 = TH1F(key2,'subruns',1,0,1)
                temp2.SetBinContent(1, subruns)
                #for i in range(subruns):
                #    temp2.Fill(0.5)
                data[key]['subruns_hist'] = temp2
                data[key]['subruns'] = subruns

                runtime = subruns*subrunTime*60.*60.
                key3 = key+'_runtime'
                temp3 = TH1F(key3,'runtime',1,0,1)
                temp3.SetBinContent(1, runtime)
                #for i in range(runtime):
                #    temp3.Fill(0.5)
                data[key]['runtime_hist'] = temp3
                data[key]['runtime'] = runtime
        
        else:
            print 'ERROR: No data files found... quitting...'
            sys.exit()
        
    return data


def buildMC40(info, mc, calib=2):

    if info['reuse']:
        if info['rootfile']:
            rootfile = info['rootfile']
        else:
            print 'WARNING: no rootfile specified in backgrounds file'
            return data
        
        if not os.path.exists(rootfile):
            print 'WARNING: rootfile not found -->', rootfile
            return data
        
        rfile = TFile(rootfile, "READ")

        for e in range(2):
            key = info['key']+'-e'+str(e)
            try:
                mc[key] = {}
                mc[key]['info'] = info
                mc[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
            except:
                print "WARNING: could not find hist -->",key

            try:
                mc[key]['generated_hist'] = deepcopy(TH1F(rfile.Get(key+'_generated')))
                mc[key]['generated'] = mc[key]['generated_hist'].GetBinContent(1)
            except:
                print "WARNING: could not find generated_hist -->",key+'_generated'
                
                #mc[key]['acti'] = info['acti']
                #mc[key]['erro'] = info['erro']
                #mc[key]['fbnd'] = info['fbnd']
            #except:
            #    del mc[key]
            #    print "WARNING: could not find hist -->",key
            #    continue
                
        return mc
    
    else:
        local = amLocal()

        chain = TChain("MC","")
        nfiles=0
        if local:
            nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/sim/newGeometry/'+info['isof']+'/set2/'+'*'+info['loca']+'*root')
        else:
            nfiles = chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'+info['isof']+'/set2/'+'*'+info['loca']+'*root')

        if nfiles > 0:
            print nfiles,'MC files found for', info['loca'], info['isof']
            for E in range(2):
                
                i = info['xstl'] - 1
                
                ### see if i can use unique names for all the histograms
                #key  = 'x'+info['xstl']
                #key += '-'+info['loca']
                #if info['chst'] == info['chsp']: key += '-'+info['chst']
                #else: key += '-'+info['chst']+'_'+info['chsp']
                #key += '-c'+info['chan']
                #key += '-e'+str(E)
                key = info['key'] + '-e'+str(E)
                
                mc[key] = {}
                mc[key]['info'] = info
                
                par = histparam(E)
                histo = TH1F(key, longNames(i), par[0], par[1], par[2])

                cut1 = TCut('edep['+str(i)+']*1000. > 0')
                
                if info['loca'] == 'internal':
                    cut2 = TCut('primVolumeName == "'+volumeNames(i)+'"')
                elif info['loca'] == 'pmt':
                    cut2 = TCut('primVolumeName == "phys_pmt"')
                elif info['loca'] == 'lsveto':
                    cut2 = TCut('primVolumeName == "lsveto"')
                elif info['loca'] == 'lsvetoair':
                    cut2 = TCut('(primVolumeName == "DetPMTCover") || (primVolumeName == "DetPMTEnvelope")')
                elif info['loca'] == 'airshield':
                    cut2 = TCut('primVolumeName == "LSVetoAirRoom"')
                elif info['loca'] == 'steel':
                    # not sure this works the way I think it should?
                    cut2 = TCut('primVolumeName != ""')
                else:
                    print "WARNING: No selection criteria for  --> ", info['loca']
                    continue
                
                ### this is needed to get the generated event numbers right!
                cut3 = TCut('primParticleName == "'+info['chst']+'"')
                
                ### this is needed? to cut out crap events
                ### (event_info.Type) is new
                ### (evt_Type) is old
                #cut4 = TCut('(event_info.Type > 10) || (evt_Type > 10)')
                cut4 = TCut('event_info.Type > 10')
                
                cut100 = 0
                if info['chst'] == 'U238' and info['chsp'] == 'Rn222':
                    cut100 = TCut('(groupNo >= 11) && (groupNo <= 14)')
                if info['chst'] == 'Pb210' and info['chsp'] == 'Pb210':
                    cut100 = TCut('groupNo == 15')
                
                key2 = key+'_generated'
                temp2 = TH1F(key2, 'generated',1,0,1)
                
                ### v31 - maybe need cut100 for generated events
                ### to get the eff/normalization right?
                #chain.Draw('primVolumeName >> '+key2, cut2)
                if cut100:
                    #chain.Draw('primVolumeName >> '+key2, cut2+cut3+cut100)
                    chain.Draw('primVolumeName >> '+key2, cut2+cut3)
                else:
                    chain.Draw('primVolumeName >> '+key2, cut2+cut3)
                
                generated = temp2.GetEntries()
                
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

                mc[key]['hist'] = histo
                mc[key]['generated_hist'] = temp2
                mc[key]['generated'] = generated
                
                #mc[key]['acti'] = info['acti']
                #mc[key]['erro'] = info['erro']
                #mc[key]['fbnd'] = info['fbnd']
                
                #print 'Eff =',detected/generated
                """
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
                """
                
        else:
            #print 'WARNING:', loca, isof, 'not found...'
            print 'ERROR: No MC files found for',info['loca'], info['isof'],'... quitting...'
            sys.exit()
            
    return mc

    
def build40(infile = 'backgrounds40.txt', freuse=0):

    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo40(line, freuse)
        if 'D' in info['type']:
            data = buildData40(info, data)
        elif 'B' in info['type']:
            bkgs = buildMC40(info, bkgs)
        elif 'S' in info['type']:
            sigs = buildMC40(info, sigs)
        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs


def dataDRU40(data):
    """
    Scale data into DRU
    """
        
    # this will eventually have to be made dynamic
    keVperBin = float(1.)
    
    for key in data:
        
        i = int(key.split('-')[0].split('x')[-1]) - 1
        days = float((data[key]['runtime'])/(60.*60.*24.))
        kgs = float(cmass(i))
        scale = float(1./(days*kgs*keVperBin))
        
        data[key]['hist'].Scale(scale)
        data[key]['druScale'] = scale
        
    return data


def scaleBkgs40(bkgs):
    
    for key in bkgs:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]
        
        #druscale = data[x+'-data-'+e]['druScale']
        #runtime  = data[x+'-data-'+e]['runtime']
        
        kev  = 1.     # keV/bin
        day  = 86400. # in seconds
        
        xkgs = cmass(int(x[-1])-1)
        pmts = 16.
        lskg = 1800.
        
        generated = float(bkgs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        if loca == 'internal':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
        elif loca == 'pmt':
            #scale = bkgs[key]['info']['acti'] * (pmts) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (pmts) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
        elif loca == 'lsveto':
            #scale = bkgs[key]['info']['acti'] * (lskg) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (lskg) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
        elif loca == 'lsvetoair':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
        elif loca == 'airshield':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
        elif loca == 'steel':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
        else:
            print "WARNING: No background scaling for  --> ", loca
            continue

    return bkgs


def scaleSigs40(sigkeys, sigs):

    for key in sigkeys:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]
        
        #druscale = data[x+'-data-'+e]['druScale']
        #runtime = data[x+'-data-'+e]['runtime']

        kev  = 1.     # keV/bin
        day  = 86400. # in seconds
        
        xkgs = cmass(int(x[-1])-1)
        pmts = 16.
        lskg = 1800.
        
        generated = float(sigs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        # verbose?
        V = 1
        E = int(e[-1])
        
        if loca == 'internal':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V and E: print '!!!!', key, sigs[key]['info']['acti'],'mBq/kg'
        if loca == 'pmt':
            #fitActivity = sigs[key]['fitscale'] * (1./pmts) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./pmts) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V and E: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        if loca == 'lsveto':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./lskg) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V and E: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        if loca == 'lsvetoair':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V and E: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        if loca == 'airshield':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V and E: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        if loca == 'steel':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V and E: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'

    return sigs

