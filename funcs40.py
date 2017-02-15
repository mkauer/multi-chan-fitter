#!/usr/bin/env python
######################################################################
# funcs40.py
# 
# Get single-hit data into the processing stream
# 
# Works with v40 and later versions
# 
# version: 2017-02-15
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + add new calib41() new calib points from Pushpa
# ~ revamped the way the data calibration is handled for lo/hi energy
# + add new calib40() with old calib points and new format
# ~ cleaned up some of the code
# + added edep and lsveto cuts to MC selection
# + added lsveto cuts to data selection
# + added ncluster cuts to data selection
# ~ revamped the sigle/multi hit TCut making
# ~ tweak MC rootfile selection code so that when selecting lsveto
#   you don't also get lsvetoair - so now the code selects files
#   for *lsvetoU238* for example
# + user _MERGED rootfiles for data
# + finally fixed issue with single/multi hit selection criteria
# + add TCut('1') as default - it works - no more if-else statements
# + add fchan force selection 
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
# /data/COSINE/NTP/phys/V00-02-03_MERGED
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


def getInfo40(line, freuse=0, fchan=0):
    
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
    
    #-----------------------------------------------------------------
    # channel
    info['chan'] = str(bits[1])
    # force a channel of data selection
    if fchan == 1: info['chan'] = 'A'
    if fchan == 2: info['chan'] = 'S'
    if fchan == 3: info['chan'] = 'M'
    #-----------------------------------------------------------------
    
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
            #nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/data/phys/'+info['build']+'/ntp_I*'+info['run']+'*root*')
            nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/data/phys/'+info['build']+'_MERGED/ntp_I*'+info['run']+'*root*')
            
        else:
            #nfiles = chain.Add('/data/COSINE/NTP/phys/'+info['build']+'/ntp_I*'+info['run']+'*root*')
            #nfiles = chain.Add('/data/COSINE/NTP/phys/'+info['build']+'/ntp_I*'+info['run']+'*root.000')
            nfiles = chain.Add('/data/COSINE/NTP/phys/'+info['build']+'_MERGED/ntp_I*'+info['run']+'*root*')
            
        if nfiles > 0:
            print nfiles,'data files found'
            for E in range(2):
                
                par = histparam(E)
                i = info['xstl'] - 1
                
                key = info['key'] + '-e'+str(E)

                run = int(info['run'])
                
                data[key] = {}
                data[key]['info'] = info
                
                histo = TH1F(key, longNames(i), par[0], par[1], par[2])
                
                
                ### old calib
                selection = calib40(i,E)
                ### new calib
                #selection = calib41(i,E)
                
                
                ### reminder - [A]All-hits, [S]Single-hits, [M]Multi-hits
                #-------------------------------------------------------------------------
                if info['chan'] == 'A':
                    chanCut = TCut('(1)')
                
                elif info['chan'] == 'S':
                    edepcuts = ''
                    nclustercuts = ''
                    for j in range(8):
                        if j != i:
                            edepcuts += '(crystal'+str(j+1)+'.'+str(edep)+' <= 0.0) && '
                            nclustercuts += '(crystal'+str(j+1)+'.'+'nc'+' < 4) && '

                    ### remove extra '&&' or '||'
                    edepcuts = edepcuts[:-4]
                    nclustercuts = nclustercuts[:-4]

                    ### liquid scint veto cut - from Pushpa
                    #lsvetocut = '(BLSveto.Charge/110. < 50.)'
                    lsvetocut = '(BLSveto.Charge/110. < 20.)'
                    
                    #chanCut = TCut('(('+edepcuts+'))')
                    #chanCut = TCut('(('+nclustercuts+'))')
                    #chanCut = TCut('(('+edepcuts+') && ('+nclustercuts+'))')
                    
                    ### my cuts
                    #chanCut = TCut('(('+edepcuts+') && ('+nclustercuts+') && ('+lsvetocut+'))')
                    ### Pushpa cuts
                    chanCut = TCut('(('+nclustercuts+') && ('+lsvetocut+'))')
                    
                    
                elif info['chan'] == 'M':
                    edepcuts = ''
                    nclustercuts = ''
                    for j in range(8):
                        if j != i:
                            edepcuts += '(crystal'+str(j+1)+'.'+str(edep)+' > 0.0) || '
                            nclustercuts += '(crystal'+str(j+1)+'.'+'nc'+' >= 4) || '
                    
                    ### remove extra '&&' or '||'
                    edepcuts = edepcuts[:-4]
                    nclustercuts = nclustercuts[:-4]
                    
                    ### liquid scint veto cut - from Pushpa
                    #lsvetocut = '(BLSveto.Charge/110. > 50.)'
                    lsvetocut = '(BLSveto.Charge/110. > 20.)'
                    
                    #chanCut = TCut('(('+edepcuts+'))')
                    #chanCut = TCut('(('+nclustercuts+'))')
                    #chanCut = TCut('(('+edepcuts+') || ('+nclustercuts+'))')
                    
                    ### my cuts
                    #chanCut = TCut('(('+edepcuts+') || ('+nclustercuts+') || ('+lsvetocut+'))')
                    ### Pushpa cuts
                    chanCut = TCut('(('+nclustercuts+') || ('+lsvetocut+'))')
                    
                else:
                    print 'ERROR: I do not know what to do with channel -->',info['chan']
                    print 'Available channels are [A]All-hits, [S]Single-hits, [M]Multi-hits'
                    sys.exit()

                    
                ###-----------------------------------------------------------------------
                chain.Draw(selection+' >> '+key, chanCut)
                ###-----------------------------------------------------------------------
                
                    
                histo.Sumw2()
                histo.SetLineColor(kBlack)
                histo.SetMarkerColor(kBlack)
                histo.SetLineWidth(1)

                data[key]['hist'] = histo

                subruns = nfiles
                key2 = key+'_subruns'
                temp2 = TH1F(key2,'subruns',1,0,1)
                temp2.SetBinContent(1, subruns)
                data[key]['subruns_hist'] = temp2
                data[key]['subruns'] = subruns

                runtime = subruns*subrunTime*60.*60.
                key3 = key+'_runtime'
                temp3 = TH1F(key3,'runtime',1,0,1)
                temp3.SetBinContent(1, runtime)
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
            nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/sim/newGeometry/'
                               +info['isof']+'/set2/'+'*'+info['loca']+info['isof']+'*root')
            
        else:
            nfiles = chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
                               +info['isof']+'/set2/'+'*'+info['loca']+info['isof']+'*root')

        if nfiles > 0:
            print nfiles,'MC files found for', info['loca'], info['isof']
            for E in range(2):
                
                i = info['xstl'] - 1
                
                key = info['key'] + '-e'+str(E)
                
                mc[key] = {}
                mc[key]['info'] = info
                
                par = histparam(E)
                histo = TH1F(key, longNames(i), par[0], par[1], par[2])
                
                energyCut = TCut('(edep['+str(i)+']*1000. > 0.0)')
                
                if info['chan'] == 'A':
                    chanCut = TCut('(1)')
                
                elif info['chan'] == 'S':
                    ### main single/multi hit cut
                    edepcuts = '((singleHitTag['+str(i)+'] > 0.0) && (multipleHitTag['+str(i)+'] < 0.0))'
                    
                    ### liquid scint veto cut - from Estella
                    lsvetocut = '(edep[8]*1000. < 20.0)'
                    
                    #chanCut = TCut('(('+edepcuts+'))')
                    chanCut = TCut('(('+edepcuts+') && ('+lsvetocut+'))')
                    
                elif info['chan'] == 'M':
                    ### main single/multi hit cut
                    edepcuts = '((singleHitTag['+str(i)+'] < 0.0) && (multipleHitTag['+str(i)+'] > 0.0))'
                    
                    ### liquid scint veto cut - from Estella
                    lsvetocut = '(edep[8]*1000. > 20.0)'
                    
                    #chanCut = TCut('(('+edepcuts+'))')
                    chanCut = TCut('(('+edepcuts+') || ('+lsvetocut+'))')
                    
                else:
                    print 'ERROR: I do not know what to do with channel -->',info['chan']
                    print 'Available channels are [A]All-hits, [S]Single-hits, [M]Multi-hits'
                    sys.exit()
                    
                
                volumeCut = TCut('(1)')
                if info['loca'] == 'internal':
                    volumeCut = TCut('(primVolumeName == "'+volumeNames(i)+'")')
                elif info['loca'] == 'pmt':
                    volumeCut = TCut('(primVolumeName == "phys_pmt")')
                elif info['loca'] == 'lsveto':
                    volumeCut = TCut('(primVolumeName == "lsveto")')
                elif info['loca'] == 'lsvetoair':
                    volumeCut = TCut('((primVolumeName == "DetPMTCover") || (primVolumeName == "DetPMTEnvelope"))')
                elif info['loca'] == 'airshield':
                    volumeCut = TCut('(primVolumeName == "LSVetoAirRoom")')
                elif info['loca'] == 'steel':
                    # Not sure this works the way I think it should?
                    # I want all volumes except for the blank volume
                    volumeCut = TCut('(primVolumeName != "")')
                else:
                    print "WARNING: No selection criteria for  --> ", info['loca']
                    continue
                
                
                ### this is needed to get the generated event numbers right!
                motherCut = TCut('(primParticleName == "'+info['chst']+'")')
                
                ### this is needed? to cut out crap events
                ### (event_info.Type) is new
                ### (evt_Type) is old
                #eventTypeCut = TCut('(event_info.Type > 10) || (evt_Type > 10)')
                eventTypeCut = TCut('(event_info.Type > 10)')
                
                #brokenChainCut = 0
                brokenChainCut = TCut('(1)')
                if info['chst'] == 'U238' and info['chsp'] == 'Rn222':
                    brokenChainCut = TCut('((groupNo >= 11) && (groupNo <= 14))')
                if info['chst'] == 'Pb210' and info['chsp'] == 'Pb210':
                    brokenChainCut = TCut('(groupNo == 15)')
                
                key2 = key+'_generated'
                temp2 = TH1F(key2, 'generated',1,0,1)
                
                totalCuts = TCut(volumeCut.GetTitle()+'&&'+motherCut.GetTitle())
                chain.Draw('primVolumeName >> '+key2, totalCuts)
                
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
                
                masterCut = TCut('('+
                                 energyCut.GetTitle()+' && '+
                                 #volumeCut.GetTitle()+' && '+
                                 eventTypeCut.GetTitle()+' && '+
                                 brokenChainCut.GetTitle()+' && '+
                                 chanCut.GetTitle()
                                 +')')
                selection = '((edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng))'

                
                ###-------------------------------------------------------------------------------
                chain.Draw(selection+' >> '+key, masterCut)
                ###-------------------------------------------------------------------------------

                
                detected = histo.GetEntries()
                                
                histo.SetLineColor(kBlack)
                histo.SetMarkerColor(kBlack)
                histo.SetLineWidth(1)

                mc[key]['hist'] = histo
                mc[key]['generated_hist'] = temp2
                mc[key]['generated'] = generated
                
        else:
            #print 'WARNING:', loca, isof, 'not found...'
            print 'ERROR: No MC files found for',info['loca'], info['isof'],'... quitting...'
            sys.exit()
            
    return mc

    
def build40(infile = 'backgrounds40.txt', freuse=0, fchan=0):

    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo40(line, freuse, fchan)
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


def calib40(i, E=0):
    """
    Return the crystal calibrations
    """
    # from Pushpa
    # this is the old calib that uses crystal1.

    # adc = crystalX.energyD
    # E = (adc+[0])*[1]
    hiEcalib = [
        [0.0, 1.000],
	[0.0, 0.992],
	[0.0, 1.000],
	[0.0, 1.000],
	[0.0, 0.966],
	[0.0, 1.014],
	[0.0, 0.984],
	[0.0, 1.000]
    ]
    
    # adc = crystalX.qc5
    # E = (adc+[0])*[1]
    loEcalib = [
        [0.0, 1.06000e-4],
	[0.0, 1.05676e-4],
	[0.0, 1.15900e-4],
	[0.0, 1.12700e-4],
	[0.0, 2.86907e-4],
	[0.0, 1.17634e-4],
	[0.0, 1.11900e-4],
	[0.0, 3.84263e-4]
    ]
    
    if E:
        edep = '(crystal'+str(i+1)+'.energyD)'
        b,m = hiEcalib[int(i)]
    else:
        edep = '(crystal'+str(i+1)+'.qc5)'
        b,m = loEcalib[int(i)]
        
    selection = '(('+edep+'+'+str(b)+')*'+str(m)+')'
    #if E: return hiEcalib[int(i)]
    #else: return loEcalib[int(i)]
    return selection


def calib41(i, E=0):
    """
    Return the crystal calibrations
    """

    # from Pushpa
    # https://cupwiki.ibs.re.kr/Kims/NaICalibration?validation_key=fae031e9908735aeff81aacfcbf83931

    # adc = pmtX1.rqcD1_5+pmtX2.rqcD1_5
    # E = (adc+[0])*[1]
    hiEcalib = [
        [3942.0, 1./232.5],
 	[4327.0, 1./241.9],
 	[3834.0, 1./262.1],
 	[3909.0, 1./231.3],
 	[   0.0, 1./125.0],
        #[   0.0, 8.401e-3],
 	[-231.0, 1./175.7],
 	[4769.0, 1./260.9],
        [   0.0, 1./ 40.2]
        #[   0.0, 2.487e-2]
    ]

    # adc = crystalX.qc5
    # E = (adc+[0])*[1]
    loEcalib = [
        [0.0, 0.000109373],
 	[0.0, 0.000108797],
 	[0.0, 0.000117414],
 	[0.0, 0.000114566],
 	[0.0, 0.000342601],
 	[0.0, 0.000118527],
 	[0.0, 0.000113608],
 	[0.0, 0.000500613]
    ]
    
    if E:
        edep = '(pmt'+str(i+1)+'1.rqcD1_5 + pmt'+str(i+1)+'2.rqcD1_5)'
        b,m = hiEcalib[int(i)]
    else:
        edep = '(crystal'+str(i+1)+'.qc5)'
        b,m = loEcalib[int(i)]
        
    selection = '(('+edep+'+'+str(b)+')*'+str(m)+')'
    #if E: return hiEcalib[int(i)]
    #else: return loEcalib[int(i)]
    return selection

