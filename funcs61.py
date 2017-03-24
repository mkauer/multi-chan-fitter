#!/usr/bin/env python
######################################################################
# funcs61.py
# 
# Tweak some funcs to get Estella's pmt_id working
# 
# Works with v60 and later versions
# 
# version: 2017-03-24
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ fixed sim path+name so not to combine lsveto and lsvetoair
# ~ think I have the primPMTid cuts working right
# ~ change number of pmts to '2' for scaleSigs61()
# ~ change number of pmts to '2' for scaleBkgs61()
# + added new primPMTid cuts to buildMC61()
# + building off funcs60.py
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
from funcs60 import *


def build61(infile = 'backgrounds60.txt', freuse=0, fchans=0):

    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo60(line, freuse, fchans)
        if 'D' in info['type']:
            data = buildData60(info, data)
        elif 'B' in info['type']:
            bkgs = buildMC61(info, bkgs)
        elif 'S' in info['type']:
            sigs = buildMC61(info, sigs)
        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs


def buildMC61(info, mc):

    for c in info['chans']:
        info['chan'] = c
    
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
                key  = info['key']
                key += '-c'+info['chan']
                key += '-e'+str(e)
                
                try:
                    mc[key] = {}
                    mc[key]['info'] = info
                    mc[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                    #mc[key]['hist'].Sumw2()
                except:
                    print "WARNING: could not find hist -->",key

                try:
                    mc[key]['generated_hist'] = deepcopy(TH1F(rfile.Get(key+'_generated')))
                    mc[key]['generated'] = mc[key]['generated_hist'].GetBinContent(1)
                except:
                    print "WARNING: could not find generated_hist -->",key+'_generated'

                
        else:
            local = amLocal()
            
            ### top level path to the MC
            path1 = 0
            if local: path1 = '/home/mkauer/COSINE/CUP/mc-fitting/sim/newGeometry/'
            else: path1 = '/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
            
            ### 2nd level path to the specific files
            #path2 = info['isof']+'/set2/'+'*'+info['loca']+'*'+info['isof']+'*root'
            path2 = info['isof']+'/set2/'+'*'+info['loca']+info['isof']+'*root'

            # but there will be a few exceptions...
            if info['loca'] == 'internalsurf':
                #path2 = info['isof']+'/set2/surf/10um/'+'*'+info['loca']+'*'+'C'+str(info['xstl'])+'*'+info['isof']+'*root'
                path2 = info['isof']+'/set2/surf/10um/'+'*'+info['loca']+'*'+info['isof']+'*'+'C'+str(info['xstl'])+'*root'

            
            chain  = TChain("MC","")
            nfiles = 0
            nfiles = chain.Add(path1+path2)
                        
            if nfiles > 0:
                print nfiles,'MC files found for', info['loca'], info['isof']
                for e in range(2):

                    i = info['xstl'] - 1

                    key  = info['key']
                    key += '-c'+info['chan']
                    key += '-e'+str(e)
                    
                    mc[key] = {}
                    mc[key]['info'] = info
                    
                    par = histparam(e)
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
                    if   info['loca'] == 'internal':
                        volumeCut = TCut('(primVolumeName == "'+volumeNames(i)+'")')

                    elif info['loca'] == 'internalsurf':
                        volumeCut = TCut('(primVolumeName == "'+volumeNames(i)+'")')
                    
                    
                    
                    ### testing the new pmt cuts
                    #=====================================================================
                    #=====================================================================
                    elif info['loca'] == 'pmt':
                        #volumeCut = TCut('(primVolumeName == "phys_pmt")')
                        volumeCut = TCut('((primPMTid[0] == '+str((int(info['xstl'])*2)-2)+') || (primPMTid[0] == '+str((int(info['xstl'])*2)-1)+'))')
                        #print '!!! Crystal',info['xstl'],'is using primPMTid[0] ==',(int(info['xstl'])*2)-2,'or',(int(info['xstl'])*2)-1
                        #volumeCut = TCut('((primPMTid['+str(i)+'] == '+str((int(info['xstl'])*2)-2)+') || (primPMTid['+str(i)+'] == '+str((int(info['xstl'])*2)-1)+'))')
                        #print '!!! Crystal',info['xstl'],'is using primPMTid['+str(i)+'] ==',(int(info['xstl'])*2)-2,'or',(int(info['xstl'])*2)-1
                    #=====================================================================
                    #=====================================================================
                    
                        
                        
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
                    
                    totalCuts = TCut(volumeCut.GetTitle()+' && '+motherCut.GetTitle())

                    
                    ### I think that has to change?
                    #=====================================================================
                    #=====================================================================
                    #chain.Draw('primVolumeName >> '+key2, totalCuts)
                    
                    if info['loca'] == 'pmt':
                        chain.Draw('primPMTid[0] >> '+key2, totalCuts)
                        #chain.Draw('primPMTid['+str(i)+'] >> '+key2, totalCuts)
                    else:
                        chain.Draw('primVolumeName >> '+key2, totalCuts)
                    #=====================================================================
                    #=====================================================================
                                        
                    
                    generated = temp2.GetEntries()
                    
                    
                    ### using Box-Muller for resolution smearing method
                    ### setting "rng" as the smearing alias
                    chain.SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                    resolFunc = resol60(i,e)
                    chain.SetAlias('sigma', resolFunc)

                    
                    ### do i need to change this too the pmt id cut?
                    ### ie add the volume cut back in? - maybe not...
                    #=====================================================================
                    #=====================================================================
                    masterCut = TCut('('+
                                     energyCut.GetTitle()+' && '+
                                     #volumeCut.GetTitle()+' && '+
                                     eventTypeCut.GetTitle()+' && '+
                                     brokenChainCut.GetTitle()+' && '+
                                     chanCut.GetTitle()
                                     +')')
                    #=====================================================================
                    #=====================================================================
                    
                    
                    selection = '((edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng))'
                    chain.Draw(selection+' >> '+key, masterCut)

                    detected = histo.GetEntries()

                    #histo.SetLineColor(kBlack)
                    #histo.SetMarkerColor(kBlack)
                    #histo.SetLineWidth(1)

                    mc[key]['hist'] = histo
                    #mc[key]['hist'].Sumw2()

                    mc[key]['generated_hist'] = temp2
                    mc[key]['generated'] = generated

            else:
                #print 'WARNING:', loca, isof, 'not found...'
                print 'ERROR: No MC files found for',info['loca'], info['isof'],'... quitting...'
                sys.exit()

    return mc


def scaleBkgs61(bkgs):
    
    for key in bkgs:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]
        
        #druscale = data[x+'-data-'+e]['druScale']
        #runtime  = data[x+'-data-'+e]['runtime']
        
        kev  = 1.     # keV/bin
        day  = 86400. # in seconds
        
        xkgs = cmass(int(x[-1])-1)
        #pmts = 16.
        pmts = 2.
        lskg = 1800.
        
        generated = float(bkgs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        if loca == 'internal':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'internalsurf':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'pmt':
            #scale = bkgs[key]['info']['acti'] * (pmts) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (pmts) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'lsveto':
            #scale = bkgs[key]['info']['acti'] * (lskg) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (lskg) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'lsvetoair':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'airshield':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'steel':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        else:
            print "WARNING: No background scaling for  --> ", loca
            continue

    return bkgs


def scaleSigs61(sigkeys, sigs):

    for key in sigkeys:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]

        kev  = 1.     # keV/bin
        day  = 86400. # in seconds
        
        xkgs = cmass(int(x[-1])-1)
        #pmts = 16.
        pmts = 2.
        lskg = 1800.
        
        generated = float(sigs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        # verbose?
        V = 0
                
        if loca == 'internal':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/kg'
        elif loca == 'internalsurf':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/kg'
        elif loca == 'pmt':
            #fitActivity = sigs[key]['fitscale'] * (1./pmts) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./pmts) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'lsveto':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./lskg) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'lsvetoair':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'airshield':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'steel':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        else:
            print "WARNING: No background scaling for  --> ", loca
            continue

    return sigs

