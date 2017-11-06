#!/usr/bin/env python
######################################################################
# funcs70.py
# 
# New funcs for extending histos
# 
# Works with v70 and later versions
# 
# version: 2017-07-26
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ modify updateBkgsFile70() to exclude lines between the """
# ~ tweaked updateBkgsFile70() to handle new fit results format
# + add uniqString() to unique 'SM' and 'S' for example
# ~ (evt_Type > 10) still works with new rootfiles
# ~ primParticalName is 'triton' for H3
# ~ tweaks to buildMC70() to get Cd109 and H3 working
# + new paths to pushpa's Cd109 and H3
# ~ fixed issue with Te names not showing the 'm'
# + added setBinError() to test different error schemes
# ~ updated updateBkgsFile70()
# ~ treat surface and bulk teflon and copper as a unit of 1
# + added fit errors to the signal scaling function
# + added new criteria for surface things, teflon, and copper
# + include the volume cut for all things!
# 
# email me: mkauer@physics.wisc.edu
######################################################################

import os,sys
import re
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs64 import *


def build70(infile = 'backgrounds640.txt', freuse=0, fchans=0):

    data = {}
    bkgs = {}
    sigs = {}
    runtime = 0
    
    for line in readFile(infile):
        info = getInfo64(line, freuse, fchans)
        
        if 'D' in info['type']:
            data, runtime = buildData70(info, data)
            
        elif 'B' in info['type']:
            bkgs = buildMC70(info, bkgs)
            
        elif 'F' in info['type']:
            sigs = buildMC70(info, sigs)
            
        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs, runtime


def buildMC70(info, mc):

    base = baseDir()

    if info['reuse']:
        for c in info['chans']:
            info['chan'] = c
            
            if not info['rootfile']:
                print 'ERROR: no rootfile specified in backgrounds file'
                sys.exit()
            
            rootfile = base+info['rootfile']
            if not os.path.exists(rootfile):
                print 'ERROR: rootfile not found -->', rootfile
                sys.exit()
            
            rfile = TFile(rootfile, "READ")

            for e in range(2):
                key  = info['key']
                key += '-c'+info['chan']
                key += '-e'+str(e)
                
                try:
                    mc[key] = {}
                    mc[key]['info'] = info
                    mc[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                    mc[key]['pars'] = getPars(mc[key]['hist'])
                    #mc[key]['hist'].Sumw2()
                except:
                    print "ERROR: could not find hist -->",key
                    sys.exit()
                    
                try:
                    mc[key]['generated_hist'] = deepcopy(TH1F(rfile.Get(key+'_generated')))
                    #mc[key]['generated'] = mc[key]['generated_hist'].GetBinContent(1)
                    mc[key]['generated'] = mc[key]['generated_hist'].GetEntries()
                except:
                    print "ERROR: could not find generated_hist -->",key+'_generated'
                    sys.exit()
        
    else:
        
        path1 = '/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
        if info['isof'] == 'Cd109' or info['isof'] == 'H3':
            path1 = '/data/COSINE/WORK/pushpa/sim/process/Crystal7/'
        if not onCup():
            path1 = '/home/mkauer/COSINE/CUP/mc-fitting'+path1
        
        location = info['floca']
        if location == 'extpmt':
            location = 'pmt'
        
        ### 2nd level path to the specific files
        path2 = info['isof']+'/set2/'+'*'+location+info['isof']+'*root'
        
        if location == 'internal-surf-10um':
            path2 = info['isof']+'/set2/surf/10um/'+'*'+location+'-C'+str(info['xstl'])+'*root'
        if location == 'cu-surf-10um':
            path2 = info['isof']+'/set2/surf/cu-10um/'+'*'+location+'-C'+str(info['xstl'])+'*root'
        if location == 'teflon-surf-10um':
            path2 = info['isof']+'/set2/surf/teflon-10um/'+'*'+location+'-C'+str(info['xstl'])+'*root'
        
        if info['isof'] == 'Cd109' or info['isof'] == 'H3':
            path2 = '*'+info['isof']+'*.root'

        
        print 'INFO: looking for files with -->',path1+path2
        
        
        chain  = TChain("MC","")
        nfiles = 0
        nfiles = chain.Add(path1+path2)
        
        if nfiles > 0:
            
            print 'INFO:',nfiles,'files found for', info['loca'], info['isof']
            
            for c in info['chans']:
                info['chan'] = c
                
                for e in range(2):
                    
                    # DEFINE HIST AND KEY
                    #-----------------------------------------------------------------------
                    i = info['xstl'] - 1

                    key  = info['key']
                    key += '-c'+info['chan']
                    key += '-e'+str(e)
                    
                    mc[key] = {}
                    mc[key]['info'] = info
                    
                    pars = histparam64(e)
                    mc[key]['pars'] = pars
                    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
                    #-----------------------------------------------------------------------
                    
                    ### totally hacking away at this bit of code
                    #=====================================================================================
                    #=====================================================================================

                    energyCut = TCut('(edep['+str(i)+']*1000. > 0.0)')
                    ### do we need an energy cut?
                    ### yes needed to get the generated events normalization right!!!
                    ### will this circumvent the primVolumeName bug? - NOPE!
                    #energyCut = TCut('(1)')
                    
                    if info['chan'] == 'A':
                        chanCut = TCut('(1)')

                    elif info['chan'] == 'S':
                        ### main single/multi hit cut
                        #hitCut = '((singleHitTag['+str(i)+'] > 0.0) && (multipleHitTag['+str(i)+'] < 0.0))'
                        hitCut = '((singleHitTag['+str(i)+'] == 1.0))'
                        ### do i have the single multi hit cuts wrong??
                        #hitCut = '((singleHitTag['+str(i)+'] > 0.0))'
                        #hitCut = '((singleHitTag['+str(i)+'] < 0.0))'
                        
                        ### liquid scint veto cut - from Estella
                        lsvetocut = '(edep[8]*1000. < 20.0)'
                        #lsvetocut = '(edep[8]*1000. < 50.0)'
                        ### do we need a lsveto cut for MC?
                        #lsvetocut = '(1)'
                        
                        chanCut = TCut('(('+hitCut+') && ('+lsvetocut+'))')
                        ### skip the lsvetocut because it's always 1
                        #chanCut = TCut('(('+hitCut+'))')
                        
                    elif info['chan'] == 'M':
                        ### main single/multi hit cut
                        #hitCut = '((singleHitTag['+str(i)+'] < 0.0) && (multipleHitTag['+str(i)+'] > 0.0))'
                        hitCut = '((multipleHitTag['+str(i)+'] == 1.0))'
                        ### do i have the single multi hit cuts wrong??
                        #hitCut = '((multipleHitTag['+str(i)+'] > 0.0))'
                        #hitCut = '((multipleHitTag['+str(i)+'] < 0.0))'
                        
                        ### liquid scint veto cut - from Estella
                        lsvetocut = '(edep[8]*1000. > 20.0)'
                        #lsvetocut = '(edep[8]*1000. > 50.0)'
                        #lsvetocut = '(1)'
                        
                        chanCut = TCut('(('+hitCut+') || ('+lsvetocut+'))')
                        #chanCut = TCut('(('+hitCut+') && ('+lsvetocut+'))')
                        ### skip the lsvetocut because it's always 1
                        #chanCut = TCut('(('+hitCut+'))')
                        
                    else:
                        print 'ERROR: I do not know what to do with channel -->',info['chan']
                        print 'Available channels are [A]All-hits, [S]Single-hits, [M]Multi-hits'
                        sys.exit()


                    ### trying to figure out lsveto weirdness
                    #if info['loca'] == 'lsveto':
                        #chanCut = TCut('('+hitCut+')')
                        #chanCut = TCut('('+lsvetocut+')')
                        #chanCut = TCut('(1)')
                    
                    
                    #=====================================================================================
                    #=====================================================================================
                    
                    
                    pmt1 = str((int(info['xstl'])*2)-2)
                    pmt2 = str((int(info['xstl'])*2)-1)
                    #print pmt1,pmt2
                    
                    volumeCut = TCut('(1)')
                    if info['loca'] == 'internal':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Crystal")')
                    
                    elif 'internalsurf' in info['loca']:
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Crystal")')
                    
                    elif 'cusurf' in info['loca']:
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Case")')
                    
                    elif 'teflonsurf' in info['loca']:
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Teflon")')
                    
                    elif info['loca'] == 'cucase':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Case")')
                    
                    elif info['loca'] == 'teflon':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Teflon")')
                    
                    elif info['loca'] == 'pmt':
                        volumeCut = TCut('((primPMTid[0] == '+pmt1+') || (primPMTid[0] == '+pmt2+'))')
                    
                    elif info['loca'] == 'extpmt':
                        volumeCut = TCut('((primPMTid[0] != '+pmt1+') && (primPMTid[0] != '+pmt2+'))')
                        """
                        tcan = TCanvas('tcan','tcan',800,600)
                        testing = TH2I('testing','primPMTid[1] : primPMTid[0]',17,0,17,17,0,17)
                        chain.Draw('primPMTid[1] : primPMTid[0] >> testing')
                        testing.Draw()
                        raw_input('enter to continue')
                        """
                    
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
                        #volumeCut = TCut('(1)')

                    else:
                        print "WARNING: No selection criteria for  --> ", info['loca']
                        continue


                    ### this is needed to get the generated event numbers right!
                    pname = info['chst']
                    if info['chst'].startswith('Te') \
                       and info['chst'].endswith('m'):
                        pname = info['chst'][:-1]
                    if info['chst'] == 'H3':
                        pname = 'triton'
                    #motherIsoCut = TCut('(primParticleName == "'+info['chst']+'")')
                    motherIsoCut = TCut('(primParticleName == "'+pname+'")')
                    
                    
                    ### This is needed to cut out crap events
                    ### (evt_Type) is old
                    ### (event_info.Type) is new
                    ### (evt_Type) still works with new files
                    #eventTypeCut = TCut('(event_info.Type > 10)')
                    #if info['isof'] == 'Cd109' or info['isof'] == 'H3':
                    #    eventTypeCut = TCut('(evt_Type > 10)')
                    eventTypeCut = TCut('(evt_Type > 10)')
                    
                    
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    #   CRITIAL POINT IN THE SIM !!!
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    brokenChainCut = TCut('(1)')
                    groupCuts = groupNum62(info)
                    
                    if groupCuts:
                        """
                        print 'INFO:',info['isof'],\
                            'chan start =',info['chst'],\
                            'chan stop =',info['chsp'],\
                            'groupNo start >=',groupCuts[0],\
                            'groupNo stop <',groupCuts[1]
                        """
                        
                        ### this is what i used before but it might be wrong??
                        ### it should probably be 'groupNo < number' instead of 'groupNo <= number' ???
                        # include last group
                        #brokenChainCut = TCut('((groupNo >= '+str(groupCuts[0])+') && (groupNo <= '+str(groupCuts[1])+'))')

                        ### I think it actually needs to exclude the last decay product groupNo - not sure??
                        ### minor change '<=' to '<' but I think it's a big difference... 
                        # do not include last group
                        brokenChainCut = TCut('((groupNo >= '+str(groupCuts[0])+') && (groupNo < '+str(groupCuts[1])+'))')

                    ### not sure how this groupNo thing works - I need to ask Eunju for details!
                    ### everything should have groupNo > 1 right?
                    ### NO! Cosmogenics are groupNo == 1
                    #else:
                    #    brokenChainCut = TCut('(groupNo > 1)')
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                    
                    
                    ### create a hist of the generated events
                    #---------------------------------------------------------------------
                    key2 = key+'_generated'
                    temp2 = TH1F(key2, 'generated', 1, 1, 2)
                    
                    #totalCuts = TCut(volumeCut.GetTitle()+' && '+motherIsoCut.GetTitle())
                    totalCuts = TCut(motherIsoCut.GetTitle()+' && '+volumeCut.GetTitle())
                    
                    ### I think that has to change?
                    #=====================================================================
                    #=====================================================================
                    #chain.Draw('primVolumeName >> '+key2, totalCuts)
                    
                    if info['loca'] == 'pmt':
                        chain.Draw('1 >> '+key2, totalCuts)
                        #chain.Draw('primPMTid >> '+key2, totalCuts)
                        #chain.Draw('primPMTid[0] >> '+key2, totalCuts)
                        #chain.Draw('primPMTid['+str(i)+'] >> '+key2, totalCuts)
                        #chain.Draw('primVolumeName >> '+key2, totalCuts)
                    else:
                        #chain.Draw('primVolumeName >> '+key2, totalCuts)
                        chain.Draw('1 >> '+key2, totalCuts)
                        
                    #=====================================================================
                    #=====================================================================
                    
                    #mc[key]['generated_hist'] = temp2
                    #mc[key]['generated'] = generated

                    mc[key]['generated_hist'] = temp2
                    #mc[key]['generated_hist'] = deepcopy(TH1F(temp2))
                    #mc[key]['generated'] = mc[key]['generated_hist'].GetBinContent(1)
                    #mc[key]['generated'] = temp2.GetBinContent(1)
                    mc[key]['generated'] = temp2.Integral()
                    #mc[key]['generated'] = temp2.GetEntries()
                    
                    entries = temp2.GetEntries()
                    generated = mc[key]['generated']
                    #print '!!!! ',key2,'',generated,' generated events  (',entries,'entries )' 
                    
                    #---------------------------------------------------------------------
                    
                    
                    ### using Box-Muller for resolution smearing method
                    ### setting "rng" as the smearing alias
                    chain.SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                    #resolFunc = resol60(i,e)
                    resolFunc = resol64(i,e)
                    chain.SetAlias('sigma', resolFunc)

                    
                    ### do i need to change this too the pmt id cut?
                    ### ie add the volume cut back in? - maybe not...
                    #=====================================================================
                    #=====================================================================
                    
                    masterCut = TCut('('+
                                     energyCut.GetTitle()+' && '+
                                     eventTypeCut.GetTitle()+' && '+
                                     brokenChainCut.GetTitle()+' && '+
                                     chanCut.GetTitle()+' && '+
                                     volumeCut.GetTitle()
                                     +')')
                    """
                    ### keep old cuts for now
                    masterCut = TCut('('+
                                     energyCut.GetTitle()+' && '+
                                     #volumeCut.GetTitle()+' && '+
                                     eventTypeCut.GetTitle()+' && '+
                                     brokenChainCut.GetTitle()+' && '+
                                     chanCut.GetTitle()
                                     +')')
                    ### need the volumeCut for the pmt selection
                    if 'pmt' in info['loca']:
                        masterCut = TCut('('+
                                         energyCut.GetTitle()+' && '+
                                         eventTypeCut.GetTitle()+' && '+
                                         brokenChainCut.GetTitle()+' && '+
                                         chanCut.GetTitle()+' && '+
                                         volumeCut.GetTitle()
                                         +')')
                    """
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

                    
        else:
            #print 'WARNING:', loca, isof, 'not found...'
            print 'WARNING: No MC files found for --> ',\
                'x'+str(info['xstl']), info['floca'], info['isof']
            #sys.exit()
    
    return mc


def scaleData70(data, dru=0):
    """
    Scale for DRU or not
    """
    for key in data:
        i = int(key.split('-')[0].split('x')[-1]) - 1
        days = 1.
        xkgs = 1.
        keVperBin = 1.
        if dru:
            days = float((data[key]['runtime'])/86400.)
        xkgs = float(cmass(i))
        keVperBin = 1./float(data[key]['pars'][3])
        scale = float(1./(days*xkgs*keVperBin))
        data[key]['hist'].Scale(scale)
        data[key]['druScale'] = scale
    return data


def scaleBkgs70(bkgs,runtime=0):
    
    for key in bkgs:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]

        ### 1 day in seconds
        day = 86400.
        xkgs = 1.
        keVperBin = 1.
        if runtime:
            day = runtime
        xkgs = cmass(int(x[-1])-1)
        keVperBin = 1./float(bkgs[key]['pars'][3])

        nmass     = cmass(int(x[-1])-1)
        pmts      = 2.
        extpmts   = 14.
        lskg      = 1800.
        
        ### treat the NaI surface and the copper and teflon as 1 unit
        surf = 1.
        
        generated = float(bkgs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        if loca == 'internal':
            scale = bkgs[key]['info']['acti'] * (nmass) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif 'surf' in loca:
            scale = bkgs[key]['info']['acti'] * (surf) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif loca == 'cucase':
            scale = bkgs[key]['info']['acti'] * (surf) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif loca == 'teflon':
            scale = bkgs[key]['info']['acti'] * (surf) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif loca == 'pmt':
            scale = bkgs[key]['info']['acti'] * (pmts) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif loca == 'extpmt':
            scale = bkgs[key]['info']['acti'] * (extpmts) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif loca == 'lsveto':
            scale = bkgs[key]['info']['acti'] * (lskg) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif loca == 'lsvetoair':
            scale = bkgs[key]['info']['acti'] * (nmass) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif loca == 'airshield':
            scale = bkgs[key]['info']['acti'] * (nmass) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif loca == 'steel':
            scale = bkgs[key]['info']['acti'] * (nmass) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        else:
            print "WARNING: No background scaling for  --> ", loca
            continue

    return bkgs


def scaleSigs70(sigkeys, sigs, runtime=0):

    for key in sigkeys:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]

        ### 1 day in seconds
        day = 86400.
        xkgs = 1.
        keVperBin = 1.
        if runtime:
            day = runtime
        xkgs = cmass(int(x[-1])-1)
        keVperBin = 1./float(sigs[key]['pars'][3])

        nmass     = cmass(int(x[-1])-1)
        pmts      = 2.
        extpmts   = 14.
        lskg      = 1800.
        
        ### treat the NaI surface and the copper and teflon as 1 unit
        surf = 1.
        
        generated = float(sigs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        if loca == 'internal':
            fitActivity = sigs[key]['fitscale'] * (1./nmass) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif 'surf' in loca:
            fitActivity = sigs[key]['fitscale'] * (1./surf) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif loca == 'cucase':
            fitActivity = sigs[key]['fitscale'] * (1./surf) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif loca == 'teflon':
            fitActivity = sigs[key]['fitscale'] * (1./surf) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif loca == 'pmt':
            fitActivity = sigs[key]['fitscale'] * (1./pmts) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif loca == 'extpmt':
            fitActivity = sigs[key]['fitscale'] * (1./extpmts) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif loca == 'lsveto':
            fitActivity = sigs[key]['fitscale'] * (1./lskg) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif loca == 'lsvetoair':
            fitActivity = sigs[key]['fitscale'] * (1./nmass) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif loca == 'airshield':
            fitActivity = sigs[key]['fitscale'] * (1./nmass) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif loca == 'steel':
            fitActivity = sigs[key]['fitscale'] * (1./nmass) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        else:
            print "WARNING: No signal scaling for  --> ", loca
            continue

    return sigs


def buildData70(info, data):
    
    base = baseDir()
    runtime = 0.
    
    if info['reuse']:
        for c in info['chans']:
            info['chan'] = c
            
            if not info['rootfile']:
                print 'ERROR: no rootfile specified in backgrounds file'
                sys.exit()

            rootfile = base+info['rootfile']
            if not os.path.exists(rootfile):
                print 'ERROR: rootfile not found -->', rootfile
                sys.exit()

            rfile = TFile(rootfile, "READ")

            for e in range(2):
                key  = info['key']
                key += '-c'+str(c)
                key += '-e'+str(e)
                #print key
                try:
                    data[key] = {}
                    data[key]['info'] = info
                    data[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                    data[key]['pars'] = getPars(data[key]['hist'])
                    #data[key]['hist'].Sumw2()
                except:
                    print "ERROR: could not find hist -->",key
                    sys.exit()
                    
                try:
                    data[key]['subruns_hist'] = deepcopy(TH1F(rfile.Get(key+'_subruns')))
                    data[key]['subruns'] = data[key]['subruns_hist'].GetBinContent(1)
                except:
                    print "ERROR: could not find subruns_hist -->",key+'_subruns'
                    sys.exit()
                    
                try:
                    data[key]['runtime_hist'] = deepcopy(TH1F(rfile.Get(key+'_runtime')))
                    data[key]['runtime'] = data[key]['runtime_hist'].GetBinContent(1)
                    runtime = data[key]['runtime']
                except:
                    print "ERROR: could not find runtime_hist -->",key+'_runtime'
                    sys.exit()

    else:
        
        runfile = base+info['file']
        #print 'looking for -->',runfile
        if not os.path.exists(runfile):
            print 'ERROR: runfile not found -->', runfile
            sys.exit()

        runtime = 0.
        nfiles = 0
        chain = TChain("ntp","")
        for line in readFile(runfile):
            if not onCup(): fpath = '/home/mkauer/COSINE/CUP/mc-fitting'+line
            else: fpath = line
            #print 'looking for -->',fpath
            if not os.path.exists(fpath):
                continue
            else:
                #print 'adding -->',fpath
                runtime += getDuration(fpath)
                nfiles += chain.Add(fpath)

        if nfiles > 0:
            print 'INFO:',nfiles,'files found for data',info['tag']
            for c in info['chans']:
                info['chan'] = c
                for e in range(2):

                    # DEFINE HIST AND KEY
                    #-----------------------------------------------------------------------
                    i = info['xstl'] - 1
                    
                    key  = info['key']
                    key += '-c'+str(c)
                    key += '-e'+str(e)
                    
                    data[key] = {}
                    data[key]['info'] = info
                    
                    pars = histparam64(e)
                    data[key]['pars'] = pars
                    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
                    #-----------------------------------------------------------------------
                    
                    
                    # CALIBRATION
                    #-----------------------------------------------------------------------
                    ### calibrate the lo/hi energy data
                    # old calib
                    #edep, selection = calib60a(i,e)
                    # new calib
                    #edep, selection = calib60b(i,e)
                    # newer calib
                    edep, selection = calib64(i,e)
                    #-----------------------------------------------------------------------
                    
                    
                    # DEFINE CUTS
                    #-----------------------------------------------------------------------
                    if (i+1 == 7) and (e == 0):
                        masterCut = cutsBDT(i,c,e)
                    else:
                        masterCut = cutsOld(i,c,e)
                    #-----------------------------------------------------------------------
                    
                    
                    # FILL HISTOS
                    #-----------------------------------------------------------------------
                    chain.Draw(selection+' >> '+key, masterCut)
                    
                    #histo.SetLineColor(kBlack)
                    #histo.SetMarkerColor(kBlack)
                    #histo.SetLineWidth(1)
                    
                    data[key]['hist'] = histo
                    #data[key]['hist'].Sumw2()
                    
                    subruns = nfiles
                    key2 = key+'_subruns'
                    temp2 = TH1F(key2,'subruns',1,0,1)
                    temp2.SetBinContent(1, subruns)
                    data[key]['subruns_hist'] = temp2
                    data[key]['subruns'] = subruns
                    
                    #runtime = subruns*subrunTime*60.*60.
                    key3 = key+'_runtime'
                    temp3 = TH1F(key3,'runtime',1,0,1)
                    temp3.SetBinContent(1, runtime)
                    data[key]['runtime_hist'] = temp3
                    data[key]['runtime'] = runtime
                    #-----------------------------------------------------------------------
                    
        else:
            print 'ERROR: No data files found... quitting...'
            sys.exit()

    return data, runtime


def updateBkgsFile70(bkgsfile, resultsfile, newbkgs, BF='BR'):

    for thisfile in [bkgsfile, resultsfile]:
        if not os.path.exists(thisfile):
            print 'WARNING: file not found -->', thisfile
            return
    
    with open(bkgsfile) as fbkgs:
        bkgslines = fbkgs.read().splitlines()
    fbkgs.close()

    with open(resultsfile) as ffits:
        fitlines = ffits.read().splitlines()
    ffits.close()

    output = open(newbkgs, 'w')

    print ''
    print 'INFO: Updating bkgsfile -->',bkgsfile
    print '      To a new bkgsfile -->',newbkgs
    print ''

    skip = 0
    for bline in bkgslines:

        #bline = bline.strip()
        if not bline:
            output.write('\n')
            continue
        
        if bline.startswith('\"\"\"'):
            if skip == 0: skip = 1
            else: skip = 0
            output.write(bline+'\n')
            continue
        if skip:
            output.write(bline+'\n')
            continue
        
        if bline.startswith('#'):
            if 'version' in bline:
                output.write('# NEW GENERATED backgrounds file from fit!\n\n')
                output.write(bline+'\n')
            else:
                output.write(bline+'\n')
            continue
        
        #bbits = bline.split()
        bbits = filter(None, re.split("[ \s\t\n\r,:]+", bline.strip()))
        if 'F' not in bbits[0]:
            for bits in bbits:
                output.write(bits+'\t')
            output.write('\n')
            continue
        
        replaced = 0
        for fline in fitlines:
            if not replaced:
                #fline = fline.strip()
                if not fline:
                    continue
                #fbits = fline.split()
                fbits = filter(None, re.split("[ \s\t\n\r,:]+", fline.strip()))

                ### ------------------------------------------------------------
                ### if output format changes, this is generally the part that fails...
                #if fbits[-1] == 'mBq' or fbits[-3] == 'mBq':
                if fbits[0] == 'fit':
                    xstal = fbits[1].split('-')[0].split('x')[1]
                    loca  = fbits[1].split('-')[1]
                    chst  = fbits[1].split('-')[2].split('_')[0]
                    chsp  = fbits[1].split('-')[2].split('_')[1]
                    acti  = str(fbits[3])

                    lenbbits = len(bbits)
                    if bbits[2] == xstal and bbits[3].replace('-','') == loca and bbits[5] == chst and bbits[6] == chsp:
                        for i in range(lenbbits):
                            if i == 0:
                                output.write(BF+'\t')
                            elif i == 7:
                                if acti != '0.0': output.write(acti+'\t')
                                else: output.write(bbits[i]+'\t')
                            elif i==9:
                                output.write('0.1\t')
                            elif i==10:
                                output.write('10\t')
                            elif i==(lenbbits-1):
                                output.write(bbits[i])
                            else:
                                output.write(bbits[i]+'\t')
                        output.write('\n')
                        replaced = 1
        
        if not replaced:
            output.write(bline+'\n')
            print '\nWARNING: Could not match line -->',fline.split('\n')[0]
            print   '         To line -->',bline.split('\n')[0]
    
    output.close()
    return


def setBinError(histo):
    for n in range(histo.GetNbinsX()+1):
        histo.SetBinError(n, sqrt(histo.GetBinContent(n)))
    return histo


def zeroBinError(histo):
    for n in range(histo.GetNbinsX()+1):
        histo.SetBinError(n, 0.0)
    return histo


def uniqString(string):
    newstring=''
    for char in string:
        if char not in newstring:
            newstring += char
    return newstring

