#!/usr/bin/env python
######################################################################
# funcs93.py
# 
# Tweaking MC locations and calibrations
# 
# Works with v93 and later versions
# 
# version: 2018-04-04
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ do to surface MC issues, switch back to just pmt and internal "others"
# + makePlots93() - "others" for everything except lsveto, innersteel, data
# ~ do "others" for everything except for lsveto, innersteel and data
# + add copper to others in build93()
# + scaleBkgs93()
# + scaleSigs93()
# + buildMC93()
# + build93()
# ~ import funcs92
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
from funcs92 import *


def build93(infile='backgrounds900.txt', others=1, vcut=1, freuse=0, fchans=0, fxstals=[]):

    data = {}
    bkgs = {}
    sigs = {}
    runtime = 0
    
    for line in readFile(infile):
        info = getInfo91(line, freuse, fchans, fxstals)
        if len(info) == 0:
            continue
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++
        ### This isn't correct because internals, teflon, copper
        ###    scale by mass or surface area etc. But this
        ###    is okay for now...
        infos=[]
        if others and ('pmt' in info['key'] or 'internal' in info['key']):
        #if others and ('lsveto' not in info['key'] and 'innersteel' not in info['key'] and 'data' not in info['key']):
            key = info['key']
            for i in range(8):
                #info['key'] = key+'-f'+str(i+1)
                newinfo = deepcopy(info)
                newinfo['key'] = key+'-f'+str(i+1)
                newinfo['from'] = str(i+1)
                infos.append(newinfo)
        elif not others and ('pmt' in info['key'] or 'internal' in info['key']):
        #elif not others and ('lsveto' not in info['key'] and 'innersteel' not in info['key'] and 'data' not in info['key']):
            key = info['key']
            newinfo = deepcopy(info)
            newinfo['key'] = key+'-f'+str(info['xstl'])
            newinfo['from'] = info['xstl']
            infos.append(newinfo)
        # force these to not grab "others" until fixed
        elif ('surf' in info['key'] or 'teflon' in info['key'] or 'copper' in info['key'] or 'case' in info['key']):
            key = info['key']
            newinfo = deepcopy(info)
            newinfo['key'] = key+'-f'+str(info['xstl'])
            newinfo['from'] = info['xstl']
            infos.append(newinfo)
        else:
            infos.append(info)
        #+++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++
        
        for info in infos:
            #print info['key']
            if 'D' in info['type']:
                data, runtime = buildData92(info, data)
            
            elif 'B' in info['type']:
                bkgs = buildMC93(info, bkgs, vcut)
            
            elif 'F' in info['type']:
                sigs = buildMC93(info, sigs, vcut)
            
            else:
                print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
                print info['line']
                continue
    
    return data, bkgs, sigs, runtime


def buildMC93(info, mc, vcut):

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
        if info['isof'] == 'H3':
            path1 = '/data/COSINE/WORK/pushpa/sim/process/Crystal7/'
        if info['isof'] == 'Cd109':
            path1 = '/data/MC/COSINE/V3.1.1/reprocessed/G4_10/'
        
        location = info['floca']
        if location == 'extpmt':
            location = 'pmt'

        ### 2nd level path to the specific files
        path2 = info['isof']+'/set2/'+'*'+location+info['isof']+'*root'
        
        ### Estella's Surface Pb210 MC
        if location == 'internal-surf-10um':
            path2 = info['isof']+'/set2/surf/10um/'+'*'+location+'-C'+str(info['xstl'])+'*root'
        if location == 'cu-surf-10um':
            path2 = info['isof']+'/set2/surf/cu-10um/'+'*'+location+'-C'+str(info['xstl'])+'*root'
        if location == 'teflon-surf-10um':
            path2 = info['isof']+'/set2/surf/teflon-10um/'+'*'+location+'-C'+str(info['xstl'])+'*root'
        
        ### Pushpa's MC has different root variables
        pushpas = ['nai-surf-10um', 'teflon-surf-2um', 'teflon-bulk', 'copper-case']
        if location in pushpas:
            path1 = '/data/COSINE/WORK/pushpa/sim/process/Crystal7/ResolA_ResolD/'
            if location == 'nai-surf-10um':
                path2 = 'anal_lsvetofull-NaI-surface_C7-10um-Pb210-*.root'
            if location == 'teflon-surf-2um':
                path2 = 'anal_lsvetofull-Teflon-surface_2um-Pb210-*.root'
            if location == 'teflon-bulk':
                path2 = 'anal_lsvetofull-Teflon-surface_all-Pb210-*.root'
            if location == 'copper-case':
                path1 = '/data/MC/COSINE/V3.1.1/reprocessed/'
        
        if (info['isof'] == 'Cd109' or info['isof'] == 'H3') and location == 'internal':
            path2 = '*'+info['isof']+'*.root'
        if info['isof'] == 'Co60' and location == 'copper-case':
            path2 = '*'+'CuCase-'+info['isof']+'*.root'
        
        if not onCup():
            path1 = '/home/mkauer/COSINE/CUP/mc-fitting'+path1
        
        print 'INFO: looking for files with -->', os.path.join(path1, path2)

        pushpasMC = 0
        if location in pushpas: #or info['isof'] in ['Cd109', 'H3']:
            pushpasMC = 1
        
        chain  = TChain("MC","")
        nfiles = 0
        nfiles = chain.Add(os.path.join(path1, path2))
        
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
                    if pushpasMC:
                        energyCut = TCut('(edep[6]*1000. > 0.0)')

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
                        if pushpasMC:
                            hitCut = '((singleHitTag[6] == 1.0))'
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
                        if pushpasMC:
                            hitCut = '((multipleHitTag[6] == 1.0))'
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
                    
                    
                    #pmt1 = str((int(info['xstl'])*2)-2)
                    #pmt2 = str((int(info['xstl'])*2)-1)
                    pmt1 = str((int(info['from'])*2)-2)
                    pmt2 = str((int(info['from'])*2)-1)
                    #print pmt1,pmt2
                    
                    volumeCut = TCut('(1)')
                    if info['loca'] == 'internal':
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Crystal")')
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal")')
                        
                    elif 'internalsurf' in info['loca']:
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Crystal")')
                    
                    elif 'naisurf' in info['loca']:
                        # Pushpas C7 Surface MC
                        volumeCut = TCut('(primVolumeName == "NaIDet07Crystal")')
                    
                    elif 'cusurf' in info['loca']:
                        # Estellas MC
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Case")')
                    
                    elif 'teflonsurf' in info['loca']:
                        # Estellas MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Teflon")')
                        # Pushpas C7 Surface MC
                        volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                    
                    elif info['loca'] == 'cucase':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Case")')

                    elif info['loca'] == 'coppercase':
                        # Pushpas MC
                        volumeCut = TCut('((primVolumeName == "NaIDet07Fringe0") || (primVolumeName == "NaIDet07Fringe1") || (primVolumeName == "NaIDet07CalibHole"))')
                    
                    elif info['loca'] == 'teflon':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Teflon")')

                    elif info['loca'] == 'teflonbulk':
                        # Pushpas C7 Surface MC
                        volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                    
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
                        volumeCut = TCut('((primVolumeName == "SteelSupport") || (primVolumeName == "SteelSupportTop"))')
                        #volumeCut = TCut('((primVolumeName == "SteelSupport"))')
                        #volumeCut = TCut('((primVolumeName == "SteelSupportTop"))')
                        #volumeCut = TCut('((primVolumeName == "Welding"))')
                        #volumeCut = TCut('((primVolumeName != ""))')

                    elif info['loca'] == 'innersteel':
                        volumeCut = TCut('(primVolumeName == "InnerSteel")')
                    
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
                    #resolFunc = resol64(i,e)
                    resolFunc = resol80(i,e)
                    if pushpasMC:
                        resolFunc = resol80(i,e,1)
                    chain.SetAlias('sigma', resolFunc)

                    
                    #=====================================================================
                    #=====================================================================
                    # use the volume cut?
                    #vcut = 1
                    if vcut:
                        masterCut = TCut('('+
                                     energyCut.GetTitle()+' && '+
                                     eventTypeCut.GetTitle()+' && '+
                                     brokenChainCut.GetTitle()+' && '+
                                     chanCut.GetTitle()+' && '+
                                     volumeCut.GetTitle()
                                     +')')
                    else:
                        masterCut = TCut('('+
                                     energyCut.GetTitle()+' && '+
                                     eventTypeCut.GetTitle()+' && '+
                                     brokenChainCut.GetTitle()+' && '+
                                     chanCut.GetTitle()
                                     +')')
                    #=====================================================================
                    #=====================================================================
                    
                    
                    selection = '((edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng))'
                    if pushpasMC:
                        selection = '((edep[6]*1000.) + (sigma*edep[6]*1000.*rng))'
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


def scaleSigs93(sigkeys, sigs, runtime=0):

    for key in sigkeys:

        try:
            test = sigs[key]['fitscale']
        except:
            print 'WARNING: no fitscale for', key
            sigs[key]['info']['fitacti'] = 0
            sigs[key]['info']['fiterro'] = 0
            continue
        
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

        nmass      = cmass(int(x[-1])-1)
        pmts       = 2.
        extpmts    = 14.
        lskg       = 1800.
        steel      = 1600.
        innersteel = 4000.
        
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
            
        elif 'teflon' in loca:
            fitActivity = sigs[key]['fitscale'] * (1./surf) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        elif loca == 'cucase' or loca == 'coppercase':
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
            fitActivity = sigs[key]['fitscale'] * (1./steel) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']

        elif loca == 'innersteel':
            fitActivity = sigs[key]['fitscale'] * (1./innersteel) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['fitacti'] = fitActivity
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
        else:
            print "WARNING: No signal scaling for  --> ", loca
            continue

    return sigs


def scaleBkgs93(bkgs,runtime=0):
    
    for key in bkgs:
        
        #print 'scaling -->',key
        
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

        nmass      = cmass(int(x[-1])-1)
        pmts       = 2.
        extpmts    = 14.
        lskg       = 1800.
        steel      = 1600.
        innersteel = 4000.
        
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
            
        elif 'teflon' in loca:
            scale = bkgs[key]['info']['acti'] * (surf) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        elif loca == 'cucase' or loca == 'coppercase':
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
            scale = bkgs[key]['info']['acti'] * (steel) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale

        elif loca == 'innersteel':
            scale = bkgs[key]['info']['acti'] * (innersteel) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
            
        else:
            print "WARNING: No background scaling for  --> ", loca
            continue

    return bkgs


def makePlots93(bkgs, combine, others, vcut):
    
    gStyle.SetPadTopMargin    (0.07)
    gStyle.SetPadBottomMargin (0.11)
    gStyle.SetPadLeftMargin   (0.12)
    gStyle.SetPadRightMargin  (0.02)
    
    combos=[]
    for key in bkgs:
        
        crystal = key.split('-')[0]
        location = key.split('-')[1]
        isotope = key.split('-')[2]
        
        combo = crystal+location+isotope
        if combo not in combos:
            
            combos.append(combo)
            
            #canvas = TCanvas('canvas','canvas',0,0,1400,900)
            canvas = TCanvas('canvas','canvas',0,0,1200,900)
            canvas.Divide(2,2)

            i=1
            temp = [0 for k in range(4)]
            for C in ['S','M']:
                for E in [0, 1]:
                    
                    canvas.cd(i)
                    
                    #==========================================================
                    #==========================================================
                    # combine all normalized pmt and internal contributions
                    #combine = 1
                    newkey  = crystal
                    newkey += '-'+location
                    newkey += '-'+isotope
                    #if location == 'pmt' or location == 'internal':
                    if location != 'lsveto' and location != 'innersteel':
                        newkey += '-f'+crystal[-1]
                    newkey += '-c'+C
                    newkey += '-e'+str(E)
                    temp[i-1] = deepcopy(bkgs[newkey]['hist'])
                    
                    if combine and (location == 'pmt' or location == 'internal'):
                    #if combine and (location != 'lsveto' and location != 'innersteel'):
                        for k in range(1,9):
                            if k != int(crystal[-1]):
                                newkey  = crystal
                                newkey += '-'+location
                                newkey += '-'+isotope
                                newkey += '-f'+str(k)
                                newkey += '-c'+C
                                newkey += '-e'+str(E)
                                #try:
                                temp[i-1].Add(deepcopy(bkgs[newkey]['hist']))
                                #except:
                                #    pass
                    #==========================================================
                    #==========================================================
                    
                    bins = temp[i-1].GetNbinsX()
                    hmin = temp[i-1].GetBinLowEdge(1)
                    hmax = temp[i-1].GetBinLowEdge(bins+1)
                    kvpb = (hmax-hmin)/bins
                    
                    if E:
                        temp[i-1].Rebin(int(10/kvpb))
                        temp[i-1].SetAxisRange(100, 3000, 'x')
                        #temp[i-1].SetAxisRange(1e-3, 1, 'y')
                    else:
                        temp[i-1].Rebin(int(1/kvpb))
                        temp[i-1].SetAxisRange(0, 100, 'x')
                        #temp[i-1].SetAxisRange(1e-2, 10, 'y')
                        
                    what=''
                    if C == 'S': what += 'Single-hit'
                    if C == 'M': what += 'Multi-hit'
                    if E: what += ' Hi-energy'
                    else: what += ' Lo-energy'
                    
                    temp[i-1].SetTitle(crystal+' '+location+' '+isotope+'   '+what)
                    temp[i-1].GetXaxis().SetTitle('energy (keV)')
                    temp[i-1].GetYaxis().SetTitle('arb. counts')
                    temp[i-1].Draw()
                    canvas.cd(i).SetLogy(1)
                    i+=1
                    
                    canvas.Update()
            
            save  = crystal+'-'+location+'-'+isotope
            save += '-othr'+str(others)
            save += '-vnct'+str(vcut)
            save += '-cmbd'+str(combine)
            canvas.Print('./plots-isotopes/'+save+'.png')
            #raw_input()
            del canvas
            del temp

