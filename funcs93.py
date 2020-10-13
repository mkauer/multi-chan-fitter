#!/usr/bin/env python
######################################################################
# funcs93.py
# 
# Tweaking MC locations and calibrations
# 
# Works with v93 and later versions
# 
# version: 2020-03-30
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ changed high energy bin range in makePlots93() to 4000

# + added special case for external Tl208 gammas to groupNum93() in funcs93.py
# + add gamma to makePlots93() exceptions in funcs93.py
# + add plastic to makePlots93() exceptions
# + add U235 groupNo 41,42 to groupNum93()
# + add U235 to groupNum93()
# + add I129 to groupNum93()
# - do not init THF1 in makePlots93()
# + add steel to makePlots93() exceptions
# ~ tweaked makePlots93()
# ~ was using groupNo 1 for K40! changed to 31 in groupNum93()
# ~ fixed the H3 normalization - it was still wrong!
# + added calib93() and buildData93()
# + add path to Sn113 MC and groupNo for Sn113
# ~ now should have cuts for generated events correct!
# + add new groupNum93() to test out efficiency cuts
# + tweaked functions for copper Co60 MC from Pushpa
# ~ need to figure out difference between evt_type cut and primParticleName cut
# ~ something is wonky with the primParticleName cut
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
#from funcs92 import *
from funcs_misc import *


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
        if others and ('pmt' in info['key'] or 'internal' in info['key'] or 'copper' in info['key']):
        #if others and ('lsveto' not in info['key'] and 'innersteel' not in info['key'] and 'data' not in info['key']):
            key = info['key']
            for i in range(8):
                #info['key'] = key+'-f'+str(i+1)
                newinfo = deepcopy(info)
                newinfo['key'] = key+'-f'+str(i+1)
                newinfo['from'] = str(i+1)
                infos.append(newinfo)
        elif not others and ('pmt' in info['key'] or 'internal' in info['key'] or 'copper' in info['key']):
        #elif not others and ('lsveto' not in info['key'] and 'innersteel' not in info['key'] and 'data' not in info['key']):
            key = info['key']
            newinfo = deepcopy(info)
            newinfo['key'] = key+'-f'+str(info['xstl'])
            newinfo['from'] = info['xstl']
            infos.append(newinfo)
        # force these to not grab "others" until fixed
        #elif ('surf' in info['key'] or 'teflon' in info['key'] or 'copper' in info['key'] or 'case' in info['key']):
        elif ('surf' in info['key'] or 'teflon' in info['key'] or 'case' in info['key']):
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
                data, runtime = buildData93(info, data)
            
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

        ### many special cases for Pushpa's MC locations
        if (info['isof'] == 'Cd109' or info['isof'] == 'H3') and location == 'internal':
            path2 = '*'+info['isof']+'*.root'
        if info['isof'] == 'Co60' and location == 'copper-case':
            path2 = '*'+'CuCase-'+info['isof']+'*.root'
        if info['isof'] == 'Co60' and location == 'copper':
            path1 = '/data/MC/COSINE/V3.1.1/reprocessed/Co60'
            path2 = '*'+'Copper-'+info['isof']+'*.root'
        if info['isof'] == 'Sn113' and location == 'internal':
            path1 = '/data/MC/COSINE/V3.1.1/reprocessed'
            path2 = '*'+'internal-'+info['isof']+'*.root'
        
        
        if not onCup():
            path1 = '/home/mkauer/COSINE/CUP/mc-fitting'+path1
        
        #print 'INFO: looking for files with -->', os.path.join(path1, path2)

        pushpasMC = 0
        if location in pushpas: #or info['isof'] in ['Cd109', 'H3']:
            pushpasMC = 1
        
        chain  = TChain("MC","")
        nfiles = 0
        nfiles = chain.Add(os.path.join(path1, path2))
        
        if nfiles > 0:
            
            #print 'INFO:',nfiles,'files found for', info['loca'], info['isof']
            
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
                    
                    elif info['loca'] == 'teflon':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Teflon")')

                    elif info['loca'] == 'teflonbulk':
                        # Pushpas C7 Surface MC
                        volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                    
                    elif 'teflonsurf' in info['loca']:
                        # Estellas MC
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Teflon")')
                        # Pushpas C7 Surface MC
                        if pushpasMC:
                            volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                    
                    elif 'cusurf' in info['loca']:
                        # Estellas MC
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Case")')
                    
                    elif info['loca'] == 'cucase':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Case")')

                    elif info['loca'] == 'coppercase':
                        # Pushpas MC
                        volumeCut = TCut('((primVolumeName == "NaIDet07Fringe0")'
                                     +' || (primVolumeName == "NaIDet07Fringe1")'
                                     +' || (primVolumeName == "NaIDet07CalibHole"))')
                    
                    elif info['loca'] == 'copper':
                        # Pushpas MC
                        volumeCut = TCut('((primVolumeName == "NaIDet0'+str(info['from'])+'Fringe0")'
                                     +' || (primVolumeName == "NaIDet0'+str(info['from'])+'Fringe1")'
                                     +' || (primVolumeName == "NaIDet0'+str(info['from'])+'Case"))')
                    
                    elif info['loca'] == 'pmt':
                        volumeCut = TCut('((primPMTid[0] == '+pmt1+') || (primPMTid[0] == '+pmt2+'))')
                    
                    elif info['loca'] == 'extpmt':
                        volumeCut = TCut('((primPMTid[0] != '+pmt1+') && (primPMTid[0] != '+pmt2+'))')
                    
                    elif info['loca'] == 'lsveto':
                        volumeCut = TCut('(primVolumeName == "lsveto")')
                    
                    elif info['loca'] == 'lsvetoair':
                        volumeCut = TCut('((primVolumeName == "DetPMTCover") || (primVolumeName == "DetPMTEnvelope"))')

                    elif info['loca'] == 'airshield':
                        volumeCut = TCut('(primVolumeName == "LSVetoAirRoom")')

                    elif info['loca'] == 'steel':
                        volumeCut = TCut('((primVolumeName == "SteelSupport") || (primVolumeName == "SteelSupportTop"))')

                    elif info['loca'] == 'innersteel':
                        volumeCut = TCut('(primVolumeName == "InnerSteel")')
                    
                    else:
                        print "WARNING: No selection criteria for  --> ", info['loca']
                        continue

                    
                    brokenChainCut = groupNum93(info)
                    eventTypeCut = TCut('(evt_Type > 10)')
                    generatedCuts = TCut('(evt_Type < 10)'+' && '+volumeCut.GetTitle())
                    
                    if info['chst'] == 'H3':
                        generatedCuts = TCut('(evt_Type > 10)'+' && '+volumeCut.GetTitle())
                    
                        
                    ### create a hist of the generated events
                    #---------------------------------------------------------------------
                    key2 = key+'_generated'
                    temp2 = TH1F(key2, 'generated', 1, 1, 2)
                    
                    chain.Draw('1 >> '+key2, generatedCuts)
                                        
                    mc[key]['generated_hist'] = temp2
                    generated = temp2.GetEntries()
                    mc[key]['generated'] = generated
                    #---------------------------------------------------------------------
                    
                    
                    ### using Box-Muller for resolution smearing method
                    ### setting "rng" as the smearing alias
                    chain.SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                    resolFunc = resol80(i,e)
                    if pushpasMC:
                        resolFunc = resol80(i,e,1)
                    chain.SetAlias('sigma', resolFunc)

                    
                    masterCut = TCut('('+
                                energyCut.GetTitle()+' && '+
                                eventTypeCut.GetTitle()+' && '+
                                brokenChainCut.GetTitle()+' && '+
                                chanCut.GetTitle()+' && '+
                                volumeCut.GetTitle()
                                     +')')
                    
                    
                    selection = '((edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng))'
                    if pushpasMC:
                        selection = '((edep[6]*1000.) + (sigma*edep[6]*1000.*rng))'
                    chain.Draw(selection+' >> '+key, masterCut)
                                        
                    detected = histo.GetEntries()

                    if (c=='S' and e==0) or (c=='M' and e==1):
                        print 'DEBUG:', key, 'generated events =', generated
                        print 'DEBUG:', key, 'detected events =', detected
                        print 'DEBUG:', key, 'efficiency =', round(100*detected/generated, 2), '%'
                    
                    
                    mc[key]['hist'] = histo
                    #mc[key]['hist'].Sumw2()
            print ''
                    
        else:
            print 'ERROR: no MC files found for -->', \
                'x'+str(info['xstl']), info['floca'], info['isof']
            sys.exit()
    
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
            
        elif loca == 'copper':
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
            print "ERROR: no signal scaling for -->", loca
            sys.exit()

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
            
        elif loca == 'copper':
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
            print "ERROR: no background scaling for -->", loca
            sys.exit()

    return bkgs


def makePlots93(bkgs, combine, others, vcut=0):
    
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
                    # do not do this 2018-10-23
                    #temp[i-1] = TH1F('temp','temp',4000,0,4000)
                    
                    #==========================================================
                    #==========================================================
                    # combine all normalized pmt and internal contributions
                    #combine = 1
                    newkey  = crystal
                    newkey += '-'+location
                    newkey += '-'+isotope
                    
                    # add 'steel' exception 2018-10-23
                    # add 'plastic' exception 2019-01-29
                    # add 'gamma' exception 2019-06-24
                    if location != 'lsveto' and location != 'innersteel' \
                      and location != 'steel' and location != 'plastic' \
                      and location != 'gamma':
                        newkey += '-f'+crystal[-1]
                    newkey += '-c'+C
                    newkey += '-e'+str(E)
                    try:
                        temp[i-1] = deepcopy(bkgs[newkey]['hist'])
                    except:
                        #print 'skipping key -->', newkey
                        pass
                    
                    #print newkey
                    
                    #if combine and (location == 'pmt' or location == 'internal'):
                    #if combine and (location != 'lsveto' and location != 'innersteel'):
                    if combine:
                        for k in range(1,10):
                            if k != int(crystal[-1]):
                                newerkey  = crystal
                                newerkey += '-'+location
                                newerkey += '-'+isotope
                                newerkey += '-f'+str(k)
                                newerkey += '-c'+C
                                newerkey += '-e'+str(E)
                                try:
                                    temp[i-1].Add(deepcopy(bkgs[newerkey]['hist']))
                                except:
                                    #print 'skipping key -->', newerkey
                                    pass
                    #==========================================================
                    #==========================================================
                    
                    bins = temp[i-1].GetNbinsX()
                    hmin = temp[i-1].GetBinLowEdge(1)
                    hmax = temp[i-1].GetBinLowEdge(bins+1)
                    kvpb = (hmax-hmin)/bins
                    
                    temp[i-1].SetLineColor(kBlack)
                    if E:
                        temp[i-1].Rebin(int(10/kvpb))
                        temp[i-1].SetAxisRange(0, 4000, 'x')
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


def groupNum93(info):
    
    # U238 group numbers
    # ------------------------
    # 11: U238  -> Th230
    # 12: Th230 -> Ra226
    # 13: Ra226 -> Rn222
    # 14: Rn222 -> Pb210
    # 15: Pb210 -> ground
    
    # Th232 group numbers
    # ------------------------
    # 21: Th232 -> Ra228
    # 22: Ra228 -> Th228
    # 23: Th228 -> ground
    
    # U235 group numbers
    # ------------------------
    #  0: U235 -> ground? (doesn't look like it gets splip up?)
    # 41: U235 -> Pa231
    # 42: Pa231 -> ground
    
    # others
    # ------------------------
    # 31: K40

    
    if info['chst'] in ['Te121',  'Te121m',
                        'Te123m', 'Te125m', 'Te127m',
                        'I125',   'I126',   'I129',
                        'Na22',   'H3',
                        'Cd109',  'Sn113',
                        'Co60']:
        return TCut('(groupNo == 0)')
    
    if info['chst'] == 'K40':
        return TCut('(groupNo == 31)')
    
    if   info['chst'] == 'U238':  start = 11
    elif info['chst'] == 'Th230': start = 12
    elif info['chst'] == 'Ra226': start = 13
    elif info['chst'] == 'Rn222': start = 14
    elif info['chst'] == 'Pb210': start = 15
    elif info['chst'] == 'Th232': start = 21
    elif info['chst'] == 'Ra228': start = 22
    elif info['chst'] == 'Th228': start = 23
    elif info['chst'] == 'U235':  start = 41
    elif info['chst'] == 'Pa231': start = 42
    ### special case for external Tl208 gammas
    elif info['chst'] == 'Tl208': start = 0
    else: start = -1
    
    #if   info['chsp'] == 'U238':  stop = 11 # should not be a stop group
    if   info['chsp'] == 'Th230': stop = 12
    elif info['chsp'] == 'Ra226': stop = 13
    elif info['chsp'] == 'Rn222': stop = 14
    elif info['chsp'] == 'Pb210': stop = 15
    #elif info['chsp'] == 'Th232': stop = 21 # should not be a stop group
    elif info['chsp'] == 'Ra228': stop = 22
    elif info['chsp'] == 'Th228': stop = 23
    #elif info['chsp'] == 'U235':  stop = 41 # should not be a stop group
    elif info['chsp'] == 'Pa231': stop = 42
    ### handle the GRND group number better...
    #elif info['chsp'] == 'GRND':  stop = 16
    #elif info['chsp'] == 'GRND':  stop = 24
    #elif info['chsp'] == 'GRND':  stop = 43
    elif info['chsp'] == 'GRND':
        if   start in [11,12,13,14,15]: stop = 16
        elif start in [21,22,23]:       stop = 24
        elif start in [41,42]:          stop = 43
        ### special case for external Tl208 gammas
        elif start in [0]:              stop = 43
        else: stop = -1
    else: stop = -1

    if start != -1 and stop != -1:
        return TCut('((groupNo >= '+str(start)+') && (groupNo < '+str(stop)+'))')
    else:
        print 'ERROR: groupNo not found for -->', info['chst']
        sys.exit()


def buildData93(info, data):
    
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
                tmpruntime = getDuration(fpath)
                if tmpruntime > 0:
                    runtime += tmpruntime
                    nfiles += chain.Add(fpath)
                else:
                    print 'WARNING: no run duration found for',fpath
            
        if nfiles > 0:
            print 'INFO:',nfiles,'files found for data',info['tag']
            for C in info['chans']:
                info['chan'] = C
                for E in range(2):

                    # DEFINE HIST AND KEY
                    #-----------------------------------------------------------------------
                    i = info['xstl'] - 1
                    
                    key  = info['key']
                    key += '-c'+str(C)
                    key += '-e'+str(E)
                    
                    data[key] = {}
                    data[key]['info'] = info
                    
                    pars = histparam64(E)
                    data[key]['pars'] = pars
                    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
                    #-----------------------------------------------------------------------
                    
                    
                    # CALIBRATION
                    #-----------------------------------------------------------------------
                    edep, selection = calib93(i, E)
                    #-----------------------------------------------------------------------
                    
                    
                    # DEFINE CUTS
                    #-----------------------------------------------------------------------
                    masterCut = cutsBDT92(i, C, E, edep, selection)
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


def calib93(i,E):
    
    hiE = [[1.0,  -6],
           [1.0,  -6],
           [1.0,  -8],
           [1.0,  -8],
           [1.0, -20],
           [1.0,  -5],
           [1.0,  -8],
           [0.7,   0]]
    
    if E:
        edep = '(crystal'+str(i+1)+'.energyD)'
        selection = '(('+edep+'*'+str(hiE[i][0])+')+'+str(hiE[i][1])+')'
    else:
        edep = '(crystal'+str(i+1)+'.energy)'
        selection = edep
    
    return edep, selection


#=====================================================================
#  Other old functions that are needed...
#=====================================================================


def getInfo91(line, freuse=0, fchans=0, fxstals=[]):
    
    info={}
    
    info['line'] = line
    
    bits = line.split()
    ### probably don't need to use 're' since split() is working fine
    #import re
    #bits = filter(None, re.split("[ \s\t]+", line.strip()))
    
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
    info['chans'] = str(bits[1])
    # force a channel(s) of data selection
    if fchans:
        info['chans'] = fchans
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # crystal
    info['xstl'] = int(bits[2])
    # force only a particular crystal(s) and skip the others
    if len(fxstals) > 0 and info['xstl'] not in fxstals:
        return []
    # what crystal is the background from?
    # set to xstl by default - can be updated later
    info['from'] = int(bits[2])
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # reuse rootfile name
    r=0
    if 'root' in str(bits[-1]):
        r=-1
        temp = str(bits[-1])
        for i in range(2):
            if '.' in temp[0] or '/' in temp[0]:
                temp = temp[1:]
        info['rootfile'] = temp
    else:
        info['rootfile'] = 0
    #-----------------------------------------------------------------
    

    if 'D' in info['type']:
        
        ### data file
        temp = str(bits[3])
        for i in range(2):
            if '.' in temp[0] or '/' in temp[0]:
                temp = temp[1:]
        info['file'] = temp
        
        ### data tag
        info['tag'] = str(bits[4])

        ### processing version
        #info['build'] = str(bits[5])

        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-data'
        key += '-'+str(info['tag'])
        #key += '-c'+info['chan']
        info['key'] = key
        
    else:
        
        ### background location help
        info['floca'] = str(bits[3])
        # without any '-' for hist key help
        info['loca'] = str(bits[3]).replace('-','')
        
        ### top level isotope file name
        info['isof'] = str(bits[4])

        ### isotope chain break start
        info['chst'] = str(bits[5])

        ### isotope chain break stop
        info['chsp'] = str(bits[6])

        ### activity
        info['acti'] = float(bits[7])

        ### error on activity
        ### do not need this anymore
        #info['erro'] = float(bits[8])

        ### fit bounds
        #info['fbnd'] = [float(bits[9]), float(bits[10])]
        info['fbnd'] = [float(bits[r-2]), float(bits[r-1])]
        
        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-'+info['loca']
        #if info['chst'] == info['chsp']: key += '-'+info['chst']
        #if info['chsp'] == 'GRND': key += '-'+info['chst']
        #else: key += '-'+info['chst']+'_'+info['chsp']
        key += '-'+info['chst']+'_'+info['chsp']
        info['key'] = key

        ### assign a plotting group to the isotope-location
        info['group'] = setGroup(info)
        
    return info


def sortDataKeys92(data):
    datkeys=[]
    for key in data:
        datkeys.append(key)
    datkeys.sort()
    return datkeys


def sortSimKeys92(sigs):
    sigkeys=[]
    """
    delete=[]
    for key in sigs:
        if sigs[key]['hist'].Integral() > 0:
            sigkeys.append(key)
        else: delete.append(key)
    for key in delete:
        print 'INFO: deleting sigs key', key
        del sigs[key]
    """
    for key in sigs:
        sigkeys.append(key)
    sigkeys.sort()
    return sigs, sigkeys


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
        if i==8:
            xkgs = 1800.
            # NO! Why did I do that?
            #xkgs = 1900. # 2019-04-16 increased mass a little 
        keVperBin = 1./float(data[key]['pars'][3])
        
        #print 'DEBUG: data key =', key
        #print 'DEBUG: time in days =', days
        #print 'DEBUG: mass in kg =', xkgs
        #print 'DEBUG: keV/bin =', keVperBin
        #print 'DEBUG: total events =', data[key]['hist'].GetEntries()
        #print ''
        
        scale = float(1./(days*xkgs*keVperBin))
        data[key]['hist'].Scale(scale)
        data[key]['druScale'] = scale

    return data


def combineOthers92(sigs, globalMC):
    donekeys=[]
    delete=[]
    #for key in sigkeys:
    for key in sigs:
        if key in donekeys: continue
        bits = key.split('-')
        if len(bits) != 6: continue
        # don't combine if combining in global fit anyway?
        if bits[1] in globalMC: continue
        for F in range(1,9):
            X = int(bits[0][1])
            if F == X: continue
            try:
                default = 'x'+str(X)+'-'+bits[1]+'-'+bits[2]+'-f'+str(X)+'-'+bits[4]+'-'+bits[5]
                newkey  = 'x'+str(X)+'-'+bits[1]+'-'+bits[2]+'-f'+str(F)+'-'+bits[4]+'-'+bits[5]
                #print '!!! adding ', newkey, ' to ', default
                sigs[default]['hist'].Add(sigs[newkey]['hist'])
                donekeys.append(default)
                donekeys.append(newkey)
                delete.append(newkey)
            except:
                pass
    
    # delete this histograms?
    for key in delete:
        #print '!!! deleting -fx key ', key
        del sigs[key]
    
    return sigs


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
                            #elif i==9:
                            elif i==(lenbbits-3):
                                output.write('0.1\t')
                            #elif i==10:
                            elif i==(lenbbits-2):
                                output.write('10\t')
                            elif i==(lenbbits-1):
                                output.write(bbits[i])
                            else:
                                output.write(bbits[i]+'\t')
                        output.write('\n')
                        replaced = 1
        
        if not replaced:
            output.write(bline+'\n')
            #print 'WARNING: Could not match -->', filter(None, re.split("[ \s\t\n\r,:]+", bline.strip()))[3:7]
            
    output.close()
    return


def resol80(i, E=0, pushpaC7Surf=0):
    """
    https://cupwiki.ibs.re.kr/Kims/NaICalibration
    """
    # res = p[0]/sqrt(x) + p[1]
    # low energy
    loEresol = [[0.2413,  0.01799],
	        [0.2951,  0.01427],
	        [0.3106,  0.007894],
	        [0.3894, -0.001437],
                [0.8, 0], # tweaking C5
	        [0.3620,  0.0006355],
	        [0.3042,  0.009784],
                [1.3, 0]] # tweaking C8
    # high energy
    hiEresol = [[0.6729, 0.009374],
	        [0.6531, 0.006627],
	        [0.5926, 0.009506],
	        [0.7227, 0.004790],
                [1.7, 0], # tweaking C5
	        [0.6498, 0.009670],
	        [0.7034, 0.007812],
                [4.0, 0]] # tweaking C8
    if E:
        p0, p1 = hiEresol[int(i)]
    else:
        p0, p1 = loEresol[int(i)]

    selection = '(('+str(float(p0))+'/sqrt(edep['+str(i)+']*1000.)) + '+str(float(p1))+')'
    if pushpaC7Surf:
        selection = '(('+str(float(p0))+'/sqrt(edep[6]*1000.)) + '+str(float(p1))+')'
    return selection

