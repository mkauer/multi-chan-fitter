#!/usr/bin/env python
######################################################################
# funcs101.py
# 
# Remove LSveto low-energy from the histograms
# 
# version: 2019-01-23
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ change back to NOT account for surface area - use mBq/kg NaI mass
# ~ change scaleBkgs101() and scaleSigs101() to account for surface area
# + add surfArea() calculation for surface components
# + add plastic to build101() exclude list
# + add plastic to buildMC101()
# ~ move smoothXbins() to funcs_misc.py
# + add smoothXbins() for testing - this works!
# + add steel to the build101() exclude list
# + lsveto needs minimum energy deposit in any crystal > 0.1 keV
# ~ set lsveto mc energy resolution back to default
# ~ tweaked lsveto mc energy resolution to 0.95
# ~ updated functions to ignore lsveto multi-hit low-enrgy including
#   build101(), buildData101(), cutsBDT101(), and buildMC101()
# 
# email me: mkauer@physics.wisc.edu
######################################################################

import os,sys,re
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs100 import *
from funcs_misc import *


def build101(infile='backgrounds1010.txt', others=1, freuse=0, fchans=0, fxstals=[]):

    debug = 0

    data = {}
    bkgs = {}
    sigs = {}
    runtime = 0
    
    for line in readFile(infile):
        info = getInfo100(line, freuse, fchans, fxstals)
        if len(info) == 0:
            continue
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # These are mandatory excludes because they are always global
        exclude = bool('lsveto' in info['key']
                    or 'plastic' in info['key']
                    or 'innersteel' in info['key']
                    or 'steel' in info['key'])
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        infos = []
        exception = 0
        if 'data' in info['key']:
            if debug: print 'DEBUG: if data in info[key]: -->', info['key']
            infos.append(info)
        elif exception:
            if debug: print 'DEBUG: elif exception: -->', info['key']
            key = info['key']
            newinfo = deepcopy(info)
            newinfo['key'] = key+'-f'+str(info['xstl'])
            newinfo['from'] = info['xstl']
            infos.append(newinfo)
        elif info['xstl']==9 and not exclude:
            if debug: print 'DEBUG: elif info[xstl]==9 and not exclude: -->', info['key']
            key = info['key']
            for i in range(8):
                #info['key'] = key+'-f'+str(i+1)
                newinfo = deepcopy(info)
                newinfo['key'] = key+'-f'+str(i+1)
                newinfo['from'] = str(i+1)
                infos.append(newinfo)
        elif others and not exclude:
            if debug: print 'DEBUG: elif others and not exclude: -->', info['key']
            key = info['key']
            for i in range(8):
                #info['key'] = key+'-f'+str(i+1)
                newinfo = deepcopy(info)
                newinfo['key'] = key+'-f'+str(i+1)
                newinfo['from'] = str(i+1)
                infos.append(newinfo)
        elif not others and not exclude:
            if debug: print 'DEBUG: elif not others and not exclude: -->', info['key']
            key = info['key']
            newinfo = deepcopy(info)
            newinfo['key'] = key+'-f'+str(info['xstl'])
            newinfo['from'] = info['xstl']
            infos.append(newinfo)
        elif exclude:
            if debug: print 'DEBUG: elif exclude: -->', info['key']
            infos.append(info)
        else:
            print 'WARNING:',info['key'],'does not fit any known criteria'
            print '         Please check out build101() '
            infos.append(info)
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for info in infos:
            #print info['key']
            if 'D' in info['type']:
                data, runtime = buildData101(info, data)
            
            elif 'B' in info['type']:
                bkgs = buildMC101(info, bkgs)
            
            elif 'F' in info['type']:
                sigs = buildMC101(info, sigs)
            
            else:
                print 'WARNING: I do not know type',info['type'],'in line:'
                print info['line']
                continue
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    return data, bkgs, sigs, runtime


def buildData101(info, data):
    
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

                # exclude lsveto low energy
                #if 'x9' in key and e == 0:
                #    continue
                
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
                    edep, selection = calib100(i, E)
                    #-----------------------------------------------------------------------
                    
                    
                    # DEFINE CUTS
                    #-----------------------------------------------------------------------
                    masterCut = cutsBDT101(i, C, E, edep, selection)
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


def buildMC101(info, mc):

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

                # exclude lsveto low energy
                #if 'x9' in key and e == 0:
                #    continue
                
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
            path2 = '*-CuCase-'+info['isof']+'*.root'
        
        if info['isof'] == 'Co60' and location == 'copper':
            path1 = '/data/MC/COSINE/V3.1.1/reprocessed/Co60'
            path2 = '*-Copper-'+info['isof']+'*.root'
        
        if info['isof'] == 'Sn113' and location == 'internal':
            path1 = '/data/MC/COSINE/V3.1.1/reprocessed'
            path2 = '*-internal-'+info['isof']+'*.root'
        
        if info['isof'] == 'I129' and location == 'internal':
            path1 = '/data/COSINE/WORK/pushpa/sim/sim_data'
            path2 = '*-internal-'+info['isof']+'*.root'
        
        if location == 'plastic':
            path1 = '/data/MC/COSINE/V3.1.1/reprocessed'
            path2 = '*-Plasssupport-*'+info['isof']+'*.root'
        

        ### append the path if running locally 
        if not onCup():
            path1 = '/home/mkauer/COSINE/CUP/mc-fitting'+path1
        
        
        ### quick hack to test my new lsveto MC
        #if location == 'lsveto':
        #    path1 = '/home/mkauer/COSINE/CUP/processed'
        #    path2 = '*'+location+info['isof']+'*root'
        
        
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
                    

                    ###  put together all the cuts you need...
                    #=====================================================================================
                    
                    ### general energy cut
                    energyCut = TCut('(edep['+str(i)+']*1000. > 0.0)')
                    if pushpasMC: energyCut = TCut('(edep[6]*1000. > 0.0)')
                    if info['xstl'] == 9:
                        #energyCut = TCut('(edepResol[8]*1000. > 0.0)')
                        energyCut = TCut('((edepResol[8]*1000. > 0.0) && '
                                         +'(edep[0]*1e3 > 0.1 || '
                                         + 'edep[1]*1e3 > 0.1 || '
                                         + 'edep[2]*1e3 > 0.1 || '
                                         + 'edep[3]*1e3 > 0.1 || '
                                         + 'edep[4]*1e3 > 0.1 || '
                                         + 'edep[5]*1e3 > 0.1 || '
                                         + 'edep[6]*1e3 > 0.1 || '
                                         + 'edep[7]*1e3 > 0.1))')
                    
                    # skip low energy lsveto
                    if info['xstl'] == 9 and (e == 0 or info['chan'] == 'S'):
                        energyCut = TCut('(0)')
                    
                    ### single hit cuts
                    if info['chan'] == 'S':
                        
                        ### single-hit cut
                        hitCut = '((singleHitTag['+str(i)+'] == 1) && (multipleHitTag['+str(i)+'] == -1))'
                        if pushpasMC: hitCut = '((singleHitTag[6] == 1) && (multipleHitTag[6] == -1))'
                        
                        ### ls veto cut
                        ### need to use smeared resolution
                        lsvetocut = '(edepResol[8]*1000. < 20.0)'
                        if info['xstl'] == 9:
                            lsvetocut = '(edepResol[8]*1000. > 0.0)'
                            
                        ### combined cuts
                        chanCut = TCut('(('+hitCut+') && ('+lsvetocut+'))')

                    ### multi hit cuts
                    elif info['chan'] == 'M':
                        
                        ### multi-hit cut
                        hitCut = '((multipleHitTag['+str(i)+'] == 1) && (singleHitTag['+str(i)+'] == -1))'
                        if pushpasMC: hitCut = '((multipleHitTag[6] == 1) && (singleHitTag[6] == -1))'
                        
                        ### ls veto cut
                        ### need to use smeared resolution
                        lsvetocut = '(edepResol[8]*1000. > 20.0)'
                        if info['xstl'] == 9:
                            lsvetocut = '(edepResol[8]*1000. > 0.0)'
                            
                        ### combined cuts
                        chanCut = TCut('(('+hitCut+') || ('+lsvetocut+'))')
                        
                    else:
                        print 'ERROR: I do not know what to do with channel -->', info['chan']
                        print 'Available channels are [S]Single-hits, [M]Multi-hits'
                        sys.exit()
                    

                    ###  primary volume cuts
                    #-------------------------------------------------------------------------------
                    
                    pmt1 = str((int(info['from'])*2)-2)
                    pmt2 = str((int(info['from'])*2)-1)
                    
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
                        if pushpasMC: volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                    
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
                    
                    elif info['loca'] == 'plastic':
                        volumeCut = TCut('(primVolumeName == "PlasSupport")')
                    
                    else:
                        print "WARNING: No selection criteria for  --> ", info['loca']
                        continue
                    
                    
                    ### broken chain / group number cut
                    brokenChainCut = groupNum93(info)
                    
                    ### event type cut
                    evType = 'evt_Type'  # old processing - still works with new processing
                    #evType = 'event_info.Type'  # new processing
                    
                    eventTypeCut = TCut('('+evType+' > 10)')
                    generatedCuts = TCut('('+evType+' < 10)'+' && '+volumeCut.GetTitle())
                    #generatedCuts = TCut('('+evType+' < 10)'+' && '+volumeCut.GetTitle()+' && '+energyCut.GetTitle())
                    
                    ### special case for H3
                    if info['chst'] == 'H3':
                        generatedCuts = TCut('('+evType+' > 10)'+' && '+volumeCut.GetTitle())
                    
                    
                    ### create a hist of the generated events
                    #---------------------------------------------------------------------
                    key2 = key+'_generated'
                    temp2 = TH1F(key2, 'generated', 1, 1, 2)
                    
                    chain.Draw('1 >> '+key2, generatedCuts)
                                        
                    mc[key]['generated_hist'] = temp2
                    generated = temp2.GetEntries()
                    mc[key]['generated'] = generated
                    #---------------------------------------------------------------------
                    

                    ### separate real crystals from lsveto
                    #=====================================================================
                    if i >= 0 and i <=7:
                        # using Box-Muller for resolution smearing method
                        # setting "rng" as the smearing alias
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

                        ### Use default resolution smeared MC?
                        #selection = '(edepResol['+str(i)+']*1000.)'
                        #if pushpasMC: selection = '(edepResol[6]*1000.)'
                        
                        ### or use a different resolution function?
                        selection = '((edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng))'
                        if pushpasMC: selection = '((edep[6]*1000.) + (sigma*edep[6]*1000.*rng))'
                        
                    ### if lsveto
                    if i == 8:
                        # XXX : not sure if these cuts are right
                        #       but they seem to be working well...
                        masterCut = TCut('('+
                                    energyCut.GetTitle()+' && '+
                                    eventTypeCut.GetTitle()+' && '+
                                    brokenChainCut.GetTitle()+' && '+
                                    chanCut.GetTitle()+' && '+
                                    volumeCut.GetTitle()
                                         +')')
                        
                        #masterCut = TCut('1')
                        
                        # XXX : I think the resolution function is smearing too much?
                        selection = '(edepResol[8] * 1000.)'
                        #selection = '(edep[8] * 1000.)'
                        #selection = '(0.95 * edepResol[8] * 1000.)'
                        #selection = '(0.90 * edepResol[8] * 1000.)'
                        
                        
                    #=====================================================================
                    
                    chain.Draw(selection+' >> '+key, masterCut)
                    detected = histo.GetEntries()
                    mc[key]['hist'] = histo
                    
                    if (c=='S' and e==0) or (c=='M' and e==1):
                        print 'DEBUG:', key, 'generated events =', generated
                        print 'DEBUG:', key, 'detected events =', detected
                        print 'DEBUG:', key, 'efficiency =', round(100*detected/generated, 2), '%'
                    
        else:
            print 'ERROR: no MC files found for -->', \
                'x'+str(info['xstl']), info['floca'], info['isof']
            sys.exit()
    
    return mc


def scaleBkgs101(bkgs, runtime=0):
    
    for key in bkgs:
        
        #print 'scaling -->',key
        
        bits = key.split('-')
        x    = int(bits[0][-1])
        loca = bits[1]
        isos = bits[2]
        f    = 0
        c    = bits[-2][-1]
        e    = bits[-1][-1]
        if len(bits)==6:
            f = int(bits[3][-1])
        
        # conversion factors for things
        keVperBin  = 1./float(bkgs[key]['pars'][3])
        if runtime:
            day    = float(runtime)
        else:
            day    = 86400.
        
        if f:
            nmass  = cmass(f-1)
            surf   = cmass(f-1)
            #surf   = surfArea(f-1)
        else:
            nmass  = cmass(x-1)
            surf   = cmass(x-1)
            #surf   = surfArea(x-1)

        xkgs       = cmass(x-1)
        if x==9: xkgs = 1800.
        
        pmts       = 2.
        extpmts    = 14.
        plastic    = 1800. # set same as lsveto for comparison
        lsveto     = 1800.
        air        = 1.
        steel      = 1600.
        innersteel = 4000.
        generated  = float(bkgs[key]['generated'])
        
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        norm = 0
        if loca == 'internal':      norm = (nmass)
        elif 'surf' in loca:        norm = (surf)
        elif 'teflon' in loca:      norm = (surf)
        elif loca == 'copper':      norm = (surf)
        elif loca == 'cucase':      norm = (surf)
        elif loca == 'coppercase':  norm = (surf)
        elif loca == 'pmt':         norm = (pmts)
        elif loca == 'extpmt':      norm = (extpmts)
        elif loca == 'lsveto':      norm = (lsveto)
        elif loca == 'plastic':     norm = (plastic)
        elif loca == 'lsvetoair':   norm = (air)
        elif loca == 'airshield':   norm = (air)
        elif loca == 'steel':       norm = (steel)
        elif loca == 'innersteel':  norm = (innersteel)
        else:
            print "ERROR: no background scaling for -->", loca
            sys.exit()
        
        scale = bkgs[key]['info']['acti']*(norm)*(1./1000.)*(1./generated)*(day)*(1./xkgs)*(1./keVperBin)
        
        bkgs[key]['hist'].Scale(scale)
        bkgs[key]['scale'] = scale
    
    return bkgs


def scaleSigs101(sigkeys, sigs, runtime=0):

    for key in sigkeys:
        
        #print 'scaling -->',key
        
        try:
            test = sigs[key]['fitscale']
        except:
            print 'WARNING: no fitscale for', key
            sigs[key]['info']['fitacti'] = 0
            sigs[key]['info']['fiterro'] = 0
            continue
        
        bits = key.split('-')
        x    = int(bits[0][-1])
        loca = bits[1]
        isos = bits[2]
        f = 0
        c    = bits[-2][-1]
        e    = bits[-1][-1]
        if len(bits)==6:
            f = int(bits[3][-1])
        
        # conversion factors for things
        keVperBin  = 1./float(sigs[key]['pars'][3])
        if runtime:
            day    = float(runtime)
        else:
            day    = 86400.
        
        if f:
            nmass  = cmass(f-1)
            surf   = cmass(f-1)
            #surf   = surfArea(f-1)
        else:
            nmass  = cmass(x-1)
            surf   = cmass(x-1)
            #surf   = surfArea(x-1)
        
        xkgs       = cmass(x-1)
        if x==9: xkgs = 1800.
        
        pmts       = 2.
        extpmts    = 14.
        plastic    = 1800. # set same as lsveto for comparison
        lsveto     = 1800.
        air        = 1.
        steel      = 1600.
        innersteel = 4000.
        generated  = float(sigs[key]['generated'])

        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        norm = 0
        if loca == 'internal':      norm = (1./nmass)
        elif 'surf' in loca:        norm = (1./surf)
        elif 'teflon' in loca:      norm = (1./surf)
        elif loca == 'copper':      norm = (1./surf)
        elif loca == 'cucase':      norm = (1./surf)
        elif loca == 'coppercase':  norm = (1./surf)
        elif loca == 'pmt':         norm = (1./pmts)
        elif loca == 'extpmt':      norm = (1./extpmts)
        elif loca == 'lsveto':      norm = (1./lsveto)
        elif loca == 'plastic':     norm = (1./plastic)
        elif loca == 'lsvetoair':   norm = (1./air)
        elif loca == 'airshield':   norm = (1./air)
        elif loca == 'steel':       norm = (1./steel)
        elif loca == 'innersteel':  norm = (1./innersteel)
        else:
            print "ERROR: no signal scaling for -->", loca
            sys.exit()
        
        fitActivity = sigs[key]['fitscale']*(norm)*(1000.)*(generated)*(1./day)*(xkgs)*(keVperBin)
        
        sigs[key]['info']['fitacti'] = fitActivity
        sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
    return sigs


def cutsBDT101(i, C, E, edep, selection):
    """
    From: https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    """
    
    ### special cuts for the LS veto
    ### i = 0-7 is crystals
    ### i =  8  is lsveto
    if i == 8:
        # skip single-hit and low-energy for lsveto
        if C == 'S' or E == 0:
            #coinc  = '(BLSVeto.isCoincident == 0)'
            #lsveto = '(BLSVeto.Charge/143.8 > 0)'
            #return TCut(coinc+' && '+lsveto)
            return TCut('0')
        elif C == 'M' and E == 1:
            coinc  = '(BLSVeto.isCoincident == 1)'
            #lsveto = '(BLSVeto.Charge/143.8 > 0.0)'
            lsveto = '('+selection+' > 0.0)'
            return TCut(coinc+' && '+lsveto)
        else:
            print 'ERROR: I do not know what to do with channel -->', C
            print 'Available channels are [S]Single-hits, [M]Multi-hits'
            sys.exit()

    
    # values for the alpha cut
    alpha = [2.660,
             2.640,
             2.660,
             2.680,
             2.650,
             2.655,
             2.630,
             2.660]

    # same as Pushpa since 2017-12-19
    alphaCut = '( ! ((crystal'+str(i+1)+'.energyD > 1000.) && ((pmt'+str(i+1)+'1.rqtD1_5+pmt'+str(i+1)+'2.rqtD1_5)/2. < '+str(alpha[i])+')) )'
    # try with no alpha cut
    #alphaCut = '(1)'
    
    # BDT cut values
    cutbdt = [-0.1,
               0.0,
               0.0,
              -0.05,
              -0.2,
               0.0,
               0.0,
              -0.2]

    bdtebit = ['('+selection+'+20*bdt['+str(i+1)+'])',
               '1',
               '1',
               '('+selection+'+40*bdt['+str(i+1)+'])',
               '1',
               '('+selection+'+20*bdt['+str(i+1)+'])',
               '('+selection+'+20*bdt['+str(i+1)+'])',
               '1']
    
    cutbdte = [0,
               0,
               0,
               0,
               0,
               2,
               2,
               0]
    
    cutbdta = [-0.02,
               -0.02,
               -0.02,
               -0.05,
               -0.05,
               -0.02,
               -0.02,
               -0.05]
    
    ### global noise cuts
    coinc  = '(BLSVeto.isCoincident == 1)'
    muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
    charge = '(crystal'+str(i+1)+'.rqcn > -1)'
    nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
    
    noiseCut = '('+coinc+' && '+muons+' && '+charge+' && '+nc+')'
    
    bdt    = '(bdt['+str(i+1)+'] > '+str(cutbdt[i])+')'
    bdte   = '('+str(bdtebit[i])+' > '+str(cutbdte[i])+')'
    bdta   = '(bdtA['+str(i+1)+'] > '+str(cutbdta[i])+')'
    
    #bdtCuts  = '('+bdt+' && '+bdte+' && '+bdta+')'
    ### I don't like that bdte energy dependent cut
    bdtCuts  = '('+bdt+' && '+bdta+')'
    

    if C == 'S':
        lsveto = '(BLSVeto.Charge/143.8 < 20)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc < 4) && '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        masterCut = TCut(bdtCuts+' && '+noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        if E: masterCut = TCut(noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        return masterCut
    
    elif C == 'M':
        lsveto = '(BLSVeto.Charge/143.8 > 20)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc > 4) || '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        masterCut = TCut(bdtCuts+' && '+noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        if E: masterCut = TCut(noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        return masterCut

    else:
        print 'ERROR: I do not know what to do with channel -->', C
        print 'Available channels are [S]Single-hits, [M]Multi-hits'
        sys.exit()

