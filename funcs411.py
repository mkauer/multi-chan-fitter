#!/usr/bin/env python
######################################################################
# funcs411.py
# 
# version: 2020-12-14
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# 
# 
# email: mkauer@physics.wisc.edu
######################################################################

import os,sys,re
from copy import copy, deepcopy
import numpy as np

from ROOT import *
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from funcs410 import *


def build411(infile='backgrounds400.txt', others=1, freuse=0, fchans='SM', fxstals=[]):

    debug = 0
    
    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo400(line, freuse, fchans, fxstals)
        #print info
        
        if len(info) == 0:
            continue
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Do not include "others" in processing [1] - Default [0]
        exception = 0
        # turn these off after testing
        if 'naisurf' in info['key']: exception = 1
        if 'teflon' in info['key']: exception = 1
        if 'teflonsurf' in info['key']: exception = 1
        if 'reflector' in info['key']: exception = 1
        # testing Na22 activity
        #if 'internal-Na22' in info['key']:
        #    print 'DEBUG: skipping others from Na22'
        #    exception = 1
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # These are mandatory excludes because they are always global
        exclude = bool('lsveto' in info['key']
                    or 'plastic' in info['key']
                    or 'innersteel' in info['key']
                    or 'steel' in info['key']
                    or 'gamma' in info['key']
                    or 'pmt' in info['key']
                    or 'cushield' in info['key'])
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        infos = []
        
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
            print 'WARNING:', info['key'], 'does not fit any known criteria'
            print '         Please check out build400() '
            infos.append(info)
            
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        for info in infos:
            #print info['key']
            if 'D' in info['type']:
                data = buildData411(info, data)
                
            elif 'B' in info['type']:
                #if info['loca'] == 'naisurf10umexpo':
                #if 'naisurf10umexpo' in info['loca']:
                if 'naisurfexpo' in info['loca'] \
                   or 'teflonsurfexpo' in info['loca']:
                    #print "!!!   SPECIAL PB210 PROFILE WORKING   !!!"
                    bkgs = buildMCSurf411(info, bkgs)
                else:
                    bkgs = buildMC411(info, bkgs)
            
            elif 'F' in info['type']:
                #if info['loca'] == 'naisurf10umexpo':
                #if 'naisurf10umexpo' in info['loca']:
                if 'naisurfexpo' in info['loca'] \
                   or 'teflonsurfexpo' in info['loca']:
                    #print "!!!   SPECIAL PB210 PROFILE WORKING   !!!"
                    sigs = buildMCSurf411(info, sigs)
                else:
                    sigs = buildMC411(info, sigs)
            else:
                print 'WARNING: I do not know type',info['type'],'in line:'
                print info['line']
                continue
            
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    return data, bkgs, sigs


def buildMC411(info, mc):
    
    if info['reuse']:
        for c in info['chans']:
            info['chan'] = c
            
            if not info['rootfile']:
                print 'ERROR: no rootfile specified in backgrounds file'
                sys.exit()
            
            rootfile = info['rootfile']
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
        pushpasMC = 0
        ### select correct paths for simulation files
        if   info['sim'] == 'v3.1.1': path1, path2, pushpasMC = mcPath101(info) # G4.9
        elif info['sim'] == 'v4.0.1': path1, path2, pushpasMC = mcPath300(info) # G4.10 SET1
        elif info['sim'] == 'v4.0.2': path1, path2, pushpasMC = mcPath301(info) # G4.10 SET2
        elif info['sim'] == 'v4.0.3': path1, path2 = mcPath411(info) # my reprocessing
        else:
            print 'ERROR: simulation version not valid -->',info['sim']
            sys.exit()
            
        ### append the path if running locally 
        if amLocal():
            path1 = os.path.join('/home/mkauer/COSINE/CUP/mc-fitting','.'+path1)
            
        ### quick hack to test my new lsveto MC
        #if location == 'lsveto':
        #    path1 = '/home/mkauer/COSINE/CUP/processed'
        #    path2 = '*'+location+info['isof']+'*root'
        
        print 'INFO: looking for files with -->', os.path.join(path1, path2)
        
        chain  = TChain("MC","")
        nfiles = 0
        nfiles = chain.Add(os.path.join(path1, path2))
        
        if nfiles > 0:
            
            print 'INFO:',nfiles,'files found for', info['loca'], info['isof']
            
            for C in info['chans']:
                info['chan'] = C
                
                for E in range(2):
                    
                    # DEFINE HIST AND KEY
                    #-----------------------------------------------------------------------
                    i = info['xstl'] - 1

                    key  = info['key']
                    key += '-c'+str(C)
                    key += '-e'+str(E)
                    
                    mc[key] = {}
                    mc[key]['info'] = info
                    
                    pars = histparam6000(E)
                    mc[key]['pars'] = pars
                    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
                    #-----------------------------------------------------------------------
                    

                    ###  put together all the cuts you need...
                    #=====================================================================================
                    
                    ### general energy cut
                    energyCut = TCut('(edep['+str(i)+']*1000. > 0.0)')
                    if pushpasMC: energyCut = TCut('(edep[6]*1000. > 0.0)')

                    ### lsveto energy cut has crystal dependent energy cuts
                    threshold = 0.125 # in keV
                    if info['xstl'] == 9 and info['chan'] == 'S':
                        energyCut = TCut('((edepResol[8]*1e3 > 0.0) && '
                                         +'(edep[0]*1e3 < '+str(threshold)+' && '
                                         + 'edep[1]*1e3 < '+str(threshold)+' && '
                                         + 'edep[2]*1e3 < '+str(threshold)+' && '
                                         + 'edep[3]*1e3 < '+str(threshold)+' && '
                                         + 'edep[4]*1e3 < '+str(threshold)+' && '
                                         + 'edep[5]*1e3 < '+str(threshold)+' && '
                                         + 'edep[6]*1e3 < '+str(threshold)+' && '
                                         + 'edep[7]*1e3 < '+str(threshold)+'))')

                    if info['xstl'] == 9 and info['chan'] == 'M':
                        energyCut = TCut('((edepResol[8]*1e3 > 0.0) && '
                                         +'(edep[0]*1e3 >= '+str(threshold)+' || '
                                         + 'edep[1]*1e3 >= '+str(threshold)+' || '
                                         + 'edep[2]*1e3 >= '+str(threshold)+' || '
                                         + 'edep[3]*1e3 >= '+str(threshold)+' || '
                                         + 'edep[4]*1e3 >= '+str(threshold)+' || '
                                         + 'edep[5]*1e3 >= '+str(threshold)+' || '
                                         + 'edep[6]*1e3 >= '+str(threshold)+' || '
                                         + 'edep[7]*1e3 >= '+str(threshold)+'))')
                    
                    # skip low energy lsveto
                    if info['xstl'] == 9 and E == 0:
                        energyCut = TCut('(0)')
                    
                    ### single hit cuts
                    if info['chan'] == 'S':
                        
                        ### single-hit cut
                        hitCut = '((singleHitTag['+str(i)+'] == 1) && (multipleHitTag['+str(i)+'] == -1))'
                        if pushpasMC: hitCut = '((singleHitTag[6] == 1) && (multipleHitTag[6] == -1))'
                        
                        ### ls veto cut for crystals
                        ### need to use smeared resolution
                        #lsvetocut = '(edepResol[8]*1000. < 20.0)'
                        lsvetocut = '(edepResol[8]*1000. < 80.0)'

                        ### cuts for lsveto
                        if info['xstl'] == 9:
                            # single hit is actually "all hits"
                            hitCut = '((singleHitTag[8] == 1) || (multipleHitTag[8] == 1))'
                            lsvetocut = '(edepResol[8]*1000. > 0.0)'
                            
                        ### combined cuts
                        chanCut = TCut('(('+hitCut+') && ('+lsvetocut+'))')

                    ### multi hit cuts
                    elif info['chan'] == 'M':
                        
                        ### multi-hit cut
                        hitCut = '((multipleHitTag['+str(i)+'] == 1) && (singleHitTag['+str(i)+'] == -1))'
                        if pushpasMC: hitCut = '((multipleHitTag[6] == 1) && (singleHitTag[6] == -1))'
                        
                        ### ls veto cut for crystals
                        ### need to use smeared resolution
                        #lsvetocut = '(edepResol[8]*1000. > 20.0)'
                        lsvetocut = '(edepResol[8]*1000. > 80.0)'

                        ### cuts for lsveto
                        if info['xstl'] == 9:
                            hitCut = '((singleHitTag[8] == -1) && (multipleHitTag[8] == 1))'
                            lsvetocut = '(edepResol[8]*1000. > 0.0)'
                        
                        ### combined cuts
                        chanCut = TCut('(('+hitCut+') || ('+lsvetocut+'))')
                                                
                    else:
                        print 'ERROR: I do not know what to do with channel -->', info['chan']
                        print 'Available channels are [S]Single-hits, [M]Multi-hits'
                        sys.exit()
                    

                    ###  primary volume cuts
                    #-------------------------------------------------------------------------------
                    
                    #pmt1 = str((int(info['from'])*2)-2)
                    #pmt2 = str((int(info['from'])*2)-1)
                    
                    volumeCut = TCut('(1)')
                    
                    if info['loca'] == 'internal':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal")')
                        ### testing alpha cut
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal" && primDecayType != "alpha")')
                        ### test cutting on primParticleName
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal" && primParticleName == "Pb208")')
                        
                    elif info['loca'] == 'reflector':
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Reflector")')
                        print '!!! NEW ALPHA CUT - reflector !!!'
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Reflector" && primDecayType != "alpha")')
                                                
                    elif 'internalsurf' in info['loca']:
                        # Estellas MC
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal")')
                        
                    elif 'naisurf' in info['loca']:
                        # Govindas MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal")')
                        print '!!! NEW ALPHA CUT - nai-surf !!!'
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal" && primDecayType != "alpha")')
                                                
                        if pushpasMC:
                            print '!!! SHOULD NOT SEE THIS ANYMORE - for nai-surf mc !!!'
                            volumeCut = TCut('(primVolumeName == "NaIDet07Crystal")')
                        
                    elif info['loca'] == 'bulkteflon':
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Teflon")')
                        print '!!! NEW ALPHA CUT - teflon-bulk !!!'
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Teflon" && primDecayType != "alpha")')
                        
                        if pushpasMC:
                            print '!!! SHOULD NOT SEE THIS ANYMORE - for teflon-bulk mc !!!'
                            volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                        
                    elif 'teflonsurf' in info['loca']:
                        # Estellas MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Teflon")')
                        print '!!! NEW ALPHA CUT - teflon-surf !!!'
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Teflon" && primDecayType != "alpha")')
                        
                        if pushpasMC:
                            print '!!! SHOULD NOT SEE THIS ANYMORE - for teflon-surf mc !!!'
                            volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                        
                    elif 'cusurf' in info['loca']:
                        # Estellas MC
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Case")')
                        
                    elif info['loca'] == 'bulkcuCase':
                        volumeCut = TCut('((primVolumeName == "NaIDet0'+str(info['from'])+'Case")'
                                         +' || (primVolumeName == "NaIDet0'+str(info['from'])+'Fringe0")'
                                         +' || (primVolumeName == "NaIDet0'+str(info['from'])+'Fringe1"))')
                        
                    elif info['loca'] == 'cuShield':
                        volumeCut = TCut('(primVolumeName == "CuShield")')
                        
                    elif info['loca'] == 'pmt':
                        volumeCut = TCut('(primVolumeName == "phys_pmt")')
                        
                    elif info['loca'] == 'lsveto':
                        volumeCut = TCut('(primVolumeName == "lsveto")')
                        #volumeCut = TCut('(primVolumeName == "lsveto" && primDecayType != "alpha")')
                        
                    elif info['loca'] == 'lsvetoair':
                        volumeCut = TCut('((primVolumeName == "DetPMTCover")'
                                         +' || (primVolumeName == "DetPMTEnvelope"))')

                    elif info['loca'] == 'airshield':
                        volumeCut = TCut('(primVolumeName == "LSVetoAirRoom")')
                        
                    elif info['loca'] == 'steel':
                        volumeCut = TCut('((primVolumeName == "SteelSupport")'
                                         +' || (primVolumeName == "SteelSupportTop"))')
                        
                    elif info['loca'] == 'innersteel':
                        volumeCut = TCut('(primVolumeName == "InnerSteel")')
                    
                    elif info['loca'] == 'plastic':
                        volumeCut = TCut('(primVolumeName == "PlasSupport")')
                    
                    elif info['loca'] == 'gamma':
                        #volumeCut = TCut('(primVolumeName == "lsvetoair")')
                        volumeCut = TCut('(1)')
                        #volumeCut = TCut('(primParticleName == "Pb208")')
                        
                    else:
                        print "WARNING: No selection criteria for  --> ", info['loca']
                        continue

                    
                    ### select a specific particle?
                    ### this is not generally used
                    ### only in special test cases for example Tl208 spectra...
                    #particleCut = TCut('(primParticleName == "Pb208")')
                    
                    ### broken chain / group number cut
                    #brokenChainCut = groupNum93(info)
                    brokenChainCut = groupNum(info)
                    
                    ### event type cut
                    evType = 'evt_Type'  # old processing - still works with new processing
                    #evType = 'event_info.Type'  # new processing
                    
                    eventTypeCut = TCut('('+evType+' > 10)')
                    ### for testing for primParticleName cut
                    #eventTypeCut = TCut('(1)')
                    
                    generatedCuts = TCut('('+evType+' < 10)'+' && '+volumeCut.GetTitle())
                    ### for testing for primParticleName cut
                    #generatedCuts = TCut(volumeCut.GetTitle())
                                        
                    ### special case for H3
                    #exceptions = ['H3', 'Na22', 'Te121m', 'Te127m', 'Cd109', 'Sn113']
                    exceptions = ['H3']
                    if info['chst'] in exceptions:
                        generatedCuts = TCut('('+evType+' > 10)'+' && '+volumeCut.GetTitle())
                    
                    ### special case for external Tl208 gammas from "*-gamma-CuShield-SET1-Tl208-*"
                    if info['loca'] == 'gamma':
                        eventTypeCut  = TCut('(1)')
                        generatedCuts = TCut('(1)')
                    
                        
                    ### create a hist of the generated events
                    #---------------------------------------------------------------------
                    key2 = key+'_generated'
                    temp2 = TH1F(key2, 'generated', 1, 1, 2)
                    
                    chain.Draw('1 >> '+key2, generatedCuts)
                                        
                    mc[key]['generated_hist'] = temp2
                    generated = temp2.GetEntries()
                    if generated <= 0:
                        if int(info['xstl']) == int(info['from']):
                            print 'WARNING: no events generated for -->', info['floca'], info['isof'],\
                                'xstl', info['xstl'], 'from', info['from']
                            #sys.exit() # don't exit - 2019-05-15
                        else:
                            print 'WARNING: no events generated for -->', info['floca'], info['isof'],\
                                'xstl', info['xstl'], 'from', info['from']
                    mc[key]['generated'] = generated
                    #---------------------------------------------------------------------
                    

                    ### separate real crystals from lsveto
                    #=====================================================================
                    if i >= 0 and i <= 7:
                        # using Box-Muller for resolution smearing method
                        # setting "rng" as the smearing alias
                        chain.SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                        resolFunc = resol80(i,E)
                        if pushpasMC:
                            resolFunc = resol80(i,E,1)
                        chain.SetAlias('sigma', resolFunc)
                        
                        masterCut = TCut('('+
                                    energyCut.GetTitle()+' && '+
                                    eventTypeCut.GetTitle()+' && '+
                                    brokenChainCut.GetTitle()+' && '+
                                    chanCut.GetTitle()+' && '+
                                    volumeCut.GetTitle()
                                         +')')

                        ### Try using default resolution smeared MC?
                        if E:
                            selection = '(edepResolD['+str(i)+']*1000.)'
                        else:
                            selection = '(edepResolA['+str(i)+']*1000.)'
                        
                        ### or use a different resolution function?
                        #selection = '((edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng))'
                        #if pushpasMC: selection = '((edep[6]*1000.) + (sigma*edep[6]*1000.*rng))'
                        
                    ### if lsveto
                    if i == 8:
                        # XXX : not sure if these cuts are right
                        #       but they seem to be working well...
                        masterCut = TCut('('+
                                    energyCut.GetTitle()+' && '+
                                    eventTypeCut.GetTitle()+' && '+
                                    brokenChainCut.GetTitle()+' && '+
                                    chanCut.GetTitle()+' && '+
                                    #particleCut.GetTitle()+' && '+
                                    volumeCut.GetTitle()
                                         +')')
                        
                        
                        ### Use default resolution smeared MC
                        #selection = '(edepResol[8] * 1000.)'
                        if E:
                            selection = '(edepResolD['+str(i)+']*1000.)'
                        else:
                            selection = '(edepResolA['+str(i)+']*1000.)'
                            
                        ### Or no resolution smearing...
                        #selection = '(edep[8] * 1000.)'
                        
                        ### Or define your own smearing function...
                        #chain.SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                        #chain.SetAlias('sigma', '(( 4. /sqrt(edep[8]*1000.)) + 0.0)')
                        #selection = '((edep[8]*1000.) + (sigma*edep[8]*1000.*rng))'
                        
                        
                    #=====================================================================
                    
                    chain.Draw(selection+' >> '+key, masterCut)
                    detected = histo.GetEntries()
                    mc[key]['hist'] = histo
                    
                    #if (C=='S' and E==0) or (C=='M' and E==1):
                    print 'DEBUG:', key, 'generated events =', generated
                    print 'DEBUG:', key, 'detected events =', detected
                    try: effic = round(100*detected/generated, 2)
                    except: effic = 0.0
                    print 'DEBUG:', key, 'efficiency =', effic, '%'
                        
        else:
            print 'ERROR: no MC files found for -->', \
                'x'+str(info['xstl']), info['floca'], info['isof']
            sys.exit()
    
    return mc


def buildData411(info, data):
    
    if info['reuse']:
        for c in info['chans']:
            info['chan'] = c
            
            if not info['rootfile']:
                print 'ERROR: no rootfile specified in backgrounds file'
                sys.exit()

            rootfile = info['rootfile']
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
                except:
                    print "ERROR: could not find runtime_hist -->",key+'_runtime'
                    sys.exit()

    else:
        
        runfile = info['file']
        #print 'looking for -->',runfile
        if not os.path.exists(runfile):
            print 'ERROR: runfile not found -->', runfile
            sys.exit()

        duration = 0.
        nfiles = 0
        chain = TChain("ntp","")
        if info['xstl'] == 9:
            chain2 = TChain("activeLS","")
            duration2 = 0.
            nfiles2 = 0
            
        for line in readFile(runfile):
            bits = line.split()
            runNum = bits[0]
            subrun = bits[1]
            build_filename = info['build']+'/mrgd_M'+runNum.zfill(6)+'.root.'+subrun.zfill(3)
            if not onCup():
                fpath = '/home/mkauer/COSINE/CUP/mc-fitting/data/COSINE/MRGD/phys/'+build_filename
            else:
                fpath = '/data/COSINE/MRGD/phys/'+build_filename

            #print 'INFO: looking for data file -->', fpath
            
            if not os.path.exists(fpath):
                #print 'WARNING: data file not found for -->', fpath
                continue
            else:
                #print 'INFO: adding file -->', fpath
                tmpruntime = 0
                tmpruntime = getDuration(fpath)
                if tmpruntime > 0:
                    duration += tmpruntime
                    nfiles += chain.Add(fpath)

                    # only runs >= 1800 use the updated trigger algorithm for activeLS
                    if info['xstl'] == 9 and int(runNum) >= 1800:
                        duration2 += tmpruntime
                        nfiles2 += chain2.Add(fpath)
                    
                else:
                    print 'WARNING: no run duration found for -->', fpath
        
        if nfiles > 0:
            print 'INFO:',nfiles,'files found for data',info['tag'],'build',info['build']
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
                    
                    pars = histparam6000(E)
                    data[key]['pars'] = pars
                    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
                    #-----------------------------------------------------------------------
                    
                    
                    # CALIBRATION
                    #-----------------------------------------------------------------------
                    if info['build'] == 'V00-04-04':
                        ### calib100 uses Pushpa's tweaked low energy calibrations
                        #edep, selection = calib100(i, E)
                        ### calib301 uses the standard production calibrations
                        ### and uses my tweaked lsveto calibration
                        ### testing differences between 04-04, 04-12, and 04-14
                        edep, selection = calib301(i, E)
                        
                    elif info['build'] == 'V00-04-12':
                        ### calib300 uses Govinda's tweaked high energy calibrations
                        ### and uses my tweaked lsveto calibration
                        #edep, selection = calib300(i, E)
                        ### calib301 uses the standard production calibrations
                        ### and uses my tweaked lsveto calibration
                        ### testing differences between 04-04, 04-12, and 04-14
                        edep, selection = calib301(i, E)
                        
                    elif info['build'] == 'V00-04-14':
                        ### calib301 uses the standard production calibrations
                        ### and uses my tweaked lsveto calibration
                        ### testing differences between 04-04, 04-12, and 04-14
                        edep, selection = calib301(i, E)
                        
                    elif info['build'] == 'V00-04-15':
                        ### calib301 uses the standard production calibrations
                        ### and uses the lsveto poly43 calibration
                        #edep, selection = calib301(i, E)
                        ### new calib for different sources of lsveto data
                        ### and tweaked high energy calibrations of the crystals
                        edep, selection = calib302(i, E, C)
                        
                    else:
                        print "ERROR: no info['build'] for the data calibration"
                        sys.exit()
                    #-----------------------------------------------------------------------
                    
                    
                    # DEFINE CUTS
                    #-----------------------------------------------------------------------
                    ### old V00-04-04 BDT cuts are in cutsBDT101()
                    #if info['build'] == 'V00-04-04':
                    #    masterCut = cutsBDT101(i, C, E, edep, selection)
                    if info['build'] == 'V00-04-04' or info['build'] == 'V00-04-12':
                        masterCut = cutsBDT300(i, C, E, edep, selection)

                    ### new BDT cuts for V00-04-14
                    elif info['build'] == 'V00-04-14':
                        masterCut = cutsBDT301(i, C, E, edep, selection)

                    ### same cuts for v00-04-15 so I'm told...
                    elif info['build'] == 'V00-04-15':
                        #masterCut = cutsBDT301(i, C, E, edep, selection)
                        ### new cuts for different lsveto data sources
                        masterCut = cutsBDT302(i, C, E, edep, selection)

                    else:
                        print "ERROR: no info['build'] for the data cuts"
                        sys.exit()
                    #-----------------------------------------------------------------------
                    
                    
                    # FILL HISTOS
                    #-----------------------------------------------------------------------
                    ### special case for lsveto data
                    ### use the activeLS chain2 for lsveto "single-hits"/"all-hits"
                    if i==8 and C == 'S':
                        if nfiles2 > 0:
                            chain2.Draw(selection+' >> '+key, masterCut)
                        subruns = nfiles2
                        runtime = duration2
                    else:
                        chain.Draw(selection+' >> '+key, masterCut)
                        subruns = nfiles
                        runtime = duration
                    
                    # add the hist to the data dict
                    data[key]['hist'] = histo

                    # save the total number of subruns
                    key2 = key+'_subruns'
                    temp2 = TH1F(key2,'subruns',1,0,1)
                    temp2.SetBinContent(1, subruns)
                    data[key]['subruns_hist'] = temp2
                    data[key]['subruns'] = subruns
                    
                    # save the total runtime
                    key3 = key+'_runtime'
                    temp3 = TH1F(key3,'runtime',1,0,1)
                    temp3.SetBinContent(1, runtime)
                    data[key]['runtime_hist'] = temp3
                    data[key]['runtime'] = runtime
                    #-----------------------------------------------------------------------
                    
        else:
            print 'ERROR: No data files found... quitting...'
            sys.exit()

    return data


def buildMCSurf411(info, mc):
    
    #============================================================================
    #  super special case function to deal with NaI surface depth profiling...
    #  and now also teflon surface depth profiling... 2020-04-29
    #============================================================================
        
    if info['reuse']:
        for c in info['chans']:
            info['chan'] = c
            
            if not info['rootfile']:
                print 'ERROR: no rootfile specified in backgrounds file'
                sys.exit()
            
            rootfile = info['rootfile']
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
        
        ### select correct paths for simulation files
        if   info['sim'] == 'v3.1.1': path1, path2, pushpasMC = mcPath101(info) # G4.9
        elif info['sim'] == 'v4.0.1': path1, path2, pushpasMC = mcPath300(info) # G4.10 SET1
        elif info['sim'] == 'v4.0.2': path1, path2, pushpasMC = mcPath301(info) # G4.10 SET2
        elif info['sim'] == 'v4.0.3': path1, path2 = mcPath411(info) # my reprocessing
        else:
            print 'ERROR: simulation version not valid -->',info['sim']
            sys.exit()
            
        ### append the path if running locally 
        if amLocal():
            path1 = os.path.join('/home/mkauer/COSINE/CUP/mc-fitting','.'+path1)
        
        print 'INFO: looking for files with -->', os.path.join(path1, path2)
        
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
                    xstal = int(info['xstl']) - 1
                    fromx = int(info['from']) - 1
                    
                    key  = info['key']
                    key += '-c'+info['chan']
                    key += '-e'+str(e)
                    
                    mc[key] = {}
                    mc[key]['info'] = info
                    
                    pars = histparam6000(e)
                    mc[key]['pars'] = pars
                    histo = TH1F(key, longNames(xstal), pars[0], pars[1], pars[2])

                    # 2020-03-19 - add this here
                    key2 = key+'_generated'
                    temp2 = TH1F(key2, 'generated', 1, 1, 2)
                    #-----------------------------------------------------------------------

                    
                    ###  2020-03-19
                    ###  adding in the event-by-event selction criteria
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                    ### try getting the expo depth from file name?
                    ### fix the depth as par[0]
                    depth = float(info['floca'].split('-')[3])
                    func = TF1('func', surfProfile, 0, 10, 1)
                    func.SetParameter(0, depth/100.)
                    #print 'using depth', depth
                    
                                        
                    entries = chain.GetEntries()
                    #print entries
                    
                    for jentry in range(entries):
                        
                        chain.GetEntry(jentry)

                        ### pull entry data our of leaf
                        #--------------------------------------------------
                        xx = chain.GetLeaf('primX0').GetValue()
                        yy = chain.GetLeaf('primY0').GetValue()
                        zz = chain.GetLeaf('primZ0').GetValue()
                        
                        #elow = chain.GetLeaf('edepResA').GetValue(xstal)
                        #ehi = chain.GetLeaf('edepResD').GetValue(xstal)
                        ### these branch names changed in gyunho's sim
                        elow = chain.GetLeaf('edepResolA').GetValue(xstal)
                        ehi = chain.GetLeaf('edepResolD').GetValue(xstal)
                        
                        shit = chain.GetLeaf('singleHitTag').GetValue(xstal)
                        mhit = chain.GetLeaf('multipleHitTag').GetValue(xstal)
                        evType = chain.GetLeaf('evt_Type').GetValue()
                        group = chain.GetLeaf('groupNo').GetValue()
                        decayType = chain.GetLeaf('primDecayType').GetValue()
                        
                        
                        ### alpha cuts?
                        #--------------------------------------------------
                        # 97 == alpha
                        # 98 == beta
                        if decayType == 97: continue
                        
                        
                        
                        ### volume cuts
                        #--------------------------------------------------
                        X0, Y0, Z0, rad, height, dep = mcDimensions(fromx)
                        dist = sqrt((X0-xx)**2 + (Y0-yy)**2)
                        zdist = abs(Z0-zz)

                        # teflon thickness (mm) - 2020-04-29
                        #teflon = 1.015
                        # i think this changed in gyunho's new sim
                        teflon = 1.016

                        ### kinda cryptic...
                        # /100 to conver to um
                        # /1000 to convert to mm
                        # *10 to go 10 expo depths in
                        #tdep = (((depth/100.)/1000.) * 10)
                        # use dep = 0.01 = 10um because that's all that was simulated
                        # actually set dep to 100um = 0.1
                        dep = 0.1
                        tdep = dep
                        
                        # for NaI surf side
                        if info['loca'].endswith('side') or info['loca'].endswith('expo'):
                            if dist > rad or dist < rad-dep: continue
                            if zdist > height: continue

                        # for NaI surf face
                        elif info['loca'].endswith('face'):
                            if zdist > height or zdist < height-dep: continue
                            if dist > rad: continue
                        
                        # teflon in-side or out-side - 2020-04-29
                        #elif info['loca'].endswith('in') or info['loca'].endswith('out'):
                        #    if dist < rad or dist > rad+teflon: continue
                        #    if zdist > height: continue
                        
                        # teflon in-side - 2020-05-01
                        elif info['loca'].endswith('in'):
                            if dist < rad or dist > rad+tdep: continue
                            if zdist > height: continue

                        # teflon out-side - 2020-05-01
                        elif info['loca'].endswith('out'):
                            if dist < rad+teflon-tdep or dist > rad+teflon: continue
                            if zdist > height: continue

                        else:
                            print 'ERROR: I do not know what to do with', info['loca']
                            sys.exit()

                        
                        ### get the expo weighting
                        #--------------------------------------------------
                        # for NaI surf side
                        if info['loca'].endswith('side') or info['loca'].endswith('expo'):
                            wtf  = func.Eval(1000.*(rad-dist))
                        
                        # for NaI surf face
                        elif info['loca'].endswith('face'):
                            wtf  = func.Eval(1000.*(height-zdist))
                        
                        # teflon in-side - 2020-04-29
                        elif info['loca'].endswith('in'):
                            wtf  = func.Eval(1000.*(dist-rad))
                        
                        # teflon out-side - 2020-04-29
                        elif info['loca'].endswith('out'):
                            wtf  = func.Eval(1000.*(rad+teflon-dist))
                        
                        else:
                            print 'ERROR: I do not know what to do with', info['loca']
                            sys.exit()
                        
                        #wtf = 1

                        
                        ### generated events - I might need to reposition this code
                        #--------------------------------------------------
                        # to stay consisent I just need volume and event type cuts?
                        # err, I think I need to add the exp weighting...
                        if evType <= 10:
                            #wtf=1 # testing
                            temp2.Fill(1, wtf)
                        
                        
                        ### energy cut
                        #--------------------------------------------------
                        # do I need separate cuts for low and high energy
                        if e==0:
                            if elow*1000. < 0.1: continue
                        if e==1:
                            if ehi*1000. < 0.1: continue
                        
                        
                        ### single/multi hit cut
                        #--------------------------------------------------
                        if c == 'S':
                            if shit != 1 or mhit != -1: continue
                        if c == 'M':
                            if shit != -1 or mhit != 1: continue
                        
                        
                        ### do I need an lsveto cut?
                        #--------------------------------------------------
                        # -- insert here...
                        # probably not, minimal impact if any.
                        # single/multi hit cut will take care of this.
                        
                        
                        ### do I need a groupNo cut?
                        #--------------------------------------------------
                        # -- insert here...
                        # not needed here for Pb210.
                        # groups 14-15 are included by default anyway.

                        
                        ### event type cut
                        #--------------------------------------------------
                        if evType <= 10: continue
                        
                        
                        ### fill the weighted histograms
                        #--------------------------------------------------
                        if e==0:
                            histo.Fill(1000.*elow, wtf)
                        if e==1:
                            histo.Fill(1000.*ehi, wtf)
                        
                        

                    #---------------------------------------------------------------------
                    mc[key]['generated_hist'] = temp2
                    generated = temp2.GetEntries()
                    if generated <= 0:
                        if int(info['xstl']) == int(info['from']):
                            print 'WARNING: no events generated for -->', info['floca'], info['isof'],\
                                  'xstl', info['xstl'], 'from', info['from']
                        else:
                            print 'WARNING: no events generated for -->', info['floca'], info['isof'],\
                                  'xstl', info['xstl'], 'from', info['from']
                    mc[key]['generated'] = generated

                    #---------------------------------------------------------------------
                    
                    detected = histo.GetEntries()
                    mc[key]['hist'] = histo
                    
                    if (c=='S' and e==0) or (c=='M' and e==1):
                        print 'DEBUG:', key, 'generated events =', generated
                        print 'DEBUG:', key, 'detected events =', detected
                        try: effic = round(100*detected/generated, 2)
                        except: effic = 0.0
                        print 'DEBUG:', key, 'efficiency =', effic, '%'

                    #---------------------------------------------------------------------
                       
        else:
            print 'ERROR: no MC files found for -->', \
                'x'+str(info['xstl']), info['floca'], info['isof']
            sys.exit()
    
    return mc


def getRuntimes(data):
    # initialize the runtimes list
    runtimes = [{} for x in range(9)]
    for i in range(9):
        runtimes[i] = {'S': [0, 1, 2], 'M': [0, 1, 2]}
    
    for key in sortDataKeys92(data):
        X, C, E = getXCEFromKey(key)
        runtimes[X][C][E] = data[key]['runtime']
        #print('INFO: {0} runtime = {1} seconds'.format(key, data[key]['runtime']))
        
    return runtimes


def getXCEFromKey(key):
    # format like x9-data-SET3-V0000415-cS-e0
    # or x2-pmt-Pb210_GRND-cM-e1_generated
    bits = key.split('-')
    X = int(bits[0][1])-1
    C = str(bits[-2][1])
    E = int(bits[-1][1])
    return X, C, E


def scaleBkgs411(bkgs, runtimes=None):
    
    for key in bkgs:
        
        #print 'scaling -->',key
        
        bits  = key.split('-')
        x     = int(bits[0][1])-1
        loca  = bits[1]
        isos  = bits[2]
        f     = None
        c     = str(bits[-2][1])
        e     = int(bits[-1][1])
        if len(bits) == 6:
            f = int(bits[3][1])-1
        if len(bits) == 7:
            f = int(bits[3][1])-1
        
        # conversion factors for things
        keVperBin  = 1./float(bkgs[key]['pars'][3])
        if runtimes is not None:
            day    = runtimes[x][c][e]
        else:
            day    = 86400.

        # get the mass of the crystal/lsveto of origin
        if f is not None:
            nmass  = cmass(f)
            surf   = cmass(f)
            #surf   = surfArea(f)
        else:
            nmass  = cmass(x)
            surf   = cmass(x)
            #surf   = surfArea(x)

        #print('{0} --> {1}'.format(key, nmass))

        # get the mass of the target detector
        xkgs = cmass(x)
            
        xpmts      = 2.
        pmts       = 16.
        plastic    = 1800. # set same as lsveto for comparison
        lsveto     = 1800.
        air        = 1.
        # cu shield is 3000 - should change this soon
        copper     = 1800. # set same as lsveto for now
        steel      = 1600.
        innersteel = 4000.
        
        generated  = float(bkgs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            bkgs[key]['scale'] = 1
            continue
        
        if   loca == 'internal':    norm = (nmass)
        elif loca == 'reflector':   norm = (surf)
        elif 'surf' in loca:        norm = (surf)
        elif 'teflon' in loca:      norm = (surf)
        elif loca == 'copper':      norm = (surf)
        elif loca == 'cucase':      norm = (surf)
        elif loca == 'coppercase':  norm = (surf)
        elif loca == 'pmt':         norm = (pmts)
        elif loca == 'xpmt':        norm = (xpmts)
        elif loca == 'lsveto':      norm = (lsveto)
        elif loca == 'plastic':     norm = (plastic)
        elif loca == 'lsvetoair':   norm = (air)
        elif loca == 'airshield':   norm = (air)
        elif loca == 'cushield':    norm = (copper)
        elif loca == 'steel':       norm = (steel)
        elif loca == 'innersteel':  norm = (innersteel)
        elif loca == 'gamma':       norm = (innersteel)
        else:
            print "ERROR: no background scaling for -->", loca
            norm = 1
            #sys.exit()
        
        scale = bkgs[key]['info']['acti']*(norm)*(1./1000.)*(1./generated)*(day)*(1./xkgs)*(1./keVperBin)
        
        bkgs[key]['hist'].Scale(scale)
        bkgs[key]['scale'] = scale
    
    return bkgs


def scaleSigs411(sigkeys, sigs, runtimes=None):

    for key in sigkeys:
        
        if 'fitscale' not in sigs[key]:
            #print 'WARNING: no fitscale for', key
            sigs[key]['info']['fitacti'] = 0
            sigs[key]['info']['fiterro'] = 0
            continue

        bits  = key.split('-')
        x     = int(bits[0][1])-1
        loca  = bits[1]
        isos  = bits[2]
        f     = None
        c     = str(bits[-2][1])
        e     = int(bits[-1][1])
        if len(bits) == 6:
            f = int(bits[3][1])-1
        if len(bits) == 7:
            f = int(bits[3][1])-1
        
        # conversion factors for things
        keVperBin  = 1./float(sigs[key]['pars'][3])
        if runtimes is not None:
            day    = runtimes[x][c][e]
        else:
            day    = 86400.
            
        # get the mass of the crystal/lsveto of origin
        if f is not None:
            nmass  = cmass(f)
            surf   = cmass(f)
            #surf   = surfArea(f)
        else:
            nmass  = cmass(x)
            surf   = cmass(x)
            #surf   = surfArea(x)
            
        # get the mass of the target detector
        xkgs = cmass(x)
        
        xpmts      = 2.
        pmts       = 16.
        plastic    = 1800. # set same as lsveto for comparison
        lsveto     = 1800.
        air        = 1.
        # cu shield is 3000 - should change this soon
        copper     = 1800. # set same as lsveto for now
        steel      = 1600.
        innersteel = 4000.
        
        generated  = float(sigs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            sigs[key]['info']['fitacti'] = 0
            sigs[key]['info']['fiterro'] = 0
            continue
        
        if   loca == 'internal':    norm = (nmass)
        elif loca == 'reflector':   norm = (surf)
        elif 'surf' in loca:        norm = (surf)
        elif 'teflon' in loca:      norm = (surf)
        elif loca == 'copper':      norm = (surf)
        elif loca == 'cucase':      norm = (surf)
        elif loca == 'coppercase':  norm = (surf)
        elif loca == 'pmt':         norm = (pmts)
        elif loca == 'xpmt':        norm = (xpmts)
        elif loca == 'lsveto':      norm = (lsveto)
        elif loca == 'plastic':     norm = (plastic)
        elif loca == 'lsvetoair':   norm = (air)
        elif loca == 'airshield':   norm = (air)
        elif loca == 'cushield':    norm = (copper)
        elif loca == 'steel':       norm = (steel)
        elif loca == 'innersteel':  norm = (innersteel)
        elif loca == 'gamma':       norm = (innersteel)
        else:
            print "ERROR: no signal scaling for -->", loca
            norm = 1
            #sys.exit()
        
        fitActivity = sigs[key]['fitscale']*(1./norm)*(1000.)*(generated)*(1./day)*(xkgs)*(keVperBin)
        
        sigs[key]['info']['fitacti'] = fitActivity
        sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
    
    return sigs


def scaleData411(data, dru=0):
    """
    Scale for DRU or not
    """
    for key in data:
        runtime = data[key]['runtime']
        if runtime == 0:
            #print('WARNING: runtime = {0} for {1}'.format(runtime, key))
            continue
        i = int(key.split('-')[0][1])-1
        days = 1.
        xkgs = 1.
        keVperBin = 1.
        if dru:
            days = float(runtime/86400.)
        xkgs = cmass(i)
        keVperBin = 1./float(data[key]['pars'][3])
        scale = float(1./(days*xkgs*keVperBin))
        data[key]['hist'].Scale(scale)
        data[key]['druScale'] = scale

    return data


def mcPath411(info):

    """ New processing for LS veto single/multi hit """
    """ Uses version alias v4.0.3 """
    
    ### use the full location-names 'floca'
    loc = info['floca']
    iso = info['isof']
    #Cx = 'C'+str(info['xstl'])
    Cx = 'C*'
    
    name = 'none'
            
    if loc == 'internal':
        if iso == 'I128' or iso == 'Na24':
            mcpath = '/data/COSINE/WORK/gyunho/BkgModel/SET2_reprocessed/internal'
        else:
            mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-internal-*'+iso+'[-|_]*root'
        
    elif loc == 'reflector':
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        filename = '*-surfaceNaI-'+Cx+'-10um-'+iso+'[-|_]*root'

    # 2020-03-23
    # changed to set2 for equal comparison to below expo
    elif loc == 'nai-surf-10um':
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceNaI-'+Cx+'-10um-'+iso+'[-|_]*root'

    # 2020-05-12
    # special case for expo weighting NaI surface Pb210
    elif 'nai-surf-expo' in loc:
        # just for testing
        Cx = 'C'+str(info['xstl'])
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        #filename = '*-surfaceNaI-'+Cx+'-10um-'+iso+'-*root'
        # new paths to mc from gyunho
        mcpath = '/data/COSINE/WORK/gyunho/BkgModel/SET2_reprocessed/surf_NaI'
        filename = 'reprocessed-surface-'+Cx+'-Pb210-*.root'
        
    # 2020-05-12
    # special case for expo weighting teflon surface Pb210
    elif 'teflon-surf-expo' in loc:
        # just for testing
        Cx = 'C'+str(info['xstl'])
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        #filename = '*-bulkTeflon-'+Cx+'-'+iso+'-*root'
        # new paths to mc from gyunho
        mcpath = '/data/COSINE/WORK/gyunho/BkgModel/SET2_reprocessed/surf_Tef'
        filename = 'reprocessed-tef-surf-'+Cx+'-10um-Pb210-*.root'
        
    elif loc == 'nai-surf-1um':
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceNaI-'+Cx+'-1um-'+iso+'[-|_]*root'
        
    elif loc == 'nai-surf-0.1um':
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceNaI-'+Cx+'-0.1um-'+iso+'[-|_]*root'
        
    elif loc == 'nai-surf-0.01um':
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceNaI-'+Cx+'-0.01um-'+iso+'[-|_]*root'
        
    elif loc == 'teflon-surf-10um':
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        filename = '*-surfaceTeflon-'+Cx+'-'+iso+'[-|_]*root'
    
    elif loc == 'teflon-surf-1um':
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceTeflon-'+Cx+'-1um-'+iso+'[-|_]*root'
    
    elif loc == 'teflon-surf-0.1um':
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceTeflon-'+Cx+'-0.1um-'+iso+'[-|_]*root'
    
    elif loc == 'bulk-teflon':
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        #filename = '*-bulkTeflon-'+Cx+'-'+iso+'[-|_]*root'
        mcpath = '/data/COSINE/WORK/mkauer/processed/bulk-teflon'
        filename = '*-bulk-teflon-*'+iso+'[-|_]*root'
        
    elif loc == 'pmt':
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        mcpath = '/data/COSINE/WORK/mkauer/processed/pmt'
        filename = '*-pmt-*'+iso+'[-|_]*root'
        
    elif loc == 'bulk-cuCase':
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        #filename = '*bulkcuCase*'+Cx+'*'+iso+'*root'
        mcpath = '/data/COSINE/WORK/mkauer/processed/bulk-cuCase'
        filename = '*-bulk-cuCase-*'+iso+'[-|_]*root'

    elif loc == 'plastic':
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        #filename = '*Acrylic*'+iso+'*root'
        mcpath = '/data/COSINE/WORK/mkauer/processed/acrylic'
        filename = '*-acrylic-*'+iso+'[-|_]*root'
                
    elif loc == 'lsveto':
        #mcpath = '/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
        #filename = info['isof']+'/set2/'+'*lsveto'+info['isof']+'*root'
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        mcpath = '/data/COSINE/WORK/mkauer/processed/lsveto'
        filename = '*-lsveto-*'+iso+'[-|_]*root'
        
    elif loc == 'cuShield':
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        #filename = '*cuShield*'+iso+'*root'
        mcpath = '/data/COSINE/WORK/mkauer/processed/cuShield'
        filename = '*-cuShield-*'+iso+'[-|_]*root'
        
    elif loc == 'gamma':
        #mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        mcpath = '/data/COSINE/WORK/mkauer/processed/gamma'
        filename = '*-single-gamma-*'+iso+'[-|_]*root' # this one now looks more recent so I'll try
        #filename = '*-gamma-CuShield-SET1-Tl208-*' # govinda said to use this one - also looks the best
        #filename = '*-gamma-CuShield-SET2-Tl208-*'
        #filename = '*-gamma-CuShield-Tl208-*'
        #filename = '*-gamma-Tl208-*root'

    elif loc == 'innersteel':
        mcpath = '/data/COSINE/WORK/mkauer/processed/innersteel'
        filename = '*-innersteel-*'+iso+'[-|_]*root'

        
    else:
        print('ERROR: no location identified for --> {0}'
              .format(loc))
        sys.exit()
    
    return mcpath, filename

