#!/usr/bin/env python
######################################################################
# funcs410.py
# 
# version: 2020-05-18
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# 
# 
# email: mkauer@physics.wisc.edu
######################################################################

import os,sys,re
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs400 import *


def build410(infile='backgrounds400.txt', others=1, freuse=0, fchans='SM', fxstals=[]):

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
                data = buildData400(info, data)
                
            elif 'B' in info['type']:
                #if info['loca'] == 'naisurf10umexpo':
                #if 'naisurf10umexpo' in info['loca']:
                if 'naisurfexpo' in info['loca'] \
                   or 'teflonsurfexpo' in info['loca']:
                    #print "!!!   SPECIAL PB210 PROFILE WORKING   !!!"
                    bkgs = buildMCSurf400(info, bkgs)
                else:
                    bkgs = buildMC410(info, bkgs)
            
            elif 'F' in info['type']:
                #if info['loca'] == 'naisurf10umexpo':
                #if 'naisurf10umexpo' in info['loca']:
                if 'naisurfexpo' in info['loca'] \
                   or 'teflonsurfexpo' in info['loca']:
                    #print "!!!   SPECIAL PB210 PROFILE WORKING   !!!"
                    sigs = buildMCSurf400(info, sigs)
                else:
                    sigs = buildMC410(info, sigs)
            else:
                print 'WARNING: I do not know type',info['type'],'in line:'
                print info['line']
                continue
            
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    return data, bkgs, sigs


def buildMC410(info, mc):

    #here = baseDir()
    
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
                        #lsvetocut = '(edepResol[8]*1000. < 20.0)'
                        lsvetocut = '(edepResol[8]*1000. < 80.0)'
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
                        #lsvetocut = '(edepResol[8]*1000. > 20.0)'
                        lsvetocut = '(edepResol[8]*1000. > 80.0)'
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
                        
                    elif info['loca'] == 'teflon':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Teflon")')
                        
                    elif info['loca'] == 'teflonbulk':
                        # Govindas MC
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
                        
                    elif info['loca'] == 'cucase':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Case")')
                        
                    elif info['loca'] == 'coppercase':
                        # Govindas MC
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Case")')
                        if pushpasMC:
                            volumeCut = TCut('((primVolumeName == "NaIDet07Fringe0")'
                                             +' || (primVolumeName == "NaIDet07Fringe1")'
                                             +' || (primVolumeName == "NaIDet07Case"))')
                        
                    elif info['loca'] == 'copper':
                        # Pushpas MC
                        volumeCut = TCut('((primVolumeName == "NaIDet0'+str(info['from'])+'Fringe0")'
                                         +' || (primVolumeName == "NaIDet0'+str(info['from'])+'Fringe1")'
                                         +' || (primVolumeName == "NaIDet0'+str(info['from'])+'Case"))')

                    elif info['loca'] == 'cushield':
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
                                    #particleCut.GetTitle()+' && '+
                                    volumeCut.GetTitle()
                                         +')')
                        
                        
                        ### Use default resolution smeared MC
                        selection = '(edepResol[8] * 1000.)'

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
                    
                    if (c=='S' and e==0) or (c=='M' and e==1):
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


def scaleBkgs410(bkgs, runtime=0):
    
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

        xkgs = cmass(x-1)
        if x==9: xkgs = 1800.
        
        pmts       = 16.
        plastic    = 1800. # set same as lsveto for comparison
        lsveto     = 1800.
        air        = 1.
        copper     = 1600. # set same as steel for now
        steel      = 1600.
        innersteel = 4000.
        generated  = float(bkgs[key]['generated'])
        
        if generated < 1:
            #print "WARNING: 0 events generated for -->", key
            continue
        
        norm = 0
        if   loca == 'internal':    norm = (nmass)
        elif loca == 'reflector':   norm = (surf)
        elif 'surf' in loca:        norm = (surf)
        elif 'teflon' in loca:      norm = (surf)
        elif loca == 'copper':      norm = (surf)
        elif loca == 'cucase':      norm = (surf)
        elif loca == 'coppercase':  norm = (surf)
        elif loca == 'pmt':         norm = (pmts)
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
            sys.exit()
        
        scale = bkgs[key]['info']['acti']*(norm)*(1./1000.)*(1./generated)*(day)*(1./xkgs)*(1./keVperBin)
        
        bkgs[key]['hist'].Scale(scale)
        bkgs[key]['scale'] = scale
    
    return bkgs


def scaleSigs410(sigkeys, sigs, runtime=0):

    for key in sigkeys:
        
        #print 'scaling -->',key
        
        try:
            test = sigs[key]['fitscale']
        except:
            #print 'WARNING: no fitscale for', key
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
        
        xkgs = cmass(x-1)
        if x==9: xkgs = 1800.
        
        pmts       = 16.
        plastic    = 1800. # set same as lsveto for comparison
        lsveto     = 1800.
        air        = 1.
        copper     = 1600. # set same as steel for now
        steel      = 1600.
        innersteel = 4000.
        generated  = float(sigs[key]['generated'])

        if generated < 1:
            #print "WARNING: 0 events generated for -->", key
            continue
        
        norm = 0
        if   loca == 'internal':    norm = (1./nmass)
        elif loca == 'reflector':   norm = (1./surf)
        elif 'surf' in loca:        norm = (1./surf)
        elif 'teflon' in loca:      norm = (1./surf)
        elif loca == 'copper':      norm = (1./surf)
        elif loca == 'cucase':      norm = (1./surf)
        elif loca == 'coppercase':  norm = (1./surf)
        elif loca == 'pmt':         norm = (1./pmts)
        elif loca == 'lsveto':      norm = (1./lsveto)
        elif loca == 'plastic':     norm = (1./plastic)
        elif loca == 'lsvetoair':   norm = (1./air)
        elif loca == 'airshield':   norm = (1./air)
        elif loca == 'cushield':    norm = (1./copper)
        elif loca == 'steel':       norm = (1./steel)
        elif loca == 'innersteel':  norm = (1./innersteel)
        elif loca == 'gamma':       norm = (1./innersteel)
        else:
            print "ERROR: no signal scaling for -->", loca
            sys.exit()
        
        fitActivity = sigs[key]['fitscale']*(norm)*(1000.)*(generated)*(1./day)*(xkgs)*(keVperBin)
        
        sigs[key]['info']['fitacti'] = fitActivity
        sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
            
    return sigs


def ScaleEnergy410(bkgs, binShift, shiftwhat):

    shifted = []
    for key in bkgs:
        if 'cS' in key and 'e0' in key:
            for this in shiftwhat:
                if this in key and key not in shifted:
                    xstal = int(key.split('-')[0][1]) - 1
                    if binShift[xstal] > 0:
                        #print 'DEBUG: energy shifting hist ', key
                        bkgs[key]['hist'] = MCBinShift400(bkgs[key]['hist'], binShift[xstal])
                    shifted.append(key)
    
    return bkgs
