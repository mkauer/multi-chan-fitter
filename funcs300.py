#!/usr/bin/env python
######################################################################
# funcs300.py
# 
# version: 2020-10-15
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + recalibrate high energy region from my notebook fits
# ~ add paths to Gyunho's Na24 and I128
# ~ new paths for Gyunho's NaI and Teflon surface Pb210 in mcPath301()
# + added lsveto mc locations to mcPath300() and mcPath301()
# + added "alpha" cut on nai-surf, reflector, teflon-bulk, and teflon-surf
# ~ sticking with poly43 lsveto calibration for V00-04-15 data
# + added cases for V00-04-15 data in buildData400()
# + special case for surf-expo-side and surf-expo-face in mcPath301()
# + special case for surf-expo in mcPath301()
# + added special event type cut for "*-gamma-CuShield-SET1-Tl208-*"
# + add location 'gamma' to build300() global excludes in funcs300.py
# + add location 'gamma' to buildMC300() volumeCut in funcs300.py
# + add external 'gamma' file location(s) to mcPath301() in funcs300.py
# + added cutsBDT301() for new V00-04-14 BDT cuts
# + added calib301() to test differences between 04-04 and 04-14
#   while using the production calibration
# ~ do not sys.exit() if there are no generated events in buildMC300()
# + add "reflector" as an MC location option
# + add mcPath301() to use SET2 simulations
# ~ lsveto cut on MC now 80keV
# ~ lsveto cut on data now 80keV
# + added new cutsBDT300()
# ~ ah, particle cut needs to be Pb208 to get Tl208 decays
# ~ trying a bunch of different lsveto calibration functions 
# ~ tried selecting just Tl208 but it failed...
# - remove Na22, Te121m, Te127m, Cd109, Sn113 from exceptions
# + add more exceptions to the generated event type cut...
#   'Cd109', and 'Sn113'
# ~ tweaked other paths to be absolute and hopefully more dynamic
# ~ make rootfile path absolute in getInfo300()
# ~ tweak updateBkgsFile300() to accept a path for the new file
# + add more exceptions to the generated event type cut...
#   'H3', 'Na22', 'Te121m', and 'Te127m'
# ~ specify the crystal number for surface components in mcPath300()
# + add updateBkgsFile300()
# + add mcPath101() and mcPath300()
# + add calib300()
# + add buildData300() and buildMC300()
# + add getInfo300() and build300()
# 
# email: mkauer@physics.wisc.edu
######################################################################

import os,sys,re
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from funcs101 import *
from funcs_misc import *


def getInfo300(line, freuse=0, fchans='SM', fxstals=[]):

    here = baseDir()
    
    info = {}
    
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
    # channels to fit
    info['chans'] = fchans
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # crystal number
    info['xstl'] = int(bits[1])
    # force only a particular crystal(s) and skip the others
    if len(fxstals) > 0 and info['xstl'] not in fxstals:
        return []
    # what crystal is the background from?
    # set to xstl by default - can be updated later
    info['from'] = int(bits[1])
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # reuse rootfile name
    if 'root' in str(bits[-1]):
        info['rootfile'] = os.path.join(here, str(bits[-1]))
        #info['rootfile'] = str(bits[-1])
    else:
        info['rootfile'] = 0
    #-----------------------------------------------------------------
    

    if 'D' in info['type']:
        
        ### data set
        info['tag'] = str(bits[2])
        
        ### processing version
        info['build'] = str(bits[3])
        
        ### data file
        info['file'] = os.path.join(here, str(bits[4]))
        
        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-data'
        key += '-'+str(info['tag'])
        key += '-'+str(info['build']).replace('-','')
        info['key'] = key
                
    else:
        
        ### background location help
        info['floca'] = str(bits[2])
        # without any '-' for hist key help
        info['loca'] = str(bits[2]).replace('-','')
        
        ### top level isotope file name
        info['isof'] = str(bits[3])

        ### isotope chain break start
        info['chst'] = str(bits[4])

        ### isotope chain break stop
        info['chsp'] = str(bits[5])

        ### activity
        info['acti'] = float(bits[6])

        ### fit bounds
        info['fbnd'] = [float(bits[7]), float(bits[8])]

        ### simulation version
        info['sim'] = str(bits[9])
        
        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-'+info['loca']
        key += '-'+info['chst']+'_'+info['chsp']
        info['key'] = key

        ### assign a plotting group to the isotope-location
        info['group'] = setGroup(info)
        
    return info


def build300(infile='backgrounds300.txt', others=1, freuse=0, fchans='SM', fxstals=[]):

    debug = 0
    
    data = {}
    bkgs = {}
    sigs = {}
    runtime = 0
    
    for line in readFile(infile):
        info = getInfo300(line, freuse, fchans, fxstals)
        if len(info) == 0:
            continue
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # These are mandatory excludes because they are always global
        exclude = bool('lsveto' in info['key']
                    or 'plastic' in info['key']
                    or 'innersteel' in info['key']
                    or 'steel' in info['key']
                    or 'gamma' in info['key'])
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
            infos.append(info)
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for info in infos:
            #print info['key']
            if 'D' in info['type']:
                data, runtime = buildData300(info, data)
            
            elif 'B' in info['type']:
                bkgs = buildMC300(info, bkgs)
            
            elif 'F' in info['type']:
                sigs = buildMC300(info, sigs)
            
            else:
                print 'WARNING: I do not know type',info['type'],'in line:'
                print info['line']
                continue
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    return data, bkgs, sigs, runtime


def buildData300(info, data):
    
    #here = baseDir()
    runtime = 0.
    
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
                    runtime = data[key]['runtime']
                except:
                    print "ERROR: could not find runtime_hist -->",key+'_runtime'
                    sys.exit()

    else:
        
        runfile = info['file']
        #print 'looking for -->',runfile
        if not os.path.exists(runfile):
            print 'ERROR: runfile not found -->', runfile
            sys.exit()

        runtime = 0.
        nfiles = 0
        chain = TChain("ntp","")
        for line in readFile(runfile):
            #fpath = os.path.join(here, line)
            if not onCup():
                fpath = '/home/mkauer/COSINE/CUP/mc-fitting'+line
            else:
                fpath = line
            #print 'INFO: looking for data file -->', fpath
            if not os.path.exists(fpath):
                #print 'WARNING: data file not found for -->', fpath
                continue
            else:
                #print 'INFO: adding file -->', fpath
                tmpruntime = -1
                tmpruntime = getDuration(fpath)
                if tmpruntime > 0:
                    runtime += tmpruntime
                    nfiles += chain.Add(fpath)
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
                    
                    pars = histparam64(E)
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
                    else:
                        print "ERROR: no info['build'] for the data cuts"
                        sys.exit()
                    #-----------------------------------------------------------------------
                    
                    
                    # FILL HISTOS
                    #-----------------------------------------------------------------------
                    #----
                    ### apply standard cuts
                    chain.Draw(selection+' >> '+key, masterCut)
                    #----
                    ### look at the cut data...
                    #chain.Draw(selection+' >> '+key, "!("+masterCut.GetTitle()+")")
                    #----
                    ### apply no cuts
                    #chain.Draw(selection+' >> '+key)
                    #----
                    
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


def buildMC300(info, mc):

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
                    
                    pmt1 = str((int(info['from'])*2)-2)
                    pmt2 = str((int(info['from'])*2)-1)
                    
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
                        
                    elif info['loca'] == 'pmt':
                        volumeCut = TCut('(((primPMTid[0] == '+pmt1+')'
                                         +' || (primPMTid[0] == '+pmt2+'))'
                                         +' && primVolumeName == "phys_pmt")')
                        
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
                    brokenChainCut = groupNum93(info)
                    
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
                        # testing 4.0
                        #chain.SetAlias('sigma', '(( 4.0 /sqrt(edep[8]*1000.)) + 0.0)')
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


def calib300(i, E):
    
    if i==8: # C9 - lsveto
        edep = '(BLSVeto.Charge)'
        selection = '('+edep+' / 143.8)'
        #return edep, selection
        
        ### from my line fit and some tweaking
        #return edep, '(1.15*'+selection+')'
        
        ### try a 2nd order poly fit again - starting to look more promissing
        ###----------------------------------------------------------
        ### force a 0 point - "poly1_data"
        #A = 6.76183307e-05
        #B = 9.18732788e-01
        #C = 2.65151485e+01
        
        ### no 0 point - "poly2_data"
        #A = 9.9173553719e-05
        #B = 0.781818181818
        #C = 140.0
        
        ### no 0 point - little tweak - "poly3_data"
        #A = 5.71428571429e-05
        #B = 0.98
        #C = 43.7142857143
        
        #return edep, '(('+str(A)+'*'+selection+'**2) + ('+str(B)+'*'+selection+') + '+str(C)+')'
        
        ### try a 3rd order poly fit just to see...
        ###----------------------------------------------------------
        ### force a 0 point - "poly4_data"
        #A = 6.01051840721e-08
        #B = -0.000187327823691
        #C = 1.17575757576
        #D = -1.95618797396e-28
        
        ### tweak the data points and force a 0 point - "poly5_data"
        #A = 1.36717005791e-08
        #B = -1.42724884995e-05
        #C = 1.03372204054
        #D = 7.01845066262
        
        ### push the last 2 MC points out a little bit more - "poly6_data"
        #A = 4.61895420918e-10
        #B = 7.49137762259e-05
        #C = 0.927133764862
        #D = 18.1545308099

        ### try again based off the poly4 data points - "poly41_data"
        A = 6.14216259037e-08
        B = -0.000190223995721
        C = 1.17695224672
        D = 7.44693291248e-28
        
        ### try again and force a 1460 point - "poly42_data"
        #A = 8.81422486824e-08
        #B = -0.000277066019752
        #C = 1.21663237155
        #D = 2.41933873005e-20
        
        return edep, '(('+str(A)+'*'+selection+'**3) + ('+str(B)+'*'+selection+'**2) + ('+str(C)+'*'+selection+') + '+str(D)+')'


        
    ### use Govindha's fits (does not have C5 or C8)
    hiE = [[1813.0, 316.0],
           [1485.0, 315.0],
           [1370.0, 323.0],
           [1213.0, 294.0],
           [   0.0,   1.0],
           [1479.0, 286.0],
           [1069.0, 295.0],
           [   0.0,   1.0]]
    """
    loE = [[1,1],
           [1,1],
           [1,1],
           [1,1],
           [1,1],
           [-1591, 8491],
           [1,1],
           [1,1]]
    """
    if E:
        ### default but wrong in V00-04-12
        #edep = '(crystal'+str(i+1)+'.energyD)'
        #selection = '('+edep+')'

        ### Govinda's calibrations
        edep = '(crystal'+str(i+1)+'.rqcD1_5)'
        selection = '('+edep+'-'+str(hiE[i][0])+')/('+str(hiE[i][1])+')'
        
    else:
        ### default
        edep = '(crystal'+str(i+1)+'.energy)'
        selection = edep
        """
        if i==5: # C6 low energy
            edep = '(crystal'+str(i+1)+'.qc5)'
            selection = '('+edep+'-'+str(loE[i][0])+')/('+str(loE[i][1])+')'
        """
    return edep, selection


def calib301(i, E):
    """
    Using the standard SET3 calibrations for the crystals
    Tweaking the lsveto calibration...
    """
    
    if i==8: # C9 - lsveto
        
        edep = '(BLSVeto.Charge)'
        selection = '('+edep+' / 143.8)'
        #return edep, selection
    
        ### try again based off the poly4 data points - "poly41_data"
        #A = 6.14216259037e-08
        #B = -0.000190223995721
        #C = 1.17695224672
        #D = 7.44693291248e-28
        
        ### try again and force a 1460 point - "poly42_data"
        #A = 8.81422486824e-08
        #B = -0.000277066019752
        #C = 1.21663237155
        #D = 2.41933873005e-20

        ### try un-saturating the high energy a little bit = "poly43"
        ### but this still seems to give the best multi-crystal fits
        A = 4.0381226152e-08
        B = -0.000143935116267
        C = 1.15785808395
        D = 8.83648917811e-29

        ### try un-saturating the high energy a little more = "poly44"
        ### this was too much!
        #A = 2.48681708029e-08
        #B = -0.000109806394499
        #C = 1.14377998622
        #D = -6.09628447918e-27
        
        ### try un-saturating the high energy just a tiny = "poly45"
        #A = 4.22053979702e-08
        #B = -0.000147948294267
        #C = 1.15951351987
        #D = -5.7576250417e-22
        
        ### try un-saturating the high energy just a tiny = "poly46"
        #A = 4.80251584646e-08
        #B = -0.000160751767355
        #C = 1.16479495252
        #D = 4.5440374256e-28
        
        ### try un-saturating the high energy just a tiny = "poly47"
        #A = 5.22150630743e-08
        #B = -0.000169969557496
        #C = 1.16859729095
        #D = -1.00574241496e-21
        
        ### try tweaking the low energy part of the fit a litle = "poly50"
        # that went the wrong direction, oops...
        #A = 4.02698077117e-08
        #B = -9.26205577368e-05
        #C = 1.04107520387
        #D = -6.35274464267e-22
        
        ### try tweaking the low energy part of the fit a litle = "poly51"
        #A = 5.7582815735e-08
        #B = -0.000245923913043
        #C = 1.30856625259
        #D = 1.37579420822e-27

        
        return edep, '(('+str(A)+'*'+selection+'**3) '\
                    '+ ('+str(B)+'*'+selection+'**2) '\
                    '+ ('+str(C)+'*'+selection+') + '+str(D)+')'
    
    
    if E:
        ### default high energy production calibration
        edep = '(crystal'+str(i+1)+'.energyD)'
        selection = '('+edep+')'

    else:
        ### default low energy production calibration
        edep = '(crystal'+str(i+1)+'.energy)'
        selection = '('+edep+')'

    return edep, selection


def calib302(i, Engy, Chan):
    """
    Tweaked high energy crystal calibrations
    Tweaked the lsveto calibration
    Now adjust for different lsveto data sources
    """

    Cx = i+1
    
    if Cx == 9: # C9 / lsveto
        # return same calib for all Engy
        # return same calib for all Chan
        
        edep = '(BLSVeto.Charge)'
        selection = '('+edep+' / 143.8)'
        #return edep, selection
    
        ### try again based off the poly4 data points - "poly41_data"
        #A = 6.14216259037e-08
        #B = -0.000190223995721
        #C = 1.17695224672
        #D = 7.44693291248e-28
        
        ### try again and force a 1460 point - "poly42_data"
        #A = 8.81422486824e-08
        #B = -0.000277066019752
        #C = 1.21663237155
        #D = 2.41933873005e-20

        ### try un-saturating the high energy a little bit = "poly43"
        ### but this still seems to give the best multi-crystal fits
        A = 4.0381226152e-08
        B = -0.000143935116267
        C = 1.15785808395
        D = 8.83648917811e-29

        ### try un-saturating the high energy a little more = "poly44"
        ### this was too much!
        #A = 2.48681708029e-08
        #B = -0.000109806394499
        #C = 1.14377998622
        #D = -6.09628447918e-27
        
        ### try un-saturating the high energy just a tiny = "poly45"
        #A = 4.22053979702e-08
        #B = -0.000147948294267
        #C = 1.15951351987
        #D = -5.7576250417e-22
        
        ### try un-saturating the high energy just a tiny = "poly46"
        #A = 4.80251584646e-08
        #B = -0.000160751767355
        #C = 1.16479495252
        #D = 4.5440374256e-28
        
        ### try un-saturating the high energy just a tiny = "poly47"
        #A = 5.22150630743e-08
        #B = -0.000169969557496
        #C = 1.16859729095
        #D = -1.00574241496e-21
        
        ### try tweaking the low energy part of the fit a litle = "poly50"
        # that went the wrong direction, oops...
        #A = 4.02698077117e-08
        #B = -9.26205577368e-05
        #C = 1.04107520387
        #D = -6.35274464267e-22
        
        ### try tweaking the low energy part of the fit a litle = "poly51"
        #A = 5.7582815735e-08
        #B = -0.000245923913043
        #C = 1.30856625259
        #D = 1.37579420822e-27

        # return same calib for all E
        selection = '(('+str(A)+'*'+selection+'**3) '\
                    '+ ('+str(B)+'*'+selection+'**2) '\
                    '+ ('+str(C)+'*'+selection+') + '+str(D)+')'
        
        return edep, selection

    if Engy == 0:
        # default low energy production calibration
        edep = '(crystal{0}.energy)'.format(Cx)
        selection = '('+edep+')'

        return edep, selection
        
    elif Engy == 1:
        # default high energy production calibration
        edep = '(crystal{0}.energyD)'.format(Cx)
        selection = '('+edep+')'

        # tweaks per crystal
        if Cx == 1: # C1
            A = -9.20702415903431e-09
            B = 3.1582891760183624e-05
            C = 0.9731283907370495
            D = 3.492470586962138
            
        if Cx == 2: # C2
            return edep, selection
        
        if Cx == 3: # C3
            A = -4.668682649387045e-09
            B = 1.5052344932980456e-05
            C = 0.9874023184977824
            D = 1.6686064222603654
            
        if Cx == 4: # C4
            A = -7.966374771112426e-09
            B = 2.544049246165629e-05
            C = 0.9801665355547313
            D = 2.2024520321640186
            
        if Cx == 5: # C5
            A = 0.813447107922342
            B = 13.863912116535856
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            return edep, selection
        
        if Cx == 6: # C6
            A = -4.623190276816897e-09
            B = 1.6693196676849473e-05
            C = 0.9855987614226132
            D = 1.789741155036786
        
        if Cx == 7: # C7
            #A = -6.584663125047853e-09
            #B = 2.4904010290725124e-05
            #C = 0.9765006123749262
            #D = 3.45348577054466
            # first try to fix
            #A = -3.5286801110126694e-09
            #B = 6.692048421004304e-06
            #C = 1.0004270013767433
            #D = -1.3444385176574656
            # second try to fix
            #A = -9.787584607766034e-09
            #B = 3.3865923406543445e-05
            #C = 0.971206025113478
            #D = 3.6683048699245857
            # third try
            A = -9.93101832337613e-09
            B = 3.6976416493554536e-05
            C = 0.9661209376849521
            D = 4.747875210154714
            
        if Cx == 8: # C8
            A = 0.8011463218736552
            B = -21.712318299469676
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            return edep, selection

        
        selection = '(('+str(A)+'*'+edep+'**3) '\
                   '+ ('+str(B)+'*'+edep+'**2) '\
                   '+ ('+str(C)+'*'+edep+') + '+str(D)+')'

        return edep, selection    

    elif Engy == 2:
        ### use default high energy calib for alphas
        edep = '(crystal{0}.energyD)'.format(Cx)
        selection = '('+edep+')'

        # except for C5 and C8
        if Cx == 5:
            A = 0.813447107922342
            B = 13.863912116535856
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            return edep, selection
        if Cx == 8:
            A = 0.8011463218736552
            B = -21.712318299469676
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            return edep, selection
        
        return edep, selection

    else:
        print('WARNING: no calibration for x{0} c{1} e{2}'.format(Cx, Chan, Engy))
        sys.exit()


def mcPath101(info):
    """
    For the old v3.1.1 G4.9 simulations
    """
    
    ### either use full location-names 'floca'
    ### OR use the hyphen removed locationnames 'loca'
    location = info['floca']

    ### I don't use extpmt anymore...
    #if location == 'extpmt':
    #    location = 'pmt'
    
    ### special paths for H3 and Cd109
    path1 = '/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
    if info['isof'] == 'H3':
        path1 = '/data/COSINE/WORK/pushpa/sim/process/Crystal7/'
    if info['isof'] == 'Cd109':
        path1 = '/data/MC/COSINE/V3.1.1/reprocessed/G4_10/'
    
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
        #path2 = '*-Copper-'+info['isof']+'*.root'
        path2 = '*-Copper-'+info['isof']+'-??.root'

    if info['isof'] == 'Sn113' and location == 'internal':
        path1 = '/data/MC/COSINE/V3.1.1/reprocessed'
        path2 = '*-internal-'+info['isof']+'*.root'

    if info['isof'] == 'I129' and location == 'internal':
        path1 = '/data/COSINE/WORK/pushpa/sim/sim_data'
        path2 = '*-internal-'+info['isof']+'*.root'

    if location == 'plastic':
        path1 = '/data/MC/COSINE/V3.1.1/reprocessed'
        path2 = '*-Plasssupport-*'+info['isof']+'*.root'
    
    pushpasMC = 0
    #if location in pushpas: pushpasMC = 1
    if info['floca'] in pushpas: pushpasMC = 1

    ### special case for the U235
    if info['isof'] == 'U235':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        if location == 'internal': path2 = '*-internal-'+iso+'-*root'
        if location == 'pmt':      path2 = '*-pmt-'+iso+'-*root'
        if location == 'plastic':  path2 = '*-Acrylic-'+iso+'-*root'
    
    return path1, path2, pushpasMC


def mcPath300(info):
    """
    For v4.0.1 G4.10 SET1 simulations
    """
    
    ### use the full location-names 'floca'
    loc = info['floca']
    iso = info['isof']
    #Cx = 'C'+str(info['xstl'])
    Cx = 'C*'
    
    name = 'none'
    if loc == 'reflector':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-surfaceNaI-'+Cx+'-10um-'+iso+'-*root'
        
    elif loc == 'nai-surf-10um':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-surfaceNaI-'+Cx+'-10um-'+iso+'-*root'
        
    elif loc == 'teflon-bulk':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-bulkTeflon-'+Cx+'-'+iso+'-*root'
        
    elif loc == 'teflon-surf':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-surfaceTeflon-'+Cx+'-'+iso+'-*root'
        
    elif loc == 'coppercase':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-bulkcuCase-'+Cx+'-'+iso+'-*root'
        
    elif loc == 'internal':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-internal-'+iso+'-*root'
        
    elif loc == 'pmt':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-pmt-'+iso+'-*root'
        
    elif loc == 'plastic':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-Acrylic-'+iso+'-*root'
        
    elif loc == 'lsveto':
        path1 = '/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
        path2 = info['isof']+'/set2/'+'*lsveto'+info['isof']+'*root'
        
    else:
        print 'ERROR: no location identified for -->',loc
        sys.exit()
        
    pushpasMC = 0
    return path1, path2, pushpasMC


def mcPath301(info):
    """
    For v4.0.1 G4.10 SET2 simulations but using alias v4.0.2
    """
    
    ### use the full location-names 'floca'
    loc = info['floca']
    iso = info['isof']
    #Cx = 'C'+str(info['xstl'])
    Cx = 'C*'
    
    name = 'none'
    if loc == 'reflector':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-surfaceNaI-'+Cx+'-10um-'+iso+'-*root'

    # 2020-03-23
    # changed to set2 for equal comparison to below expo
    elif loc == 'nai-surf-10um':
        #path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*-surfaceNaI-'+Cx+'-10um-'+iso+'-*root'

    # 2020-05-12
    # special case for expo weighting NaI surface Pb210
    elif 'nai-surf-expo' in loc:
        # just for testing
        Cx = 'C'+str(info['xstl'])
        #path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        #path2 = '*-surfaceNaI-'+Cx+'-10um-'+iso+'-*root'
        # new paths to mc from gyunho
        path1 = '/data/COSINE/WORK/gyunho/BkgModel/SET2_reprocessed/surf_NaI'
        path2 = 'reprocessed-surface-'+Cx+'-Pb210-*.root'
        
    # 2020-05-12
    # special case for expo weighting teflon surface Pb210
    elif 'teflon-surf-expo' in loc:
        # just for testing
        Cx = 'C'+str(info['xstl'])
        #path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        #path2 = '*-bulkTeflon-'+Cx+'-'+iso+'-*root'
        # new paths to mc from gyunho
        path1 = '/data/COSINE/WORK/gyunho/BkgModel/SET2_reprocessed/surf_Tef'
        path2 = 'reprocessed-tef-surf-'+Cx+'-10um-Pb210-*.root'
        
    elif loc == 'nai-surf-1um':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*-surfaceNaI-'+Cx+'-1um-'+iso+'-*root'
        
    elif loc == 'nai-surf-0.1um':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*-surfaceNaI-'+Cx+'-0.1um-'+iso+'-*root'
        
    elif loc == 'nai-surf-0.01um':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*-surfaceNaI-'+Cx+'-0.01um-'+iso+'-*root'
        
    elif loc == 'teflon-bulk':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*-bulkTeflon-'+Cx+'-'+iso+'-*root'
    
    elif loc == 'teflon-surf-10um':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-surfaceTeflon-'+Cx+'-'+iso+'-*root'
    
    elif loc == 'teflon-surf-1um':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*-surfaceTeflon-'+Cx+'-1um-'+iso+'-*root'
    
    elif loc == 'teflon-surf-0.1um':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*-surfaceTeflon-'+Cx+'-0.1um-'+iso+'-*root'
    
    elif loc == 'coppercase':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*bulkcuCase*'+Cx+'*'+iso+'*root'
        
    elif loc == 'internal':
        if iso == 'I128' or iso == 'Na24':
            path1 = '/data/COSINE/WORK/gyunho/BkgModel/SET2_reprocessed/internal'
        else:
            path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*internal*'+iso+'*root'
        
    elif loc == 'pmt':
        #path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        #path2 = '*-pmt-'+iso+'-*root'
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*pmt*'+iso+'*root'
        
    elif loc == 'plastic':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*Acrylic*'+iso+'*root'
        
    elif loc == 'lsveto':
        #path1 = '/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
        #path2 = info['isof']+'/set2/'+'*lsveto'+info['isof']+'*root'
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*lsveto*'+iso+'*root'
        
    elif loc == 'gamma':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*single-gamma*'+iso+'*root' # this one now looks more recent so I'll try
        #path2 = '*-gamma-CuShield-SET1-Tl208-*' # govinda said to use this one - also looks the best
        #path2 = '*-gamma-CuShield-SET2-Tl208-*'
        #path2 = '*-gamma-CuShield-Tl208-*'
        #path2 = '*-gamma-Tl208-*root'

    elif loc == 'cushield':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        path2 = '*cuShield*'+iso+'*root'
        
    else:
        print 'ERROR: no location identified for -->',loc
        sys.exit()
        
    pushpasMC = 0
    return path1, path2, pushpasMC


def updateBkgsFile300(xstals, bkgsfile, resultsfile, newdir, BF='BR'):
    """
    Only print out to file what was actually used in the fit
    """
    
    if os.path.exists(bkgsfile):
        with open(bkgsfile) as fbkgs:
            bkgslines = fbkgs.read().splitlines()
    else:
        print 'WARNING: file not found -->', bkgsfile
        return
    
    if os.path.exists(resultsfile):
        with open(resultsfile) as ffits:
            fitlines = ffits.read().splitlines()
    else:
        fitlines = ''
        
    newbkgs = os.path.join(newdir, bkgsfile.split('/')[-1][:-4]+'_update.txt')
    output = open(newbkgs, 'w')
    output.write('# -*- sh -*-\n')
    
    print ''
    print 'INFO: Updating bkgsfile -->', bkgsfile
    print '      To a new bkgsfile -->', newbkgs
    print ''

    skip = 0
    for bline in bkgslines:

        if not bline:
            continue
        if bline.startswith('\"\"\"'):
            if skip == 0: skip = 1
            else: skip = 0
            continue
        if skip:
            continue
        if bline.startswith('#'):
            continue
        
        bbits = filter(None, re.split("[ \s\t\n\r,:]+", bline.strip()))

        # add [0] to xstals for v500 global mc
        # this should still work with older versions
        #if int(bbits[1]) not in xstals:
        if int(bbits[1]) not in xstals+[0]:
            continue
        
        if 'F' not in bbits[0]:
            for j in range(len(bbits)-1):
                output.write(bbits[j]+'\t')
            output.write(bbits[-1]+'\n')
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
                    if bbits[1] == xstal and bbits[2].replace('-','') == loca \
                       and bbits[4] == chst and bbits[5] == chsp:
                        for i in range(lenbbits):
                            if i == 0:
                                output.write(BF+'\t')
                            elif i == 6:
                                if acti != '0.0': output.write(acti+'\t')
                                else: output.write(bbits[i]+'\t')
                            elif i == 7:
                            #elif i==(lenbbits-3):
                                output.write('0.1\t')
                            elif i == 8:
                            #elif i==(lenbbits-2):
                                output.write('10\t')
                            elif i==(lenbbits-1):
                                output.write(bbits[i])
                            else:
                                output.write(bbits[i]+'\t')
                        output.write('\n')
                        replaced = 1
        
        #if not replaced:
        #    output.write(bline+'\n')
        #    print 'WARNING: Could not match -->', filter(None, re.split("[ \s\t\n\r,:]+", bline.strip()))[3:7]
            
    output.close()
    return


def cutsBDT300(i, C, E, edep, selection):
    """
    https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    https://cupwiki.ibs.re.kr/Kims/EventSelection
    """
    
    ### special cuts for the LS veto
    ### i = 0-7 is crystals
    ### i =  8  is lsveto
    if i == 8:
        ### skip single-hit and low-energy for lsveto
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

    ### values for the alpha cut
    alpha = [2.660,
             2.640,
             2.660,
             2.680,
             2.650,
             2.655,
             2.630,
             2.660]

    ### same as Pushpa since 2017-12-19
    alphaCut = '( ! ((crystal'+str(i+1)+'.energyD > 1000.) && ((pmt'+str(i+1)+'1.rqtD1_5+pmt'+str(i+1)+'2.rqtD1_5)/2. < '+str(alpha[i])+')) )'

    ### new BDT cuts from https://cupwiki.ibs.re.kr/Kims/EventSelection
    newBDT = ['((bdt[1]>-0.1 && (crystal1.energy + 20*bdt[1])>0) && ((crystal1.energy>=10 && bdtA[1]>-0.08) || (crystal1.energy>=8 && crystal1.energy<10 && bdtA[1]>-0.065) || (crystal1.energy>=7 && crystal1.energy<8 && bdtA[1]>0.005) || (crystal1.energy>=6 && crystal1.energy<7 && bdtA[1]>0.015) || (crystal1.energy>=5 && crystal1.energy<6 && bdtA[1]>0.01) || (crystal1.energy>=4 && crystal1.energy<5 && bdtA[1]>0.03) || (crystal1.energy>=3 && crystal1.energy<4 && bdtA[1]>0.035) || (crystal1.energy>=2 && crystal1.energy<3 && bdtA[1]>0.02) || (crystal1.energy<2 && bdtA[1]>0.02)))',
        '(bdt[2]>0.0 && bdtA[2]>-0.05)',
        '(bdt[3]>0.0 && bdtA[3]>-0.07)',
        '(bdt[4]>-0.05 && (crystal4.energy + 40*bdt[4])>0 && bdtA[4]>-0.07)',
        '(bdt[5]>-0.2 && bdtA[5]>-0.07)',
        '(bdt[6]>0.0 && (crystal6.energy + 20*bdt[6])>2.0 && bdtA[6]>-0.07)',
        '(bdt[7]>0.0 && (crystal7.energy + 20*bdt[7])>2.0 && bdtA[7]>-0.07)',
        '(bdt[8]>-0.2 && bdtA[8]>-0.07)']

    bdtCut = newBDT[i]

    ### global noise cuts
    coinc  = '(BLSVeto.isCoincident == 1)'
    muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
    charge = '(crystal'+str(i+1)+'.rqcn > -1)'
    nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
    t1     = '(pmt'+str(i+1)+'1.t1 > 2 && pmt'+str(i+1)+'2.t1 > 2)'
        
    noiseCut = '('+coinc+' && '+muons+' && '+charge+' && '+nc+' && '+t1+')'
    
    if C == 'S':
        #lsveto = '(BLSVeto.Charge/143.8 < 20)'
        lsveto = '(BLSVeto.Charge/143.8 < 80)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc < 4) && '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        masterCut = TCut(bdtCut+' && '+noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        if E: masterCut = TCut(noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        return masterCut
    
    elif C == 'M':
        #lsveto = '(BLSVeto.Charge/143.8 > 20)'
        lsveto = '(BLSVeto.Charge/143.8 > 80)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc > 4) || '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        masterCut = TCut(bdtCut+' && '+noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        if E: masterCut = TCut(noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        return masterCut

    else:
        print 'ERROR: I do not know what to do with channel -->', C
        print 'Available channels are [S]Single-hits, [M]Multi-hits'
        sys.exit()


def cutsBDT301(i, C, E, edep, selection):
    """
    https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    https://cupwiki.ibs.re.kr/Kims/EventSelection
    """
    #-------------------------------------------------------
    # This is for the V00-04-14 data with new BDT values!!!
    #-------------------------------------------------------

    ### values for the alpha cut
    alpha = [2.660,
             2.640,
             2.660,
             2.680,
             2.650,
             2.655,
             2.630,
             2.660]

    ### special cuts for the LS veto
    ### i = 0-7 is crystals
    ### i =  8  is lsveto
    if i == 8:
        ### skip single-hit and low-energy for lsveto
        if C == 'S' or E == 0:
            return TCut('0')
        
        elif C == 'M' and E == 1:
            coinc  = '(BLSVeto.isCoincident == 1)'
            lsveto = '('+selection+' > 0.0)'
            return TCut(coinc+' && '+lsveto)
            
            ### 2020-08-11
            ### test adding standard crystal cuts to lsveto data
            ### short story - this doesn't really change anything
            """
            muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
            alphaCut = '('
            nc = '('
            t1 = '('
            hits = '('
            charge = '('
            for j in range(8):
                alphaCut += '( ! ((crystal'+str(j+1)+'.energyD > 1000.) && ((pmt'+str(j+1)+'1.rqtD1_5 + pmt'+str(j+1)+'2.rqtD1_5)/2. < '+str(alpha[j])+')) ) && '
                nc += '(pmt'+str(j+1)+'1.nc > 0 && pmt'+str(j+1)+'2.nc > 0) || '
                t1 += '(pmt'+str(j+1)+'1.t1 > 0 && pmt'+str(j+1)+'2.t1 > 0) || '
                hits += '(crystal'+str(j+1)+'.nc > 4) || '
                charge += '(crystal'+str(j+1)+'.rqcn > -1) || '
            # remove trailing && or ||
            alphaCut = alphaCut[:-4]+')'
            nc = nc[:-4]+')'
            t1 = t1[:-4]+')'
            hits = hits[:-4]+')'
            charge = charge[:-4]+')'
            return TCut(coinc+' && '+lsveto+' && '+muons+' && '+alphaCut+' && '+nc+' && '+t1+' && '+hits+' && '+charge)
            """
        else:
            print 'ERROR: I do not know what to do with channel -->', C
            print 'Available channels are [S]Single-hits, [M]Multi-hits'
            sys.exit()


    ### same as Pushpa since 2017-12-19
    alphaCut = '( ! ((crystal'+str(i+1)+'.energyD > 1000.) && ((pmt'+str(i+1)+'1.rqtD1_5+pmt'+str(i+1)+'2.rqtD1_5)/2. < '+str(alpha[i])+')) )'

    ### new BDT cuts for V00-04-12 from https://cupwiki.ibs.re.kr/Kims/EventSelection
    newBDT = [
        '((bdt[1]>-0.1 && (crystal1.energy + 20*bdt[1])>0) && ((crystal1.energy>=10 && bdtA[1]>-0.08) || (crystal1.energy>=8 && crystal1.energy<10 && bdtA[1]>-0.065) || (crystal1.energy>=7 && crystal1.energy<8 && bdtA[1]>0.005) || (crystal1.energy>=6 && crystal1.energy<7 && bdtA[1]>0.015) || (crystal1.energy>=5 && crystal1.energy<6 && bdtA[1]>0.01) || (crystal1.energy>=4 && crystal1.energy<5 && bdtA[1]>0.03) || (crystal1.energy>=3 && crystal1.energy<4 && bdtA[1]>0.035) || (crystal1.energy>=2 && crystal1.energy<3 && bdtA[1]>0.02) || (crystal1.energy<2 && bdtA[1]>0.02)))',
              
        #'(bdt[2]>0.0 && bdtA[2]>-0.05)',
        # new bdt cut for V00-04-14
        '(((9.20842e-07*TMath::Exp(-104.504*bdt[2])+0.170872)-(9.70874*bdt[2])) < crystal2.energy)',
              
        #'(bdt[3]>0.0 && bdtA[3]>-0.07)',
        # new bdt cut for V00-04-14
        '(((6.81804e-08*TMath::Exp(-101.856*bdt[3])+0.148344)-(4.04826*bdt[3])) < crystal3.energy)',

        #'(bdt[4]>-0.05 && (crystal4.energy + 40*bdt[4])>0 && bdtA[4]>-0.07)',
        # new bdt cut for V00-04-14
        '(((1.37726e-07*TMath::Exp(-111.842*bdt[4])+0.446818)-(2.13498*bdt[4])) < crystal4.energy)',
              
        '(bdt[5]>-0.2 && bdtA[5]>-0.07)',
              
        #'(bdt[6]>0.0 && (crystal6.energy + 20*bdt[6])>2.0 && bdtA[6]>-0.07)',
        # new bdt cut for V00-04-14
        '(((0.000315623*TMath::Exp(-75.3998*bdt[6])-0.313482)-(13.6302*bdt[6])) < crystal6.energy)',
              
        #'(bdt[7]>0.0 && (crystal7.energy + 20*bdt[7])>2.0 && bdtA[7]>-0.07)',
        # new bdt cut for V00-04-14
        '(((1.45552e-05*TMath::Exp(-88.7196*bdt[7])+0.566336)-(7.57773*bdt[7])) < crystal7.energy)',
        
        '(bdt[8]>-0.2 && bdtA[8]>-0.07)']

    bdtCut = newBDT[i]

    ### global noise cuts
    #===========================================================================
    coinc  = '(BLSVeto.isCoincident == 1)'
    muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
    charge = '(crystal'+str(i+1)+'.rqcn > -1)'
    
    ### old cut
    #nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
    ### new cut for v00-04-14 from Govinda
    nc     = '(pmt'+str(i+1)+'1.nc > 0 && pmt'+str(i+1)+'2.nc > 0)'
    
    ### old cut
    #t1     = '(pmt'+str(i+1)+'1.t1 > 2 && pmt'+str(i+1)+'2.t1 > 2)'
    ### new cut for v00-04-14 from Govinda
    t1     = '(pmt'+str(i+1)+'1.t1 > 0 && pmt'+str(i+1)+'2.t1 > 0)'
    
    noiseCut = '('+coinc+' && '+muons+' && '+charge+' && '+nc+' && '+t1+')'
    #===========================================================================
    
    if C == 'S':
        #lsveto = '(BLSVeto.Charge/143.8 < 20)'
        lsveto = '(BLSVeto.Charge/143.8 < 80)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc < 4) && '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        masterCut = TCut(bdtCut+' && '+noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        if E: masterCut = TCut(noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        return masterCut
    
    elif C == 'M':
        #lsveto = '(BLSVeto.Charge/143.8 > 20)'
        lsveto = '(BLSVeto.Charge/143.8 > 80)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc > 4) || '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        masterCut = TCut(bdtCut+' && '+noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        if E: masterCut = TCut(noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        return masterCut

    else:
        print 'ERROR: I do not know what to do with channel -->', C
        print 'Available channels are [S]Single-hits, [M]Multi-hits'
        sys.exit()


def cutsBDT302(i, C, E, edep, selection):
    """
    https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    https://cupwiki.ibs.re.kr/Kims/EventSelection
    """
    #-------------------------------------------------------
    # This is for the V00-04-15 data
    # New lsveto single-hit data
    #-------------------------------------------------------

    ### values for the alpha cut
    alpha = [
        2.660,
        2.640,
        2.660,
        2.680,
        2.650,
        2.655,
        2.630,
        2.660
    ]

    ### special cuts for the LS veto
    ### i = 0-7 is crystals
    ### i =  8  is lsveto
    if i == 8:
        ### skip low energy and alpha energy for lsveto
        #if E in [0, 2]:
        #    return TCut('0')
        
        if C == 'S':
            coinc  = '(1)'
            energy = '('+selection+' > 0.0)'
            #muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
            muons  = '(1)'
            # time cut doesn't make much difference...
            time   = '(BLSVeto.Time > 2300 && BLSVeto.Time < 2600)'
            #time   = '(1)'
            return TCut('({0} && {1} && {2} && {3})'.format(coinc, energy, muons, time))
            
        elif C == 'M':
            coinc  = '(BLSVeto.isCoincident == 1)'
            energy = '('+selection+' > 0.0)'
            muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
            #time   = '(BLSVeto.Time >= 2800 && BLSVeto.Time <= 2950)'
            time   = '(1)'
            return TCut('({0} && {1} && {2} && {3})'.format(coinc, energy, muons, time))
            
        else:
            print 'ERROR: I do not know what to do with channel -->', C
            print 'Available channels are [S]Single-hits, [M]Multi-hits'
            sys.exit()


    ### same as Pushpa since 2017-12-19
    # only alphas for E=2
    if E == 2:
        alphaCut = '((crystal'+str(i+1)+'.energyD > 1000.) && ((pmt'+str(i+1)+'1.rqtD1_5+pmt'+str(i+1)+'2.rqtD1_5)/2. < '+str(alpha[i])+'))'
    else:
        alphaCut = '( ! ((crystal'+str(i+1)+'.energyD > 1000.) && ((pmt'+str(i+1)+'1.rqtD1_5+pmt'+str(i+1)+'2.rqtD1_5)/2. < '+str(alpha[i])+')) )'

        
    ### new BDT cuts for V00-04-12 from https://cupwiki.ibs.re.kr/Kims/EventSelection
    newBDT = [
        '((bdt[1]>-0.1 && (crystal1.energy + 20*bdt[1])>0) && ((crystal1.energy>=10 && bdtA[1]>-0.08) || (crystal1.energy>=8 && crystal1.energy<10 && bdtA[1]>-0.065) || (crystal1.energy>=7 && crystal1.energy<8 && bdtA[1]>0.005) || (crystal1.energy>=6 && crystal1.energy<7 && bdtA[1]>0.015) || (crystal1.energy>=5 && crystal1.energy<6 && bdtA[1]>0.01) || (crystal1.energy>=4 && crystal1.energy<5 && bdtA[1]>0.03) || (crystal1.energy>=3 && crystal1.energy<4 && bdtA[1]>0.035) || (crystal1.energy>=2 && crystal1.energy<3 && bdtA[1]>0.02) || (crystal1.energy<2 && bdtA[1]>0.02)))',
              
        #'(bdt[2]>0.0 && bdtA[2]>-0.05)',
        # new bdt cut for V00-04-14
        '(((9.20842e-07*TMath::Exp(-104.504*bdt[2])+0.170872)-(9.70874*bdt[2])) < crystal2.energy)',
              
        #'(bdt[3]>0.0 && bdtA[3]>-0.07)',
        # new bdt cut for V00-04-14
        '(((6.81804e-08*TMath::Exp(-101.856*bdt[3])+0.148344)-(4.04826*bdt[3])) < crystal3.energy)',

        #'(bdt[4]>-0.05 && (crystal4.energy + 40*bdt[4])>0 && bdtA[4]>-0.07)',
        # new bdt cut for V00-04-14
        '(((1.37726e-07*TMath::Exp(-111.842*bdt[4])+0.446818)-(2.13498*bdt[4])) < crystal4.energy)',
              
        '(bdt[5]>-0.2 && bdtA[5]>-0.07)',
              
        #'(bdt[6]>0.0 && (crystal6.energy + 20*bdt[6])>2.0 && bdtA[6]>-0.07)',
        # new bdt cut for V00-04-14
        '(((0.000315623*TMath::Exp(-75.3998*bdt[6])-0.313482)-(13.6302*bdt[6])) < crystal6.energy)',
              
        #'(bdt[7]>0.0 && (crystal7.energy + 20*bdt[7])>2.0 && bdtA[7]>-0.07)',
        # new bdt cut for V00-04-14
        '(((1.45552e-05*TMath::Exp(-88.7196*bdt[7])+0.566336)-(7.57773*bdt[7])) < crystal7.energy)',
        
        '(bdt[8]>-0.2 && bdtA[8]>-0.07)']

    bdtCut = newBDT[i]

    ### global noise cuts
    #===========================================================================
    coinc  = '(BLSVeto.isCoincident == 1)'
    muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
    charge = '(crystal'+str(i+1)+'.rqcn > -1)'
    
    ### old cut
    #nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
    ### new cut for v00-04-14 from Govinda
    nc     = '(pmt'+str(i+1)+'1.nc > 0 && pmt'+str(i+1)+'2.nc > 0)'
    
    ### old cut
    #t1     = '(pmt'+str(i+1)+'1.t1 > 2 && pmt'+str(i+1)+'2.t1 > 2)'
    ### new cut for v00-04-14 from Govinda
    t1     = '(pmt'+str(i+1)+'1.t1 > 0 && pmt'+str(i+1)+'2.t1 > 0)'
    
    noiseCut = '('+coinc+' && '+muons+' && '+charge+' && '+nc+' && '+t1+')'
    #===========================================================================
    
    if C == 'S':
        #lsveto = '(BLSVeto.Charge/143.8 < 20)'
        lsveto = '(BLSVeto.Charge/143.8 < 80)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc < 4) && '
        ### remove extra '&&'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        if E == 0:
            masterCut = TCut(bdtCut+' && '+noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        else:
            masterCut = TCut(noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)

        # lose cuts for C5 and C8
        if int(i) in [4, 7]:
            masterCut = TCut('('+lsveto+' && '+hits+')'+' && '+alphaCut)
        
        return masterCut
    
    elif C == 'M':
        #lsveto = '(BLSVeto.Charge/143.8 > 20)'
        lsveto = '(BLSVeto.Charge/143.8 > 80)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc > 4) || '
        ### remove extra '||'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        if E == 0:
            masterCut = TCut(bdtCut+' && '+noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        else:
            masterCut = TCut(noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)

        # lose cuts for C5 and C8
        if int(i) in [4, 7]:
            masterCut = TCut('('+lsveto+' || '+hits+')'+' && '+alphaCut)

        return masterCut

    else:
        print 'ERROR: I do not know what to do with channel -->', C
        print 'Available channels are [S]Single-hits, [M]Multi-hits'
        sys.exit()

