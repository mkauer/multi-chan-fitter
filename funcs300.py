#!/usr/bin/env python
######################################################################
# funcs300.py
# 
# Remove LSveto low-energy from the histograms
# 
# version: 2019-01-28
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
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

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs101 import *
from funcs_misc import *


def getInfo300(line, freuse=0, fchans='SM', fxstals=[]):
    
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
        info['rootfile'] = str(bits[-1])
    else:
        info['rootfile'] = 0
    #-----------------------------------------------------------------
    

    if 'D' in info['type']:
        
        ### data set
        info['tag'] = str(bits[2])
        
        ### processing version
        info['build'] = str(bits[3])
        
        ### data file
        info['file'] = str(bits[4])
        
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
            if not onCup():
                fpath = '/home/mkauer/COSINE/CUP/mc-fitting'+line
            else:
                fpath = line
            #print 'INFO: looking for -->', fpath
            if not os.path.exists(fpath):
                continue
            else:
                #print 'INFO: adding file -->', fpath
                tmpruntime = getDuration(fpath)
                if tmpruntime > 0:
                    runtime += tmpruntime
                    nfiles += chain.Add(fpath)
                else:
                    print 'WARNING: no run duration found for',fpath
        
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
                        edep, selection = calib100(i, E)
                    if info['build'] == 'V00-04-12':
                        edep, selection = calib300(i, E)
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


def buildMC300(info, mc):

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
        
        ### select correct paths for simulation files
        if   info['sim'] == 'v3.1.1': path1, path2, pushpasMC = mcPath101(info)
        elif info['sim'] == 'v4.0.1': path1, path2, pushpasMC = mcPath300(info)
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
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['xstl'])+'Crystal")')
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal")')
                        
                    elif 'internalsurf' in info['loca']:
                        # Estellas MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['xstl'])+'Crystal")')
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal")')
                        
                    elif 'naisurf' in info['loca']:
                        # Govindas MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['xstl'])+'Crystal")')
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Crystal")')
                        if pushpasMC: volumeCut = TCut('(primVolumeName == "NaIDet07Crystal")')
                        
                    elif info['loca'] == 'teflon':
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['xstl'])+'Teflon")')
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Teflon")')
                        
                    elif info['loca'] == 'teflonbulk':
                        # Govindas MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['xstl'])+'Teflon")')
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Teflon")')
                        if pushpasMC: volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                        
                    elif 'teflonsurf' in info['loca']:
                        # Estellas MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['xstl'])+'Teflon")')
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Teflon")')
                        if pushpasMC: volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                        
                    elif 'cusurf' in info['loca']:
                        # Estellas MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['xstl'])+'Case")')
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Case")')
                        
                    elif info['loca'] == 'cucase':
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['xstl'])+'Case")')
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['from'])+'Case")')
                        
                    elif info['loca'] == 'coppercase':
                        # Govindas MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(info['xstl'])+'Case")')
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
                        
                    ### I do not use this extpmt tag anymore
                    #elif info['loca'] == 'extpmt':
                    #    volumeCut = TCut('((primPMTid[0] != '+pmt1+') && (primPMTid[0] != '+pmt2+'))')
                    
                    elif info['loca'] == 'lsveto':
                        volumeCut = TCut('(primVolumeName == "lsveto")')
                    
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
                    if generated <= 0:
                        if int(info['xstl']) == int(info['from']):
                            print 'ERROR: no events generated for -->', info['floca'], info['isof'],\
                                'xstl', info['xstl'], 'from', info['from']
                            sys.exit()
                        else:
                            print 'WARNING: no events generated for -->', info['floca'], info['isof'],\
                                'xstl', info['xstl'], 'from', info['from']
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
        return edep, selection
    
    ### use Govindha's fits (does not have C5 or C8)
    hiE = [[1813, 316],
           [1485, 315],
           [1370, 323],
           [1213, 294],
           [1, 1],
           [1479, 286],
           [1069, 295],
           [1, 1]]
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
        edep = '(crystal'+str(i+1)+'.rqcD1_5)'
        selection = '('+edep+'-'+str(hiE[i][0])+')/('+str(hiE[i][1])+')'
        if i==4 or i==7: # C5 and C8 high energy
            edep = '(crystal'+str(i+1)+'.energyD)'
            selection = '('+edep+')'
    else:
        edep = '(crystal'+str(i+1)+'.energy)'
        selection = edep
        """
        if i==5: # C6 low energy
            edep = '(crystal'+str(i+1)+'.qc5)'
            selection = '('+edep+'-'+str(loE[i][0])+')/('+str(loE[i][1])+')'
        """
    return edep, selection


def mcPath101(info):
    
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
    
    pushpasMC = 0
    #if location in pushpas: pushpasMC = 1
    if info['floca'] in pushpas: pushpasMC = 1
    
    return path1, path2, pushpasMC


def mcPath300(info):
    
    ### use the full location-names 'floca'
    loc = info['floca']
    iso = info['isof']
    #Cx = 'C'+str(info['xstl'])
    Cx = 'C*'
    
    name = 'none'
    if loc == 'nai-surf-10um':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        path2 = '*-surfaceNaI-'+Cx+'-10um-'+iso+'-*root'
        
    elif loc == 'teflon-bulk':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        #path2 = '*-bulkTeflon-'+iso+'-*root'
        path2 = '*-bulkTeflon-'+Cx+'-'+iso+'-*root'
        
    elif loc == 'teflon-surf':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        #path2 = '*-surfaceTeflon-'+iso+'-*root'
        path2 = '*-surfaceTeflon-'+Cx+'-'+iso+'-*root'
        
    #elif loc == 'copper':
    #    path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
    #    path2 = '*-bulkCu-'+iso+'-*root'
    
    elif loc == 'coppercase':
        path1 = '/data/MC/COSINE/V4.0.1/reprocessed_SET1'
        #path2 = '*-bulkCu-'+iso+'-*root'
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
        
    else:
        print 'ERROR: no location identified for -->',loc
        sys.exit()
        
    pushpasMC = 0
    return path1, path2, pushpasMC


def updateBkgsFile300(xstals, bkgsfile, resultsfile, BF='BR'):
    """
    Only print out to file what was actually used in the fit
    """
    
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
    
    newbkgs = './plots/'+bkgsfile[:-4]+'_update.txt'
    output = open(newbkgs, 'w')
    output.write('# -*- sh -*-\n')
    
    print ''
    print 'INFO: Updating bkgsfile -->',bkgsfile
    print '      To a new bkgsfile -->',newbkgs
    print ''

    skip = 0
    for bline in bkgslines:

        if not bline:
            #output.write('\n')
            continue
        
        if bline.startswith('\"\"\"'):
            if skip == 0: skip = 1
            else: skip = 0
            #output.write(bline+'\n')
            continue
        if skip:
            #output.write(bline+'\n')
            continue
        
        if bline.startswith('#'):
            """
            if 'version' in bline:
                output.write('# NEW GENERATED backgrounds file from fit!\n\n')
                output.write(bline+'\n')
            else:
                output.write(bline+'\n')
            """
            continue
        
        bbits = filter(None, re.split("[ \s\t\n\r,:]+", bline.strip()))
        
        if int(bbits[1]) not in xstals:
            continue
        
        if 'F' not in bbits[0]:
            #for bits in bbits:
            #    output.write(bits+'\t')
            #output.write('\n')
            output.write(bline+'\n')
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
