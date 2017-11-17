#!/usr/bin/env python
######################################################################
# funcs80.py
# 
# New calibrations for all crystals
# 
# Works with v80 and later versions
# 
# version: 2017-11-16
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ tweak getInfo80() to be compatible with no activity error
# + added getInfo80() so I can select which xstals to build
# ~ smearing C8 hi energy resolution more
# + added resol80() and buildMC80()
# ~ fix cutsBDT80() to include the alpha cut
# + add cutsBDT80()
# ~ fixed a bug with run duration
# + new calibrations in calib80()
# + add build80(), buildData80(), calib80()
# ~ import funcs71
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
from funcs71 import *


def build80(infile='backgrounds800.txt', freuse=0, fchans=0, fxstals=[]):

    data = {}
    bkgs = {}
    sigs = {}
    runtime = 0
    
    for line in readFile(infile):
        info = getInfo80(line, freuse, fchans, fxstals)
        if len(info) == 0:
            continue
        
        if 'D' in info['type']:
            data, runtime = buildData80(info, data)
            
        elif 'B' in info['type']:
            bkgs = buildMC80(info, bkgs)
            
        elif 'F' in info['type']:
            sigs = buildMC80(info, sigs)
            
        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs, runtime


def getInfo80(line, freuse=0, fchans=0, fxstals=[]):
    
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
    #-----------------------------------------------------------------
    
    # reuse rootfile name
    if 'root' in str(bits[-1]):
        temp = str(bits[-1])
        for i in range(2):
            if '.' in temp[0] or '/' in temp[0]:
                temp = temp[1:]
        info['rootfile'] = temp
    else:
        info['rootfile'] = 0


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
        info['erro'] = float(bits[8])

        ### fit bounds
        #info['fbnd'] = [float(bits[9]), float(bits[10])]
        info['fbnd'] = [float(bits[-3]), float(bits[-2])]

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


def buildData80(info, data):
    
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
                    edep, selection = calib80(i,e)
                    #-----------------------------------------------------------------------
                    
                    
                    # DEFINE CUTS
                    #-----------------------------------------------------------------------
                    masterCut = cutsBDT80(i,c,e)
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


def calib80(i, E=0):

    # from Pushpa on 2017-10-17
    
    # adc = crystalX.qc5
    # E = adc * cal
    loE = [
        1.049981e-04,
        1.046627e-04,
        1.127174e-04,
        1.099834e-04,
        2.844480e-04,
        1.149712e-04,
        1.0,
        3.851424e-04        
    ]

    ### c7 low E is special
    c7 = [566., 9255.]

    hiE = [
        1.0,
        1.0,
        1.0,
        1.0,
        0.95,
        1.0,
        1.0,
        0.73
    ]
    
    if E:
        edep = '(crystal'+str(i+1)+'.energyD)'
        selection = '('+edep+'*'+str(hiE[i])+')'
    else:
        edep = '(crystal'+str(i+1)+'.qc5)'
        selection = '('+edep+'*'+str(loE[i])+')'
        if str(i) == str(6):
            selection = '(('+edep+'-'+str(c7[0])+')/'+str(c7[1])+')'

    return edep, selection


def cutsBDT80(i,c,e):
    """
    From: https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    """

    # values for the alpha cut
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
    
    alphaCut = '( ! ((crystal'+str(i+1)+'.energyD > 1000.) && ((pmt'+str(i+1)+'1.rqtD1_5+pmt'+str(i+1)+'2.rqtD1_5)/2. < '+str(alpha[i])+')) )'
    
    # BDT cut values
    cut = [
        -0.05,
         0.02,
         0.00,
        -0.02,
        -0.10,
         0.03,
         0.03,
        -0.15,
    ]
    
    ### global noise cuts
    coinc  = '(BLSVeto.isCoincident == 1)'
    muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
    bdt    = '(bdt['+str(i+1)+'] > '+str(cut[i])+')'
    charge = '(crystal'+str(i+1)+'.rqcn > -1)'
    nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
    
    noiseCut = coinc+' && '+muons+' && '+bdt+' && '+charge+' && '+nc
    
    if c == 'A':
        return TCut(noiseCut+' && '+alphaCut)
        
    elif c == 'S':
        lsveto = '(BLSVeto.Charge/143.8 < 20)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc < 4) && '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'
        
        masterCut = TCut(noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        return masterCut
    
    elif c == 'M':
        lsveto = '(BLSVeto.Charge/143.8 > 20)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc > 4) || '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'
        
        masterCut = TCut(noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        return masterCut
    
    else:
        print 'ERROR: I do not know what to do with channel -->',c
        print 'Available channels are [S]Single-hits, [M]Multi-hits'
        sys.exit()


def buildMC80(info, mc):

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
                        volumeCut = TCut('((primVolumeName == "SteelSupport") || (primVolumeName == "SteelSupportTop"))')
                        #volumeCut = TCut('((primVolumeName == "SteelSupport"))')
                        #volumeCut = TCut('((primVolumeName == "SteelSupportTop"))')
                        #volumeCut = TCut('((primVolumeName == "Welding"))')
                        #volumeCut = TCut('((primVolumeName != ""))')
                        
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


def resol80(i, E=0):
    """
    https://cupwiki.ibs.re.kr/Kims/NaICalibration
    """
        
    # res = p[0]/sqrt(x) + p[1]

    # low energy
    loEresol = [
        [0.2413,  0.01799],
	[0.2951,  0.01427],
	[0.3106,  0.007894],
	[0.3894, -0.001437],
        [0.8, 0], # tweaking C5
	[0.3620,  0.0006355],
	[0.3042,  0.009784],
        [1.3, 0] # tweaking C8
    ]

    # high energy
    hiEresol = [
        [0.6729, 0.009374],
	[0.6531, 0.006627],
	[0.5926, 0.009506],
	[0.7227, 0.004790],
        [1.7, 0], # tweaking C5
	[0.6498, 0.009670],
	[0.7034, 0.007812],
        [4., 0] # tweaking C8
    ]
    
    
    if E:
        p0, p1 = hiEresol[int(i)]
    else:
        p0, p1 = loEresol[int(i)]

    selection = '(('+str(float(p0))+'/sqrt(edep['+str(i)+']*1000.)) + '+str(float(p1))+')'
    return selection

