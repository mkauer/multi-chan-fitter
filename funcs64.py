#!/usr/bin/env python
######################################################################
# funcs64.py
# 
# New data format and selection criteria
# 
# Works with v64 and later versions
# 
# version: 2017-06-07
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + added setGroup() to start organizing things into groups
# + added chan [A] cuts to cutsBDT()
# ~ cleaned up some print statements
# ~ itterate over channels after reading in rootfiles - this should
#   speed up the processing time
# + added onCup() and baseDir() to make paths a little easier
# ~ change data[key][runtime] to use real duration
# + added getDuration() to get duration of data rootfiles
# + added histparam64()
# + added makeTotal64(), makeResid64(), and globalParams()
# ~ add _GRND back to histo keys
# ~ scaling functions now get histogram params to determine keV per bin
# ~ tweaked histparam() to also return bins per keV
# + added getPars(), dataDRU64(), scaleBkgs64(), and scaleSigs64()
# + add BDT cuts in cutsBDT()
# ~ move old data cuts to cutsOld()
# + added set1.dat with list of run files
# ~ new data format in the backgrounds file
#
# email me: mkauer@physics.wisc.edu
######################################################################

import os,sys
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs63 import *


def histparam64(E):
    """
    Global default histogram parameters for:
    number of bins, min, max, and bins/keV
    """
    if E:
        hmin = 0
        hmax = 4000
        bpkv = 1
        bins = (hmax-hmin)*bpkv
    else:
        hmin = 0
        hmax = 200
        bpkv = 12
        bins = (hmax-hmin)*bpkv

    pars = [bins, hmin, hmax, bpkv]
    return pars


def getInfo64(line, freuse=0, fchans=0):
    
    info={}
    
    info['line'] = line
    
    bits = line.split()
    #bits = [x.strip() for x in bits]
    
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
    
    # crystal
    info['xstl'] = int(bits[2])
    
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
        
        ### background location
        info['loca'] = str(bits[3])

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
        info['fbnd'] = [float(bits[9]), float(bits[10])]

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


def setGroup(info):
    if info['loca'] == 'internal':
        if info['isof'] in ['I125','Te121m','Te123m','Te125m','Te127']:
            return 'cosmo'
        else: return 'internal'
    elif 'pmt' in info['loca']: return 'pmt'
    elif 'surf' in info['loca']: return 'surf'
    elif info['loca'] == 'lsveto': return 'lsveto'
    else: return 'none'


def build64(infile = 'backgrounds640.txt', freuse=0, fchans=0):

    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo64(line, freuse, fchans)
        
        if 'D' in info['type']:
            data = buildData64(info, data)
            
        elif 'B' in info['type']:
            bkgs = buildMC64(info, bkgs)
            
        elif 'F' in info['type']:
            sigs = buildMC64(info, sigs)
            
        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs


def buildData64(info, data):
    
    base = baseDir()
    """
    for c in info['chans']:
        info['chan'] = c
        if info['reuse']:
    """
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
                    edep, selection = calib60a(i,e)
                    # new calib
                    #edep, selection = calib60b(i,e)
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

    return data


def cutsOld(i,c,e):
    
    ### Pushpa's noise cut
    noiseCut = TCut('('+noiseCuts60(i,e)+')')
    
    if c == 'A':
        #chanCut = TCut('(1)')
        return noiseCut
    
    elif c == 'S':
        edepcuts = ''
        nclustercuts = ''
        for j in range(8):
            if j != i:
                #edepcuts += '('+edep+' <= 0.0) && '
                nclustercuts += '(crystal'+str(j+1)+'.nc < 4) && '

        ### remove extra '&&' or '||'
        #edepcuts = edepcuts[:-4]
        nclustercuts = nclustercuts[:-4]

        ### liquid scint veto cut - from Pushpa
        #lsvetocut = '(BLSveto.Charge/110. < 50.)'
        lsvetocut = '(BLSVeto_Charge/110. < 20.)'

        ### my cuts
        #chanCut = TCut('(('+edepcuts+') && ('+nclustercuts+') && ('+lsvetocut+'))')
        ### Pushpa cuts
        chanCut = TCut('(('+nclustercuts+') && ('+lsvetocut+'))')

        masterCut = TCut('('+chanCut.GetTitle()+' && '+noiseCut.GetTitle()+')')
        return masterCut
    
    elif c == 'M':
        edepcuts = ''
        nclustercuts = ''
        for j in range(8):
            if j != i:
                #edepcuts += '(crystal'+str(j+1)+'.'+str(edep)+' > 0.0) || '
                nclustercuts += '(crystal'+str(j+1)+'.nc > 4) || '

        ### remove extra '&&' or '||'
        #edepcuts = edepcuts[:-4]
        nclustercuts = nclustercuts[:-4]

        ### liquid scint veto cut - from Pushpa
        #lsvetocut = '(BLSveto.Charge/110. > 50.)'
        lsvetocut = '(BLSVeto_Charge/110. > 20.)'

        ### my cuts
        #chanCut = TCut('(('+edepcuts+') || ('+nclustercuts+') || ('+lsvetocut+'))')
        ### Pushpa cuts
        chanCut = TCut('(('+nclustercuts+') || ('+lsvetocut+'))')

        masterCut = TCut('('+chanCut.GetTitle()+' && '+noiseCut.GetTitle()+')')
        return masterCut

    else:
        print 'ERROR: I do not know what to do with channel -->',c
        print 'Available channels are [A]All-hits, [S]Single-hits, [M]Multi-hits'
        sys.exit()


def cutsBDT(i,c,e):
    """
    From: https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    """
    if c == 'A':
        coinc  = '(BLSVeto_isCoincident == 1)'
        muons  = '(BMuon_totalDeltaT0/1.e6 > 30)'
        bdt    = '(crystal'+str(i+1)+'.bdt > -0.03)'
        charge = '(crystal'+str(i+1)+'.rqcn >= 0)'
        nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
        
        #masterCut = TCut('(1)')
        masterCut = TCut(coinc+' && '+muons+' && '+bdt+' && '+charge+' && '+nc)
        return masterCut
        
    elif c == 'S':
        coinc  = '(BLSVeto_isCoincident == 1)'
        muons  = '(BMuon_totalDeltaT0/1.e6 > 30)'
        bdt    = '(crystal'+str(i+1)+'.bdt > -0.03)'
        charge = '(crystal'+str(i+1)+'.rqcn >= 0)'
        nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
        
        lsveto = '(BLSVeto_Charge/143.8 < 20)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.'+'nc'+' < 4) && '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'
        
        masterCut = TCut(coinc+' && '+muons+' && '+bdt+' && '+charge+' && '+nc+' && '+lsveto+' && '+hits)
        return masterCut
    
    elif c == 'M':
        coinc  = '(BLSVeto_isCoincident == 1)'
        muons  = '(BMuon_totalDeltaT0/1.e6 > 30)'
        bdt    = '(crystal'+str(i+1)+'.bdt > -0.03)'
        charge = '(crystal'+str(i+1)+'.rqcn >= 0)'
        nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
        
        lsveto = '(BLSVeto_Charge/143.8 > 20)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.'+'nc'+' > 4) || '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'

        masterCut = TCut(coinc+' && '+muons+' && '+bdt+' && '+charge+' && '+nc+' && '+'('+lsveto+' || '+hits+')')
        return masterCut
    
    else:
        print 'ERROR: I do not know what to do with channel -->',c
        print 'Available channels are [S]Single-hits, [M]Multi-hits'
        sys.exit()


def buildMC64(info, mc):

    base = baseDir()
    """
    for c in info['chans']:
        info['chan'] = c
        if info['reuse']:
    """
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

        if onCup():
            path1 = '/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
        else:
            path1 = '/home/mkauer/COSINE/CUP/mc-fitting/sim/newGeometry/'

        location = info['loca']
        if location == 'extpmt': location = 'pmt'

        ### 2nd level path to the specific files
        #path2 = info['isof']+'/set2/'+'*'+info['loca']+info['isof']+'*root'
        path2 = info['isof']+'/set2/'+'*'+location+info['isof']+'*root'

        # but there will be a few exceptions...
        if location == 'internalsurf':
            #path2 = info['isof']+'/set2/surf/10um/'+'*'+info['loca']+'*'+info['isof']+'*'+'C'+str(info['xstl'])+'*root'
            path2 = info['isof']+'/set2/surf/10um/'+'*'+location+'*'+info['isof']+'*'+'C'+str(info['xstl'])+'*root'

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
                        
                    #=====================================================================================
                    #=====================================================================================
                    
                    
                    pmt1 = str((int(info['xstl'])*2)-2)
                    pmt2 = str((int(info['xstl'])*2)-1)
                    #print pmt1,pmt2
                    
                    volumeCut = TCut('(1)')
                    if   info['loca'] == 'internal':
                        volumeCut = TCut('(primVolumeName == "'+volumeNames(i)+'")')

                    elif info['loca'] == 'internalsurf':
                        volumeCut = TCut('(primVolumeName == "'+volumeNames(i)+'")')
                    
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
                        #volumeCut = TCut('(primVolumeName != "")')
                        volumeCut = TCut('(1)')

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
                    else:
                        brokenChainCut = TCut('(groupNo > 1)')
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                    
                    
                    ### create a hist of the generated events
                    #---------------------------------------------------------------------
                    key2 = key+'_generated'
                    temp2 = TH1F(key2, 'generated', 1, 1, 2)
                    
                    totalCuts = TCut(volumeCut.GetTitle()+' && '+motherCut.GetTitle())
                    
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
                    ### need the volumeCut for the pmt selection
                    if 'pmt' in info['loca']:
                        masterCut = TCut('('+
                                         energyCut.GetTitle()+' && '+
                                         volumeCut.GetTitle()+' && '+
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

                    
        else:
            #print 'WARNING:', loca, isof, 'not found...'
            print 'ERROR: No MC files found for',info['loca'], info['isof'],'... quitting...'
            sys.exit()
    
    return mc


def scaleBkgs64(bkgs):
    
    for key in bkgs:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]
        
        #keVperBin = 1./float(getPars(bkgs[key]['hist'])[3])
        keVperBin = 1./float(bkgs[key]['pars'][3])
        day       = 86400. # in seconds
        
        xkgs      = cmass(int(x[-1])-1)
        pmts      = 2.
        extpmts   = 14.
        lskg      = 1800.

        generated = float(bkgs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        if loca == 'internal':
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'internalsurf':
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
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
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'airshield':
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'steel':
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./keVperBin)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        else:
            print "WARNING: No background scaling for  --> ", loca
            continue

    return bkgs


def scaleSigs64(sigkeys, sigs):

    for key in sigkeys:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]

        #keVperBin = 1./float(getPars(sigs[key]['hist'])[3])
        keVperBin = 1./float(sigs[key]['pars'][3])
        #print '!!!!!!!!!!!!!!!',float(sigs[key]['pars'][3])
        #keVperBin = float(sigs[key]['pars'][3])
        day       = 86400. # in seconds
        
        xkgs      = cmass(int(x[-1])-1)
        pmts      = 2.
        extpmts   = 14.
        lskg      = 1800.
        
        generated = float(sigs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        # verbose?
        V = 0
        
        if loca == 'internal':
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/kg'
        elif loca == 'internalsurf':
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/kg'
        elif loca == 'pmt':
            fitActivity = sigs[key]['fitscale'] * (1./pmts) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'extpmt':
            fitActivity = sigs[key]['fitscale'] * (1./extpmts) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'lsveto':
            fitActivity = sigs[key]['fitscale'] * (1./lskg) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'lsvetoair':
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'airshield':
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'steel':
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (keVperBin)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        else:
            print "WARNING: No signal scaling for  --> ", loca
            continue

    return sigs


def dataDRU64(data):
    """
    Scale data into DRU
    """
    for key in data:
        i = int(key.split('-')[0].split('x')[-1]) - 1
        days = float((data[key]['runtime'])/(60.*60.*24.))
        kgs = float(cmass(i))
        #keVperBin = 1./float(getPars(data[key]['hist'])[3])
        keVperBin = 1./float(data[key]['pars'][3])
        scale = float(1./(days*kgs*keVperBin))
        data[key]['hist'].Scale(scale)
        data[key]['druScale'] = scale
    return data


def makeTotal64(chan, E, par):
    total = []
    #par = histparam(E)
    for i in range(8):
        key  = 'x'+str(i)
        key += '-c'+chan
        key += '-e'+str(E)
        key += '-total'
        tot = TH1F(key, longNames(i), par[0], par[1], par[2])
        tot.SetLineColor(kGray+1)
        tot.SetMarkerColor(kGray+1)
        tot.SetLineWidth(1)
        total.append(tot)
    return total


def makeResid64(chan, E, par):
    resid = []
    #par = histparam(E)
    for i in range(8):
        key  = 'x'+str(i)
        key += '-c'+chan
        key += '-e'+str(E)
        key += '-resid'
        res = TH1F(key, longNames(i), par[0], par[1], par[2])
        res.SetLineColor(kBlack)
        res.SetMarkerColor(kBlack)
        res.SetLineWidth(1)
        resid.append(res)
    return resid


def getPars(hist):
    bins = hist.GetNbinsX()
    hmin = hist.GetBinLowEdge(1)
    hmax = hist.GetBinLowEdge(bins+1)
    bpkv = bins/(hmax-hmin)
    #print bins,hmin,hmax,bpkv
    pars = [bins, hmin, hmax, bpkv]
    return pars


def globalParams(data):
    par0=0
    par1=0
    for key in data:
        if '-e0' in key and not par0:
            par0 = getPars(data[key]['hist'])
        if '-e1' in key and not par1:
            par1 = getPars(data[key]['hist'])
        if par0 and par1:
            break
    params = [par0, par1]
    return params


def getDuration(rootfile):

    chain = TChain('ntp','')
    chain.Add(rootfile)
    entries = chain.GetEntries()

    ### only second precission
    #var = "eventsec"
    #norm = 1.

    ### has nano-sec precission
    var = "trgtime"
    norm = 1.e9
    
    chain.GetEntry(0)
    start = float(chain.GetLeaf(var).GetValue())
    chain.GetEntry(entries-1)
    stop = float(chain.GetLeaf(var).GetValue())
    duration = (stop-start)/norm
    #print var, start, stop, duration
    
    return duration


def baseDir():
    """
    Define your main basedir paths here
    """
    if onCup():
        return '/home/mkauer/mc-fitting/'
    else:
        return '/home/mkauer/COSINE/CUP/mc-fitting/'


def onCup():
    """
    Check to see if running on cup
    """
    import socket
    ### master.cunpa.ibs
    if 'cunpa' in str(socket.gethostname()):
        return 1
    else:
        return 0

