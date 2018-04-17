#!/usr/bin/env python
######################################################################
# funcs92.py
# 
# New calibrations and bdt cuts for v00-04-04
# 
# Works with v92 and later versions
# 
# version: 2018-03-08
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + new combineOthers92()
# + new sortSimKeys92()
# + new sortDataKeys92()
# + new scaleSigs92()
# ~ and try offset in c6,c7 = [-8, -8]
# ~ try a slightly different offset in c6,c7 = [-4, -10]
# ~ tweaked calib92() for an offset in c6,c7 = [-6, -12]
# ~ heavily modified buildMC92() to accommodate Pushpa's C7 surface MC
# + include 'other' internal components
# + new buildMC92() to include internal 'others'
# ~ combine the 'other' pmt contributions for makePlots92()
# + makePlots92() so you can see every isotope separately
# + new bdt and bdtA cuts
# + new calibration numbers
# + build92(), buildData92(), calib92(), cutsBDT92()
# ~ import funcs91
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
from funcs91 import *


def build92(infile='backgrounds900.txt', others=1, vcut=1, freuse=0, fchans=0, fxstals=[]):

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
        #others = 0
        infos=[]
        if others and ('pmt' in info['key'] or 'internal' in info['key']):
        #if others and ('pmt' in info['key']):
            key = info['key']
            for i in range(8):
                #info['key'] = key+'-f'+str(i+1)
                newinfo = deepcopy(info)
                newinfo['key'] = key+'-f'+str(i+1)
                newinfo['from'] = str(i+1)
                infos.append(newinfo)
        elif not others and ('pmt' in info['key'] or 'internal' in info['key']):
        #elif not others and ('pmt' in info['key']):
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
                bkgs = buildMC92(info, bkgs, vcut)
            
            elif 'F' in info['type']:
                sigs = buildMC92(info, sigs, vcut)
            
            else:
                print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
                print info['line']
                continue
    
    return data, bkgs, sigs, runtime


def buildData92(info, data):
    
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
                    edep, selection = calib92(i, E)
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


def calib92(i,E):

    # from Pushpa on 2017-12-12
    
    # adc = crystalX.qc5
    # E = (adc + offset) / slope
    loE = [[1192.0, 9620.0],
           [2259.0, 9675.0],
           [2264.0, 9027.0],
           [1525.0, 9100.0],
           [   0.0, 3712.7],
           [ 735.1, 8755.0],
           [ 115.2, 9150.0],
           [   0.0, 2805.9]]

    hiE = [[1.00, 0],
           [1.00, 0],
           [1.00, 0],
           [1.00, 0],
           [0.95, 0],
           [1.00, -8], # I tweaked the C6 offset
           [1.00, -8], # I tweaked the C7 offset
           [0.73, 0]]
    
    if E:
        edep = '(crystal'+str(i+1)+'.energyD)'
        selection = '(('+edep+'*'+str(hiE[i][0])+')+'+str(hiE[i][1])+')'
    else:
        edep = '(crystal'+str(i+1)+'.qc5)'
        #selection = '('+edep+'*'+str(loE[i])+')'
        #if str(i) == str(6):
        #    selection = '(('+edep+'-'+str(c7[0])+')/'+str(c7[1])+')'
        selection = '(('+edep+'-'+str(loE[i][0])+')/'+str(loE[i][1])+')'

    return edep, selection


def cutsBDT92(i,C,E,edep,selection):
    """
    From: https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    """

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
    # i kinda don't like that bdte energy dependent cut
    bdtCuts  = '('+bdt+' && '+bdta+')'
    
    if C == 'A':
        return TCut(noiseCut+' && '+alphaCut)
        
    elif C == 'S':
        lsveto = '(BLSVeto.Charge/143.8 < 20)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc < 4) && '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        masterCut = TCut(bdtCuts+' && '+noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        if E:
            masterCut = TCut(noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
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
        if E:
            masterCut = TCut(noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        return masterCut
    
    else:
        print 'ERROR: I do not know what to do with channel -->',c
        print 'Available channels are [S]Single-hits, [M]Multi-hits'
        sys.exit()


def buildMC92(info, mc, vcut):

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
        
        ### Estella's Surface Pb210 MC
        if location == 'internal-surf-10um':
            path2 = info['isof']+'/set2/surf/10um/'+'*'+location+'-C'+str(info['xstl'])+'*root'
        if location == 'cu-surf-10um':
            path2 = info['isof']+'/set2/surf/cu-10um/'+'*'+location+'-C'+str(info['xstl'])+'*root'
        if location == 'teflon-surf-10um':
            path2 = info['isof']+'/set2/surf/teflon-10um/'+'*'+location+'-C'+str(info['xstl'])+'*root'
        
        ### Pushpa's Surface Pb210 MC
        pushpaC7Surf = ['nai-surf-10um', 'teflon-surf-2um', 'teflon-bulk']
        if location in pushpaC7Surf:
            if onCup(): path1 = '/data/COSINE/WORK/pushpa/sim/process/Crystal7/ResolA_ResolD/'
            else: path1 = '/home/mkauer/COSINE/CUP/mc-fitting/data/COSINE/WORK/pushpa/sim/process/Crystal7/ResolA_ResolD/'
            if location == 'nai-surf-10um':
                path2 = 'anal_lsvetofull-NaI-surface_C7-10um-Pb210-*.root'
            if location == 'teflon-surf-2um':
                path2 = 'anal_lsvetofull-Teflon-surface_2um-Pb210-*.root'
            if location == 'teflon-bulk':
                path2 = 'anal_lsvetofull-Teflon-surface_all-Pb210-*.root'
        
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
                    if location in pushpaC7Surf:
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
                        if location in pushpaC7Surf:
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
                        if location in pushpaC7Surf:
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
                        # Pushpas C7 Surface  MC
                        volumeCut = TCut('(primVolumeName == "NaIDet07Crystal")')
                    
                    elif 'cusurf' in info['loca']:
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Case")')
                    
                    elif 'teflonsurf' in info['loca']:
                        # Estella MC
                        #volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Teflon")')
                        # Pushpas C7 Surface MC
                        volumeCut = TCut('(primVolumeName == "NaIDet07Teflon")')
                    
                    elif info['loca'] == 'cucase':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Case")')
                    
                    elif info['loca'] == 'teflon':
                        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(i+1)+'Teflon")')

                    elif info['loca'] == 'teflonbulk':
                        # Pushpas C7 Surface  MC
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
                    if location in pushpaC7Surf:
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
                    if location in pushpaC7Surf:
                        selection = '((edep[6]*1000.) + (sigma*edep[6]*1000.*rng))'
                    chain.Draw(selection+' >> '+key, masterCut)
                    
                    detected = histo.GetEntries()
                    
                    #print '!!!! ', key, 'generated events =', generated
                    #print '!!!! ', key, 'detected events =', detected
                    #print '!!!! ', key, 'efficiency =', int(100*detected/generated), '%'
                    
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


def makePlots92(bkgs, combine, others, vcut):
    
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
                    if location == 'pmt' or location == 'internal':
                    #if location == 'pmt':
                        newkey += '-f'+crystal[-1]
                    newkey += '-c'+C
                    newkey += '-e'+str(E)
                    temp[i-1] = deepcopy(bkgs[newkey]['hist'])
                    
                    if combine and (location == 'pmt' or location == 'internal'):
                    #if combine and (location == 'pmt'):
                        for k in range(1,9):
                            if k != int(crystal[-1]):
                                newkey  = crystal
                                newkey += '-'+location
                                newkey += '-'+isotope
                                newkey += '-f'+str(k)
                                newkey += '-c'+C
                                newkey += '-e'+str(E)
                                try:
                                    temp[i-1].Add(deepcopy(bkgs[newkey]['hist']))
                                except:
                                    pass
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


def scaleSigs92(sigkeys, sigs, runtime=0):

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
            
        elif loca == 'cucase':
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

