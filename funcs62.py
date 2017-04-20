#!/usr/bin/env python
######################################################################
# funcs62.py
# 
# New background chsp change to GRND and more
# 
# Works with v62 and later versions
# 
# version: 2017-04-19
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------

# ~ change build62() to use buildMC62()
# + new build62() to use getInfo62()
# + new getInfo62() for the 'chsp' 'GRND' option
# 
# email me: mkauer@physics.wisc.edu
######################################################################
# 
# Where is MC?
# /data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/
# 
# Where is the processed data?
# Right now I'm using:
# /data/COSINE/NTP/phys/V00-02-03_MERGED
# and run ntp_I001546*
# 
######################################################################

import os,sys
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs61 import *


def getInfo62(line, freuse=0, fchans=0):
    
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
        info['rootfile'] = str(bits[-1])
    else:
        info['rootfile'] = 0


    if 'D' in info['type']:
        
        ### run number
        info['run'] = str(bits[4])

        ### processing version
        info['build'] = str(bits[5])

        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-data'
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
        if info['chsp'] == 'GRND': key += '-'+info['chst']
        else: key += '-'+info['chst']+'_'+info['chsp']
        info['key'] = key
        
    return info


def build62(infile = 'backgrounds64.txt', freuse=0, fchans=0):

    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo62(line, freuse, fchans)

        if 'D' in info['type']:
            data = buildData60(info, data)
            ### buildData61() is for the new V00-03-02 processing path
            ### but the processed data is not right, yet...
            #data = buildData61(info, data)

        elif 'B' in info['type']:
            #bkgs = buildMC61(info, bkgs)
            bkgs = buildMC62(info, bkgs)

        elif 'F' in info['type']:
            #sigs = buildMC61(info, sigs)
            sigs = buildMC62(info, sigs)

        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs


def groupNum62(info):
    
    # internal U238 "groupNo"
    # -----------------------------
    # 11: U238  -> Th230
    # 12: Th230 -> Ra226
    # 13: Ra226 -> Rn222
    # 14: Rn222 -> Pb210
    # 15: Pb210 -> ground
    
    # internal Th232 "groupNo"
    # -----------------------------
    # 21: Th232 -> Ra228
    # 22: Ra228 -> Th228
    # 23: Th228 -> ground
    
    """
    brokenChainCut = TCut('(1)')
    if info['chst'] == 'U238' and info['chsp'] == 'Rn222':
        brokenChainCut = TCut('((groupNo >= 11) && (groupNo <= 14))')
    if info['chst'] == 'Pb210' and info['chsp'] == 'Pb210':
        brokenChainCut = TCut('(groupNo == 15)')
    """             

    if   info['chst'] == 'U238':  start = 11
    elif info['chst'] == 'Th230': start = 12
    elif info['chst'] == 'Ra226': start = 13
    elif info['chst'] == 'Rn222': start = 14
    elif info['chst'] == 'Pb210': start = 15
    elif info['chst'] == 'Th232': start = 21
    elif info['chst'] == 'Ra228': start = 22
    elif info['chst'] == 'Th228': start = 23
    else:
        #print 'WARNING: do not know broken chain starting with -->', info['chst']
        #return [-1, -1]
        return 0
    
    if info['chsp'] == 'GRND':
        stop = 100 # just guessing something safe
        return [start, stop]

    if   info['chsp'] == 'U238':  stop = 11
    elif info['chsp'] == 'Th230': stop = 12
    elif info['chsp'] == 'Ra226': stop = 13
    elif info['chsp'] == 'Rn222': stop = 14
    elif info['chsp'] == 'Pb210': stop = 15
    elif info['chsp'] == 'Th232': stop = 21
    elif info['chsp'] == 'Ra228': stop = 22
    elif info['chsp'] == 'Th228': stop = 23
    else:
        print 'WARNING: do not know broken chain stopping with -->', info['chsp']
        return 0

    return [start, stop]


def buildMC62(info, mc):

    for c in info['chans']:
        info['chan'] = c
    
        if info['reuse']:
            if info['rootfile']:
                rootfile = info['rootfile']
            else:
                print 'WARNING: no rootfile specified in backgrounds file'
                return mc

            if not os.path.exists(rootfile):
                print 'WARNING: rootfile not found -->', rootfile
                return mc

            rfile = TFile(rootfile, "READ")

            for e in range(2):
                key  = info['key']
                key += '-c'+info['chan']
                key += '-e'+str(e)
                
                try:
                    mc[key] = {}
                    mc[key]['info'] = info
                    mc[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                    #mc[key]['hist'].Sumw2()
                except:
                    print "WARNING: could not find hist -->",key

                try:
                    mc[key]['generated_hist'] = deepcopy(TH1F(rfile.Get(key+'_generated')))
                    #mc[key]['generated'] = mc[key]['generated_hist'].GetBinContent(1)
                    mc[key]['generated'] = mc[key]['generated_hist'].GetEntries()
                except:
                    print "WARNING: could not find generated_hist -->",key+'_generated'

                
        else:
            local = amLocal()
            
            ### top level path to the MC
            path1 = 0
            if local: path1 = '/home/mkauer/COSINE/CUP/mc-fitting/sim/newGeometry/'
            else: path1 = '/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
            
            ### 2nd level path to the specific files
            #path2 = info['isof']+'/set2/'+'*'+info['loca']+'*'+info['isof']+'*root'
            path2 = info['isof']+'/set2/'+'*'+info['loca']+info['isof']+'*root'

            # but there will be a few exceptions...
            if info['loca'] == 'internalsurf':
                #path2 = info['isof']+'/set2/surf/10um/'+'*'+info['loca']+'*'+'C'+str(info['xstl'])+'*'+info['isof']+'*root'
                path2 = info['isof']+'/set2/surf/10um/'+'*'+info['loca']+'*'+info['isof']+'*'+'C'+str(info['xstl'])+'*root'

            
            chain  = TChain("MC","")
            nfiles = 0
            nfiles = chain.Add(path1+path2)
                        
            if nfiles > 0:
                print nfiles,'MC files found for', info['loca'], info['isof']
                for e in range(2):

                    i = info['xstl'] - 1

                    key  = info['key']
                    key += '-c'+info['chan']
                    key += '-e'+str(e)
                    
                    mc[key] = {}
                    mc[key]['info'] = info
                    
                    par = histparam(e)
                    histo = TH1F(key, longNames(i), par[0], par[1], par[2])

                    energyCut = TCut('(edep['+str(i)+']*1000. > 0.0)')

                    if info['chan'] == 'A':
                        chanCut = TCut('(1)')

                    elif info['chan'] == 'S':
                        ### main single/multi hit cut
                        edepcuts = '((singleHitTag['+str(i)+'] > 0.0) && (multipleHitTag['+str(i)+'] < 0.0))'

                        ### liquid scint veto cut - from Estella
                        lsvetocut = '(edep[8]*1000. < 20.0)'

                        #chanCut = TCut('(('+edepcuts+'))')
                        chanCut = TCut('(('+edepcuts+') && ('+lsvetocut+'))')

                    elif info['chan'] == 'M':
                        ### main single/multi hit cut
                        edepcuts = '((singleHitTag['+str(i)+'] < 0.0) && (multipleHitTag['+str(i)+'] > 0.0))'

                        ### liquid scint veto cut - from Estella
                        lsvetocut = '(edep[8]*1000. > 20.0)'

                        #chanCut = TCut('(('+edepcuts+'))')
                        chanCut = TCut('(('+edepcuts+') || ('+lsvetocut+'))')

                    else:
                        print 'ERROR: I do not know what to do with channel -->',info['chan']
                        print 'Available channels are [A]All-hits, [S]Single-hits, [M]Multi-hits'
                        sys.exit()


                    volumeCut = TCut('(1)')
                    if   info['loca'] == 'internal':
                        volumeCut = TCut('(primVolumeName == "'+volumeNames(i)+'")')

                    elif info['loca'] == 'internalsurf':
                        volumeCut = TCut('(primVolumeName == "'+volumeNames(i)+'")')
                    
                    
                    ### This code selects the correct pair of PMTs corresponding to the crystal of interest
                    #=========================================================================================
                    elif info['loca'] == 'pmt':
                        volumeCut = TCut('((primPMTid[0] == '+str((int(info['xstl'])*2)-2)+') || (primPMTid[0] == '+str((int(info['xstl'])*2)-1)+'))')
                    #=========================================================================================
                    
                        
                    elif info['loca'] == 'lsveto':
                        volumeCut = TCut('(primVolumeName == "lsveto")')

                    elif info['loca'] == 'lsvetoair':
                        volumeCut = TCut('((primVolumeName == "DetPMTCover") || (primVolumeName == "DetPMTEnvelope"))')

                    elif info['loca'] == 'airshield':
                        volumeCut = TCut('(primVolumeName == "LSVetoAirRoom")')

                    elif info['loca'] == 'steel':
                        # Not sure this works the way I think it should?
                        # I want all volumes except for the blank volume
                        volumeCut = TCut('(primVolumeName != "")')

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
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    """
                    brokenChainCut = TCut('(1)')
                    if info['chst'] == 'U238' and info['chsp'] == 'Rn222':
                        brokenChainCut = TCut('((groupNo >= 11) && (groupNo <= 14))')
                    if info['chst'] == 'Pb210' and info['chsp'] == 'Pb210':
                        brokenChainCut = TCut('(groupNo == 15)')
                    """
                    brokenChainCut = TCut('(1)')
                    groupCuts = groupNum62(info)
                    if groupCuts:
                        print 'INFO:', info['isof'], 'chan start =', info['chst'], 'chan stop =', info['chsp'], 'groupNo start =', groupCuts[0], 'groupNo stop =', groupCuts[1]
                        brokenChainCut = TCut('((groupNo >= '+str(groupCuts[0])+') && (groupNo <= '+str(groupCuts[1])+'))')
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
                    print '!!!! ',key2,'',generated,' generated events  (',entries,'entries )' 
                    
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

