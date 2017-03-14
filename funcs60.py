#!/usr/bin/env python
######################################################################
# funcs60.py
# 
# Get all needed funcs into here
# 
# Works with v60 and later versions
# 
# version: 2017-03-14
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# ~ rename everything to 60 version and add other funcs as needed
# + try to move all main funcs to here funcs60.py - it's becoming too
#   confusing to find all funcs in all various versions of funcs
# + building off funcs52.py
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


def getInfo60(line, freuse=0, fchans=0):
    
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
        if info['chst'] == info['chsp']: key += '-'+info['chst']
        else: key += '-'+info['chst']+'_'+info['chsp']
        #key += '-c'+info['chan']
        info['key'] = key
        
    return info


def build60(infile = 'backgrounds60.txt', freuse=0, fchans=0):

    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo60(line, freuse, fchans)
        if 'D' in info['type']:
            data = buildData60(info, data)
        elif 'B' in info['type']:
            bkgs = buildMC60(info, bkgs)
        elif 'S' in info['type']:
            sigs = buildMC60(info, sigs)
        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs


def dataDRU60(data):
    """
    Scale data into DRU
    """
        
    # this will eventually have to be made dynamic
    keVperBin = float(1.)
    
    for key in data:
        
        i = int(key.split('-')[0].split('x')[-1]) - 1
        days = float((data[key]['runtime'])/(60.*60.*24.))
        kgs = float(cmass(i))
        scale = float(1./(days*kgs*keVperBin))
        
        data[key]['hist'].Scale(scale)
        data[key]['druScale'] = scale
        
    return data


def buildData60(info, data):

    for c in info['chans']:
        info['chan'] = c
        
        if info['reuse']:

            if info['rootfile']:
                rootfile = info['rootfile']
                #print 'INFO: using rootfile -->',rootfile
            else:
                print 'WARNING: no rootfile specified in backgrounds file'
                return data

            if not os.path.exists(rootfile):
                print 'WARNING: rootfile not found -->', rootfile
                return data

            rfile = TFile(rootfile, "READ")

            for e in range(2):
                key  = info['key']
                key += '-c'+info['chan']
                key += '-e'+str(e)
                #print key
                try:
                    data[key] = {}
                    data[key]['info'] = info
                    data[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                    #data[key]['hist'].Sumw2()
                except:
                    print "WARNING: could not find hist -->",key

                try:
                    data[key]['subruns_hist'] = deepcopy(TH1F(rfile.Get(key+'_subruns')))
                    data[key]['subruns'] = data[key]['subruns_hist'].GetBinContent(1)
                except:
                    print "WARNING: could not find subruns_hist -->",key+'_subruns'

                try:
                    data[key]['runtime_hist'] = deepcopy(TH1F(rfile.Get(key+'_runtime')))
                    data[key]['runtime'] = data[key]['runtime_hist'].GetBinContent(1)
                except:
                    print "WARNING: could not find runtime_hist -->",key+'_runtime'

                #except:
                #    del data[key]
                #    print "WARNING: could not find hist -->",key
                #    return data

            #return data

        else:
            local = amLocal()

            ### will need to do this dynamically eventually
            ### run 1324 has 1 hour subruns
            ### run 1544 has 2 hour subruns
            ### run 1546 has 2 hour subruns
            subrunTime = float(2.) # in hours
            
            path1 = 0
            if local: path1 = '/home/mkauer/COSINE/CUP/mc-fitting/data/phys/'
            else:     path1 = '/data/COSINE/NTP/phys/'
            path2 = str(info['build'])+'_MERGED/'
            path3 = 'ntp_I*'+str(info['run'])+'*root*'

            chain = TChain("ntp","")
            nfiles=0
            nfiles = chain.Add(path1+path2+path3)

            if nfiles > 0:
                print nfiles,'data files found'
                for e in range(2):
                    
                    par = histparam(e)
                    i = info['xstl'] - 1
                    
                    key  = info['key']
                    key += '-c'+info['chan']
                    key += '-e'+str(e)
                    
                    run = int(info['run'])
                    
                    data[key] = {}
                    data[key]['info'] = info
                    
                    histo = TH1F(key, longNames(i), par[0], par[1], par[2])


                    ### old calib
                    edep, selection = calib60a(i,e)
                    ### new calib
                    #edep, selection = calib60b(i,e)


                    ### reminder - [A]All-hits, [S]Single-hits, [M]Multi-hits
                    #-------------------------------------------------------------------------
                    if info['chan'] == 'A':
                        chanCut = TCut('(1)')

                    elif info['chan'] == 'S':
                        edepcuts = ''
                        nclustercuts = ''
                        for j in range(8):
                            if j != i:
                                edepcuts += '('+edep+' <= 0.0) && '
                                nclustercuts += '(crystal'+str(j+1)+'.'+'nc'+' < 4) && '

                        ### remove extra '&&' or '||'
                        edepcuts = edepcuts[:-4]
                        nclustercuts = nclustercuts[:-4]

                        ### liquid scint veto cut - from Pushpa
                        #lsvetocut = '(BLSveto.Charge/110. < 50.)'
                        lsvetocut = '(BLSveto.Charge/110. < 20.)'

                        ### my cuts
                        #chanCut = TCut('(('+edepcuts+') && ('+nclustercuts+') && ('+lsvetocut+'))')
                        ### Pushpa cuts
                        chanCut = TCut('(('+nclustercuts+') && ('+lsvetocut+'))')


                    elif info['chan'] == 'M':
                        edepcuts = ''
                        nclustercuts = ''
                        for j in range(8):
                            if j != i:
                                edepcuts += '(crystal'+str(j+1)+'.'+str(edep)+' > 0.0) || '
                                nclustercuts += '(crystal'+str(j+1)+'.'+'nc'+' >= 4) || '

                        ### remove extra '&&' or '||'
                        edepcuts = edepcuts[:-4]
                        nclustercuts = nclustercuts[:-4]

                        ### liquid scint veto cut - from Pushpa
                        #lsvetocut = '(BLSveto.Charge/110. > 50.)'
                        lsvetocut = '(BLSveto.Charge/110. > 20.)'

                        ### my cuts
                        #chanCut = TCut('(('+edepcuts+') || ('+nclustercuts+') || ('+lsvetocut+'))')
                        ### Pushpa cuts
                        chanCut = TCut('(('+nclustercuts+') || ('+lsvetocut+'))')

                    else:
                        print 'ERROR: I do not know what to do with channel -->',info['chan']
                        print 'Available channels are [A]All-hits, [S]Single-hits, [M]Multi-hits'
                        sys.exit()


                    ### Pushpa's noise cut
                    noiseCut = TCut('('+noiseCuts60(i,e)+')')


                    ### combine all cuts
                    masterCut = TCut('('
                                     +chanCut.GetTitle()+' && '
                                     +noiseCut.GetTitle()
                                     +')')
                    
                    
                    ###-----------------------------------------------------------------------
                    chain.Draw(selection+' >> '+key, masterCut)
                    ###-----------------------------------------------------------------------

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

                    runtime = subruns*subrunTime*60.*60.
                    key3 = key+'_runtime'
                    temp3 = TH1F(key3,'runtime',1,0,1)
                    temp3.SetBinContent(1, runtime)
                    data[key]['runtime_hist'] = temp3
                    data[key]['runtime'] = runtime

            else:
                print 'ERROR: No data files found... quitting...'
                sys.exit()

    return data


def buildMC60(info, mc):

    for c in info['chans']:
        info['chan'] = c
    
        if info['reuse']:
            if info['rootfile']:
                rootfile = info['rootfile']
            else:
                print 'WARNING: no rootfile specified in backgrounds file'
                return data

            if not os.path.exists(rootfile):
                print 'WARNING: rootfile not found -->', rootfile
                return data

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
                    mc[key]['generated'] = mc[key]['generated_hist'].GetBinContent(1)
                except:
                    print "WARNING: could not find generated_hist -->",key+'_generated'

                
        else:
            local = amLocal()
            
            ### top level path to the MC
            path1 = 0
            if local: path1 = '/home/mkauer/COSINE/CUP/mc-fitting/sim/newGeometry/'
            else: path1 = '/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
            
            ### 2nd level path to the specific files
            path2 = info['isof']+'/set2/'+'*'+info['loca']+'*'+info['isof']+'*root'

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

                    elif info['loca'] == 'pmt':
                        volumeCut = TCut('(primVolumeName == "phys_pmt")')

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

                    #brokenChainCut = 0
                    brokenChainCut = TCut('(1)')
                    if info['chst'] == 'U238' and info['chsp'] == 'Rn222':
                        brokenChainCut = TCut('((groupNo >= 11) && (groupNo <= 14))')
                    if info['chst'] == 'Pb210' and info['chsp'] == 'Pb210':
                        brokenChainCut = TCut('(groupNo == 15)')

                    key2 = key+'_generated'
                    temp2 = TH1F(key2, 'generated',1,0,1)

                    totalCuts = TCut(volumeCut.GetTitle()+'&&'+motherCut.GetTitle())
                    chain.Draw('primVolumeName >> '+key2, totalCuts)

                    generated = temp2.GetEntries()


                    ### test out resolution smearing
                    ###-------------------------------------------------------------------------------
                    ### using Box-Muller? method here for "rng" alias
                    chain.SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
                    resolFunc = resol60(i,e)
                    chain.SetAlias('sigma', resolFunc)

                    masterCut = TCut('('+
                                     energyCut.GetTitle()+' && '+
                                     #volumeCut.GetTitle()+' && '+
                                     eventTypeCut.GetTitle()+' && '+
                                     brokenChainCut.GetTitle()+' && '+
                                     chanCut.GetTitle()
                                     +')')
                    selection = '((edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng))'


                    ###-------------------------------------------------------------------------------
                    chain.Draw(selection+' >> '+key, masterCut)
                    ###-------------------------------------------------------------------------------


                    detected = histo.GetEntries()

                    #histo.SetLineColor(kBlack)
                    #histo.SetMarkerColor(kBlack)
                    #histo.SetLineWidth(1)

                    mc[key]['hist'] = histo
                    #mc[key]['hist'].Sumw2()

                    mc[key]['generated_hist'] = temp2
                    mc[key]['generated'] = generated

            else:
                #print 'WARNING:', loca, isof, 'not found...'
                print 'ERROR: No MC files found for',info['loca'], info['isof'],'... quitting...'
                sys.exit()

    return mc


def scaleBkgs60(bkgs):
    
    for key in bkgs:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]
        
        #druscale = data[x+'-data-'+e]['druScale']
        #runtime  = data[x+'-data-'+e]['runtime']
        
        kev  = 1.     # keV/bin
        day  = 86400. # in seconds
        
        xkgs = cmass(int(x[-1])-1)
        pmts = 16.
        lskg = 1800.
        
        generated = float(bkgs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        if loca == 'internal':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'internalsurf':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'pmt':
            #scale = bkgs[key]['info']['acti'] * (pmts) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (pmts) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'lsveto':
            #scale = bkgs[key]['info']['acti'] * (lskg) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (lskg) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'lsvetoair':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'airshield':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        elif loca == 'steel':
            #scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (runtime) * (1./generated) * (druscale)
            scale = bkgs[key]['info']['acti'] * (xkgs) * (1./1000) * (1./generated) * (day) * (1./xkgs) * (1./kev)
            bkgs[key]['hist'].Scale(scale)
            bkgs[key]['scale'] = scale
        else:
            print "WARNING: No background scaling for  --> ", loca
            continue

    return bkgs


def scaleSigs60(sigkeys, sigs):

    for key in sigkeys:
        x = key.split('-')[0]
        loca = key.split('-')[1]
        e = key.split('-')[-1]

        kev  = 1.     # keV/bin
        day  = 86400. # in seconds
        
        xkgs = cmass(int(x[-1])-1)
        pmts = 16.
        lskg = 1800.
        
        generated = float(sigs[key]['generated'])
        if generated < 1:
            print "WARNING: 0 events generated for -->", key
            continue
        
        # verbose?
        V = 0
                
        if loca == 'internal':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/kg'
        elif loca == 'internalsurf':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/kg'
        elif loca == 'pmt':
            #fitActivity = sigs[key]['fitscale'] * (1./pmts) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./pmts) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'lsveto':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./lskg) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'lsvetoair':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'airshield':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        elif loca == 'steel':
            #fitActivity = sigs[key]['fitscale'] * (1./kgs) * (1000.) * (1./runtime) * (generated) * (1./druscale)
            fitActivity = sigs[key]['fitscale'] * (1./xkgs) * (1000.) * (generated) * (1./day) * (xkgs) * (kev)
            sigs[key]['info']['acti'] = fitActivity
            if V: print '!!!!', key, sigs[key]['info']['acti'],'mBq/pmt'
        else:
            print "WARNING: No background scaling for  --> ", loca
            continue

    return sigs


def makeTotal60(chan, E):
    total = []
    par = histparam(E)
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


def makeResid60(chan, E):
    resid = []
    par = histparam(E)
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


def calib60a(i, E=0):
    """
    Return the crystal calibrations
    """
    # from Pushpa
    # this is the old calib

    # adc = crystalX.energyD
    # E = (adc+[0])*[1]
    hiEcalib = [
        [0.0, 1.000],
	[0.0, 0.992],
	[0.0, 1.000],
	[0.0, 1.000],
	[0.0, 0.966],
	[0.0, 1.014],
	[0.0, 0.984],
	[0.0, 1.000]
    ]
    
    # adc = crystalX.qc5
    # E = (adc+[0])*[1]
    loEcalib = [
        [0.0, 1.06000e-4],
	[0.0, 1.05676e-4],
	[0.0, 1.15900e-4],
	[0.0, 1.12700e-4],
	[0.0, 2.86907e-4],
	[0.0, 1.17634e-4],
	[0.0, 1.11900e-4],
	[0.0, 3.84263e-4]
    ]
    
    if E:
        edep = '(crystal'+str(i+1)+'.energyD)'
        b, m = hiEcalib[int(i)]
    else:
        edep = '(crystal'+str(i+1)+'.qc5)'
        b, m = loEcalib[int(i)]
        
    selection = '(('+edep+'+'+str(b)+')*'+str(m)+')'
    return edep, selection


def calib60b(i, E=0):
    """
    Return the crystal calibrations
    """
    # from Pushpa
    # https://cupwiki.ibs.re.kr/Kims/NaICalibration?validation_key=fae031e9908735aeff81aacfcbf83931
    
    # adc = pmtX1.rqcD1_5+pmtX2.rqcD1_5
    # E = (adc+[0])*[1]
    hiEcalib = [
        [3942.0, 1./232.5],
 	[4327.0, 1./241.9],
 	[3834.0, 1./262.1],
 	[3909.0, 1./231.3],
 	[   0.0, 1./125.0],
 	[-231.0, 1./175.7],
 	[4769.0, 1./260.9],
        [   0.0, 1./ 40.2]
    ]
    
    # adc = crystalX.qc5
    # E = (adc+[0])*[1]
    loEcalib = [
        [0.0, 0.000109373],
 	[0.0, 0.000108797],
 	[0.0, 0.000117414],
 	[0.0, 0.000114566],
 	[0.0, 0.000342601],
 	[0.0, 0.000118527],
 	[0.0, 0.000113608],
 	[0.0, 0.000500613]
    ]
    
    if E:
        edep = '(pmt'+str(i+1)+'1.rqcD1_5 + pmt'+str(i+1)+'2.rqcD1_5)'
        b, m = hiEcalib[int(i)]
    else:
        edep = '(crystal'+str(i+1)+'.qc5)'
        b, m = loEcalib[int(i)]
    
    selection = '(('+edep+'+'+str(b)+')*'+str(m)+')'
    return edep, selection


def resol60(i, E=0):
    """
    Return the crystal resolutions
    """
    # from Pushpa
    
    # res = p[0]/sqrt(x) + p[1]
    hiEresol = [
        [0.7127, 0.004879],
	[0.7127, 0.004879],
	[0.7127, 0.004879],
	[0.7127, 0.004879],
	# tweaking resolution on C5
        #[0.7127, 0.004879],
        [1.7, 0.005],
	[0.7127, 0.004879],
	[0.7127, 0.004879],
        # tweaking resolution on C8
	#[0.7127, 0.004879]
        [2.4, 0.005]
    ]
    
    # C5 = C8 == C4
    # C6 = C7 == C1
    loEresol = [
        [0.2349,  0.018600],
	[0.2699,  0.014920],
	[0.2459,  0.015020],
	[0.3797, -0.003919],
	# tweaking resolution on C5
        #[0.3797, -0.003919],
        [0.8, -0.004],
	[0.2349,  0.018600],
	[0.2349,  0.018600],
        # tweaking resolution on C8
	#[0.3797, -0.003919]
        [1.3, -0.004]
    ]
    
    if E:
        p0, p1 = hiEresol[int(i)]
    else:
        p0, p1 = loEresol[int(i)]

    selection = '(('+str(p0)+'/sqrt(edep['+str(i)+']*1000.)) + '+str(p1)+')'
    return selection


def noiseCuts60(i, E):
    """
    Quick hack at implementing Pushpa's noise cuts
    """
    
    hiEcuts = [
        "!((crystal1.energyD > 1000.) && ((pmt11.rqtD1_5+pmt12.rqtD1_5)/2. < 2.660))",
        "!((crystal2.energyD > 1000.) && ((pmt21.rqtD1_5+pmt22.rqtD1_5)/2. < 2.640))",
        "!((crystal3.energyD > 1000.) && ((pmt31.rqtD1_5+pmt32.rqtD1_5)/2. < 2.660))",
        "!((crystal4.energyD > 1000.) && ((pmt41.rqtD1_5+pmt42.rqtD1_5)/2. < 2.680))",
        "!((crystal5.energyD > 1000.) && ((pmt51.rqtD1_5+pmt52.rqtD1_5)/2. < 2.650))",
        "!((crystal6.energyD > 1000.) && ((pmt61.rqtD1_5+pmt62.rqtD1_5)/2. < 2.655))",
        "!((crystal7.energyD > 1000.) && ((pmt71.rqtD1_5+pmt72.rqtD1_5)/2. < 2.630))",
        "!((crystal8.energyD > 1000.) && ((pmt81.rqtD1_5+pmt82.rqtD1_5)/2. < 2.660))"
    ]
    
    loEcal = [
        0.0001093478,
        0.000105676,
        0.0001159,
        0.0001127,
        2.869072e-04,
        0.000117534,
        0.0001119,
        3.842627e-04
    ]
    
    X = str(i+1)
    noise_cut = '(pmt'+X+'1.nc > 1 && pmt'+X+'2.nc > 1)'
    time_cut  = '(crystal'+X+'.t0 > 2.)'
    dama_cut  = '(crystal'+X+'.x2 / crystal'+X+'.x1 < 1.25)'
    asym_cut  = '(abs((pmt'+X+'1.qc5-pmt'+X+'2.qc5)/(pmt'+X+'1.qc5+pmt'+X+'2.qc5)) < 0.5)'
    qcnc_cut  = '(crystal'+X+'.qc5*'+str(loEcal[i])+' > ((1./1500000.) * (crystal'+X+'.qc/crystal'+X+'.nc) * (crystal'+X+'.qc/crystal'+X+'.nc)+(13./100000.) * (crystal'+X+'.qc/crystal'+X+'.nc)-0.1))'
    
    loEcut = '('+noise_cut+' && '+time_cut+' && '+dama_cut+' && '+asym_cut+' && '+qcnc_cut+')'
    
    if E:
        cut = hiEcuts[i]
    else:
        #cut = loEcuts[i]
        cut = loEcut

    return cut



#===============================================================================
#   STANDARD FUNCS
#===============================================================================

def amLocal():
    """
    Check to see what machine you're on so you can get right paths
    """
    import socket
    ### master.cunpa.ibs
    local = 1
    host = str(socket.gethostname())
    #print 'hostname =',host
    if 'cunpa' in host:
        local = 0
    return local


def sortKeys(data, bkgs, sigs):
    datkeys=[]
    bakkeys=[]
    sigkeys=[]

    for key in data:
        datkeys.append(key)
    datkeys.sort()

    for key in bkgs:
        bakkeys.append(key)
    bakkeys.sort()

    for key in sigs:
        sigkeys.append(key)
    sigkeys.sort()

    return datkeys, bakkeys, sigkeys


def histparam(energy=0):
    """
    Histogram parameters for number of bins, min, and max in keV
    """
    if energy:
        hmin = 0
        hmax = 3000
        bins = (hmax-hmin)
    else:
        hmin = 0
        hmax = 200
        bins = (hmax-hmin)
    return [bins, hmin, hmax]


def names(i):
    """
    The names of the crystals
    """
    crystals = [
        'NaI-01',
        'NaI-02',
        'NaI-03',
        'NaI-04',
        'NaI-05',
        'NaI-06',
        'NaI-07',
        'NaI-08'
    ]
    return crystals[int(i)]


def cnames(i):
    """
    The C names of the crystals
    """
    crystals = [
        'C1',
        'C2',
        'C3',
        'C4',
        'C5',
        'C6',
        'C7',
        'C8'
    ]
    return crystals[int(i)]


def cmass(i):
    """
    crystal masses in kg
    """
    mass = [
         8.26,
         9.15,
         9.16,
        18.01,
        18.28,
        12.50,
        12.50,
        18.28
        ]
    return mass[int(i)]


def longNames(i):
    """
    Full name and specs for the crystal
    """
    crystals = [
        cnames(i)+'  NaI-001  Sample-B  '+str(cmass(i))+'kg',
        cnames(i)+'  NaI-002  Sample-C  '+str(cmass(i))+'kg',
        cnames(i)+'  NaI-007  WimpScint-2  '+str(cmass(i))+'kg',
        cnames(i)+'  AS-3  WimpScint-2  '+str(cmass(i))+'kg',
        cnames(i)+'  AS-1  Sample-C  '+str(cmass(i))+'kg',
        cnames(i)+'  NaI-011  WimpScint-3  '+str(cmass(i))+'kg',
        cnames(i)+'  NaI-012  WimpScint-3  '+str(cmass(i))+'kg',
        cnames(i)+'  AS-2  Sample-C  '+str(cmass(i))+'kg'
    ]
    return crystals[int(i)]


def volumeNames(i):
    """
    For the simulation primary volume names
    """
    volumeNames = [
        'NaIDet01Crystal',
        'NaIDet02Crystal',
        'NaIDet03Crystal',
        'NaIDet04Crystal',
        'NaIDet05Crystal',
        'NaIDet06Crystal',
        'NaIDet07Crystal',
        'NaIDet08Crystal'
    ]
    return volumeNames[int(i)]


def rainbow(N):
    """
    Generate an array of colors for MC plotting
    """
    from random import randint
    colors=[]
    cis=[]
    # ci cannot be too large > 10,000!
    ci = randint(1000, 5000)
    for h in range(N):
        H = float(h)/float(N)
        ci += 1
        if H <= 1/5. :
            R=1.
            G=1.*5*H
            B=0.
        elif H > 1/5. and H <= 2/5. :
            R=1.-(1*5*(H-1/5.))
            G=1.
            B=0.
        elif H > 2/5. and H <= 3/5. :
            R=0.
            G=1.
            B=1.*5*(H-2/5.)
        elif H > 3/5. and H <= 4/5. :
            R=0.
            G=1.-(1*5*(H-3/5.))
            B=1.
        elif H > 4/5. and H <= 1. :
            R=1.*5*(H-4/5.)
            G=0.
            B=1.
        elif H > 1. :
            R=1.
            G=1.
            B=1.
        
        color = TColor(ci, R, G, B)
        # must keep the color and the ci in memory
        colors.append(color)
        cis.append(ci)
        
    return colors, cis


def readFile(fileName="backgrounds60.txt"):
    lines = []
    skip = 0
    with open(fileName) as mcfile:
        for line in mcfile:
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue
            if line.startswith('\"\"\"'):
                if skip == 0: skip = 1
                else: skip = 0
                continue
            if skip: continue
            lines.append(line)
    mcfile.close()
    return lines

