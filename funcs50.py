#!/usr/bin/env python
######################################################################
# funcs50.py
# 
# Get single-hit data into the processing stream
# 
# Works with v50 and later versions
# 
# version: 2017-03-09
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + iterate through the chars in the chans string
# + new buildMC50()
# + new buildData50()
# + new build50()
# + new getInfo50()
# + building off funcs40.py
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
from funcs40 import *


def getInfo50(line, freuse=0, fchans=0):
    
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


def build50(infile = 'backgrounds41.txt', freuse=0, fchans=0):

    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo50(line, freuse, fchans)
        if 'D' in info['type']:
            data = buildData50(info, data)
        elif 'B' in info['type']:
            bkgs = buildMC50(info, bkgs)
        elif 'S' in info['type']:
            sigs = buildMC50(info, sigs)
        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs


def buildData50(info, data):

    #for i in range(len(info['chans'])):
    #    info['chan'] = info['chans'][i]
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

            return data

        else:
            local = amLocal()

            ### will need to do this dynamically eventually
            ### run 1324 has 1 hour subruns
            ### run 1544 has 2 hour subruns
            ### run 1546 has 2 hour subruns
            subrunTime = float(2.) # in hours

            chain = TChain("ntp","")
            nfiles=0
            if local:
                #nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/data/phys/'+info['build']+'/ntp_I*'+info['run']+'*root*')
                nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/data/phys/'+info['build']+'_MERGED/ntp_I*'+info['run']+'*root*')

            else:
                #nfiles = chain.Add('/data/COSINE/NTP/phys/'+info['build']+'/ntp_I*'+info['run']+'*root*')
                #nfiles = chain.Add('/data/COSINE/NTP/phys/'+info['build']+'/ntp_I*'+info['run']+'*root.000')
                nfiles = chain.Add('/data/COSINE/NTP/phys/'+info['build']+'_MERGED/ntp_I*'+info['run']+'*root*')

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
                    edep, selection = calib40(i,e)
                    ### new calib
                    #edep, selection = calib41(i,e)


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

                        #chanCut = TCut('(('+edepcuts+'))')
                        #chanCut = TCut('(('+nclustercuts+'))')
                        #chanCut = TCut('(('+edepcuts+') && ('+nclustercuts+'))')

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

                        #chanCut = TCut('(('+edepcuts+'))')
                        #chanCut = TCut('(('+nclustercuts+'))')
                        #chanCut = TCut('(('+edepcuts+') || ('+nclustercuts+'))')

                        ### my cuts
                        #chanCut = TCut('(('+edepcuts+') || ('+nclustercuts+') || ('+lsvetocut+'))')
                        ### Pushpa cuts
                        chanCut = TCut('(('+nclustercuts+') || ('+lsvetocut+'))')

                    else:
                        print 'ERROR: I do not know what to do with channel -->',info['chan']
                        print 'Available channels are [A]All-hits, [S]Single-hits, [M]Multi-hits'
                        sys.exit()


                    ### Pushpa's noise cut
                    noiseCut = TCut('('+noiseCuts40(i,e)+')')


                    ### combine all cuts
                    masterCut = TCut('('+
                                     chanCut.GetTitle()+' && '+
                                     noiseCut.GetTitle()
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


def buildMC50(info, mc):

    #for i in range(len(info['chans'])):
    #    info['chan'] = info['chans'][i]
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

                #except:
                #    del mc[key]
                #    print "WARNING: could not find hist -->",key
                #    continue

            return mc

        else:
            local = amLocal()

            chain = TChain("MC","")
            nfiles=0
            if local:
                nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/sim/newGeometry/'
                                   +info['isof']+'/set2/'+'*'+info['loca']+info['isof']+'*root')

            else:
                nfiles = chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'
                                   +info['isof']+'/set2/'+'*'+info['loca']+info['isof']+'*root')

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
                    if info['loca'] == 'internal':
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
                    """
                    ### set the resolution function
                    if resol == 1:
                        ### from Estella
                        ### assume reso = p0/sqrt(energy)
                        p0 = resol(i,e)
                        chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.)')
                    if resol == 2:
                        ### from Pushpa
                        ### assume reso = p0/sqrt(energy) + p1
                        p0, p1 = resol2(i,e)
                        chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.) + '+str(p1))
                    """
                    resolFunc = resol40(i,e)
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

