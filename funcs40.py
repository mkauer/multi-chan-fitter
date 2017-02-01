#!/usr/bin/env python
######################################################################
# Get single-hit data into the processing stream
#
# Works with v40 and later versions
#
# version: 2017-01-31
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + add dataDRU40() in funcs40.py
# + add build() in funcs40.py
# + add buildMC40() in funcs40.py
# + add buildData40() in funcs40.py
# + add getInfo() in funcs40.py
# + funcs4.py clean slate
# 
# email me: mkauer@physics.wisc.edu
######################################################################
# 
# Where is MC?
# /data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/
# 
# Where is the processed data?
# Right now I'm using:
# /data/COSINE/NTP/phys/V00-02-00
# and run ntp_I001546*
# 
######################################################################

import os,sys
import numpy as np
from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs32 import *


def getInfo(line):
    
    info={}

    info['line'] = line
    
    bits = line.split()
    #bits = [x.strip() for x in bits]
        
    # data type
    info['type'] = str(bits[0])

    # reuse a joined rootfile
    if 'R' in info['type']:
        info['reuse'] = 1
    else:
        info['reuse'] = 0
        
    # channel
    info['chan'] = str(bits[1])
    
    # crystal
    info['xstl'] = int(bits[2])
    
    # reuse rootfile
    if 'root' in str(bits[-1]):
        info['root'] = str(bits[-1])
    else:
        info['root'] = None


    if 'D' in info['type']:
        
        # run number
        info['run'] = str(bits[4])

        # processing version
        info['build'] = str(bits[5])
        
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

    return info


def buildData40(info, data):
    
    if info['reuse']:
        print 'not ready not'
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
            nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/data/phys/'+info['build']+'/ntp_I*'+info['run']+'*root*')
        else:
            nfiles = chain.Add('/data/COSINE/NTP/phys/'+info['build']+'/ntp_I*'+info['run']+'*root*')

        if nfiles > 0:
            print nfiles,'data files found'
            for E in range(2):

                ### use crystalX.energyD for high energy
                ### use crystalX.qc5 for low energy
                if E:
                    edep = 'energyD'
                else:
                    edep = 'qc5'

                par = histparam(E)
                #for i in range(8):
                i = info['xstl'] - 1
                
                key  = 'x'+str(i+1)
                key += '-data'
                key += '-c'+info['chan']
                key += '-e'+str(E)
                
                data[key] = {}

                #histo = TH1F(key, key, par[0], par[1], par[2])
                histo = TH1F(key, longNames(i), par[0], par[1], par[2])
                #chain.Draw('(pmt'+str(i+1)+'1.'+str(edep)+'*'+str(calib(i,0,E))
                #           +' + pmt'+str(i+1)+'2.'+str(edep)+'*'+str(calib(i,1,E))+')/2.'
                #           +' >> '+str(key))
                chain.Draw('crystal'+str(i+1)+'.'+str(edep)+'*'+str(calib22(i,E))
                           +' >> '+str(key))

                histo.Sumw2()
                histo.SetLineColor(kBlack)
                histo.SetMarkerColor(kBlack)
                histo.SetLineWidth(1)

                data[key]['hist'] = histo

                key2 = key+'_subruns'
                temp = TH1F(key2,'subruns',1,0,1)
                #temp.SetBinContent(0, nfiles)
                for i in range(nfiles):
                    temp.Fill(0.5)
                data[key]['subruns'] = temp

                data[key]['runtime'] = nfiles*subrunTime*60.*60.
        
        else:
            print 'ERROR: No data files found... quitting...'
            sys.exit()
        
    return data


def buildMC40(info, mc, calib=2):

    if info['reuse']:
        print 'not ready not'
        return mc
    
    else:
        local = amLocal()

        chain = TChain("MC","")
        nfiles=0
        if local:
            nfiles = chain.Add('/home/mkauer/COSINE/CUP/mc-fitting/sim/newGeometry/'+info['isof']+'/set2/'+'*'+info['loca']+'*root')
        else:
            nfiles = chain.Add('/data/MC/KIMS-NaI/user-scratch/sim/processed/newGeometry/'+info['isof']+'/set2/'+'*'+info['loca']+'*root')

        if nfiles > 0:
            print nfiles,'MC files found for', info['loca'], info['isof']
            for E in range(2):
                
                i = info['xstl'] - 1
                
                ### see if i can use unique names for all the histograms
                key  = 'x'+str(i+1)
                key += '-'+info['loca']
                if info['chst'] == info['chsp']: key += '-'+info['chst']
                else: key += '-'+info['chst']+'_'+info['chsp']
                key += '-c'+info['chan']
                key += '-e'+str(E)
                
                mc[key] = {}
                print key

                par = histparam(E)
                histo = TH1F(key, longNames(i), par[0], par[1], par[2])

                cut1 = TCut('edep['+str(i)+']*1000. > 0')
                
                if info['loca'] == 'internal':
                    cut2 = TCut('primVolumeName == "'+volumeNames(i)+'"')
                elif info['loca'] == 'pmt':
                    cut2 = TCut('primVolumeName == "phys_pmt"')
                elif info['loca'] == 'lsveto':
                    cut2 = TCut('primVolumeName == "lsveto"')
                elif info['loca'] == 'lsvetoair':
                    cut2 = TCut('(primVolumeName == "DetPMTCover") || (primVolumeName == "DetPMTEnvelope")')
                elif info['loca'] == 'airshield':
                    cut2 = TCut('primVolumeName == "LSVetoAirRoom"')
                elif info['loca'] == 'steel':
                    # not sure this works the way I think it should?
                    cut2 = TCut('primVolumeName != ""')
                else:
                    print "WARNING: No selection criteria for  --> ", info['loca']
                    continue
                
                ### this is needed to get the generated event numbers right!
                cut3 = TCut('primParticleName == "'+info['chst']+'"')
                
                ### this is needed? to cut out crap events
                ### (event_info.Type) is new
                ### (evt_Type) is old
                #cut4 = TCut('(event_info.Type > 10) || (evt_Type > 10)')
                cut4 = TCut('event_info.Type > 10')
                
                cut100 = 0
                if info['chst'] == 'U238' and info['chsp'] == 'Rn222':
                    cut100 = TCut('(groupNo >= 11) && (groupNo <= 14)')
                if info['chst'] == 'Pb210' and info['chsp'] == 'Pb210':
                    cut100 = TCut('groupNo == 15')
                
                key2 = key+'_generated'
                temp = TH1F(key2, 'generated',1,0,1)
                
                ### v31 - maybe need cut100 for generated events
                ### to get the eff/normalization right?
                #chain.Draw('primVolumeName >> '+key2, cut2)
                if cut100:
                    #chain.Draw('primVolumeName >> '+key2, cut2+cut3+cut100)
                    chain.Draw('primVolumeName >> '+key2, cut2+cut3)
                else:
                    chain.Draw('primVolumeName >> '+key2, cut2+cut3)
                
                generated = temp.GetEntries()
                
                ### test out resolution smearing
                ###-------------------------------------------------------------------------------
                ### using Box-Muller? method here for "rng" alias
                chain.SetAlias('rng','sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')

                ### set the resolution function
                if calib == 1:
                    ### from Estella
                    ### assume reso = p0/sqrt(energy)
                    p0 = resol(i,E)
                    chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.)')
                if calib == 2:
                    ### from Pushpa
                    ### assume reso = p0/sqrt(energy) + p1
                    p0, p1 = resol2(i,E)
                    chain.SetAlias('sigma', str(p0)+'/sqrt(edep['+str(i)+']*1000.) + '+str(p1))

                if cut100:
                    chain.Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> '+key, cut1+cut2+cut4+cut100)
                else:
                    chain.Draw('(edep['+str(i)+']*1000.) + (sigma*edep['+str(i)+']*1000.*rng) >> '+key, cut1+cut2+cut4)

                detected = histo.GetEntries()
                ###-------------------------------------------------------------------------------
                
                histo.SetLineColor(kBlack)
                histo.SetMarkerColor(kBlack)
                histo.SetLineWidth(1)

                mc[key]['hist'] = histo
                mc[key]['generated'] = temp
                mc[key]['acti'] = info['acti']
                mc[key]['erro'] = info['erro']
                mc[key]['fbnd'] = info['fbnd']
                
                #print 'Eff =',detected/generated
                """
                if 'B' in bsdr:
                    bkgs[key] = {}
                    bkgs[key]['hist'] = histo
                    bkgs[key]['generated'] = temp
                    bkgs[key]['acti'] = acti
                    bkgs[key]['erro'] = erro

                if 'S' in bsdr:
                    sigs[key] = {}
                    sigs[key]['hist'] = histo
                    sigs[key]['generated'] = temp
                    sigs[key]['fbnd'] = fbnd
                """
                
        else:
            #print 'WARNING:', loca, isof, 'not found...'
            print 'ERROR: No MC files found for',info['loca'], info['isof'],'... quitting...'
            sys.exit()
            
    return mc

    
def build(infile = 'backgrounds40.txt'):

    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo(line)
        if 'D' in info['type']:
            data = buildData40(info, data)
        elif 'B' in info['type']:
            bkgs = buildMC40(info, bkgs)
        elif 'S' in info['type']:
            sigs = buildMC40(info, sigs)
        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs


def dataDRU40(data):
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


