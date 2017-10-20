#!/usr/bin/env python
######################################################################
# funcs80.py
# 
# New calibrations for all crystals
# 
# Works with v80 and later versions
# 
# version: 2017-10-18
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------

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


def build80(infile='backgrounds800.txt', freuse=0, fchans=0):

    data = {}
    bkgs = {}
    sigs = {}
    runtime = 0
    
    for line in readFile(infile):
        info = getInfo64(line, freuse, fchans)
        
        if 'D' in info['type']:
            data, runtime = buildData80(info, data)
            
        elif 'B' in info['type']:
            bkgs = buildMC71(info, bkgs)
            
        elif 'F' in info['type']:
            sigs = buildMC71(info, sigs)
            
        else:
            print 'WARNING: I do not know what to do with type',info['type'], 'in line:'
            print info['line']
            continue
        
    return data, bkgs, sigs, runtime


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
                    edep, selection = calib80(i,e)
                    #-----------------------------------------------------------------------
                    
                    
                    # DEFINE CUTS
                    #-----------------------------------------------------------------------
                    #if (i+1 == 7) and (e == 0):
                    #    masterCut = cutsBDT(i,c,e)
                    #else:
                    #    masterCut = cutsOld(i,c,e)
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

    ### c7 is special
    c7 = [566., 9255.]

    if E:
        edep = '(crystal'+str(i+1)+'.energyD)'
        selection = '('+edep+')'
    else:
        edep = '(crystal'+str(i+1)+'.qc5)'
        selection = '('+edep+'*'+str(loE[i])+')'
        if str(i) == str(6):
            selection = '(('+edep+'-'+str(c7[0])+')/'+str(c7[1])+')'

    return edep, selection


def resol80(i, E=0):
        
    # res = p[0]/sqrt(x) + p[1]
    hiEresol = [
        [0.6729, 0.009374],
	[0.6531, 0.006627],
	[0.5926, 0.009506],
	[0.7227, 0.004790],
        [1.7, 0.0], # tweaking C5
	[0.6498, 0.009670],
	[0.7034, 0.007812],
        [2.4, 0.0] # tweaking C8
    ]
    
    loEresol = [
        [0.2413,  0.01799],
	[0.2951,  0.01427],
	[0.3106,  0.007894],
	[0.3894, -0.001437],
        [0.8, 0.0], # tweaking C5
	[0.3620,  0.0006355],
	[0.3042,  0.009784],
        [1.3, 0.0] # tweaking C8
    ]
    
    if E:
        p0, p1 = hiEresol[int(i)]
    else:
        p0, p1 = loEresol[int(i)]

    selection = '(('+str(p0)+'/sqrt(edep['+str(i)+']*1000.)) + '+str(p1)+')'
    return selection


def cutsBDT80(i,c,e):
    """
    From: https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    """

    ### global noise cuts
    coinc  = '(BLSVeto_isCoincident == 1)'
    muons  = '(BMuon_totalDeltaT0/1.e6 > 30)'
    bdt    = '(crystal'+str(i+1)+'.bdt > -0.03)'
    charge = '(crystal'+str(i+1)+'.rqcn > -1)'
    nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
    
    noiseCut = coinc+' && '+muons+' && '+bdt+' && '+charge+' && '+nc
    
    if c == 'A':
        return TCut(noiseCut)
        
    elif c == 'S':
        lsveto = '(BLSVeto_Charge/143.8 < 20)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.'+'nc'+' < 4) && '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'
        
        masterCut = TCut(noiseCut+' && '+'('+lsveto+' && '+hits+')')
        return masterCut
    
    elif c == 'M':
        lsveto = '(BLSVeto_Charge/143.8 > 20)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.'+'nc'+' > 4) || '
        ### remove extra '&&' or '||'
        hits = hits[:-4]+')'
        
        masterCut = TCut(noiseCut+' && '+'('+lsveto+' || '+hits+')')
        return masterCut
    
    else:
        print 'ERROR: I do not know what to do with channel -->',c
        print 'Available channels are [S]Single-hits, [M]Multi-hits'
        sys.exit()

