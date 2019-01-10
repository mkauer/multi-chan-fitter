#!/usr/bin/env python
######################################################################
# funcs-misc.py
# 
# Get all the default functions in here
# 
# version: 2019-01-09
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + add numX()
# + add I129 to setGroup()
# + add plastic to setGroup()
# + addHistKey() to test if a key already exists and create one if not
# ~ tweaked extendHist() because it had a couple bugs
# + fitPrep() to rebin and resize the histograms
# + add rootSmoothing() - uses ROOTs TH1.Smooth() method
# + add averageSmoothing()
# ~ don't overwrite bin[i] with smoothed value until all done
# ~ renamed smoothXbins() to triangleSmoothing()
# + add an extend function to easily extend histograms
# ~ updated the naming of stuff for the lsveto (crystal 9)
# ~ tweaked readFile() so the # search is after the """ search
# ~ tweaked histparam() to also return bins per keV
# ~ fixed sim path+name so not to combine lsveto and lsvetoair
# ~ tweaked longNames() in funcs60.py
# + detailNames() to funcs60.py
# + chanNames() to funcs60.py
# + energyNames() to funcs60.py
# ~ rename everything to 60 version and add other funcs as needed
# + try to move all main funcs to here funcs60.py - it's becoming too
#   confusing to find all funcs in all various versions of funcs
# + building off funcs52.py
# 
# email me: mkauer@physics.wisc.edu
######################################################################

import os,sys
import socket
from copy import deepcopy
import numpy as np
from ROOT import *
import ROOT


def numX():
    # number of crystals (including lsveto)
    return 9


def amLocal():
    if 'cunpa' in str(socket.gethostname()):
        return 0
    else:
        return 1


def onCup():
    if 'cunpa' in str(socket.gethostname()):
        return 1
    else:
        return 0


def baseDir():
    if 'cunpa' in socket.gethostname():
        return '/home/mkauer/mc-fitting/'
    else:
        return '/home/mkauer/COSINE/CUP/mc-fitting/'


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
        'NaI-08',
        'LSveto'
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
        'C8',
        'LS'
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
        18.28,
        ### average crystal mass if grouped into lsveto
        106.14/8. ]
    return mass[int(i)]


def surfArea(i):
    """
    crystal surface area in m^2
    """
    diameter = [ 5.0,
                 4.2,
                 4.2,
                 5.0,
                 5.0,
                 4.8,
                 4.8,
                 5.0 ]
    length = [  7.0,
               11.0,
               11.0,
               15.3,
               15.5,
               11.8,
               11.8,
               15.5 ]
    if i == 8:
        ### average surface area if grouped into lesveto
        return 3.1416 * (sum(np.asarray(diameter)*np.asarray(length))/8.) * 0.00064516
    else:
        return 3.1416 * (diameter[i]*length[i]) * 0.00064516


def energyNames(E=0):
    name = ''
    if E == 0: name = 'lo-energy'
    if E == 1: name = 'hi-energy'
    return name

    
def chanNames(chan='S'):
    name = ''
    if chan == 'S': name = 'single-hit'
    if chan == 'M': name = 'multi-hit'
    return name
    

def detailNames(i):
    crystals = [
        'NaI-001  Sample-B',
        'NaI-002  Sample-C',
        'NaI-007  WimpScint-2',
        'AS-3  WimpScint-2',
        'AS-1  Sample-C',
        'NaI-011  WimpScint-3',
        'NaI-012  WimpScint-3',
        'AS-2  Sample-C',
        'LS Veto'
    ]
    return crystals[int(i)]

    
def longNames(i):
    """
    Full name and specs for the crystal
    """
    return str(cnames(i)+'  '+detailNames(i)+'  '+str(cmass(i))+'kg')


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
        'NaIDet08Crystal',
        'lsveto'
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


def readFile(fileName):
    lines = []
    skip = 0
    with open(fileName) as mcfile:
        for line in mcfile:
            line = line.strip()
            if not line: continue
            #if line.startswith('#'): continue
            if line.startswith('\"\"\"'):
                if skip == 0: skip = 1
                else: skip = 0
                continue
            if skip: continue
            if line.startswith('#'): continue
            lines.append(line)
    mcfile.close()
    return lines


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
        #print key, keVperBin
        scale = float(1./(days*kgs*keVperBin))
        print key, scale
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
    duration = -1
    try:
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
    except:
        return -1


def setBinError(histo):
    for n in range(histo.GetNbinsX()+1):
        histo.SetBinError(n, sqrt(histo.GetBinContent(n)))
    return histo


def zeroBinError(histo):
    for n in range(histo.GetNbinsX()+1):
        histo.SetBinError(n, 0.0)
    return histo


def uniqString(string):
    newstring=''
    for char in string:
        if char not in newstring:
            newstring += char
    return newstring


def setGroup(info):
    if info['loca'] == 'internal':
        if info['isof'] in ['Cd109',  'Sn113',
                            'H3',     'Na22',
                            'I125',   'I126',   'I129',
                            'Te121',  'Te121m',
                            'Te123m', 'Te125m', 'Te127m']:
            return 'cosmo'
        else: return 'internal'
    if 'pmt'     in info['loca']: return 'pmts'
    if 'surf'    in info['loca']: return 'surface'
    if 'cucase'  in info['loca']: return 'copper'
    if 'copper'  in info['loca']: return 'copper'
    if 'steel'   in info['loca']: return 'steel'
    if 'teflon'  in info['loca']: return 'surface'
    if 'veto'    in info['loca']: return 'lsveto'
    if 'plastic' in info['loca']: return 'plastic'
    return 'none'


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


def triangleSmoothing(hists, s=0):

    if s <= 0:
        print 'WARNING: not smoothing!'
        return hists

    print 'INFO: triangle smoothing...'
    for key in hists:

        ### get the number of bins
        bins = hists[key]['hist'].GetNbinsX()

        ### create a different hist to save the smoothed data to
        newhist = deepcopy(hists[key]['hist'])
        
        for i in range(s, bins-s):
            ### "i" is the bin you are smoothing
            ### "s" is the +/- bin range you are smoothing over
            smoothed = 0

            ### grab previous (s) bins
            for j, k in enumerate(range(-s, 0)):
                smoothed += hists[key]['hist'].GetBinContent(i+k)*(j+1)

            ### grab current bin
            smoothed += hists[key]['hist'].GetBinContent(i)*(s+1)

            ### grab following (s) bins
            for j, k in enumerate(range(1, s+1)):
                smoothed += hists[key]['hist'].GetBinContent(i+k)*(s-j)

            ### normalize
            newhist.SetBinContent(i, smoothed/float((s+1.)**2.))

        ### save the new hist back into the dictionary
        hists[key]['hist'] = deepcopy(newhist)
        del newhist
    
    return hists


def averageSmoothing(hists, s=0):

    if s <= 0:
        print 'WARNING: not smoothing!'
        return hists
    
    print 'INFO: average smoothing...'
    for key in hists:

        ### get the number of bins
        bins = hists[key]['hist'].GetNbinsX()

        ### create a different hist to save the smoothed data to
        newhist = deepcopy(hists[key]['hist'])
        
        for i in range(s, bins-s):
            ### "i" is the bin you are smoothing
            ### "s" is the +/- bin range you are smoothing over
            smoothed = 0

            ### grab previous (s) bins
            for j, k in enumerate(range(-s, 0)):
                smoothed += hists[key]['hist'].GetBinContent(i+k)

            ### grab current bin
            smoothed += hists[key]['hist'].GetBinContent(i)

            ### grab following (s) bins
            for j, k in enumerate(range(1, s+1)):
                smoothed += hists[key]['hist'].GetBinContent(i+k)

            ### normalize
            newhist.SetBinContent(i, smoothed/float(2.*s+1.))

        ### save the new hist back into the dictionary
        hists[key]['hist'] = deepcopy(newhist)
        del newhist
    
    return hists


def rootSmoothing(hists, s=0):
    
    if s <= 0:
        print 'WARNING: not smoothing!'
        return hists
    
    print 'INFO: root TH1 smoothing...'
    for key in hists:
        ### smooth the hist s times
        hists[key]['hist'].Smooth(s)
        
    return hists


def smooth(hists, s=0):
    
    hists = averageSmoothing(hists, s)
    #hists = rootSmoothing(hists, s)
    #hists = triangleSmoothing(hists, s)
    
    return hists


def extendHist(hist1, hist2):
    
    name  = hist1.GetName()
    bins1 = hist1.GetNbinsX()
    bins2 = hist2.GetNbinsX()
    bins  = bins1 + bins2
    hist3 = TH1F(name, name, bins, 0, bins)
    for n in range(bins1):
        hist3.SetBinContent(n, hist1.GetBinContent(n))
    for n in range(bins2):
        hist3.SetBinContent(n+bins1, hist2.GetBinContent(n))
        
    return hist3


def fitPrep(hists, key, params, zeros=0):
    
    # params needs to be a list of length 3 [rebin, kevfmin, kevfmax]
    rebin   = params[0]
    kevfmin = params[1]
    kevfmax = params[2]

    # if you rebin now, you need to deepcopy the original hist
    hist = deepcopy(hists[key]['hist'])
    if rebin > 1: hist.Rebin(rebin)
    
    bins   = hist.GetXaxis().GetNbins()
    kevmin = hist.GetXaxis().GetBinUpEdge(0)
    kevmax = hist.GetXaxis().GetBinUpEdge(bins)

    # keV to bin number scaling
    ktb      = bins/(kevmax-kevmin)
    startbin = int(kevfmin*ktb)
    stopbin  = int(kevfmax*ktb)
    newbins  = int(stopbin-startbin)
    
    temp = TH1F(key, key, newbins, 0, newbins)
    if zeros:
        del hist
        return temp
    
    for i in range(newbins):
        temp.SetBinContent(i, hist.GetBinContent(startbin+i))

    del hist
    return temp


def addHistKey(hists, key):
    
    try:
        tmp = hists[key]
        #del tmp
    except:
        hist = TH1F(key, key, 1, 0, 1)
        hists[key] = {}
        hists[key]['hist'] = hist
        
    return hists

