#!/usr/bin/env python
######################################################################
# funcs_misc.py
# 
# Put all the general use functions in here
# 
# version: 2021-04-02
# 
# email me: mkauer@physics.wisc.edu
######################################################################

import os, sys
import socket
from copy import deepcopy
import numpy as np
from ROOT import *
import ROOT


def numX():
    # number of crystals (including lsveto)
    return 9


def amLocal():
    host = socket.gethostname()
    if 'cunpa' in host or 'node' in host or 'ibs' in host:
        return 0
    else:
        return 1


def onCup():
    host = socket.gethostname()
    if 'cunpa' in host or 'node' in host or 'ibs' in host:
        return 1
    else:
        return 0


def baseDir():
    """
    host = socket.gethostname()
    if 'cunpa' in host or 'node' in host or 'ibs' in host:
        return '/home/mkauer/mc-fitting/'
    else:
        return '/home/mkauer/COSINE/CUP/mc-fitting/'
    """
    return os.path.dirname(os.path.abspath(__file__))


def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return


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
        12.516,
        12.516,
        18.28,
        1800.
        ]
    return mass[int(i)]


def surfArea(i):
    """
    crystal surface area in m^2
    """
    # these don't look right
    # us mcDimensions() from below
    """
    diameter = [
        5.0,
        4.2,
        4.2,
        5.0,
        5.0,
        4.8,
        4.8,
        5.0
    ]
    length = [
        7.0,
        11.0,
        11.0,
        15.3,
        15.5,
        11.8,
        11.8,
        15.5
    ]
    if i == 8:
        ### average surface area if grouped into lesveto
        return 3.1416 * (sum(np.asarray(diameter)*np.asarray(length))/8.) * 0.00064516
    else:
        return 3.1416 * (diameter[i]*length[i]) * 0.00064516
    """

    ### do this a new way
    ### return area in cm^2
    pi = 3.1416
    
    if i < 8:
        x, y, z, rad, height, dep = mcDimensions(i)
        R = rad/10. # convert to cm
        H = (height*2.)/10. # full height and in cm
        return (2.*pi*R*H) + (2.*pi*(R**2))
    
    elif i == 8:
        ### ls cube area (x=120, y=130, z=130 cm)
        #return 96200.
        
        ### from ls volume estimate (81797)
        #return 4.*pi*(80.7**2)
        
        ### from mean radius of MC events (52739)
        return 4.*pi*(64.8**2)
        
    else:
        print('ERROR: no surface area calc for crystal {0}'.format(i))
        sys.exit()
        
    
def mcDimensions(i):

    ### needed for NaI and Teflon surface depth profiling
    ### dimensions in mm
    ### 2020-05-12 - updated "X" with Gyunho's new geometry
    
    if i==0: # C1
        #x = -207.4
        x = -217.4
        y = -20.7
        z = 0.
        rad = 63.500
        height = 177.800/2.0
        dep = 0.01

    elif i==1: # C2
        #x = -62.2
        x = -67.2
        y = -26.1
        z = 0.
        rad = 53.340
        height = 279.400/2.0
        dep = 0.01
        
    elif i==2: # C3
        #x = 62.2
        x = 67.2
        y = -30.9
        z = 0.
        rad = 53.340
        height = 279.400/2.0
        dep = 0.01
        
    elif i==3: # C4
        #x = 191.7
        x = 201.7
        y = -38.6
        z = 0.
        rad = 63.500
        #height = 378.350/2.0  ## this one is wrong
        height = 387.350/2.0
        dep = 0.01
    
    elif i==4: # C5
        #x = -201.7
        x = -211.7
        y = -283.4
        z = 0.
        rad = 63.500
        height = 393.700/2.0
        dep = 0.01
        
    elif i==5: # C6
        #x = -67.8
        x = -72.8
        y = -281.7
        z = 0.
        rad = 60.325
        height = 298.450/2.0
        dep = 0.01

    elif i==6: # C7
        #x = 67.8
        x = 72.8
        y = -281.7
        z = 0.
        rad = 60.325
        height = 298.450/2.0
        dep = 0.01
    
    elif i==7: # C8
        #x = 201.8
        x = 211.8
        y = -281.7
        z = 0.
        rad = 63.500
        height = 393.700/2.0
        dep = 0.01

    else:
        print('ERROR: no crystal dimensions for crystal {0}'.format(i))
        sys.exit()
        
    return x, y, z, rad, height, dep


def surfProfile(x, pars):
    # The exp scaling doesn't really matter because
    # num "generated" events are also weighted so
    # that the efficiency is always correct.
    # Double checked 2022-11-10
    #   - no diff in spectral amplitude
    #return np.exp(-x[0]/pars[0])
    return (1/pars[0])*np.exp(-x[0]/pars[0])


def surfProfile88(x):
    p0 = 0.88
    return (1/p0)*np.exp(-x[0]/p0)


def surfProfile13(x):
    p0 = 0.13
    return (1/p0)*np.exp(-x[0]/p0)


def energyNames(E=0):
    name = ''
    if E == 0: name = 'lo-energy'
    if E == 1: name = 'hi-energy'
    if E == 2: name = 'alphas'
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


def rainbowSix(keynames):
    """
    Generate an array of colors for MC plotting
    Tweaked for ROOT6
    """
    from random import randint
    
    colors = {}
    cis = {}
    ci = randint(2000, 7000)
    #ci = 1666
    for h, key in enumerate(keynames):
        h = h+1
        ci += 1
        H = float(h)/float(len(keynames)+1)
        if H < 1/5. :
            R=1.
            G=1.*5*H
            B=0.
        elif H >= 1/5. and H < 2/5. :
            R=1.-(1*5*(H-1/5.))
            G=1.
            B=0.
        elif H >= 2/5. and H < 3/5. :
            R=0.
            G=1.
            B=1.*5*(H-2/5.)
        elif H >= 3/5. and H < 4/5. :
            R=0.
            G=1.-(1*5*(H-3/5.))
            B=1.
        elif H >= 4/5. and H < 1. :
            R=1.*5*(H-4/5.)
            G=0.
            B=1.
        elif H >= 1. :
            R=0.
            G=0.
            B=0.

        #print('{0} [{1}, {2}, {3}] {4}'.format(H, R, G, B, key))
        # doesn't work with root < 6.07
        #cis[key] = TColor.GetFreeColorIndex()
        cis[key] = ci
        colors[key] = TColor(cis[key], R, G, B, key)
        
    return colors, cis


def alphaColors(N):
    ac = [
        kBlue    +1,
        kCyan    +1,
        kGreen   +1,
        kOrange  +1,
        kMagenta +1,
        kBlue    -7,
        kCyan    -7,
        kGreen   -7,
        kOrange  -7,
        kMagenta -7,
        kBlue    +3,
        kCyan    +3,
        kGreen   +3,
        kOrange  +3,
        kMagenta +3,
    ]
    n = N % len(ac)
    return ac[n]


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
        print(key, scale)
        data[key]['hist'].Scale(scale)
        data[key]['druScale'] = scale
    return data


def makeTotal64(chan, E, par):
    total = []
    #par = histparams(E)
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
    #par = histparams(E)
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
    pars = [bins, hmin, hmax, bpkv]
    return pars


def globalParams():
    params = []
    for E in [0, 1, 2]:
        params.append(histparams(E))
    return params


def getDuration(rootfile):
    # get run duration in seconds
    duration = 0
    try:
        chain = TChain('ntp','')
        chain.Add(rootfile)
        entries = chain.GetEntries()
        ### trgtime has nano-sec precission
        var = "trgtime"
        norm = 1.e9
        chain.GetEntry(0)
        start = float(chain.GetLeaf(var).GetValue())
        ### just checking - trgtime is NOT indexed by crystal
        #for j in range(9):
        #    print('!!! {0}'.format(chain.GetLeaf(var).GetValue(j)))
        chain.GetEntry(entries-1)
        stop = float(chain.GetLeaf(var).GetValue())
        duration = (stop-start)/norm
        #print duration
        return duration
    except:
        print('WARNING: could not get runtime for --> \"{0}\"'
              .format(os.path.basename(rootfile)))
        return 0


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
                            'I128',   'Na24',
                            'I125',   'I126',   'I129',
                            'Te121',  'Te121m',
                            'Te123m', 'Te125m', 'Te127m']:
            return 'cosmo'
        else: return 'internal'
    elif 'pmt'     in info['loca']: return 'pmts'
    elif 'surf'    in info['loca']: return 'surface'
    elif 'cucase'  in info['loca']: return 'copper'
    elif 'copper'  in info['loca']: return 'copper'
    elif 'steel'   in info['loca']: return 'steel'
    elif 'teflon'  in info['loca']: return 'surface'
    elif 'reflec'  in info['loca']: return 'surface'
    elif 'veto'    in info['loca']: return 'lsveto'
    elif 'plastic' in info['loca']: return 'plastic'
    elif 'gamma'   in info['loca']: return 'ext-gamma'
    else:
        return 'none'


def setGroup500(info):
    if info['loca'] == 'internal':
        if info['isof'] in ['Cd109',  'Sn113',
                            'H3',     'Na22',
                            'I128',   'Na24',
                            'I125',   'I126',   'I129',
                            'Te121',  'Te121m',
                            'Te123m', 'Te125m', 'Te127m']:
            return 'cosmo'
        else: return 'internal'
    elif info['loca'] == 'pmt':      return 'pmts'
    elif info['loca'] == 'xpmt':     return 'pmts'
    elif info['loca'] == 'pmtbase':  return 'pmts'
    elif info['loca'] == 'xpmtbase': return 'pmts'
    elif 'cucase'   in info['loca']: return 'cucase'
    elif 'copper'   in info['loca']: return 'cucase'
    elif 'surf'     in info['loca']: return 'surface'
    elif 'teflon'   in info['loca']: return 'surface'
    elif 'reflec'   in info['loca']: return 'surface'
    elif 'lsveto'   in info['loca']: return 'lsveto'
    elif 'vetopmt'  in info['loca']: return 'vetopmt'
    elif 'film'     in info['loca']: return 'film'
    elif 'plastic'  in info['loca']: return 'plastic'
    elif 'cushield' in info['loca']: return 'cushield'
    elif 'steel'    in info['loca']: return 'steel'
    #elif 'gamma'    in info['loca']: return 'gamma'
    elif 'gamma'    in info['loca']: return 'cushield'
    elif 'neutron'  in info['loca']: return 'neutron'
    else:
        print('WARNING: no group defined for {0}'.format(info['loca']))
        return 'none'


def groupNum(info):
    
    # U238 group numbers
    # ------------------------
    # 11: U238  -> Th230
    # 12: Th230 -> Ra226
    # 13: Ra226 -> Rn222
    # 14: Rn222 -> Pb210
    # 15: Pb210 -> ground
    
    # Th232 group numbers
    # ------------------------
    # 21: Th232 -> Ra228
    # 22: Ra228 -> Th228
    # 23: Th228 -> ground
    
    # U235 group numbers
    # ------------------------
    #  0: U235 -> ground? (doesn't look like it gets splip up?)
    # 41: U235 -> Pa231
    # 42: Pa231 -> ground
    
    # others
    # ------------------------
    # 31: K40

    
    if info['chst'] in ['Te121',  'Te121m',
                        'Te123m', 'Te125m', 'Te127m',
                        'I125',   'I126',   'I129',
                        'Na22',   'H3',
                        'I128',   'Na24',
                        'Cd109',  'Sn113',
                        'Co60']:
        return TCut('(groupNo == 0)')
    
    if info['chst'] == 'K40':
        return TCut('(groupNo == 31)')
    
    if   info['chst'] == 'U238':  start = 11
    elif info['chst'] == 'Th230': start = 12
    elif info['chst'] == 'Ra226': start = 13
    elif info['chst'] == 'Rn222': start = 14
    elif info['chst'] == 'Pb210': start = 15
    elif info['chst'] == 'Th232': start = 21
    elif info['chst'] == 'Ra228': start = 22
    elif info['chst'] == 'Th228': start = 23
    elif info['chst'] == 'U235':  start = 41
    elif info['chst'] == 'Pa231': start = 42
    ### special case for external Tl208 gammas
    elif info['chst'] == 'Tl208': start = 0
    else: start = -1
    
    #if   info['chsp'] == 'U238':  stop = 11 # should not be a stop group
    if   info['chsp'] == 'Th230': stop = 12
    elif info['chsp'] == 'Ra226': stop = 13
    elif info['chsp'] == 'Rn222': stop = 14
    elif info['chsp'] == 'Pb210': stop = 15
    #elif info['chsp'] == 'Th232': stop = 21 # should not be a stop group
    elif info['chsp'] == 'Ra228': stop = 22
    elif info['chsp'] == 'Th228': stop = 23
    #elif info['chsp'] == 'U235':  stop = 41 # should not be a stop group
    elif info['chsp'] == 'Pa231': stop = 42
    ### handle the GRND group number better...
    #elif info['chsp'] == 'GRND':  stop = 16
    #elif info['chsp'] == 'GRND':  stop = 24
    #elif info['chsp'] == 'GRND':  stop = 43
    elif info['chsp'] == 'GRND':
        if   start in [11,12,13,14,15]: stop = 16
        elif start in [21,22,23]:       stop = 24
        elif start in [41,42]:          stop = 43
        ### special case for external Tl208 gammas
        elif start in [0]:              stop = 43
        else: stop = -1
    else: stop = -1

    if start != -1 and stop != -1:
        return TCut('((groupNo >= '+str(start)+') && (groupNo < '+str(stop)+'))')
    else:
        print('ERROR: groupNo not found for -->', info['chst'])
        sys.exit()


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


def histparam12000(E):
    """
    Global default histogram parameters for:
    number of bins, min, max, and bins/keV
    """
    E = int(E)
    # changed hi energy to go to 12 MeV
    if E == 0:
        hmin = 0
        hmax = 200
        bpkv = 12
        bins = (hmax-hmin)*bpkv
    elif E == 1:
        hmin = 0
        hmax = 12000
        bpkv = 1
        bins = (hmax-hmin)*bpkv
    elif E == 2:
        hmin = 0
        hmax = 12000
        bpkv = 1
        bins = (hmax-hmin)*bpkv
    else:
        print('ERROR: histparam12000() : invalid E = {0}'.format(E))
        sys.exit()
        
    pars = [bins, hmin, hmax, bpkv]
    return pars


def histparam20000(E):
    """
    Global default histogram parameters for:
    number of bins, min, max, and bins/keV
    """
    # changed high energy to 20 MeV
    # changed low energy to 1000 keV

    E = int(E)
    if E == 0:
        hmin = 0
        #hmax = 200
        hmax = 1000
        bpkv = 12
        bins = (hmax-hmin)*bpkv
    elif E == 1:
        hmin = 0
        hmax = 20000
        bpkv = 1
        bins = (hmax-hmin)*bpkv
    elif E == 2:
        hmin = 0
        hmax = 12000
        bpkv = 1
        bins = (hmax-hmin)*bpkv
    else:
        print('ERROR: histparam20000() : invalid E = {0}'.format(E))
        sys.exit()
        
    pars = [bins, hmin, hmax, bpkv]
    return pars


def rawHistPars(E):
    # for checking/plotting charge in ADC
    E = int(E)
    if E == 0:
        hmin = 0
        hmax = 1e7
        bpkv = 1
        bins = int((hmax-hmin)*bpkv)
    elif E == 1:
        hmin = 0
        hmax = 1e7
        bpkv = 1
        bins = int((hmax-hmin)*bpkv)
    elif E == 2:
        hmin = 0
        hmax = 1e7
        bpkv = 1
        bins = int((hmax-hmin)*bpkv)
    else:
        print('ERROR: rawHistPars() : invalid E = {0}'.format(E))
        sys.exit()
        
    pars = [bins, hmin, hmax, bpkv]
    return pars

    
def histparams(E):
    E = int(E)
    if E not in [0, 1, 2]:
        print('ERROR: histparams() : invalid E = {0}'.format(E))
        sys.exit()
    #return histparam12000(E)
    #return rawHistPars(E)
    return histparam20000(E)


def triangleSmoothing(bkgs, key, s=0):

    #print('DEBUG: triangle smoothing', key)

    ### get the number of bins
    bins = bkgs[key]['hist'].GetNbinsX()

    ### create a different hist to save the smoothed data to
    newhist = deepcopy(bkgs[key]['hist'])

    for i in range(s, bins-s):
        ### "i" is the bin you are smoothing
        ### "s" is the +/- bin range you are smoothing over
        smoothed = 0

        ### grab previous (s) bins
        for j, k in enumerate(range(-s, 0)):
            smoothed += bkgs[key]['hist'].GetBinContent(i+k)*(j+1)

        ### grab current bin
        smoothed += bkgs[key]['hist'].GetBinContent(i)*(s+1)

        ### grab following (s) bins
        for j, k in enumerate(range(1, s+1)):
            smoothed += bkgs[key]['hist'].GetBinContent(i+k)*(s-j)

        ### normalize
        newhist.SetBinContent(i, smoothed/float((s+1.)**2.))

    ### save the new hist back into the dictionary
    bkgs[key]['hist'] = deepcopy(newhist)
    del newhist
    
    return bkgs


def averageSmoothing(bkgs, key, s=0):

    #print('DEBUG: average smoothing', key)
    
    ### get the number of bins
    bins = bkgs[key]['hist'].GetNbinsX()

    ### create a different hist to save the smoothed data to
    newhist = deepcopy(bkgs[key]['hist'])

    for i in range(s, bins-s):
        ### "i" is the bin you are smoothing
        ### "s" is the +/- bin range you are smoothing over
        smoothed = 0

        ### grab previous (s) bins
        for j, k in enumerate(range(-s, 0)):
            smoothed += bkgs[key]['hist'].GetBinContent(i+k)

        ### grab current bin
        smoothed += bkgs[key]['hist'].GetBinContent(i)

        ### grab following (s) bins
        for j, k in enumerate(range(1, s+1)):
            smoothed += bkgs[key]['hist'].GetBinContent(i+k)

        ### normalize
        newhist.SetBinContent(i, smoothed/float(2.*s+1.))

    ### save the new hist back into the dictionary
    bkgs[key]['hist'] = deepcopy(newhist)
    del newhist
    
    return bkgs


def rootSmoothing(bkgs, key, s=0):
    
    #print('DEBUG: root TH1 smoothing', key)
    
    bkgs[key]['hist'].Smooth(s)

    return bkgs


def smooth(bkgs, smoothwhat, s=0):

    if int(s) <= 0:
        return bkgs
    
    smoothed = []
    for key in bkgs:
        for this in smoothwhat:
            if this in key and key not in smoothed:
                bkgs = rootSmoothing(bkgs, key, s)
                #bkgs = averageSmoothing(bkgs, key, s)
                #bkgs = triangleSmoothing(bkgs, key, s)
                smoothed.append(key)
    
    return bkgs


def extendHist(hist1, hist2):
    
    name  = hist1.GetName()
    bins1 = hist1.GetNbinsX()
    if bins1 == 1:
        bins1 = 0
    bins2 = hist2.GetNbinsX()
    if bins2 == 1:
        bins2 = 0
    bins  = bins1 + bins2
    hist3 = TH1F(name, name, bins, 0, bins)
    if bins1 > 0:
        for n in range(bins1):
            hist3.SetBinContent(n, hist1.GetBinContent(n))
    if bins2 > 0:
        for n in range(bins2):
            hist3.SetBinContent(n+bins1, hist2.GetBinContent(n))
        
    return hist3


def fitPrep(hists, key, params, zeros=0):
    
    # params needs to be a list of length 3 [rebin, kevfmin, kevfmax]
    rebin   = params[0]
    kevfmin = params[1]
    kevfmax = params[2]
    if kevfmin == kevfmax:
        temp = TH1F(key, key, 1, 0, 1)
        return temp
    
    # if you rebin now, you need to deepcopy the original hist
    hist = deepcopy(hists[key]['hist'])
    if rebin > 1: hist.Rebin(rebin)
    
    bins   = hist.GetXaxis().GetNbins()
    kevmin = hist.GetXaxis().GetBinUpEdge(0)
    kevmax = hist.GetXaxis().GetBinUpEdge(bins)
    #print bins, kevmin, kevmax
    
    # keV to bin number scaling
    bpk      = bins/(kevmax-kevmin)
    startbin = int(kevfmin*bpk)
    stopbin  = int(kevfmax*bpk)
    newbins  = int(stopbin-startbin)
    #print newbins, startbin, stopbin
    
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


