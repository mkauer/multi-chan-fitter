#!/usr/bin/env python
######################################################################
# funcs_new.py
# 
# Some new funcs to help with fitter hist building
# 
# version: 2021-05-18
# 
# email me: mkauer@physics.wisc.edu
######################################################################

import sys
import os
import re
import socket
from copy import deepcopy
import numpy as np
from numpy import exp, sqrt, log, pi
from scipy.special import erf

from ROOT import *
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from funcs_misc import *



def ScaleEnergy500(bkgs, shiftkey, shiftlist):
    shifted = []
    for key in bkgs:
        if 'cS' in key and 'e0' in key:
            if shiftkey in key and key not in shifted:
                xstal = int(key.split('-')[0][1]) - 1
                if shiftlist[xstal] != 0:
                    nbins = bkgs[key]['hist'].GetNbinsX()
                    S = shiftlist[xstal]
                    #print('shifting {0} by {1} bins'.format(key, S))
                    if S < 0:
                        for i in range(0, nbins, 1):
                            j = i+S
                            if j >= 0 and j < nbins:
                                bkgs[key]['hist'].SetBinContent(j, bkgs[key]['hist'].GetBinContent(i))
                            else:
                                bkgs[key]['hist'].SetBinContent(j, 0)
                    if S > 0:
                        for i in range(nbins, 0, -1):
                            j = i+S
                            if j >= 0 and j < nbins:
                                bkgs[key]['hist'].SetBinContent(j, bkgs[key]['hist'].GetBinContent(i))
                            else:
                                bkgs[key]['hist'].SetBinContent(j, 0)
                shifted.append(key)
    return bkgs


def unifyHistBins(data):
    # works for data and MC
    for key in data:
        bits = key.split('-')
        E = -1
        for bit in bits:
            if bit[0] == 'e' and len(bit) == 2:
                E = int(bit[1])
        if E not in [0, 1, 2]:
            print('ERROR: unifyHistBins() : invalid E = {0} [{1}]'.format(E, key))
            sys.exit()
        
        # pars = [bins, hmin, hmax, bpkv]
        pars = histparams(E)
        hist = data[key]['hist']
        hbins = hist.GetXaxis().GetNbins()
        
        if hbins > pars[0]:
            #print('INFO: truncating hist {0}'.format(key))
            temp = TH1F(key, key, pars[0], pars[1], pars[2])
            for i in range(pars[0]):
                temp.SetBinContent(i+1, hist.GetBinContent(i+1))
            data[key]['hist'] = deepcopy(temp)
        elif hbins < pars[0]:
            #print('INFO: extending hist {0}'.format(key))
            temp = TH1F(key, key, pars[0], pars[1], pars[2])
            for i in range(hbins):
                temp.SetBinContent(i+1, hist.GetBinContent(i+1))
            for i in range(hbins, pars[0]):
                temp.SetBinContent(i+1, 0)
            data[key]['hist'] = deepcopy(temp)
        else:
            continue
    
    return data


def makeFitBinRanges(xstals, fitranges):
    """
    bin = 0;        underflow bin
    bin = 1;        first bin with low-edge xlow INCLUDED
    bin = nbins;    last bin with upper-edge xup EXCLUDED
    bin = nbins+1;  overflow bin
    """
    # start with hist bin 1
    
    binranges = [{} for x in range(9)]
    prev_nbins = 1
    nbins = 0
    for i in range(9):
        for C in ['S', 'M']:
            for E in ['0', '1', '2']:
                if i+1 in xstals:
                    fitpars = fitranges[i][C+E]
                    #print(i+1, C+E, fitpars)
                    histpars = histparams(E)
                    #print(i+1, C+E, histpars)
                    bins = histpars[0]/float(fitpars[0])
                    bpk = bins/float(histpars[2]-histpars[1])
                    prev_nbins += nbins
                    nbins = int(bpk*(fitpars[2]-fitpars[1]))
                    #print(i+1, C+E, bins, bpk, nbins)
                    binranges[i][C+E] = [prev_nbins, prev_nbins+nbins]
                else:
                    binranges[i][C+E] = [prev_nbins+nbins, prev_nbins+nbins]

                #print(i+1, C+E, prev_nbins, nbins)
                
    tot_nbins = prev_nbins+nbins-1
    """
    # print out sanity check
    for i in range(9):
        for C in ['S', 'M']:
            for E in ['0', '1', '2']:
                print('{0} {1} = {2}'.format(str(i+1), C+E, binranges[i][C+E]))
    """
    return binranges, tot_nbins


def cropHist(hists, key, params):

    # params needs to be a list of length 3 [rebin, kevfmin, kevfmax]
    rebin   = params[0]
    kevfmin = params[1]
    kevfmax = params[2]
    if kevfmin == kevfmax:
        temp = TH1F(key, key, 1, 0, 1)
        return temp
    
    # if you rebin now, you need to deepcopy the original hist
    hist = deepcopy(hists[key]['hist'])
    if rebin > 1:
        hist.Rebin(rebin)
    
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
        
    for i in range(newbins):
        temp.SetBinContent(i+1, hist.GetBinContent(startbin+i))

    del hist
    #temp.Draw()
    #input('enter to continue')
    return temp


def addToFitHist(fithist, hist, binrange):

    startbin = binrange[0]
    stopbin = binrange[1]
    if startbin == stopbin:
        return fithist
    
    nbins = hist.GetXaxis().GetNbins()
    #print('fithist bins = {0}  hist bins = {1}'.format(stopbin-startbin, nbins))
    for i in range(nbins):
        bincontent = fithist.GetBinContent(startbin+i)
        addcontent = hist.GetBinContent(i+1)
        #print(i, bincontent, addcontent)
        fithist.SetBinContent(startbin+i, bincontent+addcontent)

    return fithist


def getSkewFunc(i):
    # only for internal pb210 alphas
    C = i+1
    xmin = 1000
    xmax = 5000
    npx = (xmax-xmin)*10
    
    if C == 1:
        # c1 is best described with 4 gaussians
        p = [
            9.86409e+00,
            2.60850e+03,
            4.94162e+01,
            1.04216e+02,
            2.72777e+03,
            1.09594e+02,
            1.03362e+02,
            3.00968e+03,
            9.37064e+01,
            3.99333e+01,
            3.21718e+03,
            6.86090e+01
        ]
        func = TF1("func",
                   "[0]*ROOT::Math::normal_pdf(x, [2], [1]) + "
                   "[3]*ROOT::Math::normal_pdf(x, [5], [4]) + "
                   "[6]*ROOT::Math::normal_pdf(x, [8], [7]) + "
                   "[9]*ROOT::Math::normal_pdf(x, [11], [10])"
                   , xmin, xmax)
    
    elif C == 2:
        # c2 is best described with 4 gaussians
        p = [
            2.45817e+01,
            2.78070e+03,
            6.16788e+01,
            9.53228e+01,
            2.93878e+03,
            7.49503e+01,
            4.33176e+01,
            3.18164e+03,
            7.20358e+01,
            2.76205e+00,
            3.35342e+03,
            4.07312e+01
        ]
        func = TF1("func",
                   "[0]*ROOT::Math::normal_pdf(x, [2], [1]) + "
                   "[3]*ROOT::Math::normal_pdf(x, [5], [4]) + "
                   "[6]*ROOT::Math::normal_pdf(x, [8], [7]) + "
                   "[9]*ROOT::Math::normal_pdf(x, [11], [10])"
                   , xmin, xmax)
    
    elif C == 3:
        # c3 is best described with 3 gaussians
        p = [
            1.78997e+01,
            2.69859e+03,
            8.59922e+01,
            2.25988e+01,
            2.83230e+03,
            5.38046e+01,
            1.89998e+01,
            2.98436e+03,
            7.06889e+01
        ]
        func = TF1("func",
                   "[0]*ROOT::Math::normal_pdf(x, [2], [1]) + "
                   "[3]*ROOT::Math::normal_pdf(x, [5], [4]) + "
                   "[6]*ROOT::Math::normal_pdf(x, [8], [7])"
                   , xmin, xmax)
    
    elif C == 4:
        # c4 is weird because the skew sigma is negative
        # need to look at a different fit
        """
        p = [
            5.32905e+01,
            2.98261e+03,
            1.07637e+02,
            -2.45587e+01,
            6.28300e+01,
            3.10023e+03,
            7.43147e+01,
            -3.23123e+10
        ]
        func = TF1("func",
            "([0]*ROOT::Math::normal_pdf(x, [2], [1])*ROOT::Math::normal_cdf(x, [3], [1])) + "
            "([4]*ROOT::Math::normal_pdf(x, [6], [5])*ROOT::Math::normal_cdf_c(x, [7], [5]))"
            , xmin, xmax)
        """
        # trying 3 gaussians for the weird tail
        p = [
            2.06520e+01,
            2.90890e+03,
            1.91932e+02,
            1.87622e+01,
            2.93014e+03,
            5.06025e+01,
            2.22460e+01,
            3.11140e+03,
            5.96238e+01
        ]
        func = TF1("func",
                   "[0]*ROOT::Math::normal_pdf(x, [2], [1]) + "
                   "[3]*ROOT::Math::normal_pdf(x, [5], [4]) + "
                   "[6]*ROOT::Math::normal_pdf(x, [8], [7])"
                   , xmin, xmax)
    
    elif C == 6:
        """
        p = [
            1.74236e+02,
            2.45069e+03,
            1.62452e+02,
            5.61403e+01,
            1.06604e+02,
            3.02751e+03,
            7.51212e+01,
            3.40676e+01
        ]
        func = TF1("func",
            "([0]*ROOT::Math::normal_pdf(x, [2], [1])*ROOT::Math::normal_cdf(x, [3], [1])) + "
            "([4]*ROOT::Math::normal_pdf(x, [6], [5])*ROOT::Math::normal_cdf_c(x, [7], [5]))"
            , xmin, xmax)
        """
        
        # try 2 skewgaus + 1 gaus
        p = [
            1.67279e+02,
            2.45330e+03,
            1.50104e+02,
            5.74175e+01,
            7.60304e+01,
            3.02368e+03,
            5.73436e+01,
            3.13194e+01,
            2.43030e+01,
            2.87646e+03,
            1.09060e+02
        ]
        func = TF1("func",
            "([0]*ROOT::Math::normal_pdf(x, [2], [1])*ROOT::Math::normal_cdf(x, [3], [1])) + "
            "([4]*ROOT::Math::normal_pdf(x, [6], [5])*ROOT::Math::normal_cdf_c(x, [7], [5])) + "
            "([8]*ROOT::Math::normal_pdf(x, [10], [9]))"
            , xmin, xmax)

        """
        # try 2 skewgaus + 2 gaus for high energy tail - Po216 version
        # this didn't help at all
        p = [
            1.41839e+02,
            2.45330e+03,
            1.51216e+02,
            4.50000e+01,
            6.77567e+01,
            3.02288e+03,
            5.95688e+01,
            3.34970e+01,
            1.81370e+01,
            2.87047e+03,
            1.06364e+02,
            5.00000e-01,
            3.14500e+03,
            2.30000e+02
        ]
        func = TF1("func",
            "([0]*ROOT::Math::normal_pdf(x, [2], [1])*ROOT::Math::normal_cdf(x, [3], [1])) + "
            "([4]*ROOT::Math::normal_pdf(x, [6], [5])*ROOT::Math::normal_cdf_c(x, [7], [5])) + "
            "([8]*ROOT::Math::normal_pdf(x, [10], [9])) + "
            "([11]*ROOT::Math::normal_pdf(x, [13], [12]))"      
            , xmin, xmax)
        """
        
        
    elif C == 7:
        """
        p = [
            1.65282e+02,
            2.45036e+03,
            1.61662e+02,
            6.14838e+01,
            1.14664e+02,
            3.03628e+03,
            8.81275e+01,
            4.56556e+01
        ]
        func = TF1("func",
            "([0]*ROOT::Math::normal_pdf(x, [2], [1])*ROOT::Math::normal_cdf(x, [3], [1])) + "
            "([4]*ROOT::Math::normal_pdf(x, [6], [5])*ROOT::Math::normal_cdf_c(x, [7], [5]))"
            , xmin, xmax)
        """
        # try 2 skewgaus + 1 gaus
        p = [
            1.48310e+02,
            2.45900e+03,
            1.32815e+02,
            6.95863e+01,
            8.25291e+01,
            3.02487e+03,
            5.90734e+01,
            5.13633e+01,
            3.04874e+01,
            2.85116e+03,
            1.26970e+02
        ]
        func = TF1("func",
            "([0]*ROOT::Math::normal_pdf(x, [2], [1])*ROOT::Math::normal_cdf(x, [3], [1])) + "
            "([4]*ROOT::Math::normal_pdf(x, [6], [5])*ROOT::Math::normal_cdf_c(x, [7], [5])) + "
            "([8]*ROOT::Math::normal_pdf(x, [10], [9]))"
            , xmin, xmax)
    
    else:
        print('WARNING: alpha skew params not found for crystal {0}'.format(C))
        return None
    

    func.SetParameters(np.asarray(p))
    func.SetNpx(npx)
    return func


def setQis(info):
    if 'nai-surf' in info['floca']:
        qis = [1, 2]
    elif 'teflon-surf' in info['floca']:
        #qis = [2]
        qis = [1, 2]
    elif info['floca'] == 'teflon':
        #qis = [2]
        qis = [1, 2]
    elif info['floca'] == 'internal' and info['chst'] == 'Pb210':
        if info['xstl'] in [5, 8]: # special case for C5 and C8
            qis = [1, 2]
        else:
            qis = [1, 2] # for new Po210 fits - 2023-08-30
            #qis = [1]
    else:
        qis = [1, 2]

    #return [1, 2] # testing C6 same Q1/Q2 ratio
    return qis


def getQratio(i, q, info):
    i = int(i)
    q = int(q)
    C = i+1

    # testing C6 same Q1/Q2 ratio
    #if q==1: return 0.3865
    #if q==2: return 0.6135
    
    # special case for C5 and C8
    if C in [5, 8]: return 0.5

    # special case for nai surface
    if q==1 and 'nai-surf' in info['floca']: return 0.5
    if q==2 and 'nai-surf' in info['floca']: return 0.5

    # special case for teflon surface
    if q==1 and 'teflon-surf' in info['floca']: return 0.5
    if q==2 and 'teflon-surf' in info['floca']: return 0.5

    # using parameterized combined Qs for Pb210 
    #if q==1 and info['floca']=='internal' and info['chst']=='Pb210': return 1.0
    #if q==2 and info['floca']=='internal' and info['chst']=='Pb210': return 0.0
    
    # using parameterized separate Qs for Pb210 - 2023-08-30
    if C==1 and q==1 and info['floca']=='internal' and info['chst']=='Pb210': return 0.529
    if C==1 and q==2 and info['floca']=='internal' and info['chst']=='Pb210': return 0.471

    if C==2 and q==1 and info['floca']=='internal' and info['chst']=='Pb210': return 0.258
    if C==2 and q==2 and info['floca']=='internal' and info['chst']=='Pb210': return 0.742
    
    if C==3 and q==1 and info['floca']=='internal' and info['chst']=='Pb210': return 0.351
    if C==3 and q==2 and info['floca']=='internal' and info['chst']=='Pb210': return 0.649
    
    if C==4 and q==1 and info['floca']=='internal' and info['chst']=='Pb210': return 0.401
    if C==4 and q==2 and info['floca']=='internal' and info['chst']=='Pb210': return 0.599
    
    if C==6 and q==1 and info['floca']=='internal' and info['chst']=='Pb210': return 0.504
    if C==6 and q==2 and info['floca']=='internal' and info['chst']=='Pb210': return 0.496
    
    if C==7 and q==1 and info['floca']=='internal' and info['chst']=='Pb210': return 0.618
    if C==7 and q==2 and info['floca']=='internal' and info['chst']=='Pb210': return 0.382
    
    
    
    # from the dedicated Th228 fits
    # assumed Rn222 and Pa231 would have same Q fraction
    if C==1 and q==1 and info['floca']=='internal' and info['chst']=='Th228': return 0.246
    if C==1 and q==2 and info['floca']=='internal' and info['chst']=='Th228': return 0.754
    
    if C==2 and q==1 and info['floca']=='internal' and info['chst']=='Th228': return 0.461
    if C==2 and q==2 and info['floca']=='internal' and info['chst']=='Th228': return 0.539
    
    if C==3 and q==1 and info['floca']=='internal' and info['chst']=='Th228': return 0.392
    if C==3 and q==2 and info['floca']=='internal' and info['chst']=='Th228': return 0.608
    
    if C==4 and q==1 and info['floca']=='internal' and info['chst']=='Th228': return 0.481
    if C==4 and q==2 and info['floca']=='internal' and info['chst']=='Th228': return 0.519
    
    if C==6 and q==1 and info['floca']=='internal' and info['chst']=='Th228': return 0.161
    if C==6 and q==2 and info['floca']=='internal' and info['chst']=='Th228': return 0.839
    
    if C==7 and q==1 and info['floca']=='internal' and info['chst']=='Th228': return 0.193
    if C==7 and q==2 and info['floca']=='internal' and info['chst']=='Th228': return 0.807
    
    
    # split in half by default
    return 0.5


def dataCutsByEvent501(chain, i, C, E):
    # https://cupwiki.ibs.re.kr/Kims/EventSelectionforSet2
    
    i = int(i)
    Cx = i+1
    E = int(E)
    
    # LS veto selection criteria
    if Cx == 9:
        
        # energy cut
        energy = (chain.GetBranch('BLSVeto').GetLeaf('Charge').GetValue()) / 143.8
        if energy <= 0:
            return False
        
        if C == 'S':
            # muon paddle cuts
            BTop = chain.GetBranch('BTop').GetLeaf('Charge').GetValue()
            BBottom = chain.GetBranch('BBottom').GetLeaf('Charge').GetValue()
            BLeft = chain.GetBranch('BLeft').GetLeaf('Charge').GetValue()
            BRight = chain.GetBranch('BRight').GetLeaf('Charge').GetValue()
            BFront = chain.GetBranch('BFront').GetLeaf('Charge').GetValue()
            BRear = chain.GetBranch('BRear').GetLeaf('Charge').GetValue()
            if BTop > 4000 or BBottom > 4000 \
               or BLeft > 4000 or BRight > 4000 \
               or BFront > 4000 or BRear > 4000:
                return False

            if (BTop + BBottom + BLeft + BRight + BFront + BRear) > 8000:
                return False
        
        elif C == 'M':
            # muon paddle cuts
            BTop = chain.GetBranch('BTop').GetLeaf('Charge').GetValue()
            BBottom = chain.GetBranch('BBottom').GetLeaf('Charge').GetValue()
            BLeft = chain.GetBranch('BLeft').GetLeaf('Charge').GetValue()
            BRight = chain.GetBranch('BRight').GetLeaf('Charge').GetValue()
            BFront = chain.GetBranch('BFront').GetLeaf('Charge').GetValue()
            BRear = chain.GetBranch('BRear').GetLeaf('Charge').GetValue()
            if BTop > 4000 or BBottom > 4000 \
               or BLeft > 4000 or BRight > 4000 \
               or BFront > 4000 or BRear > 4000:
                return False

            if (BTop + BBottom + BLeft + BRight + BFront + BRear) > 8000:
                return False

            isCoinc = chain.GetBranch('BLSVeto').GetLeaf('isCoincident').GetValue()
            if not isCoinc:
                return False

            dt0 = chain.GetBranch('BMuon').GetLeaf('totalDeltaT0').GetValue()
            if dt0/1.e6 <= 30:
                return False
        
        else:
            print('ERROR: invalid channel [{0}]'.format(C))
            sys.exit()
    
    
    # Crystal selection criteria
    elif Cx < 9:

        if E:
            energy = chain.GetBranch('crystal{0}'.format(Cx)).GetLeaf('energyD').GetValue()
        else:
            energy = chain.GetBranch('crystal{0}'.format(Cx)).GetLeaf('energy').GetValue()
        
        isCoinc = chain.GetBranch('BLSVeto').GetLeaf('isCoincident').GetValue()
        if not isCoinc:
            return False
        
        dt0 = chain.GetBranch('BMuon').GetLeaf('totalDeltaT0').GetValue()
        if dt0/1.e6 <= 30:
            return False

        # noise cuts are too strict for C5 and C8, skip
        if Cx not in [5, 8]:
            rqcn = chain.GetBranch('crystal{0}'.format(Cx)).GetLeaf('rqcn').GetValue()
            if rqcn <= -1:
                return False

            pmt1_nc = chain.GetBranch('pmt{0}1'.format(Cx)).GetLeaf('nc').GetValue()
            pmt2_nc = chain.GetBranch('pmt{0}2'.format(Cx)).GetLeaf('nc').GetValue()
            if pmt1_nc <= 0 or pmt2_nc <= 0:
                return False

            pmt1_t1 = chain.GetBranch('pmt{0}1'.format(Cx)).GetLeaf('t1').GetValue()
            pmt2_t1 = chain.GetBranch('pmt{0}2'.format(Cx)).GetLeaf('t1').GetValue()
            if pmt1_t1 <= 0 or pmt2_t1 <= 0:
                return False
        
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

        # the alpha cut
        if E:
            pmt1_d1_5 = chain.GetBranch('pmt{0}1'.format(Cx)).GetLeaf('rqtD1_5').GetValue()
            pmt2_d1_5 = chain.GetBranch('pmt{0}2'.format(Cx)).GetLeaf('rqtD1_5').GetValue()
            if E == 2:
                if energy < 1000:
                    return False
                if (pmt1_d1_5 + pmt2_d1_5)/2. > alpha[i]:
                    return False
            if E == 1:
                if energy >= 1000 and (pmt1_d1_5 + pmt2_d1_5)/2. <= alpha[i]:
                    return False

        # single/multiple hit cuts
        LS = (chain.GetBranch('BLSVeto').GetLeaf('Charge').GetValue()) / 143.8
        LS_hit = True if LS > 80 else False
        C_hit = 0
        for j in range(8):
            if j == i:
                continue
            # why don't we use energy for this cut?
            cnc = chain.GetBranch('crystal{0}'.format(j+1)).GetLeaf('nc').GetValue()
            if cnc >= 4:
                C_hit += 1
        if C == 'S':
            if LS_hit or C_hit:
                return False
        elif C == 'M':
            if not (LS_hit or C_hit):
                return False
        else:
            print('ERROR: invalid channel [{0}]'.format(C))
            sys.exit()

        # BDT cut
        if E == 0:
            bdt = chain.GetLeaf('bdt').GetValue(Cx)
            bdtA = chain.GetLeaf('bdtA').GetValue(Cx)

            if Cx == 1:
                #if bdt <= -0.2 or bdtA <= -0.07:
                if bdt <= -0.18 or bdtA <= -0.07:
                    return False
            elif Cx == 2:
                #if 9.20842e-07*exp(-104.504*bdt) + 0.170872 - 9.70874*bdt >= energy:
                if 1.52e-06*exp(-107.840*bdt) - 1.94 - 29.288*bdt >= energy:
                    return False
            elif Cx == 3:
                #if 6.81804e-08*exp(-101.856*bdt) + 0.148344 - 4.04826*bdt >= energy:
                if 1.51e-09*exp(-176.686*bdt) + 0.34 - 5.469*bdt >= energy:
                    return False
            elif Cx == 4:
                #if 1.37726e-07*exp(-111.842*bdt) + 0.446818 - 2.13498*bdt >= energy:
                if 1.86e-11*exp(-137.413*bdt) - 10.1 - 81.519*bdt >= energy:
                    return False
            elif Cx == 5:
                #if bdt <= -0.2 or bdtA <= -0.07:
                if bdt <= -0.22 or bdtA <= -0.1:
                    return False
            elif Cx == 6:
                #if 0.000315623*exp(-75.3998*bdt) - 0.313482 - 13.6302*bdt >= energy:
                if 3.30e-07*exp(-124.845*bdt) - 2.74 - 52.57*bdt >= energy:
                    return False
            elif Cx == 7:
                #if 1.45552e-05*exp(-88.7196*bdt) + 0.566336 - 7.57773*bdt >= energy:
                if 6.01e-08*exp(-120.177*bdt) - 0.24 - 31.615*bdt >= energy:
                    return False
            elif Cx == 8:
                #if bdt <= -0.2 or bdtA <= -0.07:
                if bdt <= -0.40 or bdtA <= -0.1:
                    return False
            else:
                print('ERROR: invalid crystal [{0}]'.format(Cx))
                sys.exit()
    else:
        print('ERROR: invalid crystal [{0}]'.format(Cx))
        sys.exit()
    
    return True


def dataCalibByEvent501(i, En, energy):
    i = int(i)
    Cx = i + 1
    En = int(En)

    # LS veto calib
    if Cx == 9:
        """
        # the classic "poly43"
        A = 4.0381226152e-08
        B = -0.000143935116267
        C = 1.15785808395
        D = 8.83648917811e-29
        return A*energy**3 + B*energy**2 + C*energy + D
        """
        """
        # n01
        A = 4.9118806973988345e-08
        B = -0.00017985370456661143
        C = 1.1918514156593005
        D = -2.9975829231342
        return A*energy**3 + B*energy**2 + C*energy + D
        """
        # n02 - looks pretty good
        A = 5.8788706575037026e-08
        B = -0.00021960488721729316
        C = 1.2294719215153322
        D = -6.315023331151148
        return A*energy**3 + B*energy**2 + C*energy + D

    
    # crystal low-energy
    if En == 0:
        if Cx == 1: # low-energy
            """
            # poly3
            A = -1.3294559078502743e-05
            B = 0.0008016340469073624
            C = 0.991977667645866
            D = -0.019766403922790943
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            # xe33p3 - with 200 point
            p0 = 9.258672197971397e-06
            p1 = -0.0026193102346206442
            p2 = 0.16919709267267818
            p3 = -3.1365362995871062
            zerox = 46.390
            x = energy
            if energy < zerox:
                return energy
            else:
                return x + ((1+erf((x-33.)/(.1)))*(p0*x**3 + p1*x**2 + p2*x + p3))
            
        
        elif Cx == 2: # low-energy
            """
            # poly3
            A = -1.500717582974678e-05
            B = 0.0009887540389559842
            C = 0.9868482301908721
            D = 0.004006929027896427
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            # xe33p3 - with 200 point
            p0 = 9.514402977519678e-06
            p1 = -0.0026916720431469734
            p2 = 0.17387395675895947
            p3 = -3.2232876983848
            zerox = 46.391
            x = energy
            if energy < zerox:
                return energy
            else:
                return x + ((1+erf((x-33.)/(.1)))*(p0*x**3 + p1*x**2 + p2*x + p3))
            
        
        elif Cx == 3: # low-energy
            """
            # poly3
            A = -1.0021160815609775e-05
            B = 0.0004817996633758037
            C = 0.999829990056206
            D = -0.05175734620563088
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            # xe33p3 - with 200 point
            p0 = 9.811727917358307e-06
            p1 = -0.0027757867470467076
            p2 = 0.17930751410890067
            p3 = -3.3240153531717755
            zerox = 46.391
            x = energy
            if energy < zerox:
                return energy
            else:
                return x + ((1+erf((x-33.)/(.1)))*(p0*x**3 + p1*x**2 + p2*x + p3))
            
        
        elif Cx == 4: # low-energy
            """
            # poly3
            A = -1.500717582974678e-05
            B = 0.0009887540389559842
            C = 0.9868482301908721
            D = 0.004006929027896427
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            # xe33p3 - with 200 point
            p0 = 9.811727917358307e-06
            p1 = -0.0027757867470467076
            p2 = 0.17930751410890067
            p3 = -3.3240153531717755
            zerox = 46.391
            x = energy
            if energy < zerox:
                return energy
            else:
                return x + ((1+erf((x-33.)/(.1)))*(p0*x**3 + p1*x**2 + p2*x + p3))

            
        elif Cx == 6: # low-energy
            """
            # poly3
            A = -1.971668166600867e-05
            B = 0.0014737450489396252
            C = 0.9742724382958354
            D = 0.05881029315449602
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            """
            # p0g0
            p0 = 1.075897147543578
            p1 = -12.739162401930667
            p2 = 118.96570651557325
            p3 = 45.8795702145503
            return p0*energy + p1*exp(-0.5*((energy-p2)/(p3))**2)
            """
            """
            # p0g1
            p0 = 1.012542663092105
            p1 = 0.27651655063310504
            p2 = 140.03137079445665
            p3 = 34.52641217881225
            return p0*energy + p1*(energy-p2)*exp(-0.5*((energy-p2)/(p3))**2)
            """
            """
            # p0g0 - v2
            p0 = 1.0742729002816178
            p1 = -12.453646998186175
            p2 = 118.18310158848722
            p3 = 45.654444893054226
            zerox = 42.62
            if energy < zerox:
                return energy
            else:
                return p0*energy + p1*exp(-0.5*((energy-p2)/(p3))**2)
            """
            """
            # p1g0
            p0 = 1.035145682889945
            p1 = -1.3219509977430535
            p2 = -6.1567192802657855
            p3 = 98.55343132785725
            p4 = 21.495034641530523
            zerox = 45.966
            if energy < zerox:
                return energy
            else:
                return p0*energy + p1 + p2*exp(-0.5*((energy-p3)/(p4))**2)
            """
            """
            # poly1 - try basic line to 80keV peak
            p0 = 0.9019607843137253
            p1 = 4.705882352941181
            zerox = 48.0
            if energy < zerox:
                return energy
            else:
                return p0*energy + p1
            """
            # xe33p3 - with 200 point
            p0 = 9.041448250244587e-06
            p1 = -0.0025578287017501123
            p2 = 0.16522069386699081
            p3 = -3.0627217076692945
            zerox = 46.388
            x = energy
            if energy < zerox:
                return energy
            else:
                return x + ((1+erf((x-33.)/(.1)))*(p0*x**3 + p1*x**2 + p2*x + p3))
            
            
        elif Cx == 7: # low-energy
            """
            # poly3
            A = -1.971668166600867e-05
            B = 0.0014737450489396252
            C = 0.9742724382958354
            D = 0.05881029315449602
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            """
            # p0g0
            p0 = 1.0505969310793273
            p1 = -7.7710896665498534
            p2 = 100.75521465195206
            p3 = 37.35504113822099
            return p0*energy + p1*exp(-0.5*((energy-p2)/(p3))**2)
            """
            # xe33p3 - with 200 point
            p0 = 9.30058880215625e-06
            p1 = -0.002631154165154706
            p2 = 0.16995963401542158
            p3 = -3.150620314205207
            zerox = 46.389
            x = energy
            if energy < zerox:
                return energy
            else:
                return x + ((1+erf((x-33.)/(.1)))*(p0*x**3 + p1*x**2 + p2*x + p3))
            
        else:
            return energy

        
    # crystal high-energy
    elif En == 1:
        if Cx == 1: # high-energy
            """
            A = -9.20702415903431e-09
            B = 3.1582891760183624e-05
            C = 0.9731283907370495
            D = 3.492470586962138
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            # for neutrons
            """
            A = 3.196944827665124e-15
            B = -2.981914465100653e-11
            C = 8.925660048168864e-08
            D = -0.00010376055403377454
            E = 1.0406860506992937
            F = -2.014860418216931
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
            """
            A = 3.602138700084496e-15
            B = -3.323076072071916e-11
            C = 9.912220565571808e-08
            D = -0.00011512134124263555
            E = 1.0451372879256497
            F = -2.235404420611547
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
        
        elif Cx == 2: # high-energy
            # for neutrons
            """
            A = 1.6948778686943696e-15
            B = -1.4494311236995248e-11
            C = 4.3553422373913564e-08
            D = -5.407316154016161e-05
            E = 1.0244084989218878
            F = -2.126621693953173
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
            """
            A = 2.4136503386011915e-15
            B = -2.0614828037442594e-11
            C = 6.185052288730988e-08
            D = -7.664789370700738e-05
            E = 1.0345184753958556
            F = -2.99725128253043
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
        
        elif Cx == 3: # high-energy
            """
            A = -4.668682649387045e-09
            B = 1.5052344932980456e-05
            C = 0.9874023184977824
            D = 1.6686064222603654
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            # for neutrons
            """
            A = 2.5211596609478237e-15
            B = -2.2913088875737618e-11
            C = 6.875906074476028e-08
            D = -8.226999832656399e-05
            E = 1.033728782723561
            F = -1.7988187051363322
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
            """
            A = 2.21291799608444e-15
            B = -2.031927860732983e-11
            C = 6.120669357442806e-08
            D = -7.341038129573121e-05
            E = 1.030159311388201
            F = -1.6124986888548951
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
        
        elif Cx == 4: # high-energy
            """
            A = -7.966374771112426e-09
            B = 2.544049246165629e-05
            C = 0.9801665355547313
            D = 2.2024520321640186
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            # for neutrons
            """
            A = 2.4620030100103767e-15
            B = -2.218245227911145e-11
            C = 6.233369632829584e-08
            D = -6.629170390094409e-05
            E = 1.0226895100321736
            F = -0.7589541402007167
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
            """
            A = 2.1716056943476836e-15
            B = -1.9683676028828895e-11
            C = 5.485539912613195e-08
            D = -5.722194195028613e-05
            E = 1.0188805330134927
            F = -0.5433872378465342
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
        
        elif Cx == 5: # high-energy
            """
            A = 0.813447107922342
            B = 13.863912116535856
            return A*energy + B
            """
            # for neutrons
            A = 1.287147512440593e-10
            B = -1.5919892208931184e-05
            C = 0.8617383079362346
            D = 0.42313767279705417
            return A*energy**3 + B*energy**2 + C*energy + D
        
        elif Cx == 6: # high-energy
            """
            A = -4.623190276816897e-09
            B = 1.6693196676849473e-05
            C = 0.9855987614226132
            D = 1.789741155036786
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            # for neutrons
            A = 6.715711222916333e-16
            B = -6.7075705824691325e-12
            C = 1.913234218174146e-08
            D = -1.830689730510006e-05
            E = 1.0042085806248093
            F = 0.17540811843090862
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
        
        elif Cx == 7: # high-energy
            """
            A = -9.93101832337613e-09
            B = 3.6976416493554536e-05
            C = 0.9661209376849521
            D = 4.747875210154714
            return A*energy**3 + B*energy**2 + C*energy + D
            """
            # for neutrons
            """
            A = 1.006970304584831e-14
            B = -8.509104310584813e-11
            C = 2.448168780847727e-07
            D = -0.0002808467701494711
            E = 1.1093253143120627
            F = -5.2644472738988455
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
            """
            A = 2.3692044301608016e-14
            B = -1.9077835787944898e-10
            C = 5.284769930437812e-07
            D = -0.0005873506581351882
            E = 1.2220712167491847
            F = -10.20632848257608
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
        
        elif Cx == 8: # high-energy
            """
            A = 0.8011463218736552
            B = -21.712318299469676
            return A*energy + B
            """
            # for neutrons
            A = 2.5144270503067e-10
            B = -2.325285315009937e-05
            C = 0.8025672088693712
            D = 21.91843846381534
            return A*energy**3 + B*energy**2 + C*energy + D

        else:
            return energy

        
    # crystal alphas
    elif En == 2:
        if Cx == 7: # alphas
            A = 1.2992745062332737e-14
            B = -9.808397355150091e-11
            C = 2.588343347431707e-07
            D = -0.0002802711310597177
            E = 1.1090809943187847
            F = -7.871191739422699
            return A*energy**5 + B*energy**4 + C*energy**3 + D*energy**2 + E*energy + F
        
        elif Cx == 5: # alphas
            # use neutron calib
            A = 1.287147512440593e-10
            B = -1.5919892208931184e-05
            C = 0.8617383079362346
            D = 0.42313767279705417
            return A*energy**3 + B*energy**2 + C*energy + D

        elif Cx == 8: # alphas
            # use neutron calib
            A = 2.5144270503067e-10
            B = -2.325285315009937e-05
            C = 0.8025672088693712
            D = 21.91843846381534
            return A*energy**3 + B*energy**2 + C*energy + D

        else:
            return energy

    else:
        print('ERROR: no calibration for C{0} E{1}'.format(Cx, En))
        sys.exit()

