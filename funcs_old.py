#!/usr/bin/env python
######################################################################
# funcs_old.py
# 
# Some old funcs that are still used
# 
# version: 2021-04-02
# 
# email me: mkauer@physics.wisc.edu
######################################################################

import os
import sys
import re
import socket
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from funcs_misc import *


def getRuntimes(data):
    # initialize the runtimes list
    runtimes = [{} for x in range(9)]
    for i in range(9):
        runtimes[i] = {'S': [0, 1, 2], 'M': [0, 1, 2]}
    
    for key in sortDataKeys92(data):
        X, C, E = getXCEFromKey(key)
        runtimes[X][C][E] = data[key]['runtime']
        #print('INFO: {0} runtime = {1} seconds'.format(key, data[key]['runtime']))
        
    return runtimes


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


def getXCEFromKey(key):
    # format like x9-data-SET3-V0000415-cS-e0
    # or x2-pmt-Pb210_GRND-cM-e1_generated
    bits = key.split('-')
    X = int(bits[0][1])-1
    C = str(bits[-2][1])
    E = int(bits[-1][1])
    return X, C, E


def scaleData411(data, dru=0):
    """
    Scale for DRU or not
    """
    for key in data:
        runtime = data[key]['runtime']
        if runtime == 0:
            #print('WARNING: runtime = {0} for {1}'.format(runtime, key))
            continue
        i = int(key.split('-')[0][1])-1
        days = 1.
        xkgs = 1.
        keVperBin = 1.
        if dru:
            days = float(runtime/86400.)
        xkgs = cmass(i)
        keVperBin = 1./float(data[key]['pars'][3])
        scale = float(1./(days*xkgs*keVperBin))
        data[key]['hist'].Scale(scale)
        data[key]['druScale'] = scale

    return data


def scaleBkgs411(bkgs, runtimes=None):
    
    for key in bkgs:
        
        #print 'scaling -->',key
        
        bits  = key.split('-')
        x     = int(bits[0][1])-1
        loca  = bits[1]
        isos  = bits[2]
        f     = None
        c     = str(bits[-2][1])
        e     = int(bits[-1][1])
        if len(bits) == 6:
            f = int(bits[3][1])-1
        if len(bits) == 7:
            f = int(bits[3][1])-1
        
        # conversion factors for things
        keVperBin  = 1./float(bkgs[key]['pars'][3])
        if runtimes is not None:
            day    = runtimes[x][c][e]
        else:
            day    = 86400.

        # get the mass of the crystal/lsveto of origin
        if f is not None:
            nmass  = cmass(f)
            surf   = cmass(f)
            #surf   = surfArea(f)
        else:
            nmass  = cmass(x)
            surf   = cmass(x)
            #surf   = surfArea(x)

        #print('{0} --> {1}'.format(key, nmass))

        # get the mass of the target detector
        xkgs = cmass(x)
            
        xpmts      = 2.
        pmts       = 16.
        plastic    = 1800. # set same as lsveto for comparison
        lsveto     = 1800.
        air        = 1.
        # cu shield is 3000 - should change this soon
        copper     = 1800. # set same as lsveto for now
        steel      = 1600.
        innersteel = 4000.
        
        generated  = float(bkgs[key]['generated'])
        if generated < 1:
            print("WARNING: 0 events generated for -->", key)
            bkgs[key]['scale'] = 1
            continue
        
        if   loca == 'internal':    norm = (nmass)
        elif loca == 'reflector':   norm = (surf)
        elif 'surf' in loca:        norm = (surf)
        elif 'teflon' in loca:      norm = (surf)
        elif loca == 'copper':      norm = (surf)
        elif loca == 'cucase':      norm = (surf)
        elif loca == 'coppercase':  norm = (surf)
        elif loca == 'pmt':         norm = (pmts)
        elif loca == 'xpmt':        norm = (xpmts)
        elif loca == 'lsveto':      norm = (lsveto)
        elif loca == 'plastic':     norm = (plastic)
        elif loca == 'lsvetoair':   norm = (air)
        elif loca == 'airshield':   norm = (air)
        elif loca == 'cushield':    norm = (copper)
        elif loca == 'steel':       norm = (steel)
        elif loca == 'innersteel':  norm = (innersteel)
        elif loca == 'gamma':       norm = (innersteel)
        else:
            print("ERROR: no background scaling for -->", loca)
            norm = 1
            #sys.exit()
        
        scale = bkgs[key]['info']['acti']*(norm)*(1./1000.)*(1./generated)*(day)*(1./xkgs)*(1./keVperBin)
        
        bkgs[key]['hist'].Scale(scale)
        bkgs[key]['scale'] = scale
    
    return bkgs


def scaleSigs411(sigkeys, sigs, runtimes=None):

    for key in sigkeys:
        
        if 'fitscale' not in sigs[key]:
            #print 'WARNING: no fitscale for', key
            sigs[key]['info']['fitacti'] = 0
            sigs[key]['info']['fiterro'] = 0
            continue

        bits  = key.split('-')
        x     = int(bits[0][1])-1
        loca  = bits[1]
        isos  = bits[2]
        f     = None
        c     = str(bits[-2][1])
        e     = int(bits[-1][1])
        if len(bits) == 6:
            f = int(bits[3][1])-1
        if len(bits) == 7:
            f = int(bits[3][1])-1
        
        # conversion factors for things
        keVperBin  = 1./float(sigs[key]['pars'][3])
        if runtimes is not None:
            day    = runtimes[x][c][e]
        else:
            day    = 86400.
            
        # get the mass of the crystal/lsveto of origin
        if f is not None:
            nmass  = cmass(f)
            surf   = cmass(f)
            #surf   = surfArea(f)
        else:
            nmass  = cmass(x)
            surf   = cmass(x)
            #surf   = surfArea(x)
            
        # get the mass of the target detector
        xkgs = cmass(x)
        
        xpmts      = 2.
        pmts       = 16.
        plastic    = 1800. # set same as lsveto for comparison
        lsveto     = 1800.
        air        = 1.
        # cu shield is 3000 - should change this soon
        copper     = 1800. # set same as lsveto for now
        steel      = 1600.
        innersteel = 4000.
        
        generated  = float(sigs[key]['generated'])
        if generated < 1:
            print("WARNING: 0 events generated for -->", key)
            sigs[key]['info']['fitacti'] = 0
            sigs[key]['info']['fiterro'] = 0
            continue
        
        if   loca == 'internal':    norm = (nmass)
        elif loca == 'reflector':   norm = (surf)
        elif 'surf' in loca:        norm = (surf)
        elif 'teflon' in loca:      norm = (surf)
        elif loca == 'copper':      norm = (surf)
        elif loca == 'cucase':      norm = (surf)
        elif loca == 'coppercase':  norm = (surf)
        elif loca == 'pmt':         norm = (pmts)
        elif loca == 'xpmt':        norm = (xpmts)
        elif loca == 'lsveto':      norm = (lsveto)
        elif loca == 'plastic':     norm = (plastic)
        elif loca == 'lsvetoair':   norm = (air)
        elif loca == 'airshield':   norm = (air)
        elif loca == 'cushield':    norm = (copper)
        elif loca == 'steel':       norm = (steel)
        elif loca == 'innersteel':  norm = (innersteel)
        elif loca == 'gamma':       norm = (innersteel)
        else:
            print("ERROR: no signal scaling for -->", loca)
            norm = 1
            #sys.exit()
        
        fitActivity = sigs[key]['fitscale']*(1./norm)*(1000.)*(generated)*(1./day)*(xkgs)*(keVperBin)
        
        sigs[key]['info']['fitacti'] = fitActivity
        if 'fiterror' in sigs[key]:
            sigs[key]['info']['fiterro'] = fitActivity * sigs[key]['fiterror']
        else:
            print('WARNING: no \"fiterror\" for [{0}]'.format(key))
            sigs[key]['info']['fiterro'] = 0
            
            
    return sigs


def ScaleEnergy410(bkgs, binShift, shiftwhat):

    shifted = []
    for key in bkgs:
        if 'cS' in key and 'e0' in key:
            for this in shiftwhat:
                if this in key and key not in shifted:
                    xstal = int(key.split('-')[0][1]) - 1
                    if binShift[xstal] > 0:
                        #print 'DEBUG: energy shifting hist ', key
                        bkgs[key]['hist'] = MCBinShift400(bkgs[key]['hist'], binShift[xstal])
                    shifted.append(key)
    
    return bkgs


def MCBinShift400(histo, binShift):

    ### this currently only works with shifting hists to lower energy
    
    """
    name = hist.GetName()
    title = hist.GetTitle()
    
    newhist = TH1F(key, longNames(i), par[0], par[1], par[2])
    pars = getPars(hist)
    bins = pars[0]
    """
    
    for n in range(histo.GetNbinsX()+1):
        histo.SetBinContent(n, histo.GetBinContent(n+binShift))

    return histo


def calib302(i, Engy, Chan):
    """
    Tweaked high energy crystal calibrations
    Tweaked the lsveto calibration
    Now adjust for different lsveto data sources
    """

    Cx = i+1
    
    if Cx == 9: # C9 / lsveto
        # return same calib for all Engy
        # return same calib for all Chan
        
        edep = '(BLSVeto.Charge)'
        selection = '('+edep+' / 143.8)'
        #return edep, selection
    
        ### try again based off the poly4 data points - "poly41_data"
        #A = 6.14216259037e-08
        #B = -0.000190223995721
        #C = 1.17695224672
        #D = 7.44693291248e-28
        
        ### try again and force a 1460 point - "poly42_data"
        #A = 8.81422486824e-08
        #B = -0.000277066019752
        #C = 1.21663237155
        #D = 2.41933873005e-20

        ### try un-saturating the high energy a little bit = "poly43"
        ### but this still seems to give the best multi-crystal fits
        A = 4.0381226152e-08
        B = -0.000143935116267
        C = 1.15785808395
        D = 8.83648917811e-29

        ### try un-saturating the high energy a little more = "poly44"
        ### this was too much!
        #A = 2.48681708029e-08
        #B = -0.000109806394499
        #C = 1.14377998622
        #D = -6.09628447918e-27
        
        ### try un-saturating the high energy just a tiny = "poly45"
        #A = 4.22053979702e-08
        #B = -0.000147948294267
        #C = 1.15951351987
        #D = -5.7576250417e-22
        
        ### try un-saturating the high energy just a tiny = "poly46"
        #A = 4.80251584646e-08
        #B = -0.000160751767355
        #C = 1.16479495252
        #D = 4.5440374256e-28
        
        ### try un-saturating the high energy just a tiny = "poly47"
        #A = 5.22150630743e-08
        #B = -0.000169969557496
        #C = 1.16859729095
        #D = -1.00574241496e-21
        
        ### try tweaking the low energy part of the fit a litle = "poly50"
        # that went the wrong direction, oops...
        #A = 4.02698077117e-08
        #B = -9.26205577368e-05
        #C = 1.04107520387
        #D = -6.35274464267e-22
        
        ### try tweaking the low energy part of the fit a litle = "poly51"
        #A = 5.7582815735e-08
        #B = -0.000245923913043
        #C = 1.30856625259
        #D = 1.37579420822e-27

        # return same calib for all E
        selection = '(('+str(A)+'*'+selection+'**3) '\
                    '+ ('+str(B)+'*'+selection+'**2) '\
                    '+ ('+str(C)+'*'+selection+') + '+str(D)+')'
        
        return edep, selection

    if Engy == 0:
        # default low energy production calibration
        edep = '(crystal{0}.energy)'.format(Cx)
        selection = '('+edep+')'

        return edep, selection
        
    elif Engy == 1:
        # default high energy production calibration
        edep = '(crystal{0}.energyD)'.format(Cx)
        selection = '('+edep+')'

        # tweaks per crystal
        if Cx == 1: # C1
            A = -9.20702415903431e-09
            B = 3.1582891760183624e-05
            C = 0.9731283907370495
            D = 3.492470586962138
            
        if Cx == 2: # C2
            return edep, selection
        
        if Cx == 3: # C3
            A = -4.668682649387045e-09
            B = 1.5052344932980456e-05
            C = 0.9874023184977824
            D = 1.6686064222603654
            
        if Cx == 4: # C4
            A = -7.966374771112426e-09
            B = 2.544049246165629e-05
            C = 0.9801665355547313
            D = 2.2024520321640186
            
        if Cx == 5: # C5
            A = 0.813447107922342
            B = 13.863912116535856
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            return edep, selection
        
        if Cx == 6: # C6
            A = -4.623190276816897e-09
            B = 1.6693196676849473e-05
            C = 0.9855987614226132
            D = 1.789741155036786
        
        if Cx == 7: # C7
            #A = -6.584663125047853e-09
            #B = 2.4904010290725124e-05
            #C = 0.9765006123749262
            #D = 3.45348577054466
            # first try to fix
            #A = -3.5286801110126694e-09
            #B = 6.692048421004304e-06
            #C = 1.0004270013767433
            #D = -1.3444385176574656
            # second try to fix
            #A = -9.787584607766034e-09
            #B = 3.3865923406543445e-05
            #C = 0.971206025113478
            #D = 3.6683048699245857
            # third try
            A = -9.93101832337613e-09
            B = 3.6976416493554536e-05
            C = 0.9661209376849521
            D = 4.747875210154714
            
        if Cx == 8: # C8
            A = 0.8011463218736552
            B = -21.712318299469676
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            return edep, selection

        
        selection = '(('+str(A)+'*'+edep+'**3) '\
                   '+ ('+str(B)+'*'+edep+'**2) '\
                   '+ ('+str(C)+'*'+edep+') + '+str(D)+')'

        return edep, selection    

    elif Engy == 2:
        ### use default high energy calib for alphas
        edep = '(crystal{0}.energyD)'.format(Cx)
        selection = '('+edep+')'

        # except for C5 and C8
        if Cx == 5:
            A = 0.813447107922342
            B = 13.863912116535856
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            return edep, selection
        if Cx == 8:
            A = 0.8011463218736552
            B = -21.712318299469676
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            return edep, selection
        
        return edep, selection

    else:
        print('WARNING: no calibration for x{0} c{1} e{2}'.format(Cx, Chan, Engy))
        sys.exit()


def cutsBDT302(i, C, E, edep, selection):
    """
    https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    https://cupwiki.ibs.re.kr/Kims/EventSelection
    """
    #-------------------------------------------------------
    # This is for the V00-04-15 data
    # New lsveto single-hit data
    #-------------------------------------------------------

    ### values for the alpha cut
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

    ### special cuts for the LS veto
    ### i = 0-7 is crystals
    ### i =  8  is lsveto
    if i == 8:
        ### skip low energy and alpha energy for lsveto
        #if E in [0, 2]:
        #    return TCut('0')
        
        if C == 'S':
            coinc  = '(1)'
            energy = '('+selection+' > 0.0)'
            #muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
            muons  = '(1)'
            # time cut doesn't make much difference...
            time   = '(BLSVeto.Time > 2300 && BLSVeto.Time < 2600)'
            #time   = '(1)'
            return TCut('({0} && {1} && {2} && {3})'.format(coinc, energy, muons, time))
            
        elif C == 'M':
            coinc  = '(BLSVeto.isCoincident == 1)'
            energy = '('+selection+' > 0.0)'
            muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
            #time   = '(BLSVeto.Time >= 2800 && BLSVeto.Time <= 2950)'
            time   = '(1)'
            return TCut('({0} && {1} && {2} && {3})'.format(coinc, energy, muons, time))
            
        else:
            print('ERROR: I do not know what to do with channel -->', C)
            print('Available channels are [S]Single-hits, [M]Multi-hits')
            sys.exit()


    ### same as Pushpa since 2017-12-19
    # only alphas for E=2
    if E == 2:
        alphaCut = '((crystal'+str(i+1)+'.energyD > 1000.) && ((pmt'+str(i+1)+'1.rqtD1_5+pmt'+str(i+1)+'2.rqtD1_5)/2. < '+str(alpha[i])+'))'
    else:
        alphaCut = '( ! ((crystal'+str(i+1)+'.energyD > 1000.) && ((pmt'+str(i+1)+'1.rqtD1_5+pmt'+str(i+1)+'2.rqtD1_5)/2. < '+str(alpha[i])+')) )'

        
    ### new BDT cuts for V00-04-12 from https://cupwiki.ibs.re.kr/Kims/EventSelection
    newBDT = [
        '((bdt[1]>-0.1 && (crystal1.energy + 20*bdt[1])>0) && ((crystal1.energy>=10 && bdtA[1]>-0.08) || (crystal1.energy>=8 && crystal1.energy<10 && bdtA[1]>-0.065) || (crystal1.energy>=7 && crystal1.energy<8 && bdtA[1]>0.005) || (crystal1.energy>=6 && crystal1.energy<7 && bdtA[1]>0.015) || (crystal1.energy>=5 && crystal1.energy<6 && bdtA[1]>0.01) || (crystal1.energy>=4 && crystal1.energy<5 && bdtA[1]>0.03) || (crystal1.energy>=3 && crystal1.energy<4 && bdtA[1]>0.035) || (crystal1.energy>=2 && crystal1.energy<3 && bdtA[1]>0.02) || (crystal1.energy<2 && bdtA[1]>0.02)))',
              
        #'(bdt[2]>0.0 && bdtA[2]>-0.05)',
        # new bdt cut for V00-04-14
        '(((9.20842e-07*TMath::Exp(-104.504*bdt[2])+0.170872)-(9.70874*bdt[2])) < crystal2.energy)',
              
        #'(bdt[3]>0.0 && bdtA[3]>-0.07)',
        # new bdt cut for V00-04-14
        '(((6.81804e-08*TMath::Exp(-101.856*bdt[3])+0.148344)-(4.04826*bdt[3])) < crystal3.energy)',

        #'(bdt[4]>-0.05 && (crystal4.energy + 40*bdt[4])>0 && bdtA[4]>-0.07)',
        # new bdt cut for V00-04-14
        '(((1.37726e-07*TMath::Exp(-111.842*bdt[4])+0.446818)-(2.13498*bdt[4])) < crystal4.energy)',
              
        '(bdt[5]>-0.2 && bdtA[5]>-0.07)',
              
        #'(bdt[6]>0.0 && (crystal6.energy + 20*bdt[6])>2.0 && bdtA[6]>-0.07)',
        # new bdt cut for V00-04-14
        '(((0.000315623*TMath::Exp(-75.3998*bdt[6])-0.313482)-(13.6302*bdt[6])) < crystal6.energy)',
              
        #'(bdt[7]>0.0 && (crystal7.energy + 20*bdt[7])>2.0 && bdtA[7]>-0.07)',
        # new bdt cut for V00-04-14
        '(((1.45552e-05*TMath::Exp(-88.7196*bdt[7])+0.566336)-(7.57773*bdt[7])) < crystal7.energy)',
        
        '(bdt[8]>-0.2 && bdtA[8]>-0.07)']

    bdtCut = newBDT[i]

    ### global noise cuts
    #===========================================================================
    coinc  = '(BLSVeto.isCoincident == 1)'
    muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'
    charge = '(crystal'+str(i+1)+'.rqcn > -1)'
    
    ### old cut
    #nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
    ### new cut for v00-04-14 from Govinda
    nc     = '(pmt'+str(i+1)+'1.nc > 0 && pmt'+str(i+1)+'2.nc > 0)'
    
    ### old cut
    #t1     = '(pmt'+str(i+1)+'1.t1 > 2 && pmt'+str(i+1)+'2.t1 > 2)'
    ### new cut for v00-04-14 from Govinda
    t1     = '(pmt'+str(i+1)+'1.t1 > 0 && pmt'+str(i+1)+'2.t1 > 0)'
    
    noiseCut = '('+coinc+' && '+muons+' && '+charge+' && '+nc+' && '+t1+')'
    #===========================================================================
    
    if C == 'S':
        #lsveto = '(BLSVeto.Charge/143.8 < 20)'
        lsveto = '(BLSVeto.Charge/143.8 < 80)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc < 4) && '
        ### remove extra '&&'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        if E == 0:
            masterCut = TCut(bdtCut+' && '+noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)
        else:
            masterCut = TCut(noiseCut+' && '+'('+lsveto+' && '+hits+')'+' && '+alphaCut)

        # lose cuts for C5 and C8
        if int(i) in [4, 7]:
            masterCut = TCut('('+lsveto+' && '+hits+')'+' && '+alphaCut)
        
        return masterCut
    
    elif C == 'M':
        #lsveto = '(BLSVeto.Charge/143.8 > 20)'
        lsveto = '(BLSVeto.Charge/143.8 > 80)'
        hits   = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc > 4) || '
        ### remove extra '||'
        hits = hits[:-4]+')'

        # BDT cuts ONLY for low energy spectrum!!
        if E == 0:
            masterCut = TCut(bdtCut+' && '+noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)
        else:
            masterCut = TCut(noiseCut+' && '+'('+lsveto+' || '+hits+')'+' && '+alphaCut)

        # lose cuts for C5 and C8
        if int(i) in [4, 7]:
            masterCut = TCut('('+lsveto+' || '+hits+')'+' && '+alphaCut)

        return masterCut

    else:
        print('ERROR: I do not know what to do with channel -->', C)
        print('Available channels are [S]Single-hits, [M]Multi-hits')
        sys.exit()


def updateBkgsFile300(xstals, bkgsfile, resultsfile, newdir, BF='BR'):
    """
    Only print out to file what was actually used in the fit
    """
    
    if os.path.exists(bkgsfile):
        with open(bkgsfile) as fbkgs:
            bkgslines = fbkgs.read().splitlines()
    else:
        print('WARNING: file not found -->', bkgsfile)
        return
    
    if os.path.exists(resultsfile):
        with open(resultsfile) as ffits:
            fitlines = ffits.read().splitlines()
    else:
        fitlines = ''
        
    newbkgs = os.path.join(newdir, bkgsfile.split('/')[-1][:-4]+'_update.txt')
    output = open(newbkgs, 'w')
    output.write('# -*- sh -*-\n')
    
    print('INFO: Updating bkgsfile --> {0}'.format(bkgsfile))
    print('      To a new bkgsfile --> {0}'.format(newbkgs))
    
    skip = 0
    for bline in bkgslines:

        if not bline:
            continue
        if bline.startswith('\"\"\"'):
            if skip == 0: skip = 1
            else: skip = 0
            continue
        if skip:
            continue
        if bline.startswith('#'):
            continue
        
        #bbits = filter(None, re.split("[ \s\t\n\r,:]+", bline.strip()))
        bbits = re.split("[ \s\t\n\r,:]+", bline.strip())
        
        # add [0] to xstals for v500 global mc
        # this should still work with older versions
        #if int(bbits[1]) not in xstals:
        if int(bbits[1]) not in xstals+[0]:
            continue
        
        if 'F' not in bbits[0]:
            for j in range(len(bbits)-1):
                output.write(bbits[j]+'\t')
            output.write(bbits[-1]+'\n')
            continue
        
        replaced = 0
        for fline in fitlines:
            if not replaced:
                #fline = fline.strip()
                if not fline:
                    continue
                #fbits = fline.split()
                #fbits = filter(None, re.split("[ \s\t\n\r,:]+", fline.strip()))
                fbits = re.split("[ \s\t\n\r,:]+", fline.strip())
                
                ### ------------------------------------------------------------
                ### if output format changes, this is generally the part that fails...
                #if fbits[-1] == 'mBq' or fbits[-3] == 'mBq':
                if fbits[0] == 'fit':
                    xstal = fbits[1].split('-')[0].split('x')[1]
                    loca  = fbits[1].split('-')[1]
                    chst  = fbits[1].split('-')[2].split('_')[0]
                    chsp  = fbits[1].split('-')[2].split('_')[1]
                    acti  = str(fbits[3])

                    lenbbits = len(bbits)
                    if bbits[1] == xstal and bbits[2].replace('-','') == loca \
                       and bbits[4] == chst and bbits[5] == chsp:
                        for i in range(lenbbits):
                            if i == 0:
                                output.write(BF+'\t')
                            elif i == 6:
                                if acti != '0.0': output.write(acti+'\t')
                                else: output.write(bbits[i]+'\t')
                            elif i == 7:
                            #elif i==(lenbbits-3):
                                output.write('0.1\t')
                            elif i == 8:
                            #elif i==(lenbbits-2):
                                output.write('10\t')
                            elif i==(lenbbits-1):
                                output.write(bbits[i])
                            else:
                                output.write(bbits[i]+'\t')
                        output.write('\n')
                        replaced = 1
        
        if not replaced:
            #output.write(bline+'\n')
            print('WARNING: Could not match -->', bbits[3:7])
            
    output.close()
    return


def makeTotal100(chan, E, par):
    total = []
    #par = histparam(E)
    for i in range(numX()):
        key  = 'x'+str(i+1)
        key += '-c'+chan
        key += '-e'+str(E)
        key += '-total'
        tot = TH1F(key, longNames(i), par[0], par[1], par[2])
        tot.SetLineColor(kGray+1)
        tot.SetMarkerColor(kGray+1)
        tot.SetLineWidth(1)
        total.append(tot)
    return total


def makeResid100(chan, E, par):
    resid = []
    #par = histparam(E)
    for i in range(numX()):
        key  = 'x'+str(i+1)
        key += '-c'+chan
        key += '-e'+str(E)
        key += '-resid'
        res = TH1F(key, longNames(i), par[0], par[1], par[2])
        res.SetLineColor(kBlack)
        res.SetMarkerColor(kBlack)
        res.SetLineWidth(1)
        resid.append(res)
    return resid

