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


def scaleBkgs411(bkgs, runtimes=None, one_mBq=False):
    
    for key in bkgs:
        
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

        # get num of generated events
        generated = float(bkgs[key]['generated'])
        if generated < 1:
            print("WARNING: 0 events generated for --> {0}".format(key))
            bkgs[key]['scale'] = 1
            continue

        # get keV per bin
        keVperBin = 1./float(bkgs[key]['pars'][3])

        # get seconds in a day or runtime in seconds
        if runtimes is not None:
            day = runtimes[x][c][e]
        else:
            day = 86400.

        # to scale counts/sec to mBq
        mBq = 1000.
        
        # get the mass of the target detector
        xkgs = cmass(x)

        # get the mass of the material of origin
        if f is not None:
            nmass = cmass(f)
            surf  = cmass(f)
            #surf  = surfArea(f)
        else:
            nmass = cmass(x)
            surf  = cmass(x)
            #surf  = surfArea(x)

        # get mass/number of other components
        xpmts      = 2.
        #if x == 8: xpmts = 1.  # this seems double counted (check code)
        pmts       = 16.
        vetopmts   = 18.
        plastic    = 1800.  # set same as lsveto
        lsveto     = 1800.
        film       = 1.
        air        = 1.
        copper     = 3000.
        steel      = 1600.
        innersteel = 4000.

        xfactor = 1.
        if x == 8: xfactor = 8.
        
        if   loca == 'internal':    norm = (nmass)
        elif loca == 'mystery':     norm = (nmass)
        elif loca == 'reflector':   norm = (surf)
        elif 'surf' in loca:        norm = (surf)
        elif 'teflon' in loca:      norm = (surf)
        elif loca == 'copper':      norm = (surf)
        elif loca == 'cucase':      norm = (surf)
        elif loca == 'coppercase':  norm = (surf)
        elif loca == 'pmt':         norm = (pmts)
        elif loca == 'pmtbase':     norm = (pmts)
        elif loca == 'xpmt':        norm = (xpmts)
        elif loca == 'xpmtbase':    norm = (xpmts)
        elif loca == 'lsveto':      norm = (lsveto)
        elif loca == 'vetopmt':     norm = (vetopmts)
        elif loca == 'film':        norm = (film)
        #elif 'plastic' in loca:     norm = (plastic)
        elif loca == 'plastic':     norm = (plastic)
        elif loca == 'xplastic':    norm = (plastic/xfactor)
        elif loca == 'lsvetoair':   norm = (air)
        elif loca == 'airshield':   norm = (air)
        elif loca == 'cushield':    norm = (copper)
        elif loca == 'steel':       norm = (steel)
        elif loca == 'innersteel':  norm = (innersteel)
        #elif 'gamma' in loca:       norm = (innersteel)
        elif loca == 'gamma':       norm = (innersteel)
        elif loca == 'xgamma':      norm = (innersteel)
        elif loca == 'neutron':
            ### neutrons are in different units (cm^2 * sec)
            # needs to be the area at which the events were generated
            # from MC primX0/Y0/Z0 they were generated at radius 3500mm
            gen_radius = 350.  # in cm
            gen_area = 4.*3.1416*(gen_radius**2)
            ### Eunju said twice as many events are really generated
            generated = (2*generated)
            ### normalize to integral of A5 measurements
            ### 3.04 = sum of my relative scaling factors
            ### and *1000 to put in units of mBq
            norm = gen_area * 1000 * (1/3.04)
            
        else:
            print("ERROR: no background scaling for --> {0}".format(loca))
            norm = 1
            sys.exit()

        if one_mBq:
            activity = 1.0
        else:
            activity = bkgs[key]['info']['acti']
            
        scale = (activity)*(norm)*(1./mBq)*(1./generated)*(day)*(1./xkgs)*(1./keVperBin)
        
        bkgs[key]['hist'].Scale(scale)
        bkgs[key]['scale'] = scale
    
    return bkgs


def scaleSigs411(sigkeys, sigs, runtimes=None):

    #for key in sigkeys:
    for key in sigs:
        
        if 'fitscale' not in sigs[key]:
            print('WARNING: no fitscale for --> {0}'.format(key))
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
        
        # get num of generated events
        generated = float(sigs[key]['generated'])
        if generated < 1:
            print("WARNING: 0 events generated for --> {0}".format(key))
            sigs[key]['info']['fitacti'] = 0
            sigs[key]['info']['fiterro'] = 0
            continue
                
        # get keV per bin
        keVperBin = 1./float(sigs[key]['pars'][3])
        
        # get seconds in a day or runtime in seconds
        if runtimes is not None:
            day = runtimes[x][c][e]
        else:
            day = 86400.

        # to scale counts/sec to mBq
        mBq = 1000.

        # get the mass of the target detector
        xkgs = cmass(x)

        # get the mass of the crystal/lsveto of origin
        if f is not None:
            nmass = cmass(f)
            surf  = cmass(f)
            #surf  = surfArea(f)
        else:
            nmass = cmass(x)
            surf  = cmass(x)
            #surf  = surfArea(x)
            
        # get mass/number of other components
        xpmts      = 2.
        #if x == 8: xpmts = 1.  # this seems double counted (check code)
        pmts       = 16.
        vetopmts   = 18.
        plastic    = 1800. # set same as lsveto
        lsveto     = 1800.
        film       = 1.
        air        = 1.
        copper     = 3000.
        steel      = 1600.
        innersteel = 4000.
        
        xfactor = 1.
        if x == 8: xfactor = 8.
        
        if   loca == 'internal':    norm = (nmass)
        elif loca == 'mystery':     norm = (nmass)
        elif loca == 'reflector':   norm = (surf)
        elif 'surf' in loca:        norm = (surf)
        elif 'teflon' in loca:      norm = (surf)
        elif loca == 'copper':      norm = (surf)
        elif loca == 'cucase':      norm = (surf)
        elif loca == 'coppercase':  norm = (surf)
        elif loca == 'pmt':         norm = (pmts)
        elif loca == 'pmtbase':     norm = (pmts)
        elif loca == 'xpmt':        norm = (xpmts)
        elif loca == 'xpmtbase':    norm = (xpmts)
        elif loca == 'lsveto':      norm = (lsveto)
        elif loca == 'vetopmt':     norm = (vetopmts)
        elif loca == 'film':        norm = (film)
        #elif 'plastic' in loca:     norm = (plastic)
        elif loca == 'plastic':     norm = (plastic)
        elif loca == 'xplastic':    norm = (plastic/xfactor)
        elif loca == 'lsvetoair':   norm = (air)
        elif loca == 'airshield':   norm = (air)
        elif loca == 'cushield':    norm = (copper)
        elif loca == 'steel':       norm = (steel)
        elif loca == 'innersteel':  norm = (innersteel)
        #elif 'gamma' in loca:       norm = (innersteel)
        elif loca == 'gamma':       norm = (innersteel)
        elif loca == 'xgamma':      norm = (innersteel)
        elif loca == 'neutron':
            ### neutrons are in different units (cm^2 * sec)
            # needs to be the area at which the events were generated
            # from MC primX0/Y0/Z0 they were generated at radius 3500mm
            gen_radius = 350.  # in cm
            gen_area = 4.*3.1416*(gen_radius**2)
            ### Eunju said twice as many events are really generated
            generated = (2*generated)
            ### normalize to integral of A5 measurements
            ### 3.04 = sum of my relative scaling factors
            ### and *1000 to put in units of mBq
            norm = gen_area * 1000 * (1/3.04)
            
        else:
            print("ERROR: no signal scaling for --> {0}".format(loca))
            norm = 1
            sys.exit()
        
        fitActivity = sigs[key]['fitscale']*(1./norm)*(mBq)*(generated)*(1./day)*(xkgs)*(keVperBin)
        
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
                        #print('DEBUG: energy shifting hist {0}'.format(key))
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


def calib302(i, Engy, Chan, norecal=False):

    # or don't recalibrate - for testing
    #norecal = True
    
    i = int(i)
    Engy = int(Engy)
    Cx = i+1

    if norecal:
        #print('INFO: not recalibrating data')
        if Cx == 9:
            #selection = '(BLSVeto.Charge)'
            selection = '(BLSVeto.Charge/143.8)'
            return selection
        if Engy == 0:
            #selection = '(crystal{0}.qc1_5)'.format(Cx)
            selection = '(crystal{0}.energy)'.format(Cx)
            return selection
        else:
            #selection = '(crystal{0}.rqcD1_5)'.format(Cx)
            selection = '(crystal{0}.energyD)'.format(Cx)
            return selection
        
            
    # LS veto calib
    if Cx == 9:

        edep = '(BLSVeto.Charge/143.8)'
        selection = '('+edep+')'
    
        ### the classic "poly43"
        #A = 4.0381226152e-08
        #B = -0.000143935116267
        #C = 1.15785808395
        #D = 8.83648917811e-29

        # n01
        #A = 4.9118806973988345e-08
        #B = -0.00017985370456661143
        #C = 1.1918514156593005
        #D = -2.9975829231342
        
        # n02 - looks pretty good
        A = 5.8788706575037026e-08
        B = -0.00021960488721729316
        C = 1.2294719215153322
        D = -6.315023331151148
        
        
        selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+({4}))'
                     .format(edep, A, B, C, D))
        
        return selection
    
    
    # low energy calib
    elif Engy == 0:
        # default low energy production calibration
        edep = '(crystal{0}.energy)'.format(Cx)
        #edep = '(crystal{0}.qc1_5)'.format(Cx)
        selection = '('+edep+')'

        if Cx == 1: # low energy
            """
            # poly3
            A = -1.3294559078502743e-05
            B = 0.0008016340469073624
            C = 0.991977667645866
            D = -0.019766403922790943
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            """
            # xe33p3 - with 200 point
            # still has a discontinuity - use by entry with zero crossing point
            p0 = 9.258672197971397e-06
            p1 = -0.0026193102346206442
            p2 = 0.16919709267267818
            p3 = -3.1365362995871062
            selection = ('({0} + (1+ROOT::Math::erf(({0}-33.)/.1)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            """
            
            # poly3x35v2 - 2022-11-16
            p0 = 2.2813057393143775e-05
            p1 = -0.0025546892730003994
            zerox = 35.0
            selection = ('(({0}<{3})*{0} + ({0}>={3})*'
                         '({1}*({0}-{3})**3 + {2}*({0}-{3})**2 + {0}))'
                         .format(edep, p0, p1, zerox))
            
            
            return selection

        
        elif Cx == 2: # low energy
            """
            # poly3
            A = -1.500717582974678e-05
            B = 0.0009887540389559842
            C = 0.9868482301908721
            D = 0.004006929027896427
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            """
            # xe33p3 - with 200 point
            # still has a discontinuity - use by entry with zero crossing point
            p0 = 9.514402977519678e-06
            p1 = -0.0026916720431469734
            p2 = 0.17387395675895947
            p3 = -3.2232876983848
            selection = ('({0} + (1+ROOT::Math::erf(({0}-33.)/.1)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            """
            
            # poly3x35v2 - 2022-11-16
            p0 = 1.8593703725424738e-05
            p1 = -0.0023084529959055745
            zerox = 35.0
            selection = ('(({0}<{3})*{0} + ({0}>={3})*'
                         '({1}*({0}-{3})**3 + {2}*({0}-{3})**2 + {0}))'
                         .format(edep, p0, p1, zerox))
            
            
            return selection

        
        elif Cx == 3: # low energy
            """
            # poly3
            A = -1.0021160815609775e-05
            B = 0.0004817996633758037
            C = 0.999829990056206
            D = -0.05175734620563088
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            """
            # xe33p3 - with 200 point
            # still has a discontinuity - use by entry with zero crossing point
            p0 = 9.811727917358307e-06
            p1 = -0.0027757867470467076
            p2 = 0.17930751410890067
            p3 = -3.3240153531717755
            selection = ('({0} + (1+ROOT::Math::erf(({0}-33.)/.1)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            """
            
            # poly3x35v2 - 2022-11-16
            p0 = 1.4921828348252173e-05
            p1 = -0.0019219510403052549
            zerox = 35.0
            selection = ('(({0}<{3})*{0} + ({0}>={3})*'
                         '({1}*({0}-{3})**3 + {2}*({0}-{3})**2 + {0}))'
                         .format(edep, p0, p1, zerox))
            
            
            return selection

        
        elif Cx == 4: # low energy
            """
            # poly3
            A = -1.500717582974678e-05
            B = 0.0009887540389559842
            C = 0.9868482301908721
            D = 0.004006929027896427
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            """
            # xe33p3 - with 200 point
            # still has a discontinuity - use by entry with zero crossing point
            p0 = 9.811727917358307e-06
            p1 = -0.0027757867470467076
            p2 = 0.17930751410890067
            p3 = -3.3240153531717755
            selection = ('({0} + (1+ROOT::Math::erf(({0}-33.)/.1)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            """
            """
            # poly3x35 - 2022-10-18
            p0 = 1.501437887099597e-05
            p1 = -0.0013275691087104651
            p2 = 0.919843226912901
            zerox = 35.0
            selection = ('(({0}<{4})*{0} + ({0}>={4})*'
                         '({1}*({0}-{4})**3 + {2}*({0}-{4})**2 + {3}*({0}-{4}) + {4}))'
                         .format(edep, p0, p1, p2, zerox))
            """
            """
            # poly5 - 2022-10-19
            #this shifts the pb210 peak too much
            p0 = -4.857486126973767e-10
            p1 = 4.153345638619488e-07
            p2 = -0.00011174470882442253
            p3 = 0.013741358580952545
            p4 = 0.17752347264219503
            p5 = 16.16611963104454
            zerox = 35.095
            selection = ('(({0}<{7})*{0} + ({0}>={7})*'
                         '( ({1}*{0}**5 + {2}*{0}**4 + {3}*{0}**3 + {4}*{0}**2 + {5}*{0} + {6}) ))'
                         .format(edep, p0, p1, p2, p3, p4, p5, zerox))
            """
            """
            # poly3 - 2022-10-19
            p0 = 1.8290139833404742e-05
            p1 = -0.004260526122702743
            p2 = 1.2303917391707309
            p3 = -3.6287561160219557
            zerox = 35.0
            selection = ('(({0}<{5})*{0} + ({0}>={5})*'
                         '( {1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3, zerox))
            """
            """
            # poly3 v2 tweaked 180keV data point - 2022-10-20
            p0 = 9.475418812685678e-06
            p1 = -0.0013869270699483254
            p2 = 0.9929848579197037
            p3 = 1.53825705190317
            zerox = 35.0
            selection = ('(({0}<{5})*{0} + ({0}>={5})*'
                         '( {1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3, zerox))
            """
            """
            # xpoly3z45 - 2022-10-20
            p0 = 8.504073511752402e-06
            p1 = -0.000933308814412496
            p2 = -0.0730316541129224
            p3 = 4.401441085508372
            zerox = 45.0
            selection = ('(({0}<{5})*{0} + ({0}>={5})*'
                         '( {0} + {1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3, zerox))
            """
            
            # poly3x35v2 - 2022-11-16
            p0 = 1.9589071979858656e-05
            p1 = -0.002693799013704164
            zerox = 35.0
            selection = ('(({0}<{3})*{0} + ({0}>={3})*'
                         '({1}*({0}-{3})**3 + {2}*({0}-{3})**2 + {0}))'
                         .format(edep, p0, p1, zerox))
            
            
            return selection

        
        elif Cx == 6: # low energy
            """
            # poly3
            A = -1.971668166600867e-05
            B = 0.0014737450489396252
            C = 0.9742724382958354
            D = 0.05881029315449602
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            """
            # p0g0
            p0 = 1.075897147543578
            p1 = -12.739162401930667
            p2 = 118.96570651557325
            p3 = 45.8795702145503
            selection = ('({1}*{0} + {2}*exp(-0.5*(({0}-{3})/({4}))**2))'
                         .format(edep, p0, p1, p2, p3))
            """
            """
            # p0g1
            p0 = 1.012542663092105
            p1 = 0.27651655063310504
            p2 = 140.03137079445665
            p3 = 34.52641217881225
            selection = ('({1}*{0} + {2}*({0}-{3})*exp(-0.5*(({0}-{3})/({4}))**2))'
                         .format(edep, p0, p1, p2, p3))
            """
            """
            # p0g0 - v2
            p0 = 1.0742729002816178
            p1 = -12.453646998186175
            p2 = 118.18310158848722
            p3 = 45.654444893054226
            selection = ('({1}*{0} + {2}*exp(-0.5*(({0}-{3})/({4}))**2))'
                         .format(edep, p0, p1, p2, p3))
            """
            """
            # xe48
            p0 = -0.049019607840496845
            p1 = 2.3529411762648573
            selection = ('({0} + (1+ROOT::Math::erf(({0}-48.)/.1)) * ({1}*{0} + {2}))'
                         .format(edep, p0, p1))
            """
            """
            # xe40
            p0 = -0.038860103626942984
            p1 = 1.5544041450777168
            selection = ('({0} + (1+ROOT::Math::erf(({0}-40.)/.1)) * ({1}*{0} + {2}))'
                         .format(edep, p0, p1))
            """
            """
            # xe50
            p0 = -0.05244755244755254
            p1 = 2.6223776223776283
            selection = ('({0} + (1+ROOT::Math::erf(({0}-50.)/.1)) * ({1}*{0} + {2}))'
                         .format(edep, p0, p1))
            """
            """
            # xe33p2
            p0 = -0.0011089519545521077
            p1 = 0.09221305045708757
            p2 = -1.8897613619854816
            selection = ('({0} + (1+ROOT::Math::erf(({0}-33.)/.1)) * ({1}*{0}**2 + {2}*{0} + {3}))'
                         .format(edep, p0, p1, p2))
            """
            """
            # xe33p3 - with 200 point
            # still has a discontinuity - use by entry with zero crossing point
            p0 = 9.041448250244587e-06
            p1 = -0.0025578287017501123
            p2 = 0.16522069386699081
            p3 = -3.0627217076692945
            selection = ('({0} + (1+ROOT::Math::erf(({0}-33.)/.1)) * '
                         '({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            """
            """
            # testing multi-hit calib - 2022-10-11
            p0 = 9.08711381736282e-06
            p1 = -0.0017094456420253348
            p2 = 1.0419993176810325
            p3 = -0.17455474782700844
            selection = ('{1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}'.format(edep, p0, p1, p2, p3))
            """
            """
            # with points out to 600 keV - 2022-10-13
            p0 = 1.5661107920566318e-05
            p1 = -0.003440225941090345
            p2 = 1.1407726352480092
            p3 = -1.088332097750311
            zerox = 52.628
            # poly3 - full data process for this one
            #selection = ('{1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}'.format(edep, p0, p1, p2, p3))
            # Figured out how to do "if-else" like statement in TFormula! 2022-10-13
            # poly3x50
            selection = ('(({0}<{5})*{0} + ({0}>={5})*({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3, zerox))
            """
            """
            # try the above again a different way - 2022-10-18
            # poly3x50
            p0 = 1.614928038384384e-05
            p1 = -0.0012862900590242035
            p2 = 0.9282888886957914
            zerox = 50.0
            selection = ('(({0}<{4})*{0} + ({0}>={4})*'
                         '({1}*({0}-{4})**3 + {2}*({0}-{4})**2 + {3}*({0}-{4}) + {4}))'
                         .format(edep, p0, p1, p2, zerox))
            """
            """
            # poly3x35 - 2022-10-18 - this is pretty good
            p0 = 1.6377990673646178e-05
            p1 = -0.0021384117996877247
            p2 = 0.9972430012260727
            zerox = 35.0
            selection = ('(({0}<{4})*{0} + ({0}>={4})*'
                         '({1}*({0}-{4})**3 + {2}*({0}-{4})**2 + {3}*({0}-{4}) + {4}))'
                         .format(edep, p0, p1, p2, zerox))
            """
            """
            # poly3 - 2022-10-19
            p0 = 1.5065019955578992e-05
            p1 = -0.003494926287521847
            p2 = 1.1844209887113952
            p3 = -2.8193626332800177
            zerox = 35.0
            selection = ('(({0}<{5})*{0} + ({0}>={5})*'
                         '( {1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3, zerox))
            """
            
            # poly3x35v2 - 2022-11-16
            p0 = 1.6587709628593584e-05
            p1 = -0.002260430243687135
            zerox = 35.0
            selection = ('(({0}<{3})*{0} + ({0}>={3})*'
                         '({1}*({0}-{3})**3 + {2}*({0}-{3})**2 + {0}))'
                         .format(edep, p0, p1, zerox))
            
            """
            # try the erf method again - 2022-10-13
            # xe33p3 - with lin50()
            p0 = 8.24034239801928e-06
            p1 = -0.0019582101088150643
            p2 = 0.10772836040310248
            p3 = -1.6635102995967714
            zerox = 43.214
            selection = ('(({0}<{5})*{0} + ({0}>={5})*({0} + '
                         '(1+ROOT::Math::erf(({0}-33.)/.1)) * '
                         '({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4})))'
                         .format(edep, p0, p1, p2, p3, zerox))
            """
            """
            # raw adc refit - 2022-10-17
            # p3g2
            p0 = 3.7089887772057903e-17
            p1 = -7.242128675325004e-11
            p2 = 0.00015211027391750048
            p3 = -4.479975232869551
            p4 = 228624.11263357548
            p5 = 127953.95569168364
            p6 = 1.032059954809929
            p7 = 15052.585637414644
            p8 = 15634.815392873556
            
            edep = '(crystal{0}.qc1_5)'.format(Cx)
            selection = ('( ({1}*{0}**3 + {2}*{0}**2 + {3}*{0}) '
                         '+ ({4}*exp(-0.5*(({0}-{5})/({6}))**2)) '
                         '+ ({7}*exp(-0.5*(({0}-{8})/({9}))**2)) )'
                         .format(edep, p0, p1, p2, p3, p4, p5, p6, p7, p8))
            """
            """
            # poly5 (force tail up) - 2022-10-19
            # did not look good at all
            p0 = -7.783841010884831e-10
            p1 = 7.047259715631928e-07
            p2 = -0.00021486903893859648
            p3 = 0.029450991345908655
            p4 = -0.795061282515554
            p5 = 35.575875272785204
            zerox = 36.6185
            selection = ('(({0}<{7})*{0} + ({0}>={7})*'
                         '( ({1}*{0}**5 + {2}*{0}**4 + {3}*{0}**3 + {4}*{0}**2 + {5}*{0} + {6}) ))'
                         .format(edep, p0, p1, p2, p3, p4, p5, zerox))
            """
            
            return selection

        
        elif Cx == 7: # low energy
            """
            # poly3
            A = -1.971668166600867e-05
            B = 0.0014737450489396252
            C = 0.9742724382958354
            D = 0.05881029315449602
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            """
            # p0g0
            p0 = 1.0505969310793273
            p1 = -7.7710896665498534
            p2 = 100.75521465195206
            p3 = 37.35504113822099
            selection = ('({1}*{0} + {2}*exp(-0.5*(({0}-{3})/({4}))**2))'
                         .format(edep, p0, p1, p2, p3))
            """
            """
            # xe33p3 - with 200 point
            # still has a discontinuity - use by entry with zero crossing point
            p0 = 9.30058880215625e-06
            p1 = -0.002631154165154706
            p2 = 0.16995963401542158
            p3 = -3.150620314205207
            selection = ('({0} + (1+ROOT::Math::erf(({0}-33.)/.1)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            """
            
            # poly3x35v2 - 2022-11-16
            p0 = 1.9753483176035022e-05
            p1 = -0.0024573733349558687
            zerox = 35.0
            selection = ('(({0}<{3})*{0} + ({0}>={3})*'
                         '({1}*({0}-{3})**3 + {2}*({0}-{3})**2 + {0}))'
                         .format(edep, p0, p1, zerox))
            
            
            return selection
            
        else:
            return selection

        
    # high energy calib
    elif Engy == 1:
        # default high energy production calibration
        edep = '(crystal{0}.energyD)'.format(Cx)
        selection = '('+edep+')'

        # tweaks per crystal
        if Cx == 1: # high energy
            """
            # poly3
            A = -9.20702415903431e-09
            B = 3.1582891760183624e-05
            C = 0.9731283907370495
            D = 3.492470586962138
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            # for neutrons
            """
            # poly5
            A = 3.196944827665124e-15
            B = -2.981914465100653e-11
            C = 8.925660048168864e-08
            D = -0.00010376055403377454
            E = 1.0406860506992937
            F = -2.014860418216931
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            """
            # poly5 v2
            A = 3.602138700084496e-15
            B = -3.323076072071916e-11
            C = 9.912220565571808e-08
            D = -0.00011512134124263555
            E = 1.0451372879256497
            F = -2.235404420611547
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            
            # xe2400p3 - c1n2
            p0 = 4.662591384599011e-08
            p1 = -0.00040965076865556624
            p2 = 1.1482651150501608
            p3 = -1040.8044816712918
            selection = ('({0} + (1+ROOT::Math::erf(({0}-2400.)/1.)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            
                        
            return selection

        
        elif Cx == 2: # high energy
            #return selection
            # for neutrons
            """
            # poly5
            A = 1.6948778686943696e-15
            B = -1.4494311236995248e-11
            C = 4.3553422373913564e-08
            D = -5.407316154016161e-05
            E = 1.0244084989218878
            F = -2.126621693953173
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            """
            # poly5 v2
            A = 2.4136503386011915e-15
            B = -2.0614828037442594e-11
            C = 6.185052288730988e-08
            D = -7.664789370700738e-05
            E = 1.0345184753958556
            F = -2.99725128253043
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            # xe3600p2 - c2n2
            p0 = 0.00019857404920856444
            p1 = -1.4299317283509279
            p2 = 2574.2345443203453
            selection = ('({0} + (1+ROOT::Math::erf(({0}-3600.)/1.)) * ({1}*{0}**2 + {2}*{0} + {3}))'
                         .format(edep, p0, p1, p2))
            
            return selection

        
        elif Cx == 3: # high energy
            """
            # poly3
            A = -4.668682649387045e-09
            B = 1.5052344932980456e-05
            C = 0.9874023184977824
            D = 1.6686064222603654
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            # for neutrons
            """
            # poly5
            A = 2.5211596609478237e-15
            B = -2.2913088875737618e-11
            C = 6.875906074476028e-08
            D = -8.226999832656399e-05
            E = 1.033728782723561
            F = -1.7988187051363322
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            """
            # poly5 v2
            A = 2.21291799608444e-15
            B = -2.031927860732983e-11
            C = 6.120669357442806e-08
            D = -7.341038129573121e-05
            E = 1.030159311388201
            F = -1.6124986888548951
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            # xe2400p3 - c3n2
            p0 = 2.5877555148993577e-08
            p1 = -0.00019727046410586484
            p2 = 0.4620786786031633
            p3 = -330.4422777774975
            selection = ('({0} + (1+ROOT::Math::erf(({0}-2400.)/1.)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            
            return selection

        
        elif Cx == 4: # high energy
            """
            # poly3
            A = -7.966374771112426e-09
            B = 2.544049246165629e-05
            C = 0.9801665355547313
            D = 2.2024520321640186
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            # for neutrons
            """
            # poly5
            A = 2.4620030100103767e-15
            B = -2.218245227911145e-11
            C = 6.233369632829584e-08
            D = -6.629170390094409e-05
            E = 1.0226895100321736
            F = -0.7589541402007167
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            """
            # poly5 v2
            A = 2.1716056943476836e-15
            B = -1.9683676028828895e-11
            C = 5.485539912613195e-08
            D = -5.722194195028613e-05
            E = 1.0188805330134927
            F = -0.5433872378465342
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            # xe2400p3 - c4n3
            p0 = 3.6860459987325365e-08
            p1 = -0.00033337030262278316
            p2 = 0.9533597156670861
            p3 = -877.4093733585614
            selection = ('({0} + (1+ROOT::Math::erf(({0}-2400.)/1.)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            
            return selection

        
        elif Cx == 5: # high energy
            """
            # poly1
            A = 0.813447107922342
            B = 13.863912116535856
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            """
            # for neutrons
            # poly3 v2
            A = 1.287147512440593e-10
            B = -1.5919892208931184e-05
            C = 0.8617383079362346
            D = 0.42313767279705417
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            
            return selection

        
        elif Cx == 6: # high energy
            """
            # poly3
            A = -4.623190276816897e-09
            B = 1.6693196676849473e-05
            C = 0.9855987614226132
            D = 1.789741155036786
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            # for neutrons
            """
            # poly5
            A = 6.715711222916333e-16
            B = -6.7075705824691325e-12
            C = 1.913234218174146e-08
            D = -1.830689730510006e-05
            E = 1.0042085806248093
            F = 0.17540811843090862
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            """
            # xe2600p2 - c6n2
            p0 = 3.116962097956178e-05
            p1 = -0.19404526941906328
            p2 = 293.81106266772673
            selection = ('({0} + (1+ROOT::Math::erf(({0}-2600.)/1.)) '
                         '* ({1}*{0}**2 + {2}*{0} + {3}))'
                         .format(edep, p0, p1, p2))
            """
            # xe2650p2 - c6n3
            p0 = 3.1597760543552394e-05
            p1 = -0.19854373181791274
            p2 = 304.2456159003718
            selection = ('({0} + (1+ROOT::Math::erf(({0}-2650.)/1.)) '
                         '* ({1}*{0}**2 + {2}*{0} + {3}))'
                         .format(edep, p0, p1, p2))
            
            return selection

        
        elif Cx == 7: # high energy
            """
            # poly3
            A = -9.93101832337613e-09
            B = 3.6976416493554536e-05
            C = 0.9661209376849521
            D = 4.747875210154714
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            """
            # for neutrons
            """
            # poly5
            A = 1.006970304584831e-14
            B = -8.509104310584813e-11
            C = 2.448168780847727e-07
            D = -0.0002808467701494711
            E = 1.1093253143120627
            F = -5.2644472738988455
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            """
            # poly5 v2
            A = 2.3692044301608016e-14
            B = -1.9077835787944898e-10
            C = 5.284769930437812e-07
            D = -0.0005873506581351882
            E = 1.2220712167491847
            F = -10.20632848257608
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+({6}))'
                         .format(edep, A, B, C, D, E, F))
            """
            # xe2600p3 - c7n2
            p0 = 1.728038721149522e-07
            p1 = -0.001496848531870412
            p2 = 4.266347526349869
            p3 = -4011.0083493580746
            selection = ('({0} + (1+ROOT::Math::erf(({0}-2600.)/1.)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            
            return selection
        
        
        elif Cx == 8: # high energy
            """
            # poly1
            A = 0.8011463218736552
            B = -21.712318299469676
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            """
            # for neutrons
            # poly3 v2
            A = 2.5144270503067e-10
            B = -2.325285315009937e-05
            C = 0.8025672088693712
            D = 21.91843846381534
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            
            return selection

        else:
            return selection

        
    # alpha calib
    elif Engy == 2:
        ### use default high energy calib for alphas
        edep = '(crystal{0}.energyD)'.format(Cx)
        selection = '('+edep+')'

        # try tweaking c7 alpha spectra
        if Cx == 7: # alphas
            """
            # poly5 - separate fit just for alphas
            A = 1.2992745062332737e-14
            B = -9.808397355150091e-11
            C = 2.588343347431707e-07
            D = -0.0002802711310597177
            E = 1.1090809943187847
            F = -7.871191739422699
            selection = ('(({1}*{0}**5)+({2}*{0}**4)+({3}*{0}**3)+({4}*{0}**2)+({5}*{0})+{6})'
                         .format(edep, A, B, C, D, E, F))
            """
            """
            # try the c7 high energy calib
            # nope - shifts the Q1 pb210 peak a bit
            # xe2600p3
            p0 = 1.728038721149522e-07
            p1 = -0.001496848531870412
            p2 = 4.266347526349869
            p3 = -4011.0083493580746
            selection = ('({0} + (1+ROOT::Math::erf(({0}-2600.)/1.)) '
                         '* ({1}*{0}**3 + {2}*{0}**2 + {3}*{0} + {4}))'
                         .format(edep, p0, p1, p2, p3))
            """
            # try erf for the alpha fit
            # xe3200p2 - c7a2
            p0 = 0.0003888833014851269
            p1 = -2.3375171812665165
            p2 = 3497.8899728451534
            selection = ('({0} + (1+ROOT::Math::erf(({0}-3200.)/1.)) '
                         '* ({1}*{0}**2 + {2}*{0} + {3}))'
                         .format(edep, p0, p1, p2))
            
            return selection

        
        # tweak C5 alphas
        elif Cx == 5: # alphas
            """
            # poly1
            A = 0.813447107922342
            B = 13.863912116535856
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            """
            # use neutron calib
            # poly3 v2
            A = 1.287147512440593e-10
            B = -1.5919892208931184e-05
            C = 0.8617383079362346
            D = 0.42313767279705417
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            
            return selection

        
        # tweak C8 alphas
        elif Cx == 8: # alphas
            """
            # poly1
            A = 0.8011463218736552
            B = -21.712318299469676
            selection = '({0}*{1} + {2})'.format(edep, A, B)
            """
            # use neutron calib
            # poly3 v2
            A = 2.5144270503067e-10
            B = -2.325285315009937e-05
            C = 0.8025672088693712
            D = 21.91843846381534
            selection = ('(({1}*{0}**3)+({2}*{0}**2)+({3}*{0})+{4})'
                         .format(edep, A, B, C, D))
            
            return selection

        else:
            return selection

    else:
        print('WARNING: no calibration for x{0} c{1} e{2}'.format(Cx, Chan, Engy))
        sys.exit()


def cutsBDT302(i, C, E):
    """
    https://cupwiki.ibs.re.kr/Kims/SET1EventSelection
    https://cupwiki.ibs.re.kr/Kims/EventSelection
    """
    #-------------------------------------------------------
    # This is for the V00-04-15 data
    # New lsveto single-hit data
    #-------------------------------------------------------

    i = int(i)
    Cx = i+1
    
    ### special cuts for the LS veto
    if Cx == 9:
        ### skip low energy and alpha energy for lsveto
        #if E in [0, 2]:
        #    return TCut('0')
        
        if C == 'S':
            # remember the ls single-hit gets the activeLS chain
            # basically none of the cuts work here
            energy = '(BLSVeto.Charge/143.8 > 0.0)'
            coinc  = '(1)'
            #muons  = '(1)'
            time   = '(1)'
            
            muon1 = ('(BTop.Charge   < 4000 && BBottom.Charge < 4000 && '
                      'BLeft.Charge  < 4000 && BRight.Charge  < 4000 && '
                      'BFront.Charge < 4000 && BRear.Charge   < 4000)')

            muon2 = ('((BTop.Charge + BBottom.Charge + '
                       'BLeft.Charge + BRight.Charge + '
                       'BFront.Charge + BRear.Charge) < 8000)')

            #muons = '(!'+muon2+')'
            muons = muon1+' && '+muon2

            # this means nothing in activeLS chain
            #coinc  = '(BLSVeto.isCoincident == 1)'

            # in activeLS branch "totalDeltaT0" always 0
            # actually muons and neutrons
            #muons  = '(BMuon.totalDeltaT0/1.e6 <= 30.)'
            
            #time   = '(BLSVeto.Time > 2400 && BLSVeto.Time < 2500)'
            #time   = '(BLSVeto.Time < 2820 || BLSVeto.Time > 2920)'
            #time   = '(BLSVeto.Time < 2400 || BLSVeto.Time > 2500)'
            
            return TCut('({0} && {1} && {2} && {3})'.format(energy, coinc, muons, time))
            
        elif C == 'M':
            energy = '(BLSVeto.Charge/143.8 > 0.0)'
            coinc  = '(BLSVeto.isCoincident == 1)'
            muon0  = '(BMuon.totalDeltaT0/1.e6 > 30.)'
            # try muon charge cut too?
            muon1  = ('(BTop.Charge   < 4000 && BBottom.Charge < 4000 && '
                       'BLeft.Charge  < 4000 && BRight.Charge  < 4000 && '
                       'BFront.Charge < 4000 && BRear.Charge   < 4000)')
            muon2  = ('((BTop.Charge + BBottom.Charge + '
                        'BLeft.Charge + BRight.Charge + '
                        'BFront.Charge + BRear.Charge) < 8000)')
            #muons = muon1+' && '+muon2
            muons = '('+muon0+' && '+muon1+' && '+muon2+')'
            
            # testing a time cut
            #time   = '(BLSVeto.Time >= 2800 && BLSVeto.Time <= 2950)'
            time   = '(1)'
            
            return TCut('({0} && {1} && {2} && {3})'.format(energy, coinc, muons, time))
            
        else:
            print('ERROR: invalid LS veto channel -->', C)
            sys.exit()


    ### same as Pushpa since 2017-12-19

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

    # only alphas for E=2
    if E == 2:
        alphaCut = ('((crystal{0}.energyD > 1000.) && ((pmt{0}1.rqtD1_5+pmt{0}2.rqtD1_5)/2. < {1}))'
                    .format(Cx, alpha[i]))
    else:
        alphaCut = ('( ! ((crystal{0}.energyD > 1000.) && ((pmt{0}1.rqtD1_5+pmt{0}2.rqtD1_5)/2. < {1})) )'
                    .format(Cx, alpha[i]))

        
    ### BDT cuts
    newBDT = [
        # C1
        # from 2019, looks like crap - https://cupwiki.ibs.re.kr/Kims/EventSelection
        #'((bdt[1]>-0.1 && (crystal1.energy + 20*bdt[1])>0) && ((crystal1.energy>=10 && bdtA[1]>-0.08) || (crystal1.energy>=8 && crystal1.energy<10 && bdtA[1]>-0.065) || (crystal1.energy>=7 && crystal1.energy<8 && bdtA[1]>0.005) || (crystal1.energy>=6 && crystal1.energy<7 && bdtA[1]>0.015) || (crystal1.energy>=5 && crystal1.energy<6 && bdtA[1]>0.01) || (crystal1.energy>=4 && crystal1.energy<5 && bdtA[1]>0.03) || (crystal1.energy>=3 && crystal1.energy<4 && bdtA[1]>0.035) || (crystal1.energy>=2 && crystal1.energy<3 && bdtA[1]>0.02) || (crystal1.energy<2 && bdtA[1]>0.02)))',
        # use C5,C8 style cut for C1
        #'(bdt[1]>-0.2 && bdtA[1]>-0.07)',
        '(bdt[1]>-0.18 && bdtA[1]>-0.07)',

        # C2
        #'(bdt[2]>0.0 && bdtA[2]>-0.05)',
        # new bdt cut for V00-04-14
        #'(((9.20842e-07*TMath::Exp(-104.504*bdt[2])+0.170872)-(9.70874*bdt[2])) < crystal2.energy)',
        # new set2 1keV BDT cut - https://cupwiki.ibs.re.kr/Kims/EventSelectionforSet2
        '((1.52e-06*TMath::Exp(-107.840*bdt[2]) - 1.94 - 29.288*bdt[2]) < crystal2.energy)',

        # C3
        #'(bdt[3]>0.0 && bdtA[3]>-0.07)',
        # new bdt cut for V00-04-14
        #'(((6.81804e-08*TMath::Exp(-101.856*bdt[3])+0.148344)-(4.04826*bdt[3])) < crystal3.energy)',
        # new set2 1keV BDT cut - https://cupwiki.ibs.re.kr/Kims/EventSelectionforSet2
        '((1.51e-09*TMath::Exp(-176.686*bdt[3]) + 0.34 - 5.469*bdt[3]) < crystal3.energy)',

        # C4
        #'(bdt[4]>-0.05 && (crystal4.energy + 40*bdt[4])>0 && bdtA[4]>-0.07)',
        # new bdt cut for V00-04-14
        #'(((1.37726e-07*TMath::Exp(-111.842*bdt[4])+0.446818)-(2.13498*bdt[4])) < crystal4.energy)',
        # new set2 1keV BDT cut - https://cupwiki.ibs.re.kr/Kims/EventSelectionforSet2
        '((1.86e-11*TMath::Exp(-137.413*bdt[4]) - 10.1 - 81.519*bdt[4]) < crystal4.energy)',

        # C5
        #'(bdt[5]>-0.2 && bdtA[5]>-0.07)',
        '(bdt[5]>-0.22 && bdtA[5]>-0.1)',

        
        # C6
        #'(bdt[6]>0.0 && (crystal6.energy + 20*bdt[6])>2.0 && bdtA[6]>-0.07)',
        # new bdt cut for V00-04-14
        #'(((0.000315623*TMath::Exp(-75.3998*bdt[6])-0.313482)-(13.6302*bdt[6])) < crystal6.energy)',
        # new set2 1keV BDT cut - https://cupwiki.ibs.re.kr/Kims/EventSelectionforSet2
        '((3.30e-07*TMath::Exp(-124.845*bdt[6]) - 2.74 - 52.57*bdt[6]) < crystal6.energy)',

        # C7
        #'(bdt[7]>0.0 && (crystal7.energy + 20*bdt[7])>2.0 && bdtA[7]>-0.07)',
        # new bdt cut for V00-04-14
        #'(((1.45552e-05*TMath::Exp(-88.7196*bdt[7])+0.566336)-(7.57773*bdt[7])) < crystal7.energy)',
        # new set2 1keV BDT cut - https://cupwiki.ibs.re.kr/Kims/EventSelectionforSet2
        '((6.01e-08*TMath::Exp(-120.177*bdt[7]) - 0.24 - 31.615*bdt[7]) < crystal7.energy)',

        # C8
        #'(bdt[8]>-0.2 && bdtA[8]>-0.07)',
        '(bdt[8]>-0.40 && bdtA[8]>-0.1)',

        
    ] # end of newBDT list
    
    bdtCut = newBDT[i]

    
    ### global noise/muon cuts
    #===========================================================================
    coinc  = '(BLSVeto.isCoincident == 1)'
    muons  = '(BMuon.totalDeltaT0/1.e6 > 30)'

    charge = '(crystal'+str(Cx)+'.rqcn > -1)'
    
    #nc     = '(pmt'+str(i+1)+'1.nc > 1 && pmt'+str(i+1)+'2.nc > 1)'
    ### new cut for v00-04-14 from Govinda
    nc     = '(pmt'+str(Cx)+'1.nc > 0 && pmt'+str(Cx)+'2.nc > 0)'
    
    #t1     = '(pmt'+str(i+1)+'1.t1 > 2 && pmt'+str(i+1)+'2.t1 > 2)'
    ### new cut for v00-04-14 from Govinda
    t1     = '(pmt'+str(Cx)+'1.t1 > 0 && pmt'+str(Cx)+'2.t1 > 0)'
    
    #noiseCut = '('+coinc+' && '+muons+' && '+charge+' && '+nc+' && '+t1+')'
    #===========================================================================

    # build single/multi hit cut
    if C == 'S':
        lsveto = '(BLSVeto.Charge/143.8 < 80)'
        #lsveto = '(BLSVeto.Charge/143.8 < 50)'
        hits = '('
        #bdtOther = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc < 4) && '
                #hits += '(crystal'+str(j+1)+'.nc < 3) && '
                #bdtOther += '(!'+newBDT[i]+') && '
        ### remove extra '&&'
        hits = hits[:-4]+')'
        #bdtOther = bdtOther[:-4]+')'
        hitCut = '('+lsveto+' && '+hits+')'
        
    elif C == 'M':
        lsveto = '(BLSVeto.Charge/143.8 >= 80)'
        #lsveto = '(BLSVeto.Charge/143.8 >= 50)'
        hits = '('
        #bdtOther = '('
        for j in range(8):
            if j != i:
                hits += '(crystal'+str(j+1)+'.nc >= 4) || '
                #hits += '(crystal'+str(j+1)+'.nc >= 5) || '
                #bdtOther += '('+newBDT[i]+') || '
        ### remove extra '||'
        hits = hits[:-4]+')'
        #bdtOther = bdtOther[:-4]+')'
        hitCut = '('+lsveto+' || '+hits+')'
        
    else:
        print('ERROR: invalid crystal channel -->', C)
        sys.exit()
    
    masterCut = alphaCut+' && '+coinc+' && '+muons+' && '+hitCut
    
    # don't use bdt cut on C5 and C8
    if E == 0 and Cx not in [5, 8]:
        masterCut += ' && '+bdtCut
        #masterCut += ' && '+bdtCut+' && '+bdtOther
        
    # don't use noise cuts on C5 and C8
    if Cx not in [5, 8]:
        masterCut += ' && '+charge+' && '+nc+' && '+t1
    
    return TCut(masterCut)


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
        if bline.startswith('stop'):
            break
        
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
    # use histparams
    par = histparams(E)
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
    return deepcopy(total)


def makeResid100(chan, E, par):
    resid = []
    # use histparams
    par = histparams(E)
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
    return deepcopy(resid)


### this is old - keeping for reference
def processSurf500(chain, mc, info, C, E, fromx, Q=False):

    # expo weight the depth or use depth as a cut
    EXPO_WEIGHT = 1
    
    # get info from mc key
    bits = info['key'].split('-')
    xstal = bits[0][1]
    cx = int(xstal)-1
    fx = int(fromx)-1
    location = bits[1]
    isotope = bits[2]
    
    
    # DEFINE HIST AND KEY
    #-----------------------------------------------------------------------
    key = info['key']
    key += '-f'+str(fromx)
    if E==2 and Q:
        key += '-q'+str(Q)
    key += '-c'+str(C)
    key += '-e'+str(E)
        
    if key not in mc:
        mc[key] = {}

    if 'info' not in mc[key]:
        mc[key]['info'] = copy(info)

    pars = histparams(E)
    mc[key]['pars'] = pars
    #histo = TH1F(key, longNames(cx), pars[0], pars[1], pars[2])
    histo = TH1F(key, key, pars[0], pars[1], pars[2])
    
    key2 = key+'_generated'
    temp2 = TH1F(key2, 'generated', 1, 0, 1)
    #-----------------------------------------------------------------------
    
    
    ### try getting the expo depth from file name?
    ### fix the depth as par[0]
    depth = float(info['floca'].split('-')[3])
    func = TF1('func', surfProfile, 0, 10, 1)
    # depth/100 to convert my number to um (88 == 0.88 um)
    func.SetParameter(0, depth/100.)

    entries = chain.GetEntries()
    #print entries

    for jentry in range(entries):
        #print(jentry)
        chain.GetEntry(jentry)
        
        ### pull entry data our of leaf
        #--------------------------------------------------
        xx = chain.GetLeaf('primX0').GetValue()
        yy = chain.GetLeaf('primY0').GetValue()
        zz = chain.GetLeaf('primZ0').GetValue()

        # get all crystal/lsveto energies
        if E: evar = 'edepResolD'
        else: evar = 'edepResolA'
        edepRes = []
        for i in range(0,9):
            edepRes.append(1000.*(chain.GetLeaf(evar).GetValue(i)))
        #elow = chain.GetLeaf('edepResolA').GetValue(cx)
        #ehi = chain.GetLeaf('edepResolD').GetValue(cx)

        edep = 1000.*(chain.GetLeaf('edep').GetValue(cx))
        
        singleHit = chain.GetLeaf('singleHitTag').GetValue(cx)
        multiHit = chain.GetLeaf('multipleHitTag').GetValue(cx)
        evType = chain.GetLeaf('evt_Type').GetValue()
        group = chain.GetLeaf('groupNo').GetValue()
        #decayType = chain.GetLeaf('primDecayType').GetValue()
        # Took me forever to figure this out!
        decayType = chain.primDecayType
        #print (decayType)

        
        ### volume cuts
        #--------------------------------------------------
        X0, Y0, Z0, rad, height, _dep = mcDimensions(fx)
        dist = np.sqrt((X0-xx)**2 + (Y0-yy)**2)
        zdist = abs(Z0-zz)

        # Teflon thickness in mm - 2020-04-29
        #teflon = 1.015
        # I think this changed in gyunho's new sim
        teflon = 1.016

        ### NaI depth to cut on in mm
        if EXPO_WEIGHT:
            # my depth /100 to convert to um (100 == 1.0um)
            # /1000 to convert to mm (1.0um == 0.001mm)
            # *10 to go 10 expo depths in (0.01mm)
            dep = (((depth/100.)/1000.)*10)
            #dep = (((depth/100.)/1000.)*20)
            # just use 10um depth by default
            #dep = 0.01  # use 0.01mm for the 10um generated sim
        else:
            dep = (depth/100.)/1000.
            #dep = 0.01
            
        # for NaI surf side
        if info['loca'].endswith('side') or info['loca'].endswith('expo'):
            #print(dep, dist, rad)
            if dist > rad or dist < rad-dep: continue
            if zdist > height: continue

        # for NaI surf face
        elif info['loca'].endswith('face'):
            if zdist > height or zdist < height-dep: continue
            if dist > rad: continue

        # teflon in-side - 2020-05-01
        elif info['loca'].endswith('in'):
            if dist < rad or dist > rad+dep: continue
            if zdist > height: continue

        # teflon out-side - 2020-05-01
        elif info['loca'].endswith('out'):
            if dist < rad+teflon-dep or dist > rad+teflon: continue
            if zdist > height: continue

        else:
            print('ERROR: I do not know what to do with', info['loca'])
            sys.exit()


        ### get the expo weighting
        #--------------------------------------------------
        if EXPO_WEIGHT:
            # for NaI surf side
            if info['loca'].endswith('side') or info['loca'].endswith('expo'):
                # * 1000 to convert mm to um
                wtf = func.Eval(1000.*(rad-dist))
            # for NaI surf face
            elif info['loca'].endswith('face'):
                # * 1000 to convert mm to um
                wtf = func.Eval(1000.*(height-zdist))
            # teflon in-side - 2020-04-29
            elif info['loca'].endswith('in'):
                # * 1000 to convert mm to um
                wtf = func.Eval(1000.*(dist-rad))
            # teflon out-side - 2020-04-29
            elif info['loca'].endswith('out'):
                # * 1000 to convert mm to um
                wtf = func.Eval(1000.*(rad+teflon-dist))
            else:
                print('ERROR: I do not know what to do with', info['loca'])
                sys.exit()
        else:
            wtf = 1.0
        
        #print(rad-dist, wtf)
        
        ### now it's a valid generated event
        #--------------------------------------------------
        if evType <= 10:
            temp2.Fill(0, wtf)
        
        
        ### event type cut
        #--------------------------------------------------
        if evType <= 10: continue
        
        
        ### alpha cuts?
        #--------------------------------------------------
        #print(decayType)
        if E == 2:
            if decayType != 'alpha':
                continue
        else:
            if decayType == 'alpha':
                continue

        
        ### energy cut
        #--------------------------------------------------
        threshold = 0.125 # in keV
        if edepRes[cx] < threshold: continue

        # special handling of lsveto
        if cx == 8 and C == 'S':
            if not (edepRes[0] < threshold and
                    edepRes[1] < threshold and
                    edepRes[2] < threshold and
                    edepRes[3] < threshold and
                    edepRes[4] < threshold and
                    edepRes[5] < threshold and
                    edepRes[6] < threshold and
                    edepRes[7] < threshold):
                continue

        if cx == 8 and C == 'M':
            if not (edepRes[0] >= threshold or
                    edepRes[1] >= threshold or
                    edepRes[2] >= threshold or
                    edepRes[3] >= threshold or
                    edepRes[4] >= threshold or
                    edepRes[5] >= threshold or
                    edepRes[6] >= threshold or
                    edepRes[7] >= threshold):
                continue
        
        
        ### single/multi hit cut
        #--------------------------------------------------
        if cx == 8:
            if C == 'S':
                if not (singleHit == 1 and multiHit == -1):
                    continue
            if C == 'M':
                if not (singleHit == -1 and multiHit == 1):
                    continue
        else:
            if C == 'S':
                if not ((singleHit == 1 and multiHit == -1) and edepRes[8] < 80.0):
                    continue
            if C == 'M':
                if not ((singleHit == -1 and multiHit == 1) or edepRes[8] > 80.0):
                    continue
        
        
        ### do I need a groupNo cut?
        #--------------------------------------------------
        # -- insert here...
        # not needed now for just surface Pb210


        ### fill the weighted histograms
        #--------------------------------------------------
        # quenching for crystals and E=2
        
        if cx < 8 and E == 2 and Q:
            _A, _B, _C = quenchPars(cx, Q)
            quenched = _A*(edep**2) + (_B*edep) + _C
            if quenched < 0: continue
            #p0, p1 = resolPars(cx, E)
            #sigma = (p0*np.sqrt(quenched)) + (p1*quenched)
            # manually setting alpha sigma to 40
            sigma = 50.
            energy = random.gauss(quenched, sigma)
            #print(jentry, 'filling quenched histo')
            histo.Fill(energy, wtf)
            
        else:
            #print(jentry, 'filling edepRes histo')
            histo.Fill(edepRes[cx], wtf)
            
    #=====================================================================
    
    
    ### create the MC histogram
    #---------------------------------------------------------------------
    # rescale the number of entries to weighted integral
    histo.SetEntries(histo.Integral())
    if 'hist' in mc[key]:
        mc[key]['hist'].Add(histo)
    else:
        mc[key]['hist'] = copy(histo)
    #---------------------------------------------------------------------
    
    ### create a hist of the generated events
    #---------------------------------------------------------------------
    # rescale the number of entries to weighted integral
    temp2.SetEntries(temp2.Integral())
    if 'generated_hist' in mc[key]:
        mc[key]['generated_hist'].Add(temp2)
    else:
        mc[key]['generated_hist'] = copy(temp2)
    
    generated = mc[key]['generated_hist'].GetEntries()
    if generated <= 0:
        print('WARNING: no events generated for --> {0}'.format(key))
    mc[key]['generated'] = generated
    #---------------------------------------------------------------------

    ### print out efficiency numbers
    #---------------------------------------------------------------------
    this_detected = histo.GetEntries()
    this_generated = temp2.GetEntries()
    print(key, 'gen =', round(this_generated,1), 'det =', round(this_detected,1))
    """
    this_detected = histo.GetEntries()
    try: this_eff = round(100*this_detected/this_generated, 2)
    except: this_eff = 0.0
    #try: tot_eff = round(100*detected/generated, 2)
    #except: tot_eff = 0.0
    print('DEBUG: [{3}] This det/gen = eff --> {0} / {1} = {2}%'
          .format(this_detected, this_generated, this_eff, key))
    #print('DEBUG: [{3}] Total det/gen = eff --> {0} / {1} = {2}%'
    #      .format(detected, generated, tot_eff, key))
    """
    #---------------------------------------------------------------------

    del histo
    del temp2
    
    return mc


### this is old - keeping for reference
def processSurf501(chain, mc, info, C, E, fromx, Q=False):

    # expo weight the depth or use depth as a stict cut
    EXPO_WEIGHT = 1
    
    # get info from mc key
    bits = info['key'].split('-')
    xstal = bits[0][1]
    cx = int(xstal)-1
    fx = int(fromx)-1
    location = bits[1]
    isotope = bits[2]
    
    
    # DEFINE HIST AND KEY
    #-----------------------------------------------------------------------
    key = info['key']
    key += '-f'+str(fromx)
    if E==2 and Q:
        key += '-q'+str(Q)
    key += '-c'+str(C)
    key += '-e'+str(E)

    print('INFO: processSurf501 \"{0}\"'.format(key))
    
    if key not in mc:
        mc[key] = {}

    if 'info' not in mc[key]:
        mc[key]['info'] = copy(info)

    pars = histparams(E)
    mc[key]['pars'] = pars
    #histo = TH1F(key, longNames(cx), pars[0], pars[1], pars[2])
    histo = TH1F(key, key, pars[0], pars[1], pars[2])
    
    key2 = key+'_generated'
    temp2 = TH1F(key2, 'generated', 1, 0, 1)
    #-----------------------------------------------------------------------
    
    
    ### get the expo depth from location
    ### fix the depth as par[0]
    ### depth from location is in nm (nano-meter)
    _depth = float(info['floca'].split('-')[3])
    depth = _depth*1.0e-6 # convert nm to mm
    fmin = 0.0 # in mm
    #fmax = 0.1 # in mm
    fmax = 1.0 # in mm
    func = TF1('func', surfProfile, fmin, fmax, 1)
    func.SetParameter(0, depth)

    entries = chain.GetEntries()
    #print entries

    for jentry in range(entries):
        #print(jentry)
        chain.GetEntry(jentry)
        
        ### pull entry data our of leaf
        #--------------------------------------------------
        xx = chain.GetLeaf('primX0').GetValue()
        yy = chain.GetLeaf('primY0').GetValue()
        zz = chain.GetLeaf('primZ0').GetValue()

        # get all crystal/lsveto energies
        if E: evar = 'edepResolD'
        else: evar = 'edepResolA'
        edepRes = []
        for i in range(0,9):
            edepRes.append(1000.*(chain.GetLeaf(evar).GetValue(i)))
        #elow = chain.GetLeaf('edepResolA').GetValue(cx)
        #ehi = chain.GetLeaf('edepResolD').GetValue(cx)

        edep = 1000.*(chain.GetLeaf('edep').GetValue(cx))
        
        singleHit = chain.GetLeaf('singleHitTag').GetValue(cx)
        multiHit = chain.GetLeaf('multipleHitTag').GetValue(cx)
        evType = chain.GetLeaf('evt_Type').GetValue()
        groupNo = chain.GetLeaf('groupNo').GetValue()
        volName = chain.primVolumeName
        decayType = chain.primDecayType
        particle = chain.primParticleName
        #print (decayType)

        
        ### volume cuts
        #--------------------------------------------------
        ### dimensions here in mm
        X0, Y0, Z0, rad, height, _dep = mcDimensions(fx)
        dist = np.sqrt((X0-xx)**2 + (Y0-yy)**2)
        zdist = abs(Z0-zz)

        # Teflon thickness in mm - 2020-04-29
        #teflon = 1.015
        # I think this changed in gyunho's new sim
        teflon = 1.016  # in mm

        ### NaI depth to cut on in mm
        if EXPO_WEIGHT:
            # *10 to go 10 expo depths in
            dep = depth*10.
        else:
            dep = depth

        # for NaI surf side
        if info['loca'].endswith('side') or info['loca'].endswith('expo'):
            #print(dep, dist, rad)
            if dist > rad or dist < rad-dep: continue
            if zdist > height: continue

        # for NaI surf face
        elif info['loca'].endswith('face'):
            if zdist > height or zdist < height-dep: continue
            if dist > rad: continue

        # teflon in-side - 2020-05-01
        elif info['loca'].endswith('in'):
            if dist < rad or dist > rad+dep: continue
            if zdist > height: continue

        # teflon out-side - 2020-05-01
        elif info['loca'].endswith('out'):
            if dist < rad+teflon-dep or dist > rad+teflon: continue
            if zdist > height: continue

        else:
            print('ERROR: I do not know what to do with', info['loca'])
            sys.exit()


        ### get the expo weighting
        #--------------------------------------------------
        if EXPO_WEIGHT:
            # for NaI surf side
            if info['loca'].endswith('side') or info['loca'].endswith('expo'):
                wtf = func.Eval(rad-dist)
                
            # for NaI surf face
            elif info['loca'].endswith('face'):
                wtf = func.Eval(height-zdist)
                
            # teflon in-side - 2020-04-29
            elif info['loca'].endswith('in'):
                wtf = func.Eval(dist-rad)
                
            # teflon out-side - 2020-04-29
            elif info['loca'].endswith('out'):
                wtf = func.Eval(rad+teflon-dist)
                
            else:
                print('ERROR: I do not know what to do with', info['loca'])
                sys.exit()
        else:
            wtf = 1.0
        
        #print(rad-dist, wtf)
        
        ### now it's a valid generated event
        #--------------------------------------------------
        if evType <= 10:
            temp2.Fill(0, wtf)
            continue
        
        
        ### group number cut
        #--------------------------------------------------
        gstart, gstop = groupNum500(isotope, True)
        if not (groupNo >= gstart and groupNo < gstop): continue
                
        
        ### alpha cuts
        #--------------------------------------------------
        # alpha cuts for crystals
        if cx < 8:
            # alpha cuts only on internal/surface sources
            if isAlphaMC(info):
                if E == 2:
                    if decayType != 'alpha': continue
                elif E == 1:
                    if decayType == 'alpha': continue
            else:
                # skip other events from contributing to the crystal alpha spectrum
                if E == 2: continue
        
        # show alphas in lsveto E=2 (just to see how it looks)
        if cx == 8 and E == 2:
            if decayType != 'alpha': continue
        
        
        ### energy cut
        #--------------------------------------------------
        threshold = 0.125  # in keV
        #if edepRes[cx] <= 0.0: continue
        if edepRes[cx] < threshold: continue

        # special handling of lsveto
        if cx == 8 and C == 'S':
            if not (edepRes[0] < threshold and
                    edepRes[1] < threshold and
                    edepRes[2] < threshold and
                    edepRes[3] < threshold and
                    edepRes[4] < threshold and
                    edepRes[5] < threshold and
                    edepRes[6] < threshold and
                    edepRes[7] < threshold):
                continue

        if cx == 8 and C == 'M':
            if not (edepRes[0] >= threshold or
                    edepRes[1] >= threshold or
                    edepRes[2] >= threshold or
                    edepRes[3] >= threshold or
                    edepRes[4] >= threshold or
                    edepRes[5] >= threshold or
                    edepRes[6] >= threshold or
                    edepRes[7] >= threshold):
                continue
        
        
        ### single/multi hit cut
        #--------------------------------------------------
        if cx == 8:
            if C == 'S':
                if not (singleHit == 1 and multiHit == -1):
                    continue
            if C == 'M':
                if not (singleHit == -1 and multiHit == 1):
                    continue
        else:
            # LS veto threshold
            ls_thresh = 80.0  # in keV
            if C == 'S':
                if not ((singleHit == 1 and multiHit == -1) and edepRes[8] < ls_thresh):
                    continue
            if C == 'M':
                if not ((singleHit == -1 and multiHit == 1) or edepRes[8] >= ls_thresh):
                    continue
        
        

        ### fill the weighted histograms
        #--------------------------------------------------
        if cx < 8:
            if E == 2 and Q:
                # quenching for crystals and E=2
                #_A, _B, _C = quenchPars(cx, Q)
                #quenched = _A*(edep**2) + (_B*edep) + _C
                quenched = calcQuench(cx, Q, edep)
                if quenched < 0: continue
                #p0, p1, p2 = resolPars(cx, E)
                #sigma = (p0*np.sqrt(quenched)) + (p1*quenched)
                sigma = calcSigma(cx, E, quenched)
                #sigma = 50.
                energy = random.gauss(quenched, sigma)
                #print(jentry, 'filling quenched histo')
                histo.Fill(energy, wtf)
            
            elif cx in [4, 7]:
                # tweaked resolutions for C5 and C8
                #p0, p1, p2 = resolPars(cx, E)
                #sigma = (p0*np.sqrt(edep)) + (p1*edep)
                sigma = calcSigma(cx, E, edep)
                energy = random.gauss(edep, sigma)
                histo.Fill(energy, wtf)
            
            else:
                # need to quench pile-up beta-alpha events
                # assume quenching q1 or q2?
                _Q = 2
                # special treatment for pile-up events (Rn222)
                if particle == 'Po214' and edep > 6000:
                    #alphaE = 7833. # Qval from nudat
                    alphaE = 7686. # alpha energy from nudat
                    betaE = 0.
                    if edep > alphaE:
                        betaE = edep-alphaE
                    else:
                        alphaE = edep
                    # quench and smear the alpha
                    #_A, _B, _C = quenchPars(cx, _Q)
                    #quenched = _A*(alphaE**2) + (_B*alphaE) + _C
                    quenched = calcQuench(cx, _Q, alphaE)
                    #p0, p1, p2 = resolPars(cx, E)
                    #sigmaA = (p0*np.sqrt(quenched)) + (p1*quenched)
                    sigmaA = calcSigma(cx, E, quenched)
                    #sigmaA = 50.
                    alphaE = random.gauss(quenched, sigmaA)
                    #histo.Fill(alphaE)
                    if betaE > 0:
                        # just smear the beta
                        #p0, p1, p2 = resolPars(cx, E)
                        #sigmaB = (p0*np.sqrt(betaE)) + (p1*betaE)
                        sigmaB = calcSigma(cx, E, betaE)
                        betaE = random.gauss(betaE, sigmaB)
                        #histo.Fill(betaE)
                    histo.Fill(alphaE+betaE, wtf)
                    
                # special treatment for pile-up events (Th228)
                elif particle == 'Po212' and edep > 6000:
                    #alphaE = 8954. # from MC?
                    alphaE = 8785. # alpha energy from nudat
                    betaE = 0.
                    if edep > alphaE:
                        betaE = edep-alphaE
                    else:
                        alphaE = edep
                    # quench and smear the alpha
                    #_A, _B, _C = quenchPars(cx, _Q)
                    #quenched = _A*(alphaE**2) + (_B*alphaE) + _C
                    quenched = calcQuench(cx, _Q, alphaE)
                    #p0, p1, p2 = resolPars(cx, E)
                    #sigmaA = (p0*np.sqrt(quenched)) + (p1*quenched)
                    sigmaA = calcSigma(cx, E, quenched)
                    #sigmaA = 50.
                    alphaE = random.gauss(quenched, sigmaA)
                    #histo.Fill(alphaE)
                    if betaE > 0:
                        # just smear the beta
                        #p0, p1, p2 = resolPars(cx, E)
                        #sigmaB = (p0*np.sqrt(betaE)) + (p1*betaE)
                        sigmaB = calcSigma(cx, E, betaE)
                        betaE = random.gauss(betaE, sigmaB)
                        #histo.Fill(betaE)
                    histo.Fill(alphaE+betaE, wtf)
                    
                else:
                    histo.Fill(edepRes[cx], wtf)

        else:
            histo.Fill(edepRes[cx], wtf)
            
    #=====================================================================
    
    
    ### create the MC histogram
    #---------------------------------------------------------------------
    # rescale the number of entries to weighted integral
    histo.SetEntries(histo.Integral())
    if 'hist' in mc[key]:
        mc[key]['hist'].Add(histo)
    else:
        mc[key]['hist'] = copy(histo)
    #---------------------------------------------------------------------
    
    ### create a hist of the generated events
    #---------------------------------------------------------------------
    # rescale the number of entries to weighted integral
    temp2.SetEntries(temp2.Integral())
    if 'generated_hist' in mc[key]:
        mc[key]['generated_hist'].Add(temp2)
    else:
        mc[key]['generated_hist'] = copy(temp2)
    
    generated = mc[key]['generated_hist'].GetEntries()
    if generated <= 0:
        print('WARNING: no events generated for --> {0}'.format(key))
    mc[key]['generated'] = generated
    #---------------------------------------------------------------------

    ### print out efficiency numbers
    #---------------------------------------------------------------------
    #this_detected = histo.GetEntries()
    #this_generated = temp2.GetEntries()
    #print(key, 'gen =', round(this_generated,1), 'det =', round(this_detected,1))
    """
    this_detected = histo.GetEntries()
    try: this_eff = round(100*this_detected/this_generated, 2)
    except: this_eff = 0.0
    #try: tot_eff = round(100*detected/generated, 2)
    #except: tot_eff = 0.0
    print('DEBUG: [{3}] This det/gen = eff --> {0} / {1} = {2}%'
          .format(this_detected, this_generated, this_eff, key))
    #print('DEBUG: [{3}] Total det/gen = eff --> {0} / {1} = {2}%'
    #      .format(detected, generated, tot_eff, key))
    """
    #---------------------------------------------------------------------

    del histo
    del temp2
    
    return mc


