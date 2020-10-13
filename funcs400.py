#!/usr/bin/env python
######################################################################
# funcs400.py
# 
# version: 2020-05-18
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + ScaleEnergy400() and MCBinShift400()

# + add alpha cut to surface teflon
# + add exceptions for nai-surf, teflon, and reflector in build400()
# + create an exception for nai-surf-expo in build400()
# ~ changed C4 Z length
# ~ made the expo surface depth profile function dynamic
# ~ tweak volume cuts for positional degeneracy in buildMCSurf400()
# ~ fixed a from-crystal bug in buildMCSurf400()
# ~ build400() has option to use new buildMCSurf400()
# + added the buildMCSurf400() for special case expo profile
# - remove returning the "runtime" for data - get it in main script
# + add new getInfo400(), build400(), and buildData400()
# 
# email: mkauer@physics.wisc.edu
######################################################################

import os,sys,re
from copy import deepcopy
import numpy as np

from ROOT import *
import ROOT

sys.path.append("/home/mkauer/COSINE/CUP/mc-fitting/")
sys.path.append("/home/mkauer/mc-fitting/")
from funcs300 import *


def getInfo400(line, freuse=0, fchans='SM', fxstals=[]):

    here = baseDir()
    
    info = {}
    
    info['line'] = line
    
    bits = line.split()
    ### probably don't need to use 're' since split() is working fine
    #import re
    #bits = filter(None, re.split("[ \s\t]+", line.strip()))
    
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
    # channels to fit
    info['chans'] = fchans
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # crystal number
    info['xstl'] = int(bits[1])
    # force only a particular crystal(s) and skip the others
    if len(fxstals) > 0 and info['xstl'] not in fxstals:
        return []
    # what crystal is the background from?
    # set to xstl by default - can be updated later
    info['from'] = int(bits[1])
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # reuse rootfile name
    if 'root' in str(bits[-1]):
        info['rootfile'] = os.path.join(here, str(bits[-1]))
        #info['rootfile'] = str(bits[-1])
    else:
        info['rootfile'] = 0
    #-----------------------------------------------------------------
    

    if 'D' in info['type']:
        
        ### data set
        info['tag'] = str(bits[2])
        
        ### processing version
        info['build'] = str(bits[3])
        
        ### data file
        info['file'] = os.path.join(here, str(bits[4]))
        
        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-data'
        key += '-'+str(info['tag'])
        key += '-'+str(info['build']).replace('-','')
        info['key'] = key
                
    else:
        
        ### background location help
        info['floca'] = str(bits[2])
        # without any '-' for hist key help
        info['loca'] = str(bits[2]).replace('-','')
        
        ### top level isotope file name
        info['isof'] = str(bits[3])

        ### isotope chain break start
        info['chst'] = str(bits[4])

        ### isotope chain break stop
        info['chsp'] = str(bits[5])

        ### activity
        info['acti'] = float(bits[6])

        ### fit bounds
        info['fbnd'] = [float(bits[7]), float(bits[8])]

        ### simulation version
        info['sim'] = str(bits[9])
        
        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-'+info['loca']
        key += '-'+info['chst']+'_'+info['chsp']
        info['key'] = key

        ### assign a plotting group to the isotope-location
        info['group'] = setGroup(info)
        
    return info


def build400(infile='backgrounds400.txt', others=1, freuse=0, fchans='SM', fxstals=[]):

    debug = 0
    
    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):
        info = getInfo400(line, freuse, fchans, fxstals)
        #print info
        
        if len(info) == 0:
            continue
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Do not include "others" in processing [1] - Default [0]
        exception = 0
        # turn these off after testing
        if 'naisurf' in info['key']: exception = 1
        if 'teflon' in info['key']: exception = 1
        if 'teflonsurf' in info['key']: exception = 1
        if 'reflector' in info['key']: exception = 1
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # These are mandatory excludes because they are always global
        exclude = bool('lsveto' in info['key']
                    or 'plastic' in info['key']
                    or 'innersteel' in info['key']
                    or 'steel' in info['key']
                    or 'gamma' in info['key'])
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        infos = []
        
        if 'data' in info['key']:
            if debug: print 'DEBUG: if data in info[key]: -->', info['key']
            infos.append(info)
        elif exception:
            if debug: print 'DEBUG: elif exception: -->', info['key']
            key = info['key']
            newinfo = deepcopy(info)
            newinfo['key'] = key+'-f'+str(info['xstl'])
            newinfo['from'] = info['xstl']
            infos.append(newinfo)
        elif info['xstl']==9 and not exclude:
            if debug: print 'DEBUG: elif info[xstl]==9 and not exclude: -->', info['key']
            key = info['key']
            for i in range(8):
                #info['key'] = key+'-f'+str(i+1)
                newinfo = deepcopy(info)
                newinfo['key'] = key+'-f'+str(i+1)
                newinfo['from'] = str(i+1)
                infos.append(newinfo)
        elif others and not exclude:
            if debug: print 'DEBUG: elif others and not exclude: -->', info['key']
            key = info['key']
            for i in range(8):
                #info['key'] = key+'-f'+str(i+1)
                newinfo = deepcopy(info)
                newinfo['key'] = key+'-f'+str(i+1)
                newinfo['from'] = str(i+1)
                infos.append(newinfo)
        elif not others and not exclude:
            if debug: print 'DEBUG: elif not others and not exclude: -->', info['key']
            key = info['key']
            newinfo = deepcopy(info)
            newinfo['key'] = key+'-f'+str(info['xstl'])
            newinfo['from'] = info['xstl']
            infos.append(newinfo)
        elif exclude:
            if debug: print 'DEBUG: elif exclude: -->', info['key']
            infos.append(info)
        else:
            print 'WARNING:', info['key'], 'does not fit any known criteria'
            print '         Please check out build400() '
            infos.append(info)
            
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        for info in infos:
            #print info['key']
            if 'D' in info['type']:
                data = buildData400(info, data)
                
            elif 'B' in info['type']:
                #if info['loca'] == 'naisurf10umexpo':
                #if 'naisurf10umexpo' in info['loca']:
                if 'naisurfexpo' in info['loca'] \
                   or 'teflonsurfexpo' in info['loca']:
                    #print "!!!   SPECIAL PB210 PROFILE WORKING   !!!"
                    bkgs = buildMCSurf400(info, bkgs)
                else:
                    bkgs = buildMC300(info, bkgs)
            
            elif 'F' in info['type']:
                #if info['loca'] == 'naisurf10umexpo':
                #if 'naisurf10umexpo' in info['loca']:
                if 'naisurfexpo' in info['loca'] \
                   or 'teflonsurfexpo' in info['loca']:
                    #print "!!!   SPECIAL PB210 PROFILE WORKING   !!!"
                    sigs = buildMCSurf400(info, sigs)
                else:
                    sigs = buildMC300(info, sigs)
            else:
                print 'WARNING: I do not know type',info['type'],'in line:'
                print info['line']
                continue
            
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    return data, bkgs, sigs


def buildData400(info, data):
    
    if info['reuse']:
        for c in info['chans']:
            info['chan'] = c
            
            if not info['rootfile']:
                print 'ERROR: no rootfile specified in backgrounds file'
                sys.exit()

            rootfile = info['rootfile']
            if not os.path.exists(rootfile):
                print 'ERROR: rootfile not found -->', rootfile
                sys.exit()

            rfile = TFile(rootfile, "READ")

            for e in range(2):

                key  = info['key']
                key += '-c'+str(c)
                key += '-e'+str(e)
                #print key

                # exclude lsveto low energy
                #if 'x9' in key and e == 0:
                #    continue
                
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
                except:
                    print "ERROR: could not find runtime_hist -->",key+'_runtime'
                    sys.exit()

    else:
        
        runfile = info['file']
        #print 'looking for -->',runfile
        if not os.path.exists(runfile):
            print 'ERROR: runfile not found -->', runfile
            sys.exit()

        runtime = 0.
        nfiles = 0
        chain = TChain("ntp","")
        for line in readFile(runfile):
            bits = line.split()
            build_filename = info['build']+'/mrgd_M'+bits[0].zfill(6)+'.root.'+bits[1].zfill(3)
            if not onCup():
                fpath = '/home/mkauer/COSINE/CUP/mc-fitting/data/COSINE/MRGD/phys/'+build_filename
            else:
                fpath = '/data/COSINE/MRGD/phys/'+build_filename

            #print 'INFO: looking for data file -->', fpath
            
            if not os.path.exists(fpath):
                #print 'WARNING: data file not found for -->', fpath
                continue
            else:
                #print 'INFO: adding file -->', fpath
                tmpruntime = -1
                tmpruntime = getDuration(fpath)
                #print tmpruntime, fpath
                if tmpruntime > 0:
                    runtime += tmpruntime
                    #print runtime
                    nfiles += chain.Add(fpath)
                else:
                    print 'WARNING: no run duration found for -->', fpath
        
        if nfiles > 0:
            print 'INFO:',nfiles,'files found for data',info['tag'],'build',info['build']
            for C in info['chans']:
                info['chan'] = C
                for E in range(2):

                    # DEFINE HIST AND KEY
                    #-----------------------------------------------------------------------
                    i = info['xstl'] - 1
                    
                    key  = info['key']
                    key += '-c'+str(C)
                    key += '-e'+str(E)
                    
                    data[key] = {}
                    data[key]['info'] = info
                    
                    pars = histparam64(E)
                    data[key]['pars'] = pars
                    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
                    #-----------------------------------------------------------------------
                    
                    
                    # CALIBRATION
                    #-----------------------------------------------------------------------
                    if info['build'] == 'V00-04-04':
                        ### calib100 uses Pushpa's tweaked low energy calibrations
                        #edep, selection = calib100(i, E)
                        ### calib301 uses the standard production calibrations
                        ### and uses my tweaked lsveto calibration
                        ### testing differences between 04-04, 04-12, and 04-14
                        edep, selection = calib301(i, E)
                        
                    elif info['build'] == 'V00-04-12':
                        ### calib300 uses Govinda's tweaked high energy calibrations
                        ### and uses my tweaked lsveto calibration
                        #edep, selection = calib300(i, E)
                        ### calib301 uses the standard production calibrations
                        ### and uses my tweaked lsveto calibration
                        ### testing differences between 04-04, 04-12, and 04-14
                        edep, selection = calib301(i, E)
                        
                    elif info['build'] == 'V00-04-14':
                        ### calib301 uses the standard production calibrations
                        ### and uses my tweaked lsveto calibration
                        ### testing differences between 04-04, 04-12, and 04-14
                        edep, selection = calib301(i, E)
                        
                    elif info['build'] == 'V00-04-15':
                        ### calib301 uses the standard production calibrations
                        ### and uses the lsveto poly43 calibration
                        edep, selection = calib301(i, E)
                    
                    else:
                        print "ERROR: no info['build'] for the data calibration"
                        sys.exit()
                    #-----------------------------------------------------------------------
                    
                    
                    # DEFINE CUTS
                    #-----------------------------------------------------------------------
                    ### old V00-04-04 BDT cuts are in cutsBDT101()
                    #if info['build'] == 'V00-04-04':
                    #    masterCut = cutsBDT101(i, C, E, edep, selection)
                    if info['build'] == 'V00-04-04' or info['build'] == 'V00-04-12':
                        masterCut = cutsBDT300(i, C, E, edep, selection)
                    ### new BDT cuts for V00-04-14
                    elif info['build'] == 'V00-04-14':
                        masterCut = cutsBDT301(i, C, E, edep, selection)
                    ### same cuts for v00-04-15 so I'm told...
                    elif info['build'] == 'V00-04-15':
                        masterCut = cutsBDT301(i, C, E, edep, selection)
                    else:
                        print "ERROR: no info['build'] for the data cuts"
                        sys.exit()
                    #-----------------------------------------------------------------------
                    
                    
                    # FILL HISTOS
                    #-----------------------------------------------------------------------
                    #----
                    ### apply standard cuts
                    chain.Draw(selection+' >> '+key, masterCut)
                    #----
                    ### look at the cut data...
                    #chain.Draw(selection+' >> '+key, "!("+masterCut.GetTitle()+")")
                    #----
                    ### apply no cuts
                    #chain.Draw(selection+' >> '+key)
                    #----
                    
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

    return data


def buildMCSurf400(info, mc):
    
    #============================================================================
    #  super special case function to deal with NaI surface depth profiling...
    #  and now also teflon surface depth profiling... 2020-04-29
    #============================================================================
        
    if info['reuse']:
        for c in info['chans']:
            info['chan'] = c
            
            if not info['rootfile']:
                print 'ERROR: no rootfile specified in backgrounds file'
                sys.exit()
            
            rootfile = info['rootfile']
            if not os.path.exists(rootfile):
                print 'ERROR: rootfile not found -->', rootfile
                sys.exit()
            
            rfile = TFile(rootfile, "READ")

            for e in range(2):
                key  = info['key']
                key += '-c'+info['chan']
                key += '-e'+str(e)

                # exclude lsveto low energy
                #if 'x9' in key and e == 0:
                #    continue
                
                try:
                    mc[key] = {}
                    mc[key]['info'] = info
                    mc[key]['hist'] = deepcopy(TH1F(rfile.Get(key)))
                    mc[key]['pars'] = getPars(mc[key]['hist'])
                    #mc[key]['hist'].Sumw2()
                except:
                    print "ERROR: could not find hist -->",key
                    sys.exit()
                    
                try:
                    mc[key]['generated_hist'] = deepcopy(TH1F(rfile.Get(key+'_generated')))
                    #mc[key]['generated'] = mc[key]['generated_hist'].GetBinContent(1)
                    mc[key]['generated'] = mc[key]['generated_hist'].GetEntries()
                except:
                    print "ERROR: could not find generated_hist -->",key+'_generated'
                    sys.exit()
        
    else:
        
        ### select correct paths for simulation files
        if   info['sim'] == 'v3.1.1': path1, path2, pushpasMC = mcPath101(info) # G4.9
        elif info['sim'] == 'v4.0.1': path1, path2, pushpasMC = mcPath300(info) # G4.10 SET1
        elif info['sim'] == 'v4.0.2': path1, path2, pushpasMC = mcPath301(info) # G4.10 SET2
        else:
            print 'ERROR: simulation version not valid -->',info['sim']
            sys.exit()
            
        ### append the path if running locally 
        if amLocal():
            path1 = os.path.join('/home/mkauer/COSINE/CUP/mc-fitting','.'+path1)
        
        print 'INFO: looking for files with -->', os.path.join(path1, path2)
        
        chain  = TChain("MC","")
        nfiles = 0
        nfiles = chain.Add(os.path.join(path1, path2))
        
        if nfiles > 0:
            
            print 'INFO:',nfiles,'files found for', info['loca'], info['isof']
            
            for c in info['chans']:
                info['chan'] = c
                
                for e in range(2):
                    
                    # DEFINE HIST AND KEY
                    #-----------------------------------------------------------------------
                    xstal = int(info['xstl']) - 1
                    fromx = int(info['from']) - 1
                    
                    key  = info['key']
                    key += '-c'+info['chan']
                    key += '-e'+str(e)
                    
                    mc[key] = {}
                    mc[key]['info'] = info
                    
                    pars = histparam64(e)
                    mc[key]['pars'] = pars
                    histo = TH1F(key, longNames(xstal), pars[0], pars[1], pars[2])

                    # 2020-03-19 - add this here
                    key2 = key+'_generated'
                    temp2 = TH1F(key2, 'generated', 1, 1, 2)
                    #-----------------------------------------------------------------------

                    
                    ###  2020-03-19
                    ###  adding in the event-by-event selction criteria
                    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                    ### try getting the expo depth from file name?
                    ### fix the depth as par[0]
                    depth = float(info['floca'].split('-')[3])
                    func = TF1('func', surfProfile, 0, 10, 1)
                    func.SetParameter(0, depth/100.)
                    #print 'using depth', depth
                    
                                        
                    entries = chain.GetEntries()
                    #print entries
                    
                    for jentry in range(entries):
                        
                        chain.GetEntry(jentry)

                        ### pull entry data our of leaf
                        #--------------------------------------------------
                        xx = chain.GetLeaf('primX0').GetValue()
                        yy = chain.GetLeaf('primY0').GetValue()
                        zz = chain.GetLeaf('primZ0').GetValue()
                        
                        #elow = chain.GetLeaf('edepResA').GetValue(xstal)
                        #ehi = chain.GetLeaf('edepResD').GetValue(xstal)
                        ### these branch names changed in gyunho's sim
                        elow = chain.GetLeaf('edepResolA').GetValue(xstal)
                        ehi = chain.GetLeaf('edepResolD').GetValue(xstal)
                        
                        shit = chain.GetLeaf('singleHitTag').GetValue(xstal)
                        mhit = chain.GetLeaf('multipleHitTag').GetValue(xstal)
                        evType = chain.GetLeaf('evt_Type').GetValue()
                        group = chain.GetLeaf('groupNo').GetValue()
                        decayType = chain.GetLeaf('primDecayType').GetValue()
                        
                        
                        ### alpha cuts?
                        #--------------------------------------------------
                        # 97 == alpha
                        # 98 == beta
                        if decayType == 97: continue
                        
                        
                        
                        ### volume cuts
                        #--------------------------------------------------
                        X0, Y0, Z0, rad, height, dep = mcDimensions(fromx)
                        dist = sqrt((X0-xx)**2 + (Y0-yy)**2)
                        zdist = abs(Z0-zz)

                        # teflon thickness (mm) - 2020-04-29
                        #teflon = 1.015
                        # i think this changed in gyunho's new sim
                        teflon = 1.016

                        ### kinda cryptic...
                        # /100 to conver to um
                        # /1000 to convert to mm
                        # *10 to go 10 expo depths in
                        #tdep = (((depth/100.)/1000.) * 10)
                        # use dep = 0.01 = 10um because that's all that was simulated
                        # actually set dep to 100um = 0.1
                        dep = 0.1
                        tdep = dep
                        
                        # for NaI surf side
                        if info['loca'].endswith('side') or info['loca'].endswith('expo'):
                            if dist > rad or dist < rad-dep: continue
                            if zdist > height: continue

                        # for NaI surf face
                        elif info['loca'].endswith('face'):
                            if zdist > height or zdist < height-dep: continue
                            if dist > rad: continue
                        
                        # teflon in-side or out-side - 2020-04-29
                        #elif info['loca'].endswith('in') or info['loca'].endswith('out'):
                        #    if dist < rad or dist > rad+teflon: continue
                        #    if zdist > height: continue
                        
                        # teflon in-side - 2020-05-01
                        elif info['loca'].endswith('in'):
                            if dist < rad or dist > rad+tdep: continue
                            if zdist > height: continue

                        # teflon out-side - 2020-05-01
                        elif info['loca'].endswith('out'):
                            if dist < rad+teflon-tdep or dist > rad+teflon: continue
                            if zdist > height: continue

                        else:
                            print 'ERROR: I do not know what to do with', info['loca']
                            sys.exit()

                        
                        ### get the expo weighting
                        #--------------------------------------------------
                        # for NaI surf side
                        if info['loca'].endswith('side') or info['loca'].endswith('expo'):
                            wtf  = func.Eval(1000.*(rad-dist))
                        
                        # for NaI surf face
                        elif info['loca'].endswith('face'):
                            wtf  = func.Eval(1000.*(height-zdist))
                        
                        # teflon in-side - 2020-04-29
                        elif info['loca'].endswith('in'):
                            wtf  = func.Eval(1000.*(dist-rad))
                        
                        # teflon out-side - 2020-04-29
                        elif info['loca'].endswith('out'):
                            wtf  = func.Eval(1000.*(rad+teflon-dist))
                        
                        else:
                            print 'ERROR: I do not know what to do with', info['loca']
                            sys.exit()
                        
                        #wtf = 1

                        
                        ### generated events - I might need to reposition this code
                        #--------------------------------------------------
                        # to stay consisent I just need volume and event type cuts?
                        # err, I think I need to add the exp weighting...
                        if evType <= 10:
                            #wtf=1 # testing
                            temp2.Fill(1, wtf)
                        
                        
                        ### energy cut
                        #--------------------------------------------------
                        # do I need separate cuts for low and high energy
                        if e==0:
                            if elow*1000. < 0.1: continue
                        if e==1:
                            if ehi*1000. < 0.1: continue
                        
                        
                        ### single/multi hit cut
                        #--------------------------------------------------
                        if c == 'S':
                            if shit != 1 or mhit != -1: continue
                        if c == 'M':
                            if shit != -1 or mhit != 1: continue
                        
                        
                        ### do I need an lsveto cut?
                        #--------------------------------------------------
                        # -- insert here...
                        # probably not, minimal impact if any.
                        # single/multi hit cut will take care of this.
                        
                        
                        ### do I need a groupNo cut?
                        #--------------------------------------------------
                        # -- insert here...
                        # not needed here for Pb210.
                        # groups 14-15 are included by default anyway.

                        
                        ### event type cut
                        #--------------------------------------------------
                        if evType <= 10: continue
                        
                        
                        ### fill the weighted histograms
                        #--------------------------------------------------
                        if e==0:
                            histo.Fill(1000.*elow, wtf)
                        if e==1:
                            histo.Fill(1000.*ehi, wtf)
                        
                        

                    #---------------------------------------------------------------------
                    mc[key]['generated_hist'] = temp2
                    generated = temp2.GetEntries()
                    if generated <= 0:
                        if int(info['xstl']) == int(info['from']):
                            print 'WARNING: no events generated for -->', info['floca'], info['isof'],\
                                  'xstl', info['xstl'], 'from', info['from']
                        else:
                            print 'WARNING: no events generated for -->', info['floca'], info['isof'],\
                                  'xstl', info['xstl'], 'from', info['from']
                    mc[key]['generated'] = generated

                    #---------------------------------------------------------------------
                    
                    detected = histo.GetEntries()
                    mc[key]['hist'] = histo
                    
                    if (c=='S' and e==0) or (c=='M' and e==1):
                        print 'DEBUG:', key, 'generated events =', generated
                        print 'DEBUG:', key, 'detected events =', detected
                        try: effic = round(100*detected/generated, 2)
                        except: effic = 0.0
                        print 'DEBUG:', key, 'efficiency =', effic, '%'

                    #---------------------------------------------------------------------
                       
        else:
            print 'ERROR: no MC files found for -->', \
                'x'+str(info['xstl']), info['floca'], info['isof']
            sys.exit()
    
    return mc


def ScaleEnergy400(bkgs, binShift):

    for key in bkgs:
        if 'cS' in key and 'e0' in key:
            if 'teflon' in key:
                xstal = int(key.split('-')[0][1]) - 1
                if binShift[xstal] > 0:
                    bkgs[key]['hist'] = MCBinShift400(bkgs[key]['hist'], binShift[xstal])

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

