#!/usr/bin/env python
######################################################################
# funcs400.py
# 
# version: 2019-09-03
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
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
        if len(info) == 0:
            continue
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # These are mandatory excludes because they are always global
        exclude = bool('lsveto' in info['key']
                    or 'plastic' in info['key']
                    or 'innersteel' in info['key']
                    or 'steel' in info['key']
                    or 'gamma' in info['key'])
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        infos = []
        exception = 0
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
            print 'WARNING:',info['key'],'does not fit any known criteria'
            print '         Please check out build101() '
            infos.append(info)
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for info in infos:
            #print info['key']
            if 'D' in info['type']:
                data = buildData400(info, data)
                
            elif 'B' in info['type']:
                bkgs = buildMC300(info, bkgs)
            
            elif 'F' in info['type']:
                sigs = buildMC300(info, sigs)
            
            else:
                print 'WARNING: I do not know type',info['type'],'in line:'
                print info['line']
                continue
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
                if tmpruntime > 0:
                    runtime += tmpruntime
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


