#!/usr/bin/env python
######################################################################
# funcs500.py
# 
# version: 2021-04-02
# 
# Change Log (key == [+] added, [-] removed, [~] changed)
#---------------------------------------------------------------------
# + 'x0' now special for external/global backgrounds
# + 'x99' now flags to process rootfile for all crystals
# 
# email: mkauer@physics.wisc.edu
######################################################################

import sys
import os
import re
from copy import copy, deepcopy
import numpy as np
from numpy import sqrt, exp
from scipy.special import erf
import random
import json
import glob
import subprocess
import datetime
import time

from ROOT import *
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from funcs_old import *
from funcs_misc import *
from funcs_new import *

# skip processing rootfiles that already exist
SKIP_EXISTING = 0
FILE_TOO_SMALL = 10  # in kB


def main(argv):
    if len(sys.argv) == 2:
        mcfile = sys.argv[1]
        if not os.path.exists(mcfile):
            print('ERROR: file does not exists --> {0}'.format(mcfile))
            sys.exit()
    else:
        print('Specify an MC backgrounds file')
        sys.exit()
        
    # submit jobs to build on the cluster
    data, bkgs, sigs = build500(mcfile, 1, 0, ['S','M'], [], True)

    
def getInfo500(line, freuse=0, fchans=['S','M'], fxstals=[]):

    info = {}
    info['line'] = line
    bits = line.split()
    
    #-----------------------------------------------------------------
    # data/sim type
    info['type'] = str(bits[0])
    if info['type'][0] not in ['D', 'B', 'F']:
        print('ERROR: unknown data type [{0}]'.format(info['type'][0]))
        sys.exit()
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # channels to process
    info['chans'] = fchans
    # energies to process
    info['energies'] = [0, 1, 2]
    # crystals to process
    info['xstals'] = sorted(fxstals)
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    # reuse a joined rootfile
    if 'R' in info['type']:
        info['reuse'] = 1
        test = info['type'].split('R')
        if test[1]:
            tmp = []
            for e in test[1]:
                tmp.append(int(e))
            info['energies'] = tmp
    else:
        info['reuse'] = 0
    # force reuse of everything - nice for debugging
    if freuse == 1:
        info['reuse'] = 1
    # force not reuse of everything - nice for debugging
    if freuse == 2:
        info['reuse'] = 0
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # crystal number
    info['xstl'] = int(bits[1])
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    # skip non-selected crystals
    # 0 and 99 are not included in this check
    if len(fxstals) > 0:
        if info['xstl'] in list(range(1, 10)):
            if info['xstl'] not in fxstals:
                return []
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # from crystal number
    # set same as xstl by default
    info['fromx'] = int(bits[1])
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # reuse rootfile name
    if 'root' in str(bits[-1]):
        info['rootfile'] = os.path.join(HERE, str(bits[-1]))
    else:
        info['rootfile'] = False
    #-----------------------------------------------------------------
    

    if 'D' in info['type']:
        # location is always 'data'
        info['loca'] = str(bits[2])

        # data set
        info['tag'] = str(bits[3])
        
        # production version
        info['build'] = str(bits[4])
        
        # if an energy(s) is specified, set it
        # this only works when reusing rootfiles
        """
        if 'R' in info['type']:
            test = info['type'].split('R')
            if test[1]:
                tmp = []
                for e in test[1]:
                    tmp.append(int(e))
                info['energies'] = tmp
        """
        # build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-'+str(info['loca'])
        key += '-'+str(info['tag'])
        key += '-'+str(info['build']).replace('-','')
        info['key'] = key
        
    else:
        # background location as displayed in config file
        info['floca'] = str(bits[2])
        
        # background location without '-' for hist keys
        info['loca'] = str(bits[2]).replace('-','')
        
        # top level isotope file name
        info['isof'] = str(bits[3])

        # isotope chain break start
        info['chst'] = str(bits[4])

        # isotope chain break stop
        info['chsp'] = str(bits[5])

        # chain start_stop name
        info['chain'] = info['chst']+'_'+info['chsp']
        
        # activity
        info['acti'] = float(bits[6])

        # fit bounds
        info['fbnd'] = [float(bits[7]), float(bits[8])]

        # simulation version
        info['sim'] = str(bits[9])
        
        # build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-'+info['loca']
        key += '-'+info['chain']
        info['key'] = key

        # assign a plotting group to the isotope-location
        info['group'] = setGroup500(info)
        
    return info

    
def build500(infile, others=1, freuse=0, fchans=['S','M'], fxstals=[], cluster=False):

    debug = 0
    
    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):

        #print('processing line --> {0}'.format(line))

        if line.startswith('stop'):
            break
        
        info = getInfo500(line, freuse, fchans, fxstals)
        
        if not info:
            continue

        if info['reuse']:
            if 'D' in info['type']:
                data = reuseData500(info, data)

            elif 'B' in info['type']:
                bkgs = reuseMC500(info, bkgs)

            elif 'F' in info['type']:
                sigs = reuseMC500(info, sigs)
            
            else:
                print('WARNING: Unknown data/sim type \"{0}\"'.format(info['type']))
                continue
            
        else:
            if 'D' in info['type']:
                data = buildHists500(info, data, cluster)
                
            elif 'B' in info['type']:
                bkgs = buildHists500(info, bkgs, cluster)
            
            elif 'F' in info['type']:
                sigs = buildHists500(info, sigs, cluster)
            
            else:
                print('WARNING: Unknown data/sim type \"{0}\"'.format(info['type']))
                continue
            
    return data, bkgs, sigs


def reuseData500(info, data):
    
    if os.path.exists(info['rootfile']):
        tfile = TFile(info['rootfile'])
    else:
        print('ERROR: Rootfile not found --> {0}'.format(info['rootfile']))
        #return data
        sys.exit()
        
    for C in ['S', 'M']:
        #for E in [0, 1, 2]:
        for E in info['energies']:
            
            key  = info['key']
            key += '-c'+str(C)
            key += '-e'+str(E)

            try:
                _dump = TH1F(tfile.Get(key))
            except:
                print('WARNING: Hist not found --> {0}'.format(key))
                continue
                #return data
            
            data[key] = {}
            data[key]['info'] = copy(info)

            histo = tfile.Get(key)

            data[key]['hist'] = deepcopy(TH1F(histo))
            data[key]['pars'] = getPars(data[key]['hist'])
            data[key]['runtime_hist'] = deepcopy(TH1F(tfile.Get(key+'_runtime')))
            data[key]['runtime'] = data[key]['runtime_hist'].GetBinContent(1)
            
            """
            try:
                histo = tfile.Get(key)
                data[key]['hist'] = deepcopy(TH1F(histo))
                data[key]['pars'] = getPars(data[key]['hist'])
                #data[key]['hist'].Sumw2()
            except:
                print "WARNING: could not find hist -->", key
                pars = histparams(E)
                data[key]['pars'] = pars
                hist = TH1F(key, longNames(int(info['xstl'])-1), pars[0], pars[1], pars[2])
                data[key]['hist'] = deepcopy(hist)
                pass
            """
            
    return data


def reuseMC500(info, mc):
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++
    # some speed ups for hist loading
    # turn off for the real deal
    skip_alpha_fromx  = 1  # should really always skip these
    skip_lsveto       = 0
    only_alphas       = 0
    skip_alphas       = 0
    #if skip_alpha_fromx: print('WARNING: skip loading alphas fromx')
    if skip_lsveto:      print('WARNING: skip loading lsveto')
    if only_alphas:      print('WARNING: only loading alphas')
    if skip_alphas:      print('WARNING: skip loading alphas')
    #+++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    if os.path.exists(info['rootfile']):
        tfile = TFile(info['rootfile'])
    else:
        print('WARNING: Rootfile not found --> {0}'.format(info['rootfile']))
        return mc

    # see if certain crystals were defined, else set to all
    if info['xstals']:
        pcxs = info['xstals']
    else:
        pcxs = list(range(1, 10))

    xstls = []

    # do global mc
    if int(info['xstl']) == 0:
        # skip globals if only loading alphas
        if only_alphas:
            return mc
        for xstl in pcxs:
            for C in ['S', 'M']:
                for E in [0, 1, 2]:
                    
                    # only process alphas for lsveto
                    if xstl != 9 and E == 2:
                        continue
                    # only alphas from key locations
                    if xstl == 9 and E == 2 and \
                       info['floca'] not in ['lsveto', 'plastic']:
                        continue

                    newinfo = copy(info)
                    newinfo['xstl'] = xstl
                    newinfo['fromx'] = 0
                    
                    # scale neutron single-hit
                    # multi-hit fit = 3.83e-5
                    # single-hit fit = 2.22e-4
                    nsh_scale = 5.79
                    if info['floca'] == 'neutron' and C == 'S':
                        newinfo['acti'] = newinfo['acti'] * nsh_scale
                    
                    newinfo['key'] = newinfo['key'].replace('x0-', 'x'+str(xstl)+'-', 1)
                    mc = getMC500(tfile, mc, newinfo, C, E)

    # if special case '99' then process all crystals at once
    elif int(info['xstl']) == 99:
        xstls = [x for x in pcxs if x != 9]
        is99 = True
        
    elif int(info['xstl']) in pcxs and int(info['xstl']) != 9:
        xstls = [int(info['xstl'])]
        is99 = False
        
    elif int(info['xstl']) == 9:
        print('INFO: lsveto is being done in parellel with xstals')
    
    else:
        print('WARNING: unknown xstl number [{0}]'.format(info['xstl']))

    
    for xstl in xstls:
        for fromx in range(1, 9):
            
            # do some initial selection cuts if processing expo surface sim
            if isExpoSurfMC(info):
                if xstl != fromx:
                    continue
            
            for C in ['S', 'M']:
                #for E in [0, 1, 2]:
                for E in info['energies']:

                    # skip "mystery" from other channels
                    if info['floca'] == 'mystery' and (xstl != fromx or C != 'S' or E != 0):
                        continue
                    
                    # skip others if loading only alphas
                    if only_alphas and E != 2: continue

                    # skip alphas
                    if skip_alphas and E == 2: continue
                    
                    # skip alphas for special locations
                    if (info['floca'] in ['xgamma', 'xplastic', 'xpmt', 'xpmtbase']
                        and E == 2):
                        continue
                    
                    # apply 2 quenchings to alphas
                    if E == 2 and isAlphaMC(info):
                        
                        # skip alphas fromx
                        if skip_alpha_fromx and xstl != fromx: continue
                        """
                        #qis = [1, 2]
                        # or just one Q for testing
                        #qis = [1]
                        # try Q2 for just surface components
                        if 'nai-surf' in info['floca']:
                            qis = [1,2]
                        elif 'teflon-surf' in info['floca']:
                            qis = [2]
                        elif info['floca'] == 'teflon':
                            qis = [2]
                        elif info['floca'] == 'internal' and info['chst'] == 'Pb210':
                            qis = [1]
                        else:
                            qis = [1,2]
                        """
                        qis = setQis(info)
                        
                    # else false if not alphas
                    else:
                        qis = [False]
                        
                    for qi in qis:
                        # process crystal hits into crystal
                        #--------------------------------------------
                        newinfo = copy(info)
                        newinfo['xstl'] = xstl
                        newinfo['fromx'] = fromx
                        #if qi: newinfo['acti'] = newinfo['acti']/2.
                        #if len(qis) == 2: newinfo['acti'] = newinfo['acti']/2.
                        if qi: newinfo['acti'] = newinfo['acti'] * getQratio(int(xstl)-1, qi, newinfo)
                        if is99: newinfo['key'] = newinfo['key'].replace('x99-', 'x'+str(xstl)+'-', 1)
                        mc = getMC500(tfile, mc, newinfo, C, E, fromx, qi)
                        
                        if skip_lsveto:
                            continue
                        
                        # process crystal hits into lsveto
                        #--------------------------------------------
                        # skip mystery
                        if info['floca'] == 'mystery': continue
                        
                        newinfo = copy(info)
                        newinfo['xstl'] = 9
                        newinfo['fromx'] = fromx
                        #if qi: newinfo['acti'] = newinfo['acti']/2.
                        #if len(qis) == 2: newinfo['acti'] = newinfo['acti']/2.
                        if qi: newinfo['acti'] = newinfo['acti'] * getQratio(int(xstl)-1, qi, newinfo)
                        if is99: newinfo['key'] = newinfo['key'].replace('x99-', 'x9-', 1)
                        else: newinfo['key'] = newinfo['key'].replace('x'+str(xstl)+'-', 'x9-', 1)
                        mc = getMC500(tfile, mc, newinfo, C, E, fromx, qi)
                        
    tfile.Close()
    return mc

    
def getMC500(tfile, mc, info, C, E, fromx=False, Q=False):

    # create the main key
    key = info['key']
    if fromx:
        key += '-f'+str(fromx)
    if E==2 and Q:
        key += '-q'+str(Q)
    key += '-c'+str(C)
    key += '-e'+str(E)

    #print('DEBUG: loading hist --> {0}'.format(key))
    
    try:
        _dump = TH1F(tfile.Get(key))
    except:
        print('WARNING: Hist not found --> {0}'.format(key))
        return mc
    
    mc[key] = {}
    mc[key]['info'] = copy(info)
    
    histo = tfile.Get(key)
    
    mc[key]['hist'] = deepcopy(TH1F(histo))
    mc[key]['pars'] = getPars(mc[key]['hist'])
    mc[key]['generated_hist'] = deepcopy(TH1F(tfile.Get(key+'_generated')))
    mc[key]['generated'] = mc[key]['generated_hist'].GetEntries()

    return mc

    
def buildHists500(info, hists, cluster=False):

    filetype = None
    if 'B' in info['type'] or 'F' in info['type']:
        filetype = 'sim'
    elif 'D' in info['type']:
        filetype = 'data'
    else:
        print('ERROR: unknown data/sim type --> \"{0}\"'.format(info['type']))
        sys.exit()
        
    if filetype == 'sim':
        if info['sim'] == 'v3.1.1':
            directory, match, pushpasMC = mcPath101(info) # G4.9
        elif info['sim'] == 'v4.0.1':
            directory, match, pushpasMC = mcPath300(info) # G4.10 SET1
        elif info['sim'] == 'v4.0.2':
            directory, match, pushpasMC = mcPath301(info) # G4.10 SET2
        elif info['sim'] == 'v4.0.3':
            directory, match = mcPath500(info) # my reprocessed files
        else:
            print('ERROR: unknown simulation version [{0}]'.format(info['sim']))
            sys.exit()
    elif filetype == 'data':
        if info['build'] in ['V00-04-15']:
            directory = '/data/COSINE/MRGD/phys/'+info['build']
        elif info['build'] in ['V00-04-19']:
            directory = '/data/COSINE/MRGD/'+info['build']+'/*'
        elif info['build'] in ['V00-04-20']:
            directory = '/data/COSINE/MRGD/'+info['build']+'/*'
        else:
            print('ERROR: unknown data production [{0}]'.format(info['build']))
            sys.exit()
        #match = '*root*'
        match = 'mrgd_M??????.root.???'
    else:
        print('ERROR: unknown data/sim type --> \"{0}\"'.format(filetype))
        sys.exit()
    
    ### append the path if running locally 
    #if amLocal():
    #    directory = '/home/mkauer/COSINE/CUP/mc-fitting'+directory
        
    ### find all the matching files
    filenames = sorted(glob.glob(directory+'/'+match))
    if len(filenames) == 0:
        print('WARNING: no files found for --> {0}/{1}'.format(directory, match))
        return hists
    else:
        print('INFO: found {0} files'.format(len(filenames)))

    ### make sensible cuts on run numbers
    if filetype == 'data':
        newfiles = []
        for rootfile in filenames:
            # get the run number
            # like mrgd_M001858.root.010
            runnum = int(os.path.basename(rootfile)[6:12])
            if info['tag'] == 'SET1' and runnum > 1546:
                continue
            elif info['tag'] == 'SET2' and runnum > 1858:
                continue
            elif info['tag'] == 'SET3' and runnum > 1862:
                continue
            # skip runs > 1939 because not in the GRL yet
            elif info['tag'] == 'SET4' and runnum > 1939:
                continue
            # skip runs >= 1970 until C5,C8 trigger threshold is understood
            elif runnum >= 1970:
                continue
            else:
                newfiles.append(rootfile)
                
        filenames = newfiles
    
            
    if cluster:
        # do not reprocess existing files
        if SKIP_EXISTING:
            newfiles = []
            for rootfile in filenames:
                if 'naisurfexpo' in info['loca']:
                    savedir = 'nai-surf-expo'
                elif 'naisurfscrexpo' in info['loca']:
                    savedir = 'nai-surfscr-expo'
                elif 'teflonsurfexpo' in info['loca']:
                    savedir = 'teflon-surf-expo'
                else:
                    savedir = info['loca']
                if filetype == 'data':
                    savedir = info['loca']+'/'+info['build']
                outdir = os.path.join('/data/COSINE/WORK/mkauer/histograms', savedir)
                mkdir(outdir)
                outfile = info['key']+'__'+os.path.basename(rootfile)
                filename = os.path.join(outdir, outfile)
                if not os.path.exists(filename):
                    newfiles.append(rootfile)
                else:
                    if FILE_TOO_SMALL and \
                       os.path.getsize(filename)/1024 < FILE_TOO_SMALL:
                        newfiles.append(rootfile)
            filenames = newfiles

        print('INFO: submitting [{0}] jobs to cluster'.format(len(filenames)))
        
        # submit each rootfile as a separate job
        for rootfile in filenames:
            
            # wait until there's room in the queue
            while True:
                qlen = checkQueue()
                if qlen <= 0 or qlen > 990:
                    #print('sleeping...')
                    time.sleep(60)
                else:
                    break
            
            # write info to json
            infojson = '{0}-info.json'.format(info['key'])
            scratchdir = '/data/COSINE/WORK/mkauer/scratch'
            mkdir(scratchdir)
            jsonfile = os.path.join(scratchdir, infojson)
            # Don't overwrite json if it already exists
            # This was causing errors when the cluster
            # was trying to read the file while it was
            # being recreated.
            if not os.path.exists(jsonfile):
                with open(jsonfile, 'w') as jfile:
                    json.dump(info, jfile)

            # submit job to the cluster
            submitJob(filetype, info['key'], rootfile, jsonfile)

            # break after one file (for testing)
            #break
        
    else:
        # process localy
        #for rootfile in filenames[:3]:
        for rootfile in filenames:
            print('Processing file --> {0}'.format(rootfile))
            if filetype == 'data':
                hists = dataToCluster(hists, rootfile, info)
            else:
                hists = simToCluster(hists, rootfile, info)
    return hists


def submitJob(filetype, basekey, rootfile, jsonfile):
    
    scratchdir = '/data/COSINE/WORK/mkauer/scratch'
    mkdir(scratchdir)

    jobname = basekey+'__'+os.path.basename(rootfile)
    
    jobfile = os.path.join(scratchdir, jobname+'.job')
    
    myscript = os.path.join(
        '/home/mkauer/mc-fitting/root-join-read',
        'process500.py')
    
    logdir = '/data/COSINE/WORK/mkauer/logs'
    mkdir(logdir)
    
    logfile = os.path.join(logdir, jobname+'.log')

    date = str(datetime.datetime.now().date())
    
    with open(jobfile, 'w') as sfile:
        sfile.write('#!/bin/bash \n\n')
        sfile.write('source /home/mkauer/py2root.sh \n')
        #sfile.write('source /home/mkauer/py3root.sh \n')

        # log versions of things
        sfile.write('logfile='+logfile+'\n')
        sfile.write('echo \"Versions of stuff \" &> $logfile \n')
        sfile.write('hostname &>> $logfile \n')
        sfile.write('which root &>> $logfile \n')
        sfile.write('root --version &>> $logfile \n')
        sfile.write('uname -sr &>> $logfile \n')
        sfile.write('python --version &>> $logfile \n')
        
        sfile.write('export nodename=`hostname` \n')
        sfile.write('export pbsuser=$USER \n')
        sfile.write('export pbsdate='+date+' \n')
        sfile.write('export workdir='+scratchdir+' \n')
        sfile.write(myscript
                    +' '+filetype
                    +' '+rootfile
                    +' '+jsonfile
                    +' 1>> '+logfile
                    +' 2>> '+logfile
                    +' \n')
        
            
    # send data jobs to the short queue
    if filetype == 'data':
        queue = 'short'
    else:
        #queue = 'long'
        # I think medium queue should be fine for MC?
        queue = 'medium'
    
    # make the submit process smarter
    # sometimes the server does not respond
    # qsub: submit error (Pbs Server is currently too busy...
    while True:
        try:
            out = subprocess.check_output(
                ['qsub',
                 #'-q', 'short',
                 #'-q', 'medium',
                 #'-q', 'long',
                 '-q', queue,
                 '-o', logdir,
                 '-e', logdir,
                 jobfile],
                stderr=subprocess.STDOUT
            ).decode(encoding='utf-8', errors='ignore')
            print('{0}  {1}'.format(out.strip(), jobname))
            break
        except Exception as e:
            if hasattr(e, 'output'): print((e.output).strip())
            else: print(e)
            time.sleep(10)
        
    return


def checkQueue():
    # make this check smarter
    # sometimes the server does not respond
    # qstat: Pbs Server is currently too busy...
    try:
        out = subprocess.check_output(['qstat', '-u', 'mkauer'],
                                      stderr=subprocess.STDOUT
                                      ).decode(encoding='utf-8',
                                               errors='ignore').split('\n')
        return len(out)
    except Exception as e:
        if hasattr(e, 'output'): print((e.output).strip())
        else: print(e)
        return -1


def simToCluster(mc, rootfile, info, jsonfile=None, cluster=False):

    #+++++++++++++++++++++++++++++++++++++++++++++++++
    # some speed ups for local testing
    # turn off for the real processing
    skip_fromx        = 0
    skip_alphas       = 0
    skip_alphas_fromx = 1  # should really always skip these
    skip_multi_alphas = 0
    only_loE          = 0
    only_hiE          = 0
    only_alphas       = 0
    only_alphaQ1      = 0
    only_alphaQ2      = 0
    use_my_Qs         = 0
    skip_lsveto       = 0
    skip_byEvent      = 0
    if skip_fromx:        print('WARNING: skip processing fromx')
    if skip_alphas:       print('WARNING: skip processing alphas')
    #if skip_alphas_fromx: print('WARNING: skip processing alphas fromx')
    if skip_multi_alphas: print('WARNING: skip processing multi-hit alphas')
    if only_loE:          print('WARNING: only processes low-energy')
    if only_hiE:          print('WARNING: only processes high-energy')
    if only_alphas:       print('WARNING: only processes alphas')
    if only_alphaQ1:      print('WARNING: only processes alpha Q1')
    if only_alphaQ2:      print('WARNING: only processes alpha Q2')
    if use_my_Qs:         print('WARNING: only processing my defined alpha Qs')
    if skip_lsveto:       print('WARNING: skip processing lsveto')
    if skip_byEvent:      print('WARNING: skip processing by event')
    #+++++++++++++++++++++++++++++++++++++++++++++++++
    
    if info is None and not os.path.exists(jsonfile):
        print('WARNING: you must specify an info json filename')
        sys.exit()
        
    if info is None:
        with open(jsonfile) as jf:
            info = json.load(jf)
    
    if mc is None:
        mc = {}
    
    chain = TChain("MC", "")
    chain.Add(rootfile)
    
    # see if certain crystals were defined, else set to all
    if info['xstals']:
        pcxs = info['xstals']
    else:
        pcxs = list(range(1, 10))
    
    xstls = []
    
    # if external/global mc, get energy deposit in all crystals and lsveto
    if int(info['xstl']) == 0:
        for xstl in pcxs:
            for C in ['S', 'M']:
                for E in [0, 1, 2]:
                    
                    # only process external alphas for lsveto
                    if xstl != 9 and E == 2: continue
                    
                    # or skip alphas
                    if skip_alphas and E == 2: continue
                    
                    # or just high energy
                    if only_hiE and E != 1: continue

                    # or just low energy
                    if only_loE and E != 0: continue
                    
                    newinfo = copy(info)
                    newinfo['xstl'] = xstl
                    newinfo['fromx'] = 0
                    newinfo['key'] = newinfo['key'].replace('x0-', 'x'+str(xstl)+'-', 1)
                    # process lsveto data by event to quench alphas
                    if xstl == 9 and info['floca'] in ['lsveto', 'plastic'] \
                       and info['isof'] in ['U238', 'Th232', 'U235']:
                        if skip_byEvent:
                            mc = processMC500(chain, mc, newinfo, C, E)
                        elif E in [0, 2]:
                            # don't process by event if not lsveto high energy
                            mc = processMC500(chain, mc, newinfo, C, E)
                        else:
                            mc = processMCbyEvent500(chain, mc, newinfo, C, E)
                    else:
                        mc = processMC500(chain, mc, newinfo, C, E)
                        # by-event for testing
                        #if C == 'M' and E ==0:
                        #    mc = processMCbyEvent500(chain, mc, newinfo, C, E)
                    
    # if special case '99' then process all crystals at once
    elif int(info['xstl']) == 99:
        xstls = [x for x in pcxs if x != 9]
        is99 = True
        
    elif int(info['xstl']) in pcxs and int(info['xstl']) != 9:
        xstls = [int(info['xstl'])]
        is99 = False
        
    elif int(info['xstl']) == 9:
        print('INFO: lsveto is being done in parellel with xstals')
        
    else:
        print('WARNING: unknown xstl number [{0}]'.format(info['xstl']))

        
    for xstl in xstls:
        for fromx in range(1, 9):
            
            # skip_fromx for all?
            if skip_fromx and xstl != fromx: continue

            # do some initial selection cuts if processing expo surface sim
            if isExpoSurfMC(info):
                # only process hits into crystal from same crystal
                if xstl != fromx:
                    continue
                # skip processing xstl if rootfile -Cx- doesn't match
                mc_cx = (os.path.basename(rootfile)).split('-')[4]
                mc_cx = int(mc_cx[1])
                if mc_cx != xstl:
                    continue
            
            for C in ['S', 'M']:
                for E in [0, 1, 2]:

                    # skip_alphas?
                    if skip_alphas and E == 2: continue

                    # only alphas?
                    if only_alphas and E != 2: continue

                    # only high energy?
                    if only_hiE and E != 1: continue

                    # only low energy?
                    if only_loE and E != 0: continue

                    # skip alphas for special locations
                    if (info['floca'] in ['xgamma', 'xplastic', 'xpmt', 'xpmtbase']
                        and E == 2):
                        continue
                    
                    # apply 2 quenchings to alphas or false
                    if E == 2 and isAlphaMC(info):
                        # skip alphas fromx
                        if skip_alphas_fromx and xstl != fromx: continue

                        # skip multi-hit alphas
                        if skip_multi_alphas and C == 'M': continue

                        # get my set Qs
                        if use_my_Qs:
                            qis = setQis(info)
                        else:
                            qis = [1, 2]
                        
                        # just process Q1 or Q2
                        if only_alphaQ1: qis = [1]
                        if only_alphaQ2: qis = [2]
                        
                        # special case for parameterized internal Pb210
                        if info['floca'] == 'internal' and info['chst'] == 'Pb210' \
                           and xstl not in [5, 8]:  # skip C5 and C8 until fixed
                            qis = [1]
                        
                    else:
                        qis = [False]
                    
                    for qi in qis:
                        # process crystal hits into crystal
                        newinfo = copy(info)
                        newinfo['xstl'] = xstl
                        newinfo['fromx'] = fromx
                        #if qi: newinfo['acti'] = newinfo['acti']/2.
                        #if len(qis) == 2: newinfo['acti'] = newinfo['acti']/2.
                        #if qi: newinfo['acti'] = newinfo['acti'] * getQratio(int(xstl)-1, qi, newinfo)
                        if qi:
                            if use_my_Qs:
                                newinfo['acti'] = newinfo['acti'] * getQratio(int(xstl)-1, qi, newinfo)
                            elif len(qis) == 2:
                                newinfo['acti'] = newinfo['acti']/2.
                            else:
                                pass
                        if is99: newinfo['key'] = newinfo['key'].replace('x99-', 'x'+str(xstl)+'-', 1)
                        if isExpoSurfMC(info):
                            #mc = processSurf500(chain, mc, newinfo, C, E, fromx, qi)
                            #mc = processSurf501(chain, mc, newinfo, C, E, fromx, qi)
                            mc = processSurf502(chain, mc, newinfo, C, E, fromx, qi)
                        # some special cases for processing by event
                        elif (
                                # special case for high-energy Bi-Po pile-up
                                (info['floca'] == 'internal' and E == 1 and 
                                 info['chst'] in ['Rn222', 'Th228', 'Pa231']) 
                                or
                                # special case for Pb210 alpha parameterization
                                (info['floca'] == 'internal' and E == 2 and
                                 info['chst'] == 'Pb210')
                                or
                                # special case for alpha-alpha pile-up
                                (info['floca'] == 'internal' and E == 2 and
                                 info['chst'] in ['Th228', 'Pa231']) ):
                            if skip_byEvent:
                                mc = processMC500(chain, mc, newinfo, C, E, fromx, qi)
                            elif xstl != fromx:
                                # don't process by event if fromx
                                mc = processMC500(chain, mc, newinfo, C, E, fromx, qi)
                            else:
                                mc = processMCbyEvent500(chain, mc, newinfo, C, E, fromx, qi)
                                #mc = processMC500(chain, mc, newinfo, C, E, fromx, qi)
                        else:
                            mc = processMC500(chain, mc, newinfo, C, E, fromx, qi)

                        # skip_lsveto?
                        if skip_lsveto: continue

                        # skip special cases for globals as x99
                        # don't need to process all fromx into x9
                        # just gets overwritten multiple times
                        if (info['floca'] in ['xgamma', 'xplastic']
                            and xstl != fromx):
                            continue
                        # process crystal hits into lsveto
                        newinfo = copy(info)
                        newinfo['xstl'] = 9
                        newinfo['fromx'] = fromx
                        #if qi: newinfo['acti'] = newinfo['acti']/2.
                        #if len(qis) == 2: newinfo['acti'] = newinfo['acti']/2.
                        #if qi: newinfo['acti'] = newinfo['acti'] * getQratio(int(xstl)-1, qi, newinfo)
                        if qi and len(qis) == 2:
                            newinfo['acti'] = newinfo['acti']/2.
                        if is99: newinfo['key'] = newinfo['key'].replace('x99-', 'x9-', 1)
                        else: newinfo['key'] = newinfo['key'].replace('x'+str(xstl)+'-', 'x9-', 1)
                        if isExpoSurfMC(info):
                            #mc = processSurf500(chain, mc, newinfo, C, E, fromx, qi)
                            #mc = processSurf501(chain, mc, newinfo, C, E, fromx, qi)
                            mc = processSurf502(chain, mc, newinfo, C, E, fromx, qi)
                        else:
                            mc = processMC500(chain, mc, newinfo, C, E, fromx, qi)
    
    if cluster:
        if 'naisurfexpo' in info['loca']:
            savedir = 'nai-surf-expo'
        elif 'naisurfscrexpo' in info['loca']:
            savedir = 'nai-surfscr-expo'
        elif 'teflonsurfexpo' in info['loca']:
            savedir = 'teflon-surf-expo'
        else:
            savedir = info['loca']
        
        outdir = os.path.join('/data/COSINE/WORK/mkauer/histograms', savedir)
        mkdir(outdir)
        outfile = info['key']+'__'+os.path.basename(rootfile)
        filename = os.path.join(outdir, outfile)

        print('Writing rootfile --> {0}'.format(filename))
        tfile = TFile(filename, 'RECREATE')
        writeHists500(mc)
        tfile.Write()
        tfile.Close()

        # get file size
        fsize = os.path.getsize(filename) # in Bytes
        fsize = round(fsize/1024.0, 1) # in KB
        print('Rootfile size = {0} KB'.format(fsize))
        
    return mc


def isExpoSurfMC(info):
    if 'expo' in info['loca']:
        return 1
    else:
        return 0


def isAlphaMC(info):
    if (('internal' == info['loca'] or
         'teflon' == info['loca'] or
         'naisurf' in info['loca'] or
         'teflonsurf' in info['loca']) and
        (info['isof'] in ['U238', 'Th232', 'Pb210', 'U235'])):
        return 1
    else:
        return 0

    
def processMC500(chain, mc, info, C, E, fromx=False, Q=False):

    # get info from mc key
    bits = info['key'].split('-')
    xstal = bits[0][1]
    i = int(xstal) - 1
    location = bits[1]
    isotope = bits[2]
    
    
    # DEFINE HIST AND KEY
    #-----------------------------------------------------------------------
    key = info['key']
    if fromx:
        key += '-f'+str(fromx)
    if E==2 and Q:
        key += '-q'+str(Q)
    key += '-c'+str(C)
    key += '-e'+str(E)

    print('INFO: processMC500 \"{0}\"'.format(key))
    
    if key not in mc:
        mc[key] = {}

    if 'info' not in mc[key]:
        mc[key]['info'] = copy(info)

    pars = histparams(E)
    mc[key]['pars'] = pars
    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
    #-----------------------------------------------------------------------


    ###  put together all the cuts you need...
    #=====================================================================================

    ### define the edep resolution to use
    if E: edepRes = 'edepResolD'
    else: edepRes = 'edepResolA'
    
    ### define the energy cut
    threshold = 0.125  # in keV
    if i < 8:
        energyCut = TCut('('+edepRes+'['+str(i)+']*1000. >= '+str(threshold)+')')
    
    ### lsveto energy cut has crystal dependent energy cuts
    if i == 8 and C == 'S':
        energyCut = TCut('(('+edepRes+'[8]*1000 > 0.0) && ('
                         + edepRes+'[0]*1000 < '+str(threshold)+' && '
                         + edepRes+'[1]*1000 < '+str(threshold)+' && '
                         + edepRes+'[2]*1000 < '+str(threshold)+' && '
                         + edepRes+'[3]*1000 < '+str(threshold)+' && '
                         + edepRes+'[4]*1000 < '+str(threshold)+' && '
                         + edepRes+'[5]*1000 < '+str(threshold)+' && '
                         + edepRes+'[6]*1000 < '+str(threshold)+' && '
                         + edepRes+'[7]*1000 < '+str(threshold)+'))')

    if i == 8 and C == 'M':
        energyCut = TCut('(('+edepRes+'[8]*1000 > 0.0) && ('
                         + edepRes+'[0]*1000 >= '+str(threshold)+' || '
                         + edepRes+'[1]*1000 >= '+str(threshold)+' || '
                         + edepRes+'[2]*1000 >= '+str(threshold)+' || '
                         + edepRes+'[3]*1000 >= '+str(threshold)+' || '
                         + edepRes+'[4]*1000 >= '+str(threshold)+' || '
                         + edepRes+'[5]*1000 >= '+str(threshold)+' || '
                         + edepRes+'[6]*1000 >= '+str(threshold)+' || '
                         + edepRes+'[7]*1000 >= '+str(threshold)+'))')

    # just testing here
    #if i == 8: energyCut = TCut('(1)')
    
    # skip low energy lsveto
    #if i == 8 and E == 0:
    #    energyCut = TCut('(0)')

    # LS veto threshold
    ls_thresh = 80.0  # in keV
    
    ### single hit cuts
    if C == 'S':

        ### for xstals
        hitCut = '((singleHitTag['+str(i)+'] == 1) && (multipleHitTag['+str(i)+'] == -1))'
        lsvetocut = '('+edepRes+'[8]*1000. < '+str(ls_thresh)+')'
        chanCut = TCut('(('+hitCut+') && ('+lsvetocut+'))')
        
        ### for lsveto
        if i == 8:
            hitCut = '((singleHitTag[8] == 1) && (multipleHitTag[8] == -1))'
            lsvetocut = '('+edepRes+'[8]*1000. > 0.0)'
            chanCut = TCut('(('+hitCut+') && ('+lsvetocut+'))')

    ### multi hit cuts
    elif C == 'M':
        
        ### for xstals
        hitCut = '((multipleHitTag['+str(i)+'] == 1) && (singleHitTag['+str(i)+'] == -1))'
        lsvetocut = '('+edepRes+'[8]*1000. >= '+str(ls_thresh)+')'
        # multi-hit crystal can have multi-hit tag or charge in lsveto
        chanCut = TCut('(('+hitCut+') || ('+lsvetocut+'))')
        
        ### for lsveto
        if i == 8:
            hitCut = '((multipleHitTag[8] == 1) && (singleHitTag[8] == -1))'
            lsvetocut = '('+edepRes+'[8]*1000. > 0.0)'
            # this is a little confusing but
            # multi-hit lsveto must have multi-hit tag
            chanCut = TCut('(('+hitCut+') && ('+lsvetocut+'))')

    else:
        print('ERROR: I do not know what to do with channel -->', C)
        print('Available channels are [S]Single-hits, [M]Multi-hits')
        sys.exit()


    ###  alpha cuts
    #-------------------------------------------------------------------------------
    primDecayCut = TCut('(1)')
    # alpha cuts for crystals
    if i < 8:
        # alpha cuts only on internal/surface sources
        if isAlphaMC(info):
            if E == 2:
                primDecayCut = TCut('(primDecayType == "alpha")')
            else:
                primDecayCut = TCut('(primDecayType != "alpha")')
        else:
            # skip other events from contributing to the crystal alpha spectrum
            if E == 2:
                primDecayCut = TCut('(0)')
                
    # show alphas in lsveto E=2 (just to see how it looks)
    if i == 8 and E == 2:
        primDecayCut = TCut('(primDecayType == "alpha")')
    
    
    ###  primary volume cuts
    #-------------------------------------------------------------------------------
    
    pmt1 = str((int(info['fromx'])*2)-2)
    pmt2 = str((int(info['fromx'])*2)-1)
    
    
    if location == 'internal':
        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(fromx)+'Crystal")')
            
    elif 'naisurf' in location:
        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(fromx)+'Crystal")')
        
    elif 'teflonsurf' in location:
        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(fromx)+'Teflon")')
        
    elif location == 'teflon':
        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(fromx)+'Teflon")')
        
    elif location == 'xpmt':
        volumeCut = TCut('(((primPMTid[0] == '+pmt1+')'
                         +' || (primPMTid[0] == '+pmt2+'))'
                         +' && primVolumeName == "phys_pmt")')
    
    elif location == 'xpmtbase':
        volumeCut = TCut('(((primPMTid[0] == '+pmt1+')'
                         +' || (primPMTid[0] == '+pmt2+'))'
                         +' && primVolumeName == "phys_pmtbase")')
    
    elif location == 'pmt':
        volumeCut = TCut('(primVolumeName == "phys_pmt")')

    elif location == 'pmtbase':
        volumeCut = TCut('(primVolumeName == "phys_pmtbase")')

    elif location == 'vetopmt':
        volumeCut = TCut('(primVolumeName == "phys_pmt")')

    elif location == 'cucase':
        volumeCut = TCut('((primVolumeName == "NaIDet0'+str(fromx)+'Case")'
                         +' || (primVolumeName == "NaIDet0'+str(fromx)+'Fringe0")'
                         +' || (primVolumeName == "NaIDet0'+str(fromx)+'Fringe1"))')

    elif 'cucasesurf' in location:
        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(fromx)+'Case")')

    elif location == 'plastic' or location == 'xplastic':
        volumeCut = TCut('(primVolumeName == "PlasSupport")')

    elif location == 'lsveto':
        volumeCut = TCut('(primVolumeName == "lsveto")')

    elif location == 'lsvetoair':
        volumeCut = TCut('((primVolumeName == "DetPMTCover")'
                         +' || (primVolumeName == "DetPMTEnvelope"))')

    elif location == 'film':
        volumeCut = TCut('(primVolumeName == "Reflector")')

    elif location == 'airshield':
        volumeCut = TCut('(primVolumeName == "LSVetoAirRoom")')

    elif location == 'cushield':
        volumeCut = TCut('(primVolumeName == "CuShield")')
 
    elif location == 'gamma' or location == 'xgamma':
        #volumeCut = TCut('(1)')
        volumeCut = TCut('(primVolumeName == "CuShield")')

    elif location == 'innersteel':
        volumeCut = TCut('(primVolumeName == "InnerSteel")')

    elif location == 'steel':
        volumeCut = TCut('((primVolumeName == "SteelSupport")'
                         +' || (primVolumeName == "SteelSupportTop"))')

    elif location == 'neutron':
        volumeCut = TCut('(1)')
        
    else:
        print('ERROR: No volume criteria for --> {0}'.format(location))
        sys.exit()


    ### select a specific particle?
    ### this is not generally used
    ### only in special test cases for example Tl208 spectra...
    #particleCut = TCut('(primParticleName == "Pb208")')

    ### broken chain / group number cut
    #brokenChainCut = groupNum93(info)
    brokenChainCut = groupNum500(isotope)

    ### event type cut
    evType = 'evt_Type'  # old processing - still works with new processing
    #evType = 'event_info.Type'  # new processing - doesn't work with my processing

    eventTypeCut = TCut('('+evType+' > 10)')
    ### for testing for primParticleName cut
    #eventTypeCut = TCut('(1)')

    generatedCuts = TCut('('+evType+' < 10)'+' && '+volumeCut.GetTitle())
    ### for testing for primParticleName cut
    #generatedCuts = TCut(volumeCut.GetTitle())

    ### special case for H3
    ### I think this works as expected now
    #exceptions = ['H3']
    #if isotope.split('_')[0] in exceptions:
    #    generatedCuts = TCut('('+evType+' > 10)'+' && '+volumeCut.GetTitle())

    ### special case for external Tl208 gammas
    ### and now also for neutrons
    ### and also for vetopmt
    if location in ['gamma', 'xgamma', 'neutron', 'vetopmt']:
        eventTypeCut  = TCut('(1)')
        generatedCuts = TCut('(1)'+' && '+volumeCut.GetTitle())
    
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ### JUST FOR TESTING
    primParticleCut = TCut('(1)')
    #primParticleCut = TCut('(primParticleName == "Pb211")')
    #brokenChainCut = TCut('(1)')
    #primDecayCut = TCut('(1)')
    #primDecayCut = TCut('(primDecayType == "beta")')
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    ### separate real crystals from lsveto
    #=====================================================================
    if i < 8:

        masterCut = TCut('('+
                    energyCut.GetTitle()+' && '+
                    eventTypeCut.GetTitle()+' && '+
                    brokenChainCut.GetTitle()+' && '+
                    chanCut.GetTitle()+' && '+
                    primDecayCut.GetTitle()+' && '+
                    primParticleCut.GetTitle()+' && '+
                    volumeCut.GetTitle()
                         +')')
        
        # special quenching and resolution smearing for alphas
        if E == 2 and Q:

            ### default smearing method
            # alias strings are defined as TFormulas
            chain.SetAlias('rng', 'sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
            quenchFunc = quench500(i, Q, '(edep[{0}]*1000.)'.format(i))
            #quenchFunc = simpleQuench(i, Q, info['chst'], '(edep[{0}]*1000.)'.format(i))
            chain.SetAlias('quench', quenchFunc)
            sigmaFunc = sigma500(i, E, quenchFunc)
            chain.SetAlias('sigma', sigmaFunc)
            selection = '((quench) + (sigma*rng))'
            
            ### or no quench and no smear for testing
            #selection = '(edep['+str(i)+']*1000.)'

            
        # different resolution smearing for C5 and C8
        elif i in [4, 7]:
            chain.SetAlias('rng', 'sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
            sigmaFunc = sigma500(i, E, '(edep[{0}]*1000.)'.format(i))
            chain.SetAlias('sigma', sigmaFunc)
            selection = '((edep[{0}]*1000.) + (sigma*rng))'.format(i)
        
        # tweak Na22 simulated energy
        elif E == 1 and isotope == 'Na22_GRND':
            edep = '(edep['+str(i)+']*1000.)'
            chain.SetAlias('rng', 'sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
            
            # 01 - I think this one looks better than 02
            _A = 1.0137560718825872e-08
            _B = -1.6924372865296272e-05
            _C = 1.0060126807346137
            _D = -0.1293348751221094
            
            # 02
            #_A = -4.97461464903768e-08
            #_B = 0.00021916603720066027
            #_C = 0.7997267665534115
            #_D = 30.61334935395221
            
                        
            ecal = ('(({1}*({0}**3)) + ({2}*({0}**2)) + ({3}*({0})) + ({4}))'
                    .format(edep, _A, _B, _C, _D))
            sigmaFunc = sigma500(i, E, ecal)
            chain.SetAlias('sigma', sigmaFunc)
            selection = '(({0}) + (sigma*rng))'.format(ecal)
            
            #_A = 0.20313968633685495
            #_B = 2.6081599657745107
            #_C = 0.4204609009760968
            #_D = 0.11051443449822619
            #_E = 0.004690334465459173
            #chain.SetAlias('rng', 'sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
            #edep = '(edep['+str(i)+'])'
            #ecal = ('({0}+({1}*((exp(-(({0}-{2})/{3}))+exp(({0}-{2})/{4}))**-1))+{5})'
            #        .format(edep,_A,_B,_C,_D,_E))
            ## ecal is in MeV
            #sigmaFunc = sigma500(i, E, ecal*1000.)
            #chain.SetAlias('sigma', sigmaFunc)
            #selection = '(({0}*1000.) + (sigma*rng))'.format(ecal)
        
        # use default resolution
        else:
            selection = '('+edepRes+'['+str(i)+']*1000.)'
            #selection = '(edep['+str(i)+']*1000.)'
        
    ### if lsveto
    if i == 8:

        masterCut = TCut('('+
                    energyCut.GetTitle()+' && '+
                    eventTypeCut.GetTitle()+' && '+
                    brokenChainCut.GetTitle()+' && '+
                    chanCut.GetTitle()+' && '+
                    primDecayCut.GetTitle()+' && '+
                    #particleCut.GetTitle()+' && '+
                    volumeCut.GetTitle()
                         +')')
        
        selection = '('+edepRes+'['+str(i)+']*1000.)'
        # for no smearing test
        #selection = '(edep[{0}]*1000.)'.format(i)
        
        ### or use my resolution
        #-----------------------------------------------------------
        chain.SetAlias('rng', 'sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
        sigmaFunc = sigma500(i, E, '(edep[{0}]*1000.)'.format(i))
        chain.SetAlias('sigma', sigmaFunc)
        selection = '((edep[{0}]*1000.) + (sigma*rng))'.format(i)


        
    #=====================================================================

    
    ### create the MC histogram
    #---------------------------------------------------------------------
    chain.Draw(selection+' >> '+key, masterCut)
    this_detected = histo.GetEntries()
    #print('{0} has events = {1}'.format(key, this_detected))
    if 'hist' in mc[key]:
        mc[key]['hist'].Add(histo)
    else:
        mc[key]['hist'] = copy(histo)
    #---------------------------------------------------------------------
    
    ### create a hist of the generated events
    #---------------------------------------------------------------------
    key2 = key+'_generated'
    temp2 = TH1F(key2, 'generated', 1, 0, 1)

    chain.Draw('0 >> '+key2, generatedCuts)
    this_generated = temp2.GetEntries()
    if this_generated <= 0:
        print('WARNING: no events generated for --> {0}'.format(key))
        
    if 'generated_hist' in mc[key]:
        mc[key]['generated_hist'].Add(temp2)
    else:
        mc[key]['generated_hist'] = copy(temp2)
        
    generated = mc[key]['generated_hist'].GetEntries()
    mc[key]['generated'] = generated
    #---------------------------------------------------------------------

    ### print out efficiency numbers
    #---------------------------------------------------------------------
    """
    this_detected = histo.GetEntries()
    try: this_eff = round(100*this_detected/this_generated, 2)
    except: this_eff = 0.0
    #try: tot_eff = round(100*detected/generated, 2)
    #except: tot_eff = 0.0
    print('DEBUG: [{3}] This  det/gen = eff --> {0} / {1} = {2}%'
          .format(this_detected, this_generated, this_eff, key))
    #print('DEBUG: [{3}] Total det/gen = eff --> {0} / {1} = {2}%'
    #      .format(detected, generated, tot_eff, key))
    """
    #---------------------------------------------------------------------

    return mc


def processMCbyEvent500(chain, mc, info, C, E, fromx=False, Q=False):

    # get info from mc key
    bits = info['key'].split('-')
    xstal = bits[0][1]
    cx = int(xstal)-1
    location = bits[1]
    isotope = bits[2]
    
    # create new key and hist
    #-----------------------------------------------------------------------
    key = info['key']
    if fromx:
        fx = int(fromx)
        key += '-f'+str(fromx)
    if E==2 and Q:
        key += '-q'+str(Q)
    key += '-c'+str(C)
    key += '-e'+str(E)
    
    print('INFO: processMCbyEvent500 \"{0}\"'.format(key))
    
    if key not in mc:
        mc[key] = {}

    if 'info' not in mc[key]:
        mc[key]['info'] = copy(info)

    pars = histparams(E)
    mc[key]['pars'] = pars
    histo = TH1F(key, longNames(cx), pars[0], pars[1], pars[2])
    
    key2 = key+'_generated'
    temp2 = TH1F(key2, 'generated', 1, 0, 1)
    #-----------------------------------------------------------------------

    # get the pb210 alpha skew func
    if E==2 and Q and info['floca'] == 'internal' and info['chst'] == 'Pb210':
        func = getSkewFunc(cx)
    
    entries = chain.GetEntries()
    #print entries

    Qcount = 1
    for jentry in range(entries):
    #for jentry in range(100):
        
        chain.GetEntry(jentry)
        
        # get all crystal/lsveto energies in keV
        if E: evar = 'edepResolD'
        else: evar = 'edepResolA'
        edepRes = []
        for i in range(9):
            edepRes.append(1000.*(chain.GetLeaf(evar).GetValue(i)))
        #elow = chain.GetLeaf('edepResolA').GetValue(cx)
        #ehi = chain.GetLeaf('edepResolD').GetValue(cx)
        
        edep = 1000.*(chain.GetLeaf('edep').GetValue(cx))
        #equench = 1000.*(chain.GetLeaf('edepQuenched').GetValue(cx))
        
        singleHit = chain.GetLeaf('singleHitTag').GetValue(cx)
        multiHit = chain.GetLeaf('multipleHitTag').GetValue(cx)
        evType = chain.GetLeaf('evt_Type').GetValue()
        groupNo = chain.GetLeaf('groupNo').GetValue()
        volName = chain.primVolumeName
        decayType = chain.primDecayType
        particle = chain.primParticleName
        #pdg = chain.GetLeaf('primPDGcode').GetValue()

        #if edep > 0:
        #    print('{0}  {1}  {2}'.format(particle, decayType, edep))
        
        ### look at events before cuts
        #--------------------------------------------
        #if particle == 'Po210' and edep > 2000:
        #if edep > 8000 and int(edep) != 14472:
        #if int(edep) == 0: continue
        """
        if particle == 'At215':
            print('[{0}] {1} {2} {3} keV'
                  .format(jentry,
                          particle,
                          decayType,
                          int(edep)))
        #if int(edep) > 8000: sys.exit()
        #if jentry > 6845: sys.exit()
        continue
        """
        #--------------------------------------------

        
        ###  primary volume cuts
        #--------------------------------------------------
        if location == 'internal':
            if volName != 'NaIDet0{0}Crystal'.format(fx): continue

        elif 'naisurf' in location:
            if volName != 'NaIDet0{0}Crystal'.format(fx): continue

        elif location == 'teflon':
            if volName != 'NaIDet0{0}Teflon'.format(fx): continue

        elif 'teflonsurf' in location:
            if volName != 'NaIDet0{0}Teflon'.format(fx): continue

        elif location == 'xpmt':
            pmt1 = int(fx)*2 - 2
            pmt2 = int(fx)*2 - 1
            pmtid = chain.GetLeaf('primPMTid').GetValue(0)
            if not (volName == 'phys_pmt'
                    and (pmtid == pmt1 or pmtid == pmt2)):
                continue
            
        elif location == 'pmt':
            if volName != 'phys_pmt': continue

        elif location == 'cucase':
            if not (volName == 'NaIDet0{0}Case'.format(fx)
                    or volName == 'NaIDet0{0}Fringe0'.format(fx)
                    or volName == 'NaIDet0{0}Fringe1'.format(fx)):
                continue

        elif 'cucasesurf' in location:
            if volName != 'NaIDet0{0}Case'.format(fx): continue

        elif location == 'plastic':
            if volName != 'PlasSupport': continue

        elif location == 'lsveto':
            if volName != 'lsveto': continue

        elif location == 'lsvetoair':
            if not (volName == 'DetPMTCover'
                    or volName == 'DetPMTEnvelope'):
                continue

        elif location == 'airshield':
            if volName != 'LSVetoAirRoom': continue

        elif location == 'cushield':
            if volName != 'CuShield': continue

        elif location == 'gamma':
            if volName != 'CuShield': continue

        elif location == 'innersteel':
            if volName != 'InnerSteel': continue

        elif location == 'steel':
            if not (volName == 'SteelSupport'
                    or volName == 'SteelSupportTop'):
                continue

        else:
            print('ERROR: No primary volume cut for --> {0}'.format(location))
            sys.exit()
            

        ### it's a valid generated event after the prim volume cut
        #--------------------------------------------------
        if evType <= 10:
            # fill generated events hist
            temp2.Fill(0)
            # the Tl208 gammas always have evt_Type == 3
            if location != 'gamma':
                continue
            
        
        ### group number cut
        #--------------------------------------------------
        gstart, gstop = groupNum500(isotope, True)
        if not (groupNo >= gstart and groupNo < gstop): continue
            
        
        ### alpha cut
        #--------------------------------------------------
        # alpha cuts for crystals
        if cx < 8:
            # alpha cuts only on internal/surface sources
            if isAlphaMC(info):
                if E == 2:
                    if decayType != 'alpha': continue
                else:
                    if decayType == 'alpha': continue
            else:
                # skip other events from contributing to the crystal alpha spectrum
                if E == 2: continue
        
        # show alphas in lsveto E=2 (just to see how it looks)
        if cx == 8 and E == 2:
            if decayType != 'alpha': continue
        
        
        ### energy cut
        #--------------------------------------------------
        threshold = 0.125 # in keV
        if edepRes[cx] < threshold: continue
        
        
        ### special crystal energy cuts for lsveto
        #--------------------------------------------------
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


        ### fill the histogram - processMCbyEvent500
        #--------------------------------------------------
        if cx < 8:
            if E == 2 and Q:
                # Th232 alpha-alpha pile-up Rn220 -> Po216 -> Pb212
                #  U235 alpha-alpha pile-up Rn219 -> Po215 -> Pb211
                if (location == 'internal' and particle in ['Po216', 'Po215']
                    and edep > getAlphaE(particle)+500):
                    a1 = getAlphaE(particle)
                    q1 = calcQuench(cx, Q, a1)
                    s1 = calcSigma(cx, E, q1)
                    e1 = random.gauss(q1, s1)
                    
                    a2 = edep - a1
                    q2 = calcQuench(cx, Q, a2)
                    s2 = calcSigma(cx, E, q2)
                    e2 = random.gauss(q2, s2)
                    
                    histo.Fill(e1)
                    histo.Fill(e2)
                    
                elif (cx in [0, 1, 2, 3, 5, 6] and location == 'internal'
                       and info['chst'] == 'Pb210'):
                    #print('INFO: applying skewed gaussian')
                    # skewed gaussian for internal Pb210 alphas
                    # func defined at the top of this function
                    histo.Fill(func.GetRandom())
                    
                else:
                    quenched = calcQuench(cx, Q, edep)
                    if quenched < 0: continue
                    sigma = calcSigma(cx, E, quenched)
                    energy = random.gauss(quenched, sigma)
                    histo.Fill(energy)
                    
            elif cx in [4, 7]:
                # tweaked resolutions for C5 and C8
                sigma = calcSigma(cx, E, edep)
                energy = random.gauss(edep, sigma)
                histo.Fill(energy)

            else:
                #  U238 beta-alpha pile-up Bi214 -> Po214 -> Pb210
                # Th232 beta-alpha pile-up Bi212 -> Po212 -> Pb208
                if location == 'internal' and particle in ['Po214', 'Po212']:
                    # assume some alpha energy loss due to
                    # daq waveform truncation
                    Etr = 500  # keV
                    alphaE = getAlphaE(particle)
                    if edep > alphaE:
                        betaE = edep - alphaE
                    else:
                        alphaE = 0
                        betaE = edep
                        
                    # quench and smear the alpha
                    #_Q = 2
                    _Q = Qcount%2 + 1
                    Qcount += 1
                    quenched = calcQuench(cx, _Q, alphaE) - Etr
                    if quenched > 0:
                        sigmaA = calcSigma(cx, E, quenched)
                        alphaE = random.gauss(quenched, sigmaA)
                    if betaE > 0:
                        # just smear the beta
                        sigmaB = calcSigma(cx, E, betaE)
                        betaE = random.gauss(betaE, sigmaB)

                    histo.Fill(alphaE+betaE)
                    
                else:
                    histo.Fill(edepRes[cx])

        # for LSveto 
        if cx == 8:
            if E == 2:
                # don't quench, just smear for E=2
                # just for reference
                histo.Fill(edepRes[cx])
                
            elif decayType == 'alpha':
                # quench and smear
                quenched = lsquench(edep)
                if quenched <= 0: continue
                sigma = calcSigma(cx, E, quenched)
                energy = random.gauss(quenched, sigma)
                histo.Fill(energy)
                
            else:
                # check for known pile-up events
                locas = ['lsveto', 'plastic']
                
                #  U238 beta-alpha pile-up Bi214 -> Po214 -> Pb210
                # Th232 beta-alpha pile-up Bi212 -> Po212 -> Pb208
                if location in locas and particle in ['Po214', 'Po212']:
                    alphaE = getAlphaE(particle)
                    if edep > alphaE:
                        betaE = edep - alphaE
                    else:
                        alphaE = 0
                        betaE = edep
                                            
                    # quench and smear the alpha
                    quenched = lsquench(alphaE)
                    if quenched > 0:
                        sigmaA = calcSigma(cx, E, quenched)
                        alphaE = random.gauss(quenched, sigmaA)
                    if betaE > 0:
                        # just smear the beta
                        sigmaB = calcSigma(cx, E, betaE)
                        betaE = random.gauss(betaE, sigmaB)

                    histo.Fill(alphaE+betaE)
                    
                else:
                    histo.Fill(edepRes[cx])

    #=====================================================================
    
    
    ### create the MC histogram
    #---------------------------------------------------------------------
    this_detected = histo.GetEntries()
    if 'hist' in mc[key]:
        mc[key]['hist'].Add(histo)
    else:
        mc[key]['hist'] = copy(histo)
    #---------------------------------------------------------------------
    
    ### create a hist of the generated events
    #---------------------------------------------------------------------
    this_generated = temp2.GetEntries()
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

    return mc


def processSurf502(chain, mc, info, C, E, fromx, Q=False):
    """
    Quench the alpha by Qval and then subtract energy loss
    and apply energy resolution.
    """
    
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

    print('INFO: processSurf502 \"{0}\"'.format(key))
    
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

    Qcount = 1
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
                else:
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
        
        

        ### fill the weighted histograms - processSurf502
        #--------------------------------------------------
        if cx < 8:
            if E == 2 and Q:
                # Th232 alpha-alpha pile-up Rn220 -> Po216 -> Pb212
                #  U235 alpha-alpha pile-up Rn219 -> Po215 -> Pb211
                if particle in ['Po216', 'Po215'] and edep > getAlphaE(particle)+500:
                    a1 = getAlphaE(particle)
                    q1 = calcQuench(cx, Q, a1)
                    s1 = calcSigma(cx, E, q1)
                    e1 = random.gauss(q1, s1)
                    
                    a2 = edep - a1
                    q2 = calcQuench(cx, Q, a2)
                    s2 = calcSigma(cx, E, q2)
                    e2 = random.gauss(q2, s2)
                    
                    histo.Fill(e1)
                    histo.Fill(e2)
                    
                else:
                    alpha = getAlphaE(particle)
                    if edep > alpha:
                        print('WARNING: [{0} {1}]  edep [{2}] > alpha [{3}]'
                              .format(particle, decayType, round(edep,1), round(alpha,1)))
                    quenched = calcQuench(cx, Q, alpha)
                    eloss = alpha - quenched
                    if eloss > edep:
                        #print('WARNING: [{0} {1}]  eloss [{2}] > edep [{3}]'
                        #      .format(particle, decayType, round(eloss,1), round(edep,1)))
                        continue
                    energy = edep - eloss
                    sigma = calcSigma(cx, E, energy)
                    energy = random.gauss(energy, sigma)
                    histo.Fill(energy, wtf)

            elif cx in [4, 7]:
                # tweaked resolutions for C5 and C8
                sigma = calcSigma(cx, E, edep)
                energy = random.gauss(edep, sigma)
                histo.Fill(energy, wtf)
            
            else:
                #  U238 beta-alpha pile-up Bi214 -> Po214 -> Pb210
                # Th232 beta-alpha pile-up Bi212 -> Po212 -> Pb208
                if particle in ['Po214', 'Po212']:
                    alpha = getAlphaE(particle)
                    #_Q = 2
                    _Q = Qcount%2 + 1
                    Qcount += 1
                    quenched = calcQuench(cx, _Q, alpha)
                    eloss = alpha - quenched
                    if eloss > edep:
                        #print('WARNING: [{0} {1}]  eloss [{2}] > edep [{3}]'
                        #      .format(particle, decayType, round(eloss,1), round(edep,1)))
                        energy = edep
                    else:
                        energy = edep - eloss
                    sigma = calcSigma(cx, E, energy)
                    energy = random.gauss(energy, sigma)
                    histo.Fill(energy, wtf)

                else:
                    histo.Fill(edepRes[cx], wtf)

        # for LSveto
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

    return mc


def getAlphaE(name):
    """
    return the alpha energy in keV
    """
    # make name lower case
    name = name.lower()

    # FIX - some of these are alpha energies
    # and some are Q-values
    
    ### U238 chain
    u238  = 4198
    u234  = 4774
    th230 = 4687
    ra226 = 4870
    rn222 = 5591
    po218 = 6115
    po214 = 7833
    bi210 = 5037
    po210 = 5408
    
    ### Th232 chain
    th232 = 4012
    th228 = 5521
    ra224 = 5790
    rn220 = 6405
    po216 = 6907
    bi212 = 6208
    po212 = 8784
    
    ### U235 chain
    u235  = 4395
    pa231 = 5150
    th227 = 6146
    ra223 = 5716
    rn219 = 6819
    po215 = 7386
    po211 = 7450


    # particle name is actually the daughter particle
    # using a dict here seems best
    a = {}
        
    # U238 chain
    a['th232'] = u238
    a['th230'] = u234
    a['ra226'] = th230
    a['rn222'] = ra226
    a['po218'] = rn220
    a['pb214'] = po218
    a['po214'] = po214 # special case for bi-po
    #a['bi214'] = at218
    #a['tl210'] = bi214
    a['pb210'] = po214
    #a['hg206'] = pb210
    a['tl206'] = bi210
    a['pb206'] = po210
    
    # Th232 chain
    a['ra228'] = th232
    a['ra224'] = th228
    a['rn220'] = ra224
    a['po216'] = rn220
    a['pb212'] = po216
    a['po212'] = po212 # special case for bi-po
    a['tl208'] = bi212
    a['pb208'] = po212

    # U235 chain
    a['th231'] = u235
    a['ac227'] = pa231
    #a['fr223'] = ac227
    #a['at219'] = fr223
    #a['bi215'] = at219
    a['ra223'] = th227
    a['rn219'] = ra223
    a['po215'] = rn219
    a['pb211'] = po215
    #a['bi211'] = at215
    #a['tl207'] = bi211
    a['pb207'] = po211
    
    
    # particle name is actually the daughter particle
    if name in a:
        return a[name]
    else:
        print('ERROR: ['+__name__+'/getAlphaE()]: \"{0}\" is not defined'.format(name))
        return 0

    
def groupNum500(isotope, numbers=False):

    # example isotope format = "Th228_GRND"
    
    # old U238 group numbers - I had this wrong!
    # ------------------------
    # 11: U238  -> Th230
    # 12: Th230 -> Ra226
    # 13: Ra226 -> Rn222
    # 14: Rn222 -> Pb210
    # 15: Pb210 -> ground
    # My U238 groups are wrong!!!
    # It doesn't effect anything downstream.
    # I should eventually fix my MC hist names.
    # So it's really the following:
    # ------------------------
    # 11: U238  -> U234    - 1 alpha
    # 12: U234  -> Th230   - 1 alpha
    # 13: Th230 -> Ra226   - 1 alpha
    # 14: Ra226 -> Pb210   - 4 alpha
    # 15: Pb210 -> ground  - 1 alpha
    
    # Th232 group numbers
    # ------------------------
    # 21: Th232 -> Ra228
    # 22: Ra228 -> Th228
    # 23: Th228 -> ground
    
    # U235 group numbers
    # ------------------------
    # 41: U235 -> Pa231
    # 42: Pa231 -> ground
    
    # others
    # ------------------------
    # 31: K40


    chstart = isotope.split('_')[0]
    chstop = isotope.split('_')[1]
    
    if chstart in ['Te121',  'Te121m',
                   'Te123m', 'Te125m', 'Te127m',
                   'I125',   'I126',   'I129',
                   'Na22',   'H3',
                   'I128',   'Na24',
                   'Cd109',  'Sn113',
                   'Co60',   'Fe55']:
        if numbers: return 0, 1
        else: return TCut('(groupNo == 0)')
    
    if chstart == 'K40':
        if numbers: return 31, 32
        else: return TCut('(groupNo == 31)')

    if chstart == 'Tl208':  # for the gamma MC
        if numbers: return 0, 1
        else: return TCut('(groupNo == 0)')
    
    if chstart == 'neutron':
        if numbers: return 0, 1
        else: return TCut('(groupNo == 0)')
    
    if   chstart == 'U238':  start = 11
    elif chstart == 'Th230': start = 12
    elif chstart == 'Ra226': start = 13
    elif chstart == 'Rn222': start = 14
    elif chstart == 'Pb210': start = 15
    elif chstart == 'Th232': start = 21
    elif chstart == 'Ra228': start = 22
    elif chstart == 'Th228': start = 23
    elif chstart == 'U235':  start = 41
    elif chstart == 'Pa231': start = 42
    else: start = -1
    
    #if   chstop == 'U238':  stop = 11 # should not be a stop group
    if   chstop == 'Th230': stop = 12
    elif chstop == 'Ra226': stop = 13
    elif chstop == 'Rn222': stop = 14
    elif chstop == 'Pb210': stop = 15
    #elif chstop == 'Th232': stop = 21 # should not be a stop group
    elif chstop == 'Ra228': stop = 22
    elif chstop == 'Th228': stop = 23
    #elif chstop == 'U235':  stop = 41 # should not be a stop group
    elif chstop == 'Pa231': stop = 42
    elif chstop == 'GRND':
        if   start in [11,12,13,14,15]: stop = 16
        elif start in [21,22,23]:       stop = 24
        elif start in [41,42]:          stop = 43
        else: stop = -1
    else: stop = -1

    if start != -1 and stop != -1:
        if numbers: return start, stop
        else: return TCut('((groupNo >= '+str(start)+') && (groupNo < '+str(stop)+'))')
    else:
        print('ERROR: group numbers not found for -->', isotope)
        sys.exit()


def mcPath500(info):

    """ New processing for LS veto single/multi hit """
    """ Uses version alias v4.0.3 """
    
    ### use the full location-names 'floca'
    loc = info['floca']
    iso = info['isof']
    
    if loc == 'internal':
        mcpath = '/data/COSINE/WORK/mkauer/processed/internal'
        filename = '*-internal-*'+iso+'[-|_]*root'
        if iso == 'Na22':
            filename = 'newGeant_*-internal-*'+iso+'[-|_]*root'

    elif 'nai-surf-expo' in loc or loc == 'nai-surface':
        Cx = 'C*'
        if int(info['xstl']) != 99:
            Cx = 'C'+str(info['xstl'])
        mcpath = '/data/COSINE/WORK/mkauer/processed/surface-nai'
        filename = 'kimssim-*-surface-nai-*'+Cx+'*10um*'+iso+'*evt100000.root'
    
    # special case for the "screened" surface Pb210
    elif 'nai-surfscr-expo' in loc:
        Cx = 'C*'
        if int(info['xstl']) != 99:
            Cx = 'C'+str(info['xstl'])
        mcpath = '/data/COSINE/WORK/mkauer/processed/surface-nai'
        filename = 'Screen_*-surface-nai-*'+Cx+'*10um*Pb210*evt100000.root'
    
    elif 'teflon-surf-expo' in loc or loc == 'teflon-surface':
        Cx = 'C*'
        if int(info['xstl']) != 99:
            Cx = 'C'+str(info['xstl'])
        mcpath = '/data/COSINE/WORK/mkauer/processed/surface-teflon'
        filename = '*-surface-tef-*'+Cx+'*10um*'+iso+'*evt100000.root'
        
    elif loc == 'teflon':
        Cx = 'C*'
        if int(info['xstl']) != 99:
            Cx = 'C'+str(info['xstl'])
        mcpath = '/data/COSINE/WORK/mkauer/processed/bulk-teflon'
        filename = '*-bulk-teflon-*'+Cx+'*'+iso+'[-|_]*root'
    
    elif loc == 'pmt' or loc == 'xpmt':
        mcpath = '/data/COSINE/WORK/mkauer/processed/pmt'
        filename = '*-pmt-*'+iso+'[-|_]*root'
    
    elif loc == 'pmtbase' or loc == 'xpmtbase':
        mcpath = '/data/COSINE/WORK/mkauer/processed/pmtbase'
        filename = '*-pmt_wBase-*'+iso+'[-|_]*root'
    
    elif loc == 'cucase':
        mcpath = '/data/COSINE/WORK/mkauer/processed/bulk-cuCase'
        #filename = '*-bulk-cuCase-*'+iso+'[-|_]*root'
        #filename = '*-bulk-[c|C]u[c|C]ase-*'+iso+'[-|_]*root'
        filename = '*-bulk[-|_][c|C]u[c|C]ase-*'+iso+'[-|_]*root'
    
    elif loc == 'cucase-surface':
        #Cx = 'C*'
        #if int(info['xstl']) != 99:
        #    Cx = 'C'+str(info['xstl'])
        mcpath = '/data/COSINE/WORK/mkauer/processed/surface-cuCase'
        #filename = '*-surf_cuCase-*'+Cx+'*'+iso+'[-|_]*root'
        filename = '*-surf_cuCase-*'+iso+'[-|_]*root'
    
    elif loc == 'plastic' or loc == 'xplastic':
        mcpath = '/data/COSINE/WORK/mkauer/processed/acrylic'
        filename = '*-acrylic-*'+iso+'[-|_]*root'
    
    elif loc == 'film':
        mcpath = '/data/COSINE/WORK/mkauer/processed/film'
        filename = '*-film-*'+iso+'[-|_]*root'
    
    elif loc == 'lsveto':
        mcpath = '/data/COSINE/WORK/mkauer/processed/lsveto'
        filename = '*-lsveto-*'+iso+'[-|_]*root'
    
    elif loc == 'vetopmt':
        mcpath = '/data/COSINE/WORK/mkauer/processed/vetopmt'
        filename = '*-vetopmt-*'+iso+'[-|_]*root'
    
    elif loc == 'cushield':
        mcpath = '/data/COSINE/WORK/mkauer/processed/cuShield'
        #filename = '*-cuShield-*'+iso+'[-|_]*root'
        filename = '*-[c|C]u[s|S]hield-*'+iso+'[-|_]*root'
    
    elif loc == 'gamma' or loc == 'xgamma':
        mcpath = '/data/COSINE/WORK/mkauer/processed/gamma'
        filename = '*-external-single-gamma-*Tl208*root'
    
    elif loc == 'innersteel':
        mcpath = '/data/COSINE/WORK/mkauer/processed/innersteel'
        filename = '*-[i|I]nner[s|S]teel-*'+iso+'[-|_]*root'
    
    elif loc == 'neutron':
        mcpath = '/data/COSINE/WORK/mkauer/processed/neutron'
        filename = '*-neutron*.root'
    
    
    else:
        print('ERROR: no file path identified for --> {0}'.format(loc))
        sys.exit()
    
    return mcpath, filename


def lsquench(energy_kev):
    # these pars are in MeV
    # so convert to MeV
    # then back to keV
    E = energy_kev/1000.
    # model 1
    #c0 = 0.030
    #c1 = 1.640
    #c2 = -0.518
    #c3 = 0.179
    # model 2 - use this one - 3 g/L PPO
    c0 = 0.031
    c1 = 1.689
    c2 = -0.505
    c3 = 0.190
    if E < 6.76:
        return (c0*(E**c1))*1000.
    else:
        return (c2+(c3*E))*1000.
    
    
def quenchPars(i, qi):
    # make sure qi is 1 or 2
    if qi not in [1, 2]:
        print('ERROR: ['+__name__+'/quenchPars()]: qi not 1 or 2')
        sys.exit()
        
    # pars in keV
    # pars for quench = A*x^2 + B*x + C
    A = 0.0
    B = 0.0
    C = 0.0

    cx = int(i)+1
    
    if cx == 1:
        ### 10
        #A = -2.0874404194288666e-05
        #B = 0.6677859857890859
        #if qi == 2:
        #    A = -1.1159676581564245e-05
        #    B = 0.5597742726432215

        ### 20
        #A = -2.0404821426955023e-05
        #B = 0.6679376941832087
        #if qi == 2: C = -120

        # using Po216 as end peak and Rn220
        ### 30
        #A = -4.7270556388344075e-05
        #B = 0.9056203621910291
        #if qi == 2: C = -350
        
        # Po216/Rn220 peaks from dt
        ### 40 - changed real energy of po216 and rn220
        #A = -1.0715729728966323e-05
        #B = 0.6386581087242886
        #if qi == 2: C = -300
        ### 40b - with po210
        #A = 6.840537176576379e-06
        #B = 0.5206123931686026
        #if qi == 2: C = -300

        # just use po216/po210 and get a q2
        ### 50
        A = 6.466093785181101e-06
        B = 0.5199832537120723
        if qi == 2:
            A = 1.4500536476867272e-05
            B = 0.42105516084934075


        
    elif cx == 2:
        ### 10
        #A = -1.3762044269888658e-05
        #B = 0.6663324984854302
        #if qi == 2:
        #    A = -2.007603340039448e-06
        #    B = 0.5519185864531456
        
        ### 20
        #A = -1.3953887014114313e-05
        #B = 0.6646435775254104
        #if qi == 2: C = -200
            
        # using Po216 as end peak and Rn220
        ### 30
        #A = -4.509296519986373e-05
        #B = 0.948492002773869
        #if qi == 2: C = -500

        # Po216/Rn220 peaks from dt
        ### 40 - changed real energy of po216 and rn220
        #A = -2.0054139371352584e-06
        #B = 0.6219301547413668
        #if qi == 2: C = -200
        ### 40b - with po210
        #A = 1.2292740478513289e-05
        #B = 0.5257915002304252
        #if qi == 2: C = -200
        
        # just use po216/po210 and get a q2
        ### 50
        A = 1.0755412068728351e-05
        B = 0.5337911295188668
        if qi == 2:
            A = 2.5435086463926825e-05
            B = 0.40633810015647776

            
        
    elif cx == 3:
        ### 10
        #A = -1.824430516211777e-05
        #B = 0.65356767042117
        #if qi == 2:
        #    A = -1.580774171910664e-05
        #    B = 0.6089490675677487

        ### 20
        #A = -1.1393453789159176e-05
        #B = 0.5999049516167854
        #if qi == 2: C = -200
        
        # using Po216 as end peak and Rn220
        ### 30
        #A = -3.5296303999967136e-05
        #B = 0.820018587798426
        #if qi == 2: C = -400
            
        # Po216/Rn220 peaks from dt
        ### 40 - changed real energy of po216 and rn220
        #A = 5.9085685798822485e-06
        #B = 0.5165950474398623
        #if qi == 2: C = -150
        ### 40b - with po210
        #A = 1.4426960023557496e-05
        #B = 0.4786785054076938
        #if qi == 2: C = -200
        
        # just use po216/po210 and get a q2
        ### 50
        A = 1.418258451425156e-05
        B = 0.4782679048307174
        if qi == 2:
            A = 1.584175134947229e-05
            B = 0.43785190803901103


            

    elif cx == 4:
        ### 10
        #A = -2.063107641965955e-05
        #B = 0.6901736697414937
        #if qi == 2:
        #    A = -1.5394252786912397e-05
        #    B = 0.6261366098853233

        ### 20
        #A = -2.9246524808955656e-05
        #B = 0.7550678577574402
        #if qi == 2: C = -180
        
        # using Po216 as end peak and Rn220
        ### 30
        #A = -4.086340667570415e-05
        #B = 0.8787398580023273
        #if qi == 2: C = -320
        
        # Po216/Rn220 peaks from dt
        ### 40 - changed real energy of po216 and rn220
        #A = -3.0956054063165277e-07
        #B = 0.5993583460339025
        #if qi == 2: C = -190
        ### 40b - with po210
        #A = 1.619701715637812e-05
        #B = 0.4883705834079194
        #if qi == 2: C = -160
        
        # just use po216/po210 and get a q2
        ### 50
        A = 1.1933368778710977e-05
        B = 0.5126247233656402
        if qi == 2:
            A = 1.964771841712402e-05
            B = 0.4376246118174937
        
            
    elif cx == 5:
        # just some hacky numbers to get a feel
        B = 1
        C = -2500
        if qi == 2:
            C = -3000

            
    elif cx == 6:
        ### 10
        #A = -3.431476487555885e-05
        #B = 0.73859477910965
        #if qi == 2:
        #    A = -2.6028055516483695e-05
        #    B = 0.6124057813296581

        ### 20
        #A = 1.2883305560047147e-05
        #B = 0.41485169591531523
        #if qi == 2: C = -370

        # using Po216 as end peak and Rn220
        ### 30
        #A = 1.026816717261261e-05
        #B = 0.5082003985555014
        #if qi == 2: C = -550
        
        # Po216/Rn220 peaks from dt
        ### 40 - changed real energy of po216 and rn220
        #A = 1.7445306699219593e-05
        #B = 0.4622474122778315
        #if qi == 2: C = -550
        ### 40b - with po210
        #A = 1.980583219026452e-05
        #B = 0.4463755914979328
        #if qi == 2: C = -550

        # just use po216/po210 and get a q2
        ### 50
        A = 1.9755486646877018e-05
        B = 0.4462909993793617
        if qi == 2:
            A = 2.0929161175544806e-05
            B = 0.358555067892552

            
        
    elif cx == 7:
        ### 10
        #A = -3.37887450691247e-05
        #B = 0.7357511160360669
        #if qi == 2:
        #    A = -3.128825358082553e-05
        #    B = 0.6408424120654901
        
        ### 20
        #A = -1.7986806374364305e-05
        #B = 0.6106978801697931
        #if qi == 2: C = -190
        
        # using Po216 as end peak and Rn220
        ### 30
        #A = -3.290178674206056e-05
        #B = 0.7600454599068099
        #if qi == 2: C = -320
        
        # Po216/Rn220 peaks from dt
        ### 40 - changed real energy of po216 and rn220
        #A = -2.8868763606224143e-05
        #B = 0.731465465821067
        #if qi == 2: C = -320
        
        # just use po216/po210 and get a q2
        ### 50
        A = 1.5414960611774907e-05
        B = 0.4697558831251236
        if qi == 2:
            A = 1.6588635140442662e-05
            B = 0.38201995163831415

            
    elif cx == 8:
        # just some hacky numbers to get a feel
        B = 1
        C = -2500
        if qi == 2:
            C = -2750

    
    else:
        # no quenching
        A = 0.0
        B = 1.0
        C = 0.0
        
    return A, B, C
        

def quench500(i, qi, EkeV):

    ### no quench
    #rootstring = '({0})'.format(EkeV)
    
    ### quenched = A*x^2 + B*x + C
    A, B, C = quenchPars(i, qi)
    rootstring = '(({1}*({0}**2)) + ({2}*{0}) + {3})'.format(EkeV, A, B, C)
        
    return rootstring


def calcQuench(i, qi, EkeV):
    A, B, C = quenchPars(i, qi)
    quench = A*(EkeV**2) + B*EkeV + C
    return quench


def resolPars(i, E):
    # parameters in keV
    i = int(i)
    E = int(E)
    
    # add p2 for alphas
    # res = p0*sqrt(x) + p1*x + p2

    # low energy
    loEresol = [
        [0.2413,  0.01799],   # C1
	    [0.2951,  0.01427],   # C2
	    [0.3106,  0.007894],  # C3
	    [0.3894, -0.001437],  # C4
        #[0.2413,  0.018],     # C5
        [0, 9.5/50.],         # C5 my tweak
	    [0.3620,  0.0006355], # C6
	    [0.3042,  0.009784],  # C7
        #[0.2413,  0.018],     # C8
        [0, 14.8/50.],        # C8 my tweak
        [0, 0.10]             # C9
    ]
    
    # high energy
    hiEresol = [
        [0.6729,  0.009374],  # C1
	    [0.6531,  0.006627],  # C2
	    [0.5926,  0.009506],  # C3
	    [0.7227,  0.004790],  # C4
        #[0.6729,  0.009374],  # C5
        [14.73/np.sqrt(50.), 0], # C5 my tweak
	    [0.6498,  0.009670],  # C6
	    [0.7034,  0.007812],  # C7
        #[0.6729,  0.009374],  # C8
        [24.76/np.sqrt(50.), 0], # C8 my tweak
        [0.0,  0.10]          # C9
        #[0.0,  0.08]          # C9 - this looks better for K40?
        
    ]

    # this is in keV
    alphaRes = [
        90,  # C1
        60,  # C2
        65,  # C3
        55,  # C4
        50,  # C5
        60,  # C6
        60,  # C7
        50.  # C8
    ]        
    
    p0 = 0
    p1 = 0
    p2 = 0
    if E == 2 and i < 8:
        p2 = alphaRes[i]
    elif E == 2 and i == 8:
        p0, p1 = hiEresol[i]
    elif E == 1:
        p0, p1 = hiEresol[i]
    else:
        p0, p1 = loEresol[i]
    pars = [p0, p1, p2]

    
    # apply to high-E and Alphas
    if i == 4 and E in [1, 2]:  # C5
        p0 = 2.7069054846510676e-05
        p1 = 0.0021016264539907713
        p2 = 8.824835006951153
        p3 = "poly2"
        pars = [p0, p1, p2, p3]

        
    # apply to high-E and Alphas    
    if i == 7 and E in [1, 2]:  # C8
        p0 = 3.7465588771707855e-05
        p1 = -0.014361284301020598
        p2 = 20.62969823943649
        p3 = "poly2"
        pars = [p0, p1, p2, p3]

        
    # lsveto new resolution
    # applies to all energy ranges
    if i == 8:  # LS
        """
        # r1 test3
        p0 = 3.1278054713068997
        p1 = 0.07998457360806746
        p2 = 1.2621210390538606e-05
        p3 = -2.0842113351621527e-09
        p4 = -0.38024822482886694
        p5 = "test3"
        pars = [p0, p1, p2, p3, p4, p5]
        """
        """
        # r1 erfgaus
        p0 = 358.7781853584133
        p1 = 3555.935705262654
        p2 = 131.8230686084652
        p3 = 3420.6505551866603
        p4 = 712.7126372026862
        p5 = "erfgaus"
        pars = [p0, p1, p2, p3, p4, p5]
        """
        """
        # r2 erfgaus
        p0 = 358.5906491238849
        p1 = 3543.016119392098
        p2 = 67.45148913460137
        p3 = 2599.9999999999995
        p4 = 456.1824475712886
        p5 = "erfgaus"
        pars = [p0, p1, p2, p3, p4, p5]
        """
        # r3 erfgaus
        p0 = 353.5132152826663
        p1 = 3102.796150076346
        p2 = 48.61156415597338
        p3 = 2599.999999999936
        p4 = 517.8819636291101
        p5 = "erfgaus"
        pars = [p0, p1, p2, p3, p4, p5]

        
    return pars


def sigma500(i, E, EkeV):
    i = int(i)
    E = int(E)
    # add p2 offset for alphas
    # res = p0*sqrt(x) + p1*x + p2
    p = resolPars(i, E)

    if isinstance(p[-1], str):
        if p[-1] == 'poly2':
            rootstring = ('( {1}*{0}**2 + {2}*{0} + {3} )'
                          .format(EkeV, p[0], p[1], p[2]))
        elif p[-1] == 'erfgaus':
            rootstring = ('( ({1}*(1+ROOT::Math::erf({0}/({2}*sqrt(2))))-{1}) + '
                          '({3}*exp(-0.5*(({0}-{4})/({5}))**2)) )'
                          .format(EkeV, p[0], p[1], p[2], p[3], p[4]))
        else:
            print('ERROR: sigma500() : unknown func name \"{0}\"'.format(p[-1]))
            sys.exit()
            
    else:
        rootstring = ('((({1}*sqrt({0})) + ({2}*({0})) + ({3})))'
                      .format(EkeV, p[0], p[1], p[2]))
        
    return rootstring


def calcSigma(i, E, EkeV):
    i = int(i)
    E = int(E)
    # return calculated value instead of rootstring
    p = resolPars(i, E)

    if isinstance(p[-1], str):
        if p[-1] == 'poly2':
            sigma = p[0]*EkeV**2 + p[1]*EkeV + p[2]
        elif p[-1] == 'erfgaus':
            sigma = ( (p[0]*(1+erf(EkeV/(p[1]*sqrt(2))))-p[0])
                      + (p[2]*exp(-0.5*((EkeV-p[3])/(p[4]))**2)) )
        else:
            print('ERROR: calcSigma() : unknown func name \"{0}\"'.format(p[-1]))
            sys.exit()
            
    else:
        sigma = p[0]*np.sqrt(EkeV) + p[1]*EkeV + p[2]
        
    return sigma

    
def dataToCluster(data, rootfile, info, jsonfile=None, cluster=False):

    #+++++++++++++++++++++++++++++++++++++++++++++++++
    # some speed ups for local testing
    # turn off for the real processing
    skip_alphas       = 0
    skip_hiEnergy     = 0
    skip_loEnergy     = 0
    skip_byEvent      = 1  # no longer need to recalib low-energy by event
    if skip_alphas:     print('WARNING: skip processing alphas')
    if skip_hiEnergy:   print('WARNING: skip processing high energy')
    if skip_loEnergy:   print('WARNING: skip processing low energy')
    if skip_byEvent:    print('WARNING: skip processing by event')
    #+++++++++++++++++++++++++++++++++++++++++++++++++
    
    if info is None and not os.path.exists(jsonfile):
        print('WARNING: you must specify an info json filename')
        sys.exit()

    if info is None:
        with open(jsonfile) as jf:
            info = json.load(jf)

    if data is None:
        data = {}

    if cluster:
        location = os.path.join('data', info['build'], info['tag'], 'x'+str(info['xstl']))
        outdir = os.path.join('/data/COSINE/WORK/mkauer/histograms', location)
        mkdir(outdir)
        outfile = info['key']+'__'+os.path.basename(rootfile)
        filename = os.path.join(outdir, outfile)
        tfile = TFile(filename, 'RECREATE')

    chain = TChain("ntp", "")
    chain.Add(rootfile)

    runtime = getDuration(rootfile)
    info['runtime'] = runtime
    info['rootfile'] = rootfile

    # get the run number
    # like mrgd_M001858.root.010
    runnum = int(os.path.basename(rootfile)[6:12])
    #print(runnum)

    for C in ['S', 'M']:
        for E in [0, 1, 2]:

            # skip low-energy?
            if E == 0 and skip_loEnergy:
                continue
            
            # skip high-energy?
            if E == 1 and skip_hiEnergy:
                continue
            
            # skip alphas?
            if E == 2 and skip_alphas:
                continue
            
            # if special case '99' then process all crystals at once
            if int(info['xstl']) == 99:
                for xstl in range(1, 10):
                    basekey = info['key'].replace('x99-', 'x'+str(xstl)+'-', 1)
                    if xstl == 9 and C == 'S':
                        # runs < 1800 don't have LS single-hit
                        if runnum < 1800: continue
                        chain2 = TChain("activeLS","")
                        chain2.Add(rootfile)
                        data = processData500(chain2, data, info, basekey, C, E)
                    elif xstl == 9 and C == 'M':
                        data = processData500(chain, data, info, basekey, C, E)
                    else:
                        if E == 0 and xstl not in [5, 8]:
                            if skip_byEvent:
                                data = processData500(chain, data, info, basekey, C, E)
                            else:
                                data = processDataByEvent501(chain, data, info, basekey, C, E)
                        else:
                            data = processData500(chain, data, info, basekey, C, E)
            else:
                basekey = info['key']
                if int(info['xstl']) == 9 and C == 'S':
                    # runs < 1800 don't have LS single-hit
                    if runnum < 1800: continue
                    chain2 = TChain("activeLS","")
                    chain2.Add(rootfile)
                    data = processData500(chain2, data, info, basekey, C, E)
                elif int(info['xstl']) == 9 and C == 'M':
                    data = processData500(chain, data, info, basekey, C, E)
                else:
                    if E == 0 and int(info['xstl']) not in [5, 8]:
                        if skip_byEvent:
                            data = processData500(chain, data, info, basekey, C, E)
                        else:
                            data = processDataByEvent501(chain, data, info, basekey, C, E)
                    else:
                        data = processData500(chain, data, info, basekey, C, E)
            
    if cluster:
        print('Writing rootfile --> {0}'.format(filename))
        writeHists500(data)
        tfile.Write()
        tfile.Close()

    #del chain
    #del chain2
    return data


def processData500(chain, data, info, basekey, C, E):

    # get raw data in adc
    rawdata = False
    # don't recalibrate
    norecal = False
    
    bits = basekey.split('-')
    Cx = int(bits[0][1])
    i = int(bits[0][1]) - 1
    rootfile = info['rootfile']
    runtime = float(info['runtime'])
    tag = info['tag']
    build = info['build']
    
    # DEFINE HIST AND KEY
    #-----------------------------------------------------------------------
    key  = basekey
    key += '-c'+str(C)
    key += '-e'+str(E)

    print('INFO: processData500 \"{0}\"'.format(key))
    
    if key not in data:
        data[key] = {}
        
    if 'info' not in data[key]:
        data[key]['info'] = copy(info)
    
    pars = histparams(E)
    data[key]['pars'] = pars
    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
    #-----------------------------------------------------------------------


    # Calibrations and Cuts
    #-----------------------------------------------------------------------
    if build == 'V00-04-04':
        ### calib100 uses Pushpa's tweaked low energy calibrations
        #edep, selection = calib100(i, E)
        ### calib301 uses the standard production calibrations
        ### and uses my tweaked lsveto calibration
        edep, selection = calib301(i, E)
        masterCut = cutsBDT300(i, C, E, edep, selection)

    elif build == 'V00-04-12':
        ### calib300 uses Govinda's tweaked high energy calibrations
        ### and uses my tweaked lsveto calibration
        #edep, selection = calib300(i, E)
        ### calib301 uses the standard production calibrations
        ### and uses my tweaked lsveto calibration
        edep, selection = calib301(i, E)
        masterCut = cutsBDT300(i, C, E, edep, selection)

    elif build == 'V00-04-14':
        ### calib301 uses the standard production calibrations
        ### and uses my tweaked lsveto calibration
        edep, selection = calib301(i, E)
        masterCut = cutsBDT301(i, C, E, edep, selection)

    #elif build == 'V00-04-15':
    elif build in ['V00-04-15', 'V00-04-19', 'V00-04-20']:
        ### calib301 uses the standard production calibrations
        ### and uses the lsveto poly43 calibration
        #edep, selection = calib301(i, E)
        ### calib302 can have different calibs for lsveto "S" and "M"
        ### and tweaked high energy calibrations of the crystals
        selection = calib302(i, E, C, (rawdata or norecal))
        masterCut = cutsBDT302(i, C, E)
        if rawdata:
            if E == 0:
                selection = '(crystal{0}.qc1_5)'.format(Cx)
            else:
                selection = '(crystal{0}.rqcD1_5)'.format(Cx)
            if Cx == 9:
                selection = '(BLSVeto.Charge)'
    else:
        print('WARNING: no selection for build \"{0}\"'.format(build))
        return data
    #-----------------------------------------------------------------------


    # FILL HISTOS
    #-----------------------------------------------------------------------
    chain.Draw(selection+' >> '+key, masterCut)

    # sanity check on the data
    if histo.GetEntries() == 0:
        # Don't warn about crystal multi-hit alpha channel
        if not (Cx < 9 and C == 'M' and E == 2):
            print('WARNING: no events found for [{0}] --> \"{1}\"'
                  .format(key, os.path.basename(rootfile)))
        # the new production 04-19 doesn't have activeLS events for *root.000
        # so skip saving data hist and runtime for these
        if Cx == 9 and C == 'S':
            return data

    # save histo to data
    if 'hist' in data[key]:
        data[key]['hist'].Add(histo)
    else:
        data[key]['hist'] = copy(histo)

    # save the total runtime
    key3 = key+'_runtime'
    temp3 = TH1F(key3, 'runtime', 1, 0, 1)
    temp3.SetBinContent(1, runtime)
    temp3.SetEntries(runtime)
    if 'runtime_hist' in data[key]:
        data[key]['runtime_hist'].Add(temp3)
    else:
        data[key]['runtime_hist'] = copy(temp3)
    data[key]['runtime'] = data[key]['runtime_hist'].GetEntries()
    #-----------------------------------------------------------------------

    return data


def processDataByEvent501(chain, data, info, basekey, C, E):
    
    bits = basekey.split('-')
    Cx = int(bits[0][1])
    i = int(bits[0][1]) - 1
    rootfile = info['rootfile']
    runtime = float(info['runtime'])
    tag = info['tag']
    build = info['build']
    
    # DEFINE HIST AND KEY
    #-----------------------------------------------------------------------
    key  = basekey
    key += '-c'+str(C)
    key += '-e'+str(E)

    print('INFO: processDataByEvent501 \"{0}\"'.format(key))
    
    if key not in data:
        data[key] = {}
        
    if 'info' not in data[key]:
        data[key]['info'] = copy(info)
    
    pars = histparams(E)
    data[key]['pars'] = pars
    
    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
    #-----------------------------------------------------------------------

    
    entries = chain.GetEntries()

    for jentry in range(entries):

        chain.GetEntry(jentry)

        if Cx == 9:
            energy = (chain.GetBranch('BLSVeto').GetLeaf('Charge').GetValue()) / 143.8
            #print(energy)
        elif E:
            energy = chain.GetBranch('crystal{0}'.format(Cx)).GetLeaf('energyD').GetValue()
        else:
            energy = chain.GetBranch('crystal{0}'.format(Cx)).GetLeaf('energy').GetValue()

        # event selection
        if not dataCutsByEvent501(chain, i, C, E):
            continue

        # energy calibration
        energy = dataCalibByEvent501(i, E, energy)
        
        #print(energy)
        histo.Fill(energy)
        


    # save histo and runtime to data
    #-----------------------------------------------------------------------
    # sanity check on the data
    if histo.GetEntries() == 0:
        # Don't warn about crystal multi-hit alpha channel
        if not (Cx < 9 and C == 'M' and E == 2):
            print('WARNING: no events found for [{0}] --> \"{1}\"'
                  .format(key, os.path.basename(rootfile)))
        # the new production 04-19 doesn't have activeLS events for *root.000
        # so skip saving data hist and runtime for these
        if Cx == 9 and C == 'S':
            return data
    
    # save the histo to data
    if 'hist' in data[key]:
        data[key]['hist'].Add(histo)
    else:
        data[key]['hist'] = copy(histo)
        
    # save the total runtime
    key3 = key+'_runtime'
    temp3 = TH1F(key3, 'runtime', 1, 0, 1)
    temp3.SetBinContent(1, runtime)
    temp3.SetEntries(runtime)
    if 'runtime_hist' in data[key]:
        data[key]['runtime_hist'].Add(temp3)
    else:
        data[key]['runtime_hist'] = copy(temp3)
    data[key]['runtime'] = data[key]['runtime_hist'].GetEntries()
    #-----------------------------------------------------------------------

    return data


def writeHists500(hist_dict):
    for key in hist_dict:
        for key2 in hist_dict[key]:
            ## key is like x1-teflon-Pb210_GRND-f1-cS-e0
            ## for data, key is like x1-data-SET3-V000419-cS-e0
            ## key2 is either info, hist, pars, generated_hist, generated, scale
            ## for data, key2 is either info, hist, pars, runtime_hist, runtime, druScale
            
            #if '-e0' not in key and '-e1' not in key: continue
            #if '-data-' in key: continue
            #if 'x9-' in key: continue
            #if key2 not in ['hist']: continue
            
            try:
                hist_dict[key][key2].Write()
            except:
                pass
            
    return


def cleanSmallBins(hist):
    # very small counts in bins causes nan/inf
    # when doing data/total so just set to 0
    for n in range(hist.GetNbinsX()+1):
        if hist.GetBinContent(n) == 0.0:
            continue
        elif hist.GetBinContent(n) < 1.e-30:
            #print('WARNING: found small bin of \"{0}\" at bin {1}'.format(hist.GetBinContent(n), n))
            hist.SetBinContent(n, 0)
        else:
            continue
        
    return hist


def combineOthers500(bkgs, globalMC, skip_locas=[], keep_isos=[], debug=0):

    #debug = 1
    
    if debug: print('DEBUG: combining histograms')

            
    donekeys = []
    delete = []
    lsveto = {}
    
    for key in sorted(bkgs):

        if key in donekeys:
            if debug: print('hist already done [{0}]'.format(key))
            continue

        if debug: print('working on [{0}]'.format(key))

        bits = key.split('-')

        # don't combine globals
        if bits[1] in globalMC:
            if debug: print('global mc, not combining [{0}]'.format(key))
            continue
        
        X = int(bits[0][1])
        loca = bits[1]
        iso = bits[2]
        
        for F in range(1, numX()+1):
            if F == X:
                continue
            
            if len(bits) == 6:
                default = 'x'+str(X)+'-'+bits[1]+'-'+bits[2]+'-f'+str(X)+'-'+bits[4]+'-'+bits[5]
                newkey  = 'x'+str(X)+'-'+bits[1]+'-'+bits[2]+'-f'+str(F)+'-'+bits[4]+'-'+bits[5]
            elif len(bits) == 7:
                default = 'x'+str(X)+'-'+bits[1]+'-'+bits[2]+'-f'+str(X)+'-'+bits[4]+'-'+bits[5]+'-'+bits[6]
                newkey  = 'x'+str(X)+'-'+bits[1]+'-'+bits[2]+'-f'+str(F)+'-'+bits[4]+'-'+bits[5]+'-'+bits[6]
            else:
                print('WARNING: do not know how to combine key [{0}]'.format(key))
                
            if newkey not in bkgs:
                if debug: print('key not found [{0}]'.format(newkey))
                continue
            
            # Don't combine others for some locations
            # but keep others for some isotopes
            if loca in skip_locas and iso not in keep_isos and X != 9:
                if debug:
                    print('skipping [{0}]'.format(newkey))
                pass
            else:
                if X == 9:
                    if debug: print('adding [{0}] --> [{1}]'.format(newkey, default))
                    if default not in lsveto:
                        lsveto[default] = deepcopy(bkgs[newkey])
                    else:
                        lsveto[default]['hist'].Add(bkgs[newkey]['hist'])
                else:
                    if debug: print('adding [{0}] --> [{1}]'.format(newkey, default))
                    bkgs[default]['hist'].Add(bkgs[newkey]['hist'])
            
            donekeys.append(default)
            donekeys.append(newkey)
            delete.append(newkey)

    # transfer new lsveto keys into bkgs
    for lskey in sorted(lsveto):
        if debug: print('copy new lsveto key/hist to backgrounds [{0}]'.format(lskey))
        bkgs[lskey] = deepcopy(lsveto[lskey])
        
    # delete the other histograms
    for dkey in sorted(delete):
        if debug: print('deleting [{0}]'.format(dkey))
        del bkgs[dkey]

    return bkgs


def updateBkgsFile500(xstals, bkgsfile, resultsfile, newdir, BF='BR'):
    
    if os.path.exists(bkgsfile):
        with open(bkgsfile) as fbkgs:
            bkgslines = fbkgs.read().splitlines()
    else:
        print('WARNING: file not found -->', bkgsfile)
        return

    results = []
    
    if not os.path.exists(resultsfile):
        return
    else:
        with open(resultsfile) as ffits:
            fitlines = ffits.read().splitlines()

        for fline in fitlines:
            if not fline: continue
            fbits = re.split("[ \s\t\n\r]+", fline.strip())
            if fbits[0] == 'fit':
                xstal = fbits[1].split('-')[0].split('x')[1]
                loca  = fbits[1].split('-')[1]
                chst  = fbits[1].split('-')[2].split('_')[0]
                chsp  = fbits[1].split('-')[2].split('_')[1]
                acti  = str(fbits[3])
                results.append({
                    'xstal': xstal,
                    'loca': loca,
                    'chst': chst,
                    'chsp': chsp,
                    'acti': acti
                })
    
    #newbkgs = os.path.join(newdir, bkgsfile.split('/')[-1][:-4]+'_update.txt')
    newbkgs = os.path.join(newdir, bkgsfile.replace('.txt', '_update.txt'))
    output = open(newbkgs, 'w')
    
    print('INFO: Updating bkgsfile --> {0}'.format(bkgsfile))
    print('      To a new bkgsfile --> {0}'.format(newbkgs))
    
    skip = 0
    for bline in bkgslines:
        bline = bline.strip()
        if not bline:
            output.write('\n')
            continue
        if bline.startswith('\"\"\"'):
            if skip == 0: skip = 1
            else: skip = 0
            output.write(bline+'\n')
            continue
        if skip:
            output.write(bline+'\n')
            continue
        if bline.startswith('#'):
            output.write(bline+'\n')
            continue

        bbits = re.split("[ \s\t\n\r]+", bline)
        
        if 'F' not in bbits[0]:
            output.write(bline+'\n')
            continue
        
        xstal = bbits[1]
        loca = bbits[2].replace('-','')
        chst = bbits[4]
        chsp = bbits[5]
        acti = bbits[6]

        found = False
        for result in results:
            if found: break
            if (result['xstal'] == xstal and
                result['loca'] == loca and
                result['chst'] == chst and
                result['chsp'] == chsp):
                bline = bline.replace('F', 'B', 1)
                bline = bline.replace(acti, result['acti'], 1)
                output.write(bline+'\n')
                print('UPDATED: {0}'.format(bline))
                found = True

        if not found:
            #print('WARNING: no fit result found for x{0} {1} {2}_{3}'
            #      .format(xstal, loca, chst, chsp))
            output.write(bline+'\n')
    
    
    output.close()
    return


    
if __name__ == "__main__":
    main(sys.argv[1:])

