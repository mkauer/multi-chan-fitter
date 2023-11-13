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
#import math
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
#from funcs411 import *
from funcs_old import *
from funcs_misc import *

# skip processing rootfiles that already exist
SKIP_EXISTING = False

DEBUG = 1


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
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    # reuse a joined rootfile
    if 'R' in info['type']:
        info['reuse'] = 1
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
    # channels to process
    info['chans'] = fchans
    # energies to process
    info['energies'] = [0, 1, 2]
    # crystals to process
    info['xstals'] = sorted(fxstals)
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
        ### location is 'data'
        info['loca'] = str(bits[2])

        ### data set
        info['tag'] = str(bits[3])
        
        ### production version
        info['build'] = str(bits[4])

        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-'+str(info['loca'])
        key += '-'+str(info['tag'])
        key += '-'+str(info['build']).replace('-','')
        info['key'] = key
        
    else:
        ### background location as displayed in config file
        info['floca'] = str(bits[2])
        
        # background location without '-' for hist keys
        info['loca'] = str(bits[2]).replace('-','')
        
        ### top level isotope file name
        info['isof'] = str(bits[3])

        ### isotope chain break start
        info['chst'] = str(bits[4])

        ### isotope chain break stop
        info['chsp'] = str(bits[5])

        ### chain start_stop name
        info['chain'] = info['chst']+'_'+info['chsp']
        
        ### activity
        info['acti'] = float(bits[6])

        ### fit bounds
        info['fbnd'] = [float(bits[7]), float(bits[8])]

        ### simulation version
        info['sim'] = str(bits[9])
        
        ### build the histo key
        key  = 'x'+str(info['xstl'])
        key += '-'+info['loca']
        key += '-'+info['chain']
        info['key'] = key

        ### assign a plotting group to the isotope-location
        info['group'] = setGroup500(info)
        
    return info

    
def build500(infile, others=1, freuse=0, fchans=['S','M'], fxstals=[], cluster=False):

    debug = 0
    
    data = {}
    bkgs = {}
    sigs = {}
    
    for line in readFile(infile):

        #print('processing line --> {0}'.format(line))

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
        print('WARNING: Rootfile not found --> {0}'.format(info['rootfile']))
        return data
    
    for C in ['S', 'M']:
        for E in [0, 1, 2]:
            
            key  = info['key']
            key += '-c'+str(C)
            key += '-e'+str(E)

            try:
                _dump = TH1F(tfile.Get(key))
            except:
                print('WARNING: Hist not found --> {0}'.format(key))
                return data
            
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
                pars = histparam6000(E)
                data[key]['pars'] = pars
                hist = TH1F(key, longNames(int(info['xstl'])-1), pars[0], pars[1], pars[2])
                data[key]['hist'] = deepcopy(hist)
                pass
            """
            
    return data


def reuseMC500(info, mc):

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
        for xstl in pcxs:
            for C in ['S', 'M']:
                # only process alphas for lsveto
                if xstl == 9:
                    penergies = [0, 1, 2]
                else:
                    penergies = [0, 1]
                for E in penergies:
                    newinfo = copy(info)
                    newinfo['xstl'] = xstl
                    newinfo['fromx'] = 0
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
                for E in [0, 1, 2]:
                    
                    # apply 2 quenchings to alphas
                    if E == 2 and isAlphaMC(info):
                        qis = [1, 2]
                    # or set to false
                    else:
                        qis = [False]
                            
                    for qi in qis:
                        # process crystal hits into crystal
                        newinfo = copy(info)
                        newinfo['xstl'] = xstl
                        newinfo['fromx'] = fromx
                        if qi: newinfo['acti'] = newinfo['acti']/2.
                        if is99: newinfo['key'] = newinfo['key'].replace('x99-', 'x'+str(xstl)+'-', 1)
                        mc = getMC500(tfile, mc, newinfo, C, E, fromx, qi)
                        
                        # process crystal hits into lsveto
                        newinfo = copy(info)
                        newinfo['xstl'] = 9
                        newinfo['fromx'] = fromx
                        if qi: newinfo['acti'] = newinfo['acti']/2.
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
        if info['build'] == 'V00-04-15':
            directory = '/data/COSINE/MRGD/phys/'+info['build']
        elif info['build'] in ['V00-04-16', 'V00-04-17']:
            directory = '/data/COSINE/MRGD/'+info['build']
        else:
            print('ERROR: unknown data production [{0}]'.format(info['build']))
            sys.exit()
        #match = '*root*'
        match = 'mrgd_M??????.root.???'
    else:
        print('ERROR: unknown data/sim type --> \"{0}\"'.format(filetype))
        sys.exit()
    
    ### append the path if running locally 
    if amLocal():
        directory = '/home/mkauer/COSINE/CUP/mc-fitting'+directory
        
    ### find all the matching files
    filenames = sorted(glob.glob(directory+'/'+match))
    if len(filenames) == 0:
        print('WARNING: no files found for --> {0}/{1}'.format(directory, match))
        return hists
    else:
        print('INFO: found {0} files'.format(len(filenames)))

        
    if cluster:

        # do not reprocess existing files
        if SKIP_EXISTING:
            newfiles = []
            for rootfile in filenames:
                if 'naisurfexpo' in info['loca']:
                    savedir = 'nai-surf-expo'
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
            filenames = newfiles

        print('INFO: submitting [{0}] jobs to cluster'.format(len(filenames)))
        
        # submit each rootfile as a separate job
        for rootfile in filenames:
            
            # wait until there's room in the queue
            while checkQueue() > 990:
                print('sleeping...')
                time.sleep(60)
            
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
                    
    else:
        # process localy
        for rootfile in filenames:
            print('Processing file --> {0}'.format(rootfile))
            if filetype == 'data':
                hists = dataToCluster(hists, rootfile, info)
            else:
                hists = simToCluster(hists, rootfile, info)

    return hists


def simToCluster(mc, rootfile, info, jsonfile=None, cluster=False):

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
                # only process alphas for lsveto
                if xstl == 9:
                    penergies = [0, 1, 2]
                else:
                    penergies = [0, 1]
                for E in penergies:
                    newinfo = copy(info)
                    newinfo['xstl'] = xstl
                    newinfo['fromx'] = 0
                    newinfo['key'] = newinfo['key'].replace('x0-', 'x'+str(xstl)+'-', 1)
                    # process lsveto data by event to quench all alphas
                    # and only process mc with alphas this way
                    if xstl == 9 and info['isof'] in ['U238', 'Th232', 'U235']:
                        mc = processMCbyEvent500(chain, mc, newinfo, C, E)
                    else:
                        mc = processMC500(chain, mc, newinfo, C, E)

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
            
            ## just for quicker laptop testing ##
            #====================================
            #if xstl != fromx: continue
            #====================================
            
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
                    
                    # apply 2 quenchings to alphas or false
                    if E == 2 and isAlphaMC(info):
                        qis = [1, 2]
                    else:
                        qis = [False]
                            
                    for qi in qis:
                        # process crystal hits into crystal
                        newinfo = copy(info)
                        newinfo['xstl'] = xstl
                        newinfo['fromx'] = fromx
                        if qi: newinfo['acti'] = newinfo['acti']/2.
                        if is99: newinfo['key'] = newinfo['key'].replace('x99-', 'x'+str(xstl)+'-', 1)
                        if isExpoSurfMC(info):
                            mc = processSurf500(chain, mc, newinfo, C, E, fromx, qi)
                        else:
                            mc = processMC500(chain, mc, newinfo, C, E, fromx, qi)
                            #mc = processMCbyEvent500(chain, mc, newinfo, C, E, fromx, qi)
                            
                        # process crystal hits into lsveto
                        newinfo = copy(info)
                        newinfo['xstl'] = 9
                        newinfo['fromx'] = fromx
                        if qi: newinfo['acti'] = newinfo['acti']/2.
                        if is99: newinfo['key'] = newinfo['key'].replace('x99-', 'x9-', 1)
                        else: newinfo['key'] = newinfo['key'].replace('x'+str(xstl)+'-', 'x9-', 1)
                        if isExpoSurfMC(info):
                            mc = processSurf500(chain, mc, newinfo, C, E, fromx, qi)
                        else:
                            mc = processMC500(chain, mc, newinfo, C, E, fromx, qi)
    
    if cluster:
        if 'naisurfexpo' in info['loca']:
            savedir = 'nai-surf-expo'
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
        
    #del chain
    return mc


def isExpoSurfMC(info):
    if 'naisurfexpo' in info['loca'] \
       or 'teflonsurfexpo' in info['loca']:
        return 1
    else:
        return 0

    
def isAlphaMC(info):
    if ('internal' == info['loca'] \
       or 'teflon' == info['loca'] \
       or 'naisurf' in info['loca'] \
       or 'teflonsurf' in info['loca']) \
       and \
       (info['isof'] in ['U238', 'Th232', 'Pb210', 'U235']):
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

    print('INFO: processing \"{0}\"'.format(key))
    
    if key not in mc:
        mc[key] = {}

    if 'info' not in mc[key]:
        mc[key]['info'] = copy(info)

    pars = histparam6000(E)
    mc[key]['pars'] = pars
    histo = TH1F(key, longNames(i), pars[0], pars[1], pars[2])
    #-----------------------------------------------------------------------


    ###  put together all the cuts you need...
    #=====================================================================================

    ### define the edep resolution to use
    if E: edepRes = 'edepResolD'
    else: edepRes = 'edepResolA'

    ### define the energy cut
    energyCut = TCut('('+edepRes+'['+str(i)+']*1000. > 0.0)')
    
    ### lsveto energy cut has crystal dependent energy cuts
    threshold = 0.125 # in keV
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

    # skip low energy lsveto
    #if i == 8 and E == 0:
    #    energyCut = TCut('(0)')

    ### single hit cuts
    if C == 'S':

        ### single-hit cut
        hitCut = '((singleHitTag['+str(i)+'] == 1) && (multipleHitTag['+str(i)+'] == -1))'

        ### ls veto cut for crystals
        ### need to use smeared resolution
        #lsvetocut = '(edepResol[8]*1000. < 20.0)'
        lsvetocut = '('+edepRes+'[8]*1000. < 80.0)'

        ### cuts for lsveto
        if i == 8:
            # single hit is actually "all hits"
            #hitCut = '((singleHitTag[8] == 1) || (multipleHitTag[8] == 1))'
            hitCut = '((singleHitTag[8] == 1) && (multipleHitTag[8] == -1))'
            lsvetocut = '('+edepRes+'[8]*1000. > 0.0)'

        ### combined cuts
        chanCut = TCut('(('+hitCut+') && ('+lsvetocut+'))')

    ### multi hit cuts
    elif C == 'M':

        ### multi-hit cut
        hitCut = '((multipleHitTag['+str(i)+'] == 1) && (singleHitTag['+str(i)+'] == -1))'

        ### ls veto cut for crystals
        ### need to use smeared resolution
        #lsvetocut = '(edepResol[8]*1000. > 20.0)'
        lsvetocut = '('+edepRes+'[8]*1000. > 80.0)'

        ### cuts for lsveto
        if i == 8:
            hitCut = '((singleHitTag[8] == -1) && (multipleHitTag[8] == 1))'
            lsvetocut = '('+edepRes+'[8]*1000. > 0.0)'

        ### combined cuts
        chanCut = TCut('(('+hitCut+') || ('+lsvetocut+'))')

    else:
        print('ERROR: I do not know what to do with channel -->', C)
        print('Available channels are [S]Single-hits, [M]Multi-hits')
        sys.exit()


    ###  alpha cuts
    #-------------------------------------------------------------------------------
    primDecayCut = TCut('(1)')
    # only apply alpha cuts to crystals and NaI/teflon sources
    if i < 8:
        if isAlphaMC(info):
            if E == 2:
                primDecayCut = TCut('(primDecayType == "alpha")')
            # and remove alphas for the E=0 and E=1 channels
            else:
                primDecayCut = TCut('(primDecayType != "alpha")')
        else:
            if E == 2:
                primDecayCut = TCut('(0)')
                
    # show alphas in lsveto (just to see how it looks)
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
    
    elif location == 'pmt':
        volumeCut = TCut('(primVolumeName == "phys_pmt")')

    elif location == 'cucase':
        volumeCut = TCut('((primVolumeName == "NaIDet0'+str(fromx)+'Case")'
                         +' || (primVolumeName == "NaIDet0'+str(fromx)+'Fringe0")'
                         +' || (primVolumeName == "NaIDet0'+str(fromx)+'Fringe1"))')

    elif 'cucasesurf' in location:
        volumeCut = TCut('(primVolumeName == "NaIDet0'+str(fromx)+'Case")')

    elif location == 'plastic':
        volumeCut = TCut('(primVolumeName == "PlasSupport")')

    elif location == 'lsveto':
        volumeCut = TCut('(primVolumeName == "lsveto")')

    elif location == 'lsvetoair':
        volumeCut = TCut('((primVolumeName == "DetPMTCover")'
                         +' || (primVolumeName == "DetPMTEnvelope"))')

    elif location == 'airshield':
        volumeCut = TCut('(primVolumeName == "LSVetoAirRoom")')

    elif location == 'cushield':
        volumeCut = TCut('(primVolumeName == "CuShield")')
 
    elif location == 'gamma':
        #volumeCut = TCut('(1)')
        volumeCut = TCut('(primVolumeName == "CuShield")')

    elif location == 'innersteel':
        volumeCut = TCut('(primVolumeName == "InnerSteel")')

    elif location == 'steel':
        volumeCut = TCut('((primVolumeName == "SteelSupport")'
                         +' || (primVolumeName == "SteelSupportTop"))')

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
    ### evt_Type is always 3
    if location == 'gamma':
        eventTypeCut  = TCut('(1)')
        generatedCuts = TCut('(1)'+' && '+volumeCut.GetTitle())
    
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ### JUST FOR TESTING
    primParticleCut = TCut('(1)')
    #primParticleCut = TCut('(primParticleName == "Rn222")')
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

            chain.SetAlias('rng', 'sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
            quenchFunc = quench500(i, 'edep[{0}]'.format(i), Q)
            chain.SetAlias('quench', quenchFunc)
            sigmaFunc = sigma500(i, E, quenchFunc)
            chain.SetAlias('sigma', sigmaFunc)
            selection = '((quench*1000.) + (sigma*rng))'

            # or no smear for testing
            #selection = '(edep['+str(i)+']*1000.)'

        # different resolution smearing for C5 and C8
        elif i in [4, 7]:
            chain.SetAlias('rng', 'sin(2.*pi*rndm)*sqrt(-2.*log(rndm))')
            sigmaFunc = sigma500(i, E, 'edep[{0}]'.format(i))
            chain.SetAlias('sigma', sigmaFunc)
            selection = '((edep[{0}]*1000.) + (sigma*rng))'.format(i)

        # or use resolution from processing
        else:
            selection = '('+edepRes+'['+str(i)+']*1000.)'
            
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
        
        # no quenching for alphas in lsveto yet
        if E == 2:
            # show alphas without energy smearing
            selection = '(edep['+str(i)+']*1000.)'
        else:
            selection = '('+edepRes+'['+str(i)+']*1000.)'
            
    #=====================================================================

    
    ### create the MC histogram
    #---------------------------------------------------------------------
    chain.Draw(selection+' >> '+key, masterCut)
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
    
    if 'generated_hist' in mc[key]:
        mc[key]['generated_hist'].Add(temp2)
    else:
        mc[key]['generated_hist'] = copy(temp2)
        
    generated = mc[key]['generated_hist'].GetEntries()
    if generated <= 0:
        print('WARNING: no events generated for --> {0}'.format(key))
    mc[key]['generated'] = generated
    #---------------------------------------------------------------------

    ### create a hist of the detected events
    #---------------------------------------------------------------------
    """
    this_detected = histo.GetEntries()
    key3 = key+'_detected'
    temp3 = TH1F(key3, 'detected', 1, 0, 1)
    temp3.SetBinContent(1, this_detected)
    temp3.SetEntries(this_detected)
    
    if 'detected_hist' in mc[key]:
        mc[key]['detected_hist'].Add(temp3)
    else:
        mc[key]['detected_hist'] = copy(temp3)
    
    detected = mc[key]['detected_hist'].GetEntries()
    mc[key]['detected'] = detected
    """
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

    #del histo
    #del temp2
    #del temp3
    return mc


def processSurf500(chain, mc, info, C, E, fromx, Q=False):

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

    pars = histparam6000(E)
    mc[key]['pars'] = pars
    histo = TH1F(key, longNames(cx), pars[0], pars[1], pars[2])
    
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
        X0, Y0, Z0, rad, height, dep = mcDimensions(fx)
        dist = np.sqrt((X0-xx)**2 + (Y0-yy)**2)
        zdist = abs(Z0-zz)

        # Teflon thickness (mm) - 2020-04-29
        #teflon = 1.015
        # I think this changed in gyunho's new sim
        teflon = 1.016

        ### kinda cryptic...
        # my number /100 to conver to um (88 == 0.88um)
        # /1000 to convert to mm
        # *10 to go 10 expo depths in
        #dep = (((depth/100.)/1000.)*10)
        # use old dep for consistency check
        dep = 0.1

        # for NaI surf side
        if info['loca'].endswith('side') or info['loca'].endswith('expo'):
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


        ### now it's a valid generated event
        #--------------------------------------------------
        if evType <= 10:
            # I don't think I need to weight generated events
            # do here just as a sanity check
            temp2.Fill(0, wtf)
            #temp2.Fill(0)

            
        ### event type cut
        #--------------------------------------------------
        if evType <= 10: continue

        
        ### alpha cuts?
        #--------------------------------------------------
        # FIXME: need to figure out how to make this a string var
        # 97 == alpha
        # 98 == beta
        # skip alphas in crystals for E=0 and E=1
        #print decayType
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
            if edepRes[0] > threshold: continue
            if edepRes[1] > threshold: continue
            if edepRes[2] > threshold: continue
            if edepRes[3] > threshold: continue
            if edepRes[4] > threshold: continue
            if edepRes[5] > threshold: continue
            if edepRes[6] > threshold: continue
            if edepRes[7] > threshold: continue

        if cx == 8 and C == 'M':
            if edepRes[0] <= threshold: continue
            if edepRes[1] <= threshold: continue
            if edepRes[2] <= threshold: continue
            if edepRes[3] <= threshold: continue
            if edepRes[4] <= threshold: continue
            if edepRes[5] <= threshold: continue
            if edepRes[6] <= threshold: continue
            if edepRes[7] <= threshold: continue

        # skip lsveto low energy?
        #if cx == 8 and E == 0:
        #    continue

        
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
                if not (singleHit == 1 and multiHit == -1 and edepRes[8] < 80.0):
                    continue
            if C == 'M':
                if not (singleHit == -1 and multiHit == 1 and edepRes[8] > 80.0):
                    continue


        ### do I need a groupNo cut?
        #--------------------------------------------------
        # -- insert here...
        # not needed now for just surface Pb210


        ### fill the weighted histograms
        #--------------------------------------------------
        # quenching for crystals and E=2
        
        if cx < 8 and E == 2 and Q:
            A, B = quenchPars(cx, Q)
            quenched = A*(edep**2) + (B*edep)
            if quenched < 0: continue
            #p0, p1 = resolPars(cx, E)
            #sigma = (p0*np.sqrt(quenched)) + (p1*quenched)
            # manually setting alpha sigma to 40
            sigma = 40.
            energy = random.gauss(quenched, sigma)
            histo.Fill(energy, wtf)
            
        else:
            histo.Fill(edepRes[cx], wtf)
            
    #=====================================================================
    
    
    ### create the MC histogram
    #---------------------------------------------------------------------
    if 'hist' in mc[key]:
        mc[key]['hist'].Add(histo)
    else:
        mc[key]['hist'] = copy(histo)
    #---------------------------------------------------------------------
    
    ### create a hist of the generated events
    #---------------------------------------------------------------------
    this_generated = temp2.GetEntries()
    print(key, this_generated)
    if 'generated_hist' in mc[key]:
        mc[key]['generated_hist'].Add(temp2)
    else:
        mc[key]['generated_hist'] = copy(temp2)
        
    generated = mc[key]['generated_hist'].GetEntries()
    if generated <= 0:
        print('WARNING: no events generated for --> {0}'.format(key))
    mc[key]['generated'] = generated
    #---------------------------------------------------------------------

    ### create a hist of the detected events
    #---------------------------------------------------------------------
    """
    this_detected = histo.GetEntries()
    key3 = key+'_detected'
    temp3 = TH1F(key3, 'detected', 1, 0, 1)
    temp3.SetBinContent(1, this_detected)
    temp3.SetEntries(this_detected)
    
    if 'detected_hist' in mc[key]:
        mc[key]['detected_hist'].Add(temp3)
    else:
        mc[key]['detected_hist'] = copy(temp3)
    
    detected = mc[key]['detected_hist'].GetEntries()
    mc[key]['detected'] = detected
    """
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

    #del histo
    #del temp2
    #del temp3
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
        
    if key not in mc:
        mc[key] = {}

    if 'info' not in mc[key]:
        mc[key]['info'] = copy(info)

    pars = histparam6000(E)
    mc[key]['pars'] = pars
    histo = TH1F(key, longNames(cx), pars[0], pars[1], pars[2])
    
    key2 = key+'_generated'
    temp2 = TH1F(key2, 'generated', 1, 0, 1)
    #-----------------------------------------------------------------------
    
    
    entries = chain.GetEntries()
    #print entries

    for jentry in range(entries):
        
        chain.GetEntry(jentry)
        
        # get all crystal/lsveto energies in keV
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
        if cx < 8:
            if isAlphaMC(info):
                if E == 2:
                    if decayType != 'alpha': continue
                else:
                    if decayType == 'alpha': continue

        # just show alphas for lsveto E=2
        # just to see how they look
        if cx == 8 and E == 2:
            if decayType != 'alpha': continue
            
        
        ### energy cut
        #--------------------------------------------------
        threshold = 0.125 # in keV
        if edepRes[cx] < threshold: continue
        
        
        ### special crystal energy cuts for lsveto
        #--------------------------------------------------
        if cx == 8:
            if C == 'S':
                if edepRes[0] > threshold: continue
                if edepRes[1] > threshold: continue
                if edepRes[2] > threshold: continue
                if edepRes[3] > threshold: continue
                if edepRes[4] > threshold: continue
                if edepRes[5] > threshold: continue
                if edepRes[6] > threshold: continue
                if edepRes[7] > threshold: continue
            if C == 'M':
                if edepRes[0] <= threshold: continue
                if edepRes[1] <= threshold: continue
                if edepRes[2] <= threshold: continue
                if edepRes[3] <= threshold: continue
                if edepRes[4] <= threshold: continue
                if edepRes[5] <= threshold: continue
                if edepRes[6] <= threshold: continue
                if edepRes[7] <= threshold: continue
        
        
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
            # FIXME - I think we don't handle multi hit tagging the same in data and MC
            if C == 'S':
                if not (singleHit == 1 and multiHit == -1 and edepRes[8] < 80.0):
                    continue
            if C == 'M':
                if not (singleHit == -1 and multiHit == 1 and edepRes[8] > 80.0):
                    continue


        ### fill the histogram
        #--------------------------------------------------
        if cx < 8:
            if E == 2 and Q:
                A, B = quenchPars(cx, Q)
                quenched = A*(edep**2) + (B*edep)
                if quenched < 0: continue
                #p0, p1 = resolPars(cx, E)
                #sigma = (p0*np.sqrt(quenched)) + (p1*quenched)
                # manually setting alpha sigma to 40
                sigma = 40.
                energy = random.gauss(quenched, sigma)
                histo.Fill(energy)
            elif cx in [4, 7]:
                # tweaked resolutions for C5 and C8
                p0, p1 = resolPars(cx, E)
                sigma = (p0*np.sqrt(edep)) + (p1*edep)
                energy = random.gauss(edep, sigma)
                histo.Fill(energy)
            else:
                histo.Fill(edepRes[cx])
                
        if cx == 8:
            # alphas always included and quenched for lsveto
            if decayType == 'alpha':
                
                if DEBUG: print('{0} {1} {2}'.format(jentry, decayType, particle))
                if particle == 'Po216':
                    sys.exit()
                    
                quenched = lsquench(edep)
                if quenched < 0: continue
                # use default lsveto resolution
                sigma = 0.10*quenched
                energy = random.gauss(quenched, sigma)
            else:
                energy = edepRes[cx]
            histo.Fill(energy)
                
    #=====================================================================
    
    
    ### create the MC histogram
    #---------------------------------------------------------------------
    if 'hist' in mc[key]:
        mc[key]['hist'].Add(histo)
    else:
        mc[key]['hist'] = copy(histo)
    #---------------------------------------------------------------------
    
    ### create a hist of the generated events
    #---------------------------------------------------------------------
    this_generated = temp2.GetEntries()
    print(key, this_generated)
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

    #del histo
    #del temp2
    #del temp3
    return mc
    
    
def groupNum500(isotope, numbers=False):

    # example isotope format = "Th228_GRND"
    
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
        if numbers:
            return 0, 1
        else:
            return TCut('(groupNo == 0)')
    
    if chstart == 'K40':
        if numbers:
            return 31, 32
        else:
            return TCut('(groupNo == 31)')

    if chstart == 'Tl208':
        if numbers:
            return 0, 1
        else:
            return TCut('(groupNo == 0)')
        
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
        if numbers:
            return start, stop
        else:
            return TCut('((groupNo >= '+str(start)+') && (groupNo < '+str(stop)+'))')
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
        
    # special case for expo weighting NaI surface Pb210
    elif 'nai-surf-expo' in loc:
        Cx = 'C*'
        if int(info['xstl']) != 99:
            Cx = 'C'+str(info['xstl'])
        """
        mcpath = '/data/COSINE/WORK/gyunho/BkgModel/SET2_reprocessed/surf_NaI'
        filename = 'reprocessed-surface-'+Cx+'-Pb210-*.root'
        """
        """
        searchdir = '/data/COSINE/WORK/mkauer/processed/surface-nai'
        depth = getDepth(info, searchdir, 'surface-nai')
        if depth == 'bulk':
            mcpath = '/data/COSINE/WORK/mkauer/processed/internal'
            filename = '*-internal-*Pb210*root'
        else:
            um = str(depth)+'um'
            mcpath = '/data/COSINE/WORK/mkauer/processed/surface-nai'
            filename = '*-surface-nai-*'+Cx+'*'+um+'*Pb210*root'
        """
        # using other depths from above wasn't working right
        # so just use the 10um for everything
        mcpath = '/data/COSINE/WORK/mkauer/processed/surface-nai'
        filename = '*-surface-nai-*'+Cx+'*10um*Pb210*evt100000.root'
        
    # special case for expo weighting teflon surface Pb210
    elif 'teflon-surf-expo' in loc:
        Cx = 'C*'
        if int(info['xstl']) != 99:
            Cx = 'C'+str(info['xstl'])
        """
        mcpath = '/data/COSINE/WORK/gyunho/BkgModel/SET2_reprocessed/surf_Tef'
        filename = 'reprocessed-tef-surf-'+Cx+'-10um-Pb210-*.root'
        """
        """
        searchdir = '/data/COSINE/WORK/mkauer/processed/surface-teflon'
        depth = getDepth(info, searchdir, 'surface-teflon')
        if depth == 'bulk':
            mcpath = '/data/COSINE/WORK/mkauer/processed/bulk-teflon'
            filename = '*-bulk-teflon-*'+Cx+'*Pb210*root'
        else:
            um = str(depth)+'um'
            mcpath = '/data/COSINE/WORK/mkauer/processed/surface-teflon'
            filename = '*-surface-teflon-*'+Cx+'*'+um+'*Pb210*root'
        """
        # using other depths from above wasn't working right
        # so just use the 10um for everything
        mcpath = '/data/COSINE/WORK/mkauer/processed/surface-teflon'
        filename = '*-surface-tef-*'+Cx+'*10um*Pb210*evt100000.root'
        
    elif loc == 'nai-surf-10um':
        Cx = 'C*'
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceNaI-'+Cx+'-10um-'+iso+'[-|_]*root'
    
    elif loc == 'nai-surf-1um':
        Cx = 'C*'
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceNaI-'+Cx+'-1um-'+iso+'[-|_]*root'
        
    elif loc == 'nai-surf-0.1um':
        Cx = 'C*'
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceNaI-'+Cx+'-0.1um-'+iso+'[-|_]*root'
        
    elif loc == 'nai-surf-0.01um':
        Cx = 'C*'
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceNaI-'+Cx+'-0.01um-'+iso+'[-|_]*root'
        
    elif loc == 'teflon-surf-10um':
        Cx = 'C*'
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceTeflon-'+Cx+'-10um-'+iso+'[-|_]*root'
    
    elif loc == 'teflon-surf-1um':
        Cx = 'C*'
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceTeflon-'+Cx+'-1um-'+iso+'[-|_]*root'
    
    elif loc == 'teflon-surf-0.1um':
        Cx = 'C*'
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceTeflon-'+Cx+'-0.1um-'+iso+'[-|_]*root'
    
    elif loc == 'teflon-surf-0.01um':
        Cx = 'C*'
        mcpath = '/data/MC/COSINE/V4.0.1/reprocessed_SET2'
        filename = '*-surfaceTeflon-'+Cx+'-0.01um-'+iso+'[-|_]*root'
    
    elif loc == 'teflon':
        Cx = 'C*'
        mcpath = '/data/COSINE/WORK/mkauer/processed/bulk-teflon'
        filename = '*-bulk-teflon-*'+Cx+'*'+iso+'[-|_]*root'
        
    elif loc == 'pmt':
        mcpath = '/data/COSINE/WORK/mkauer/processed/pmt'
        filename = '*-pmt-*'+iso+'[-|_]*root'
        
    elif loc == 'xpmt':
        mcpath = '/data/COSINE/WORK/mkauer/processed/pmt'
        filename = '*-pmt-*'+iso+'[-|_]*root'
        
    elif loc == 'cucase':
        mcpath = '/data/COSINE/WORK/mkauer/processed/bulk-cuCase'
        #filename = '*-bulk-cuCase-*'+iso+'[-|_]*root'
        filename = '*-bulk-[c|C]u[c|C]ase-*'+iso+'[-|_]*root'
        
    elif loc == 'plastic':
        mcpath = '/data/COSINE/WORK/mkauer/processed/acrylic'
        filename = '*-acrylic-*'+iso+'[-|_]*root'
                
    elif loc == 'lsveto':
        mcpath = '/data/COSINE/WORK/mkauer/processed/lsveto'
        filename = '*-lsveto-*'+iso+'[-|_]*root'
        
    elif loc == 'cushield':
        mcpath = '/data/COSINE/WORK/mkauer/processed/cuShield'
        #filename = '*-cuShield-*'+iso+'[-|_]*root'
        filename = '*-[c|C]u[s|S]hield-*'+iso+'[-|_]*root'
        
    elif loc == 'gamma':
        mcpath = '/data/COSINE/WORK/mkauer/processed/gamma'
        filename = '*-external-single-gamma-*root'
        
    elif loc == 'innersteel':
        mcpath = '/data/COSINE/WORK/mkauer/processed/innersteel'
        filename = '*-innersteel-*'+iso+'[-|_]*root'

        
    else:
        print('ERROR: no file path identified for --> {0}'.format(loc))
        sys.exit()
    
    return mcpath, filename


def getDepth(info, path, match):
    debug=1
    
    if amLocal():
        path = '/home/mkauer/COSINE/CUP/mc-fitting'+path

    my_loca = info['floca']
    Cx = 'C'+str(info['xstl'])
    
    # location is like "nai-surf-expo-100-side" or "teflon-surf-expo-10-in"
    bits = my_loca.split('-')
    depth = int(bits[-2])
    um = (depth/100.)*10. # covert to um and go 10 depths in

    search = path+'/'+'*'+match+'*'+Cx+'*Pb210*run0*.root'
    if debug: print('searching --> {0}'.format(search))
    # find all available depths generated
    files = glob.glob(path+'/'+'*'+match+'*'+Cx+'*Pb210*run0*.root')
    #if debug: print files
    depths = []
    for ifile in files:
        # get the rootfile name
        ifile = os.path.basename(ifile)
        # split the name by -
        bits = ifile.split('-')
        # get the depth bit - this might change for different mc
        this_depth = bits[5]
        # truncate the "um" at the end of the depth
        this_depth = this_depth[:-2]
        # make sure this is a number
        try:
            this_depth = float(this_depth)
        except:
            print('WARNING: could not convert \"{0}\" to float from {1}'
                  .format(this_depth, ifile))
            continue
        depths.append(this_depth)
    
    # get unique and sorted list of depths
    depths = sorted(list(set(depths)), reverse=True)
    if debug: print('Found MC depths of {0}'.format(depths))
    # find closed mc depth
    # set to largest by default
    use_dep = max(depths)
    for dep in depths:
        if um <= dep:
            use_dep = dep
    # convert to int if it's an int
    if use_dep == int(use_dep):
        use_dep = int(use_dep)
    # set to bulk if too large
    if um > max(depths):
        use_dep = 'bulk'

    if debug: print('My depth: [{0}] - using mc depth: [{1}]'.format(um, use_dep))
    return use_dep


def submitJob(filetype, basekey, rootfile, jsonfile):
    
    scratchdir = '/data/COSINE/WORK/mkauer/scratch'
    mkdir(scratchdir)
        
    jobfile = os.path.join(
        scratchdir,
        basekey+'__'+os.path.basename(rootfile)+'.job'
    )
    
    myscript = os.path.join(
        '/home/mkauer/mc-fitting/root-join-read',
        'process500.py'
    )
    
    logdir = '/data/COSINE/WORK/mkauer/logs'
    mkdir(logdir)
    
    logfile = os.path.join(
        logdir,
        basekey+'__'+os.path.basename(rootfile)+'.log'
    )

    date = str(datetime.datetime.now().date())
    
    with open(jobfile, 'w') as sfile:
        sfile.write('#!/bin/bash \n\n')
        sfile.write('source /home/mkauer/soft/build/bin/thisroot.sh \n')
        sfile.write('export nodename=`hostname` \n')
        sfile.write('export pbsuser=$USER \n')
        sfile.write('export pbsdate='+date+' \n')
        sfile.write('export workdir='+scratchdir+' \n')
        sfile.write(myscript
                    +' '+filetype
                    +' '+rootfile
                    +' '+jsonfile
                    +' 1> '+logfile
                    +' 2> '+logfile
                    +' \n')
        
    subprocess.call(['qsub',
                     #'-q', 'short',
                     '-q', 'medium',
                     #'-q', 'long',
                     '-o', logdir,
                     '-e', logdir,
                     jobfile])


def checkQueue():
    #output = subprocess.check_output(['qstat', '-u', 'mkauer'])
    #output = output.decode('UTF-8').split('\n')
    #output = output[5:-3]
    #return len(output)
    ### this is much faster
    output = subprocess.check_output(['qstat -u mkauer | wc -l'], shell=True)
    return int(output)
    

def lsquench(energy_kev):
    E = energy_kev/1000.
    # model 1
    #c0 = 0.030
    #c1 = 1.640
    #c2 = -0.518
    #c3 = 0.179
    # model 2
    c0 = 0.031
    c1 = 1.689
    c2 = -0.505
    c3 = 0.190
    if E < 6.76:
        return (c0*(E**c1))*1000.
    else:
        return (c2+(c3*E))*1000.
    
    
def quenchParsTest(i, qi):
    cx = int(i)+1
    #qi = int(qi)-1

    # quen = A*x^2 + B*x
    
    if cx == 1:
        A = 4.265778143859101e-05
        B = 0.35207793809347643
    elif cx == 2:
        A = 1.5996668039471557e-05
        B = 0.5054568997410289
    elif cx == 3:
        A = 1.8190956356125193e-05
        B = 0.4519741553475922
    elif cx == 4:
        A = 1.5335156867061897e-05
        B = 0.48128608130333744
    elif cx == 6:
        # linear fit
        A = 0.6555629580279808
        B = -575.9733510992637

    elif cx == 7:
        A = 2.472212333372875e-05
        B = 0.41666666666666735
    else:
        # no quenching
        A = 0.0
        B = 1.0

    if qi == 1:
        return A, B
    else:
        return A, B-475.0
    

def quenchPars(i, qi):
    cx = int(i)+1
    qi = int(qi)-1
    
    po210 = 5406
    qs = [0, 0]

    if cx == 1:
        qs[0] = 3150
        qs[1] = 2700
    elif cx == 2:    
        qs[0] = 3225
        qs[1] = 2825
    elif cx == 3:    
        qs[0] = 2975
        qs[1] = 2650
    elif cx == 4:    
        qs[0] = 3075
        qs[1] = 2875
    elif cx == 6:    
        qs[0] = 2969
        qs[1] = 2540
    elif cx == 7:    
        qs[0] = 2975
        qs[1] = 2550
    else:
        qs[0] = po210
        qs[1] = po210
        
    return 0, float(qs[qi])/float(po210)


def quench500(i, evar_mev, qi):

    ### no quench
    #rootstring = '({0})'.format(evar_mev)
    
    ### quenched = A*x + B
    A, B = quenchParsTest(i, qi)
    rootstring = '(({0}*1000.*{1}) + {2})'.format(evar_mev, A, B)
    rootstring = '({0}/1000.)'.format(rootstring)
    
    ### quenched = A*x^2 + B*x
    #A, B = quenchPars(i, qi)
    #rootstring = '(((({0}*1000.)**2)*{1}) + (({0}*1000.)*{2}))'.format(evar_mev, A, B)
    # convert quenched E back to MeV
    #rootstring = '({0}/1000.)'.format(rootstring)

    
    return rootstring


def resolPars(i, E):
    # res = p0*sqrt(x) + p1*x

    # low energy
    loEresol = [
        [0.2413,  0.01799],   # C1
	[0.2951,  0.01427],   # C2
	[0.3106,  0.007894],  # C3
	[0.3894, -0.001437],  # C4
        #[0.2413,  0.018],    # C5
        [0, 9.5/50.],         # C5 my tweak
	[0.3620,  0.0006355], # C6
	[0.3042,  0.009784],  # C7
        #[0.2413,  0.018],    # C8
        [0, 14.8/50.],        # C8 my tweak
        [0, 0.10]             # C9
    ]
    
    # high energy
    hiEresol = [
        [0.6729,  0.009374],  # C1
	[0.6531,  0.006627],  # C2
	[0.5926,  0.009506],  # C3
	[0.7227,  0.004790],  # C4
        #[0.6729,  0.009374], # C5
        [14.73/np.sqrt(50.), 0], # C5 my tweak
	[0.6498,  0.009670],  # C6
	[0.7034,  0.007812],  # C7
        #[0.6729,  0.009374], # C8
        [24.76/np.sqrt(50.), 0], # C8 my tweak
        [0.0,  0.10]          # C9
    ]
    
    if E: p0, p1 = hiEresol[int(i)]
    else: p0, p1 = loEresol[int(i)]

    return p0, p1


def sigma500(i, E, evar_mev):
    # these are equivalent:
    # res = (p0/sqrt(x) + p1)*x
    # res = p0*sqrt(x) + p1*x
    p0, p1 = resolPars(i, E)
    #rootstring = '(({1}/sqrt({0}*1000.)) + {2})'.format(evar_mev, p0, p1)
    rootstring = '(({1}*sqrt({0}*1000.)) + ({2}*{0}*1000.))'.format(evar_mev, p0, p1)

    # different sigma for alphas
    if E == 2:
        if   i == 0: # C1
            rootstring = '(40)'
        elif i == 1: # C2
            rootstring = '(40)'
        elif i == 2: # C3
            rootstring = '(40)'
        elif i == 3: # C4
            rootstring = '(40)'
        elif i == 5: # C6
            rootstring = '(40)'
        elif i == 6: # C7
            rootstring = '(40)'
        else:
            rootstring = '(40)'
            
    # convert back to MeV?
    #rootstring = '({0}/1000.)'.format(rootstring)
    return rootstring


    
def dataToCluster(data, rootfile, info, jsonfile=None, cluster=False):

    if info is None and not os.path.exists(jsonfile):
        print('WARNING: you must specify an info json filename')
        sys.exit()
        
    if info is None:
        with open(jsonfile) as jf:
            info = json.load(jf)
    
    if data is None:
        data = {}

    if cluster:
        location = info['loca']+'/'+info['build']
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
    
    for C in ['S', 'M']:
        for E in [0, 1, 2]:

            # if special case '99' then process all crystals at once
            if int(info['xstl']) == 99:
                for xstl in range(1, 10):
                    basekey = info['key'].replace('x99-', 'x'+str(xstl)+'-', 1)
                    if xstl == 9 and C == 'S':
                        chain2 = TChain("activeLS","")
                        chain2.Add(rootfile)
                        data = processData500(chain2, data, info, basekey, C, E)
                    else:
                        data = processData500(chain, data, info, basekey, C, E)

            else:
                basekey = info['key']
                if int(info['xstl']) == 9 and C == 'S':
                    chain2 = TChain("activeLS","")
                    chain2.Add(rootfile)
                    data = processData500(chain2, data, info, basekey, C, E)
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
    
    #i = int(info['xstl']) - 1
    bits = basekey.split('-')
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

    if key not in data:
        data[key] = {}
        
    if 'info' not in data[key]:
        data[key]['info'] = copy(info)
    
    pars = histparam6000(E)
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
    elif build in ['V00-04-15', 'V00-04-16', 'V00-04-17']:
        ### calib301 uses the standard production calibrations
        ### and uses the lsveto poly43 calibration
        #edep, selection = calib301(i, E)
        ### calib302 can have different calibs for lsveto "S" and "M"
        ### and tweaked high energy calibrations of the crystals
        edep, selection = calib302(i, E, C)
        masterCut = cutsBDT302(i, C, E, edep, selection)

    else:
        print('WARNING: no selection for build \"{0}\"'.format(build))
        return data
    #-----------------------------------------------------------------------


    # FILL HISTOS
    #-----------------------------------------------------------------------
    chain.Draw(selection+' >> '+key, masterCut)
    if 'hist' in data[key]:
        data[key]['hist'].Add(histo)
    else:
        data[key]['hist'] = copy(histo)

    # sanity check on the data
    if histo.GetEntries() == 0:
        # Don't warn about multi-hit alpha channel
        if not (C == 'M' and E == 2):
            print('WARNING: no events found for [{0}] --> \"{1}\"'
                  .format(key, os.path.basename(rootfile)))
        # Don't do anything about this here
        # Check for 0 events when selecting runs from
        # the good runs list
        
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

    #del histo
    #del temp2
    #del temp3
    return data


def writeHists500(hist_dict):
    for key in hist_dict:
        for key2 in hist_dict[key]:
            try: hist_dict[key][key2].Write()
            except: pass
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


def combineOthers500(sigs, globalMC):

    print('INFO: combining histograms...')
    
    debug = 0
    
    donekeys = []
    delete = []
    lsveto = {}
    
    for key in sorted(sigs):

        if key in donekeys:
            continue

        if debug: print('working on [{0}]'.format(key))

        bits = key.split('-')

        # don't combine globals
        if bits[1] in globalMC:
            if debug: print('global mc, not combining [{0}]'.format(key))
            continue
        
        X = int(bits[0][1])
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
                
            if newkey not in sigs:
                #print('WARNING: key not found in backgrounds [{0}]'.format(newkey))
                continue

            # TESTING !!!
            # don't combine xpmt from others
            if bits[1] == 'xpmt':
            #if bits[1] == 'kaakaapoopoo':
                print('TEST: skipping xpmt [{0}]'.format(key))
            else:
                if debug: print('adding [{0}] --> [{1}]'.format(newkey, default))
                if X==9:
                    if default not in lsveto:
                        lsveto[default] = deepcopy(sigs[newkey])
                    else:
                        lsveto[default]['hist'].Add(sigs[newkey]['hist'])
                else:
                    sigs[default]['hist'].Add(sigs[newkey]['hist'])

            donekeys.append(default)
            donekeys.append(newkey)
            delete.append(newkey)

    # transfer new lsveto keys into sigs
    for lskey in sorted(lsveto):
        if debug: print('copy new lsveto key/hist to backgrounds [{0}]'.format(lskey))
        sigs[lskey] = deepcopy(lsveto[lskey])
        
    # delete the other histograms
    for dkey in sorted(delete):
        if debug: print('deleting old key [{0}]'.format(dkey))
        del sigs[dkey]

    return sigs




if __name__ == "__main__":
    main(sys.argv[1:])

