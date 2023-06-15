#! /usr/bin/env python
# -*- coding: utf-8 -*-
#    2016-May-02  shaw@noao.edu

import sys
import copy
from pyraf import iraf
from pyraf.iraf import gemini, gemtools, gmos, rv
import fileSelect as fs

class ReductParam(object):
    '''Container for data reduction file lists, MasterCal names, etc.
    '''
    def __init__(self, expType, expParams, dbFile, bias=None, trace=None, gaps=None, flat=None, wvtran=None, sensfunc=None, target=None):
        '''Construct a data reduction container for a set of exposures.
        '''
        self.type = str(expType)
        self.qd = expParams
        self.dbFile = str(dbFile)
        self.MCbias = str(bias)
        self.MCtrace = str(trace)
        self.MCgaps = str(gaps)
        self.MCresponse = str(flat)
        self.MCwavtran = str(wvtran)
        self.MCsens = str(sensfunc)
        self.targName = str(target)
        # Construct a list of files appropriate for the exposure type. 
        qd = self.qd
        if qd is not None and self.dbFile is not None:
            self.fileList = fs.fileListQuery(self.dbFile, fs.createQuery(self.type, qd), qd)
        else:
            self.fileList = []

def gmos_ifu_proc():
    '''
    GMOS Data Reduction Cookbook (DRC) companion script to the chapter:
        "Reduction of IFU Spectral Images with PyRAF"

    PyRAF script to:
    Process IFU exposures for J224024.1–092748, in programs GS-2012B-Q-26 and 
    GN-2012B-Q-226.

    The names for the relevant header keywords and their expected values are
    described in the DRC chapter entitled "Supplementary Material"

    Perform the following in the parent work directory:
        cd /path/to/work_directory

    If you have not already done so, define the LACosmic script. 
        iraf.task(lacos_spec='./lacos_spec.cl')
        
    Place the following (both are available from the DRC) in your work directory: 
        - the fileSelect.py module
        - the static BPM MasterCal: bpm_gmos-s_EEV_v2_1x1_spec_MEF.fits

    This whole example depends upon first having built an sqlite3 database of 
    metadata from the raw exposures:
       cd /path/to/work_directory/raw
       python obslog.py obsLog.sqlite3
       
    Finally, execute this script in your work directory, from the unix prompt:
        python gmos_ifu_proc.py
    '''

    print ("### Begin Processing GMOS/IFU Spectral Images ###")
    print (" ")
    dbFile='raw/obsLog.sqlite3'
    bpm = 'bpm_gmos-s_EEV_v2_1x1_spec_MEF.fits'
    iraf.set.stdimage = 'imtgmos'
 
    
    print ("=== Creating MasterCals ===")

    # Create the query dictionary of essential exposure parameter=value pairs.
    # Restrict bias exposures to within ~2 months of the target observations:
    qds1 = {'use_me':1,
           'CcdBin':'1 1',
           'DateObs':'2012-08-20:2012-09-01',
           'Instrument':'GMOS-S',
           'Disperser':'B600+_%',
           'AperMask':'IFU-2',
           'CentWave':625.0,
           'Object':'J2240-0927',
           'RoI':'Full'
           }
    # Also make dictionaries for late-epoch GMOS-S...
    qds2 = copy.deepcopy(qds1)
    qds2['DateObs'] = '2012-10-26:2012-11-03'
    qds2['Object'] = 'LTT9239'

    # ...and also for GMOS-N.
    qdn1 = copy.deepcopy(qds1)
    qdn1['DateObs'] = '2012-10-26:2012-11-03'
    qdn1['Instrument'] = 'GMOS-N'
    qdn1['CentWave'] = 499.0
    qdn2 = copy.deepcopy(qdn1)
    qdn2['DateObs'] = '2012-10-26:2012-11-03'
    qdn2['Disperser'] = 'R831+_%'
    qdn2['CentWave'] = 853.0
    
    print (" --Creating Bias MasterCals--")

    # Use primarily the default task parameters.
    gemtools.gemextn.unlearn()    # Disarm a bug in gbias
    gmos.gbias.unlearn()
    gmos.gbias.logfile = 'biasLog.txt'
    gmos.gbias.rawpath = './raw/'
    gmos.gbias.fl_vardq = 'yes'
    gmos.gbias.verbose = 'no'

    bias = {}
    # The following SQL generates the list of bias files to process.
    SQL = fs.createQuery('bias', qds1)
    bias['S1'] = fs.fileListQuery(dbFile, SQL, qds1)

    # All in one statement for the second epoch.
    bias['S2'] = fs.fileListQuery(dbFile, fs.createQuery('bias', qds2), qds2)

    # Process bias files for each epoch of observations.
    # for epoch in ['S1', 'S2']:
    for epoch in ['S1']:
        # The str.join() function is needed to transform a python list into 
        # a string of comman-separated files that IRAF can understand.
        gmos.gbias(','.join(str(x) for x in bias[epoch]), 'MCbias' + epoch)
    
    # Clean up
    iraf.imdel('gS2012*.fits')
    
    print (" --Creating GCAL Spectral Flat-Field MasterCals--")
    # Set the task parameters.
    gmos.gireduce.unlearn()
    gmos.gfreduce.unlearn()
    gmos.gfreduce.verbose = 'no'
    traceFlags = {
        'fl_addmdf':'yes', 'fl_over':'yes', 'fl_trim':'yes', 'fl_bias':'yes', 
        'fl_extract':'yes', 'fl_gsappwave':'no', 'fl_wavtran':'no',
        'fl_gscrrej':'no', 'fl_skysub':'no', 'fl_fluxcal':'no', 'fl_vardq':'yes',
        'fl_inter':'no', 'rawpath':'./raw/', 'logfile':'gfreduceLog.txt'
        }
    flatFlags = {
        'fl_addmdf':'no', 'fl_over':'no', 'fl_trim':'no', 'fl_bias':'no', 
        'fl_extract':'yes', 'fl_gsappwave':'yes', 'fl_wavtran':'no', 
        'fl_gscrrej':'no', 'fl_skysub':'no', 'fl_fluxcal':'no', 'fl_vardq':'yes', 
        'fl_inter':'no', 'rawpath':'./', 'logfile':'gfreduceLog.txt'
        }
    scatsubFlags = {
        'prefix':'b', 'xorder':'5,9,5', 'yorder':'5,7,5', 'cross':'yes'
        }
    gmos.gfscatsub.unlearn()
    gemtools.gemfix.unlearn()
    gemtools.gemfix.logfile='gemfixLog.txt' 

    # Basic reductions for each epoch of GCAL flat-fields.
    
    print ("  --Basic Flat-field reductions--")
    flat = {}
    flat['S1'] = fs.fileListQuery(dbFile, fs.createQuery('gcalFlat', qds1), qds1)
    # flat['S2'] = fs.fileListQuery(dbFile, fs.createQuery('gcalFlat', qds2), qds2)
    
    for epoch in ['S1']:
        for f in flat[epoch]:
            print "  -Processing file: %s" % (f)
            gmos.gfreduce (f, bias='MCbias'+epoch, **traceFlags)
            # Insert the BPM as the initial DQ extensions.
            rgf = 'rg' + f
            for i in ['1','2','3']:
                dq = '[DQ,'+i+']'
                dqow = '[DQ,'+i+',overwrite+]'
                iraf.imcopy(bpm+dq, rgf+'.fits'+dqow)
            print "  -Interpolating over bad pixels"
            gemtools.gemfix (rgf, 'p'+rgf, method='fit1d', bitmask=1, order=32)

            prgf = 'prg' + f
            print "  -Extracting from file: %s" % (prgf)
            outFile = 'e' + prgf + '_trace'
            gmos.gfreduce (prgf, outimages=outFile, **flatFlags)

            print "  -Subtracting scattered light from file: %s" % (prgf)
            traceFile = 'e' + prgf + '_trace'
            gapsFile = 'e' + prgf + '_gaps'
            gmos.gffindblocks ('prg'+f, traceFile, gapsFile)
            gmos.gfscatsub (prgf, gapsFile, **scatsubFlags)
     
    # Apply QE correction to old GMOS-S CCDs here.
    # Re-extract GCAL flat-fields on scattered-light corrected exposures, and normalize.
    print (" -- Normalize Flat-fields --")
    gmos.gfresponse.unlearn()
    gmos.gfresponse.verbose = 'no'
    responseFlags = {
        'skyimage':'', 'fl_inter':'no', 'fl_fit':'no', 'sample':'1:1,30:2855,2879:2879', 
        'logfile':'gfresponseLog.txt', 'function':'spline3', 'order':47
        }
    flatFlags = {
        'fl_addmdf':'no', 'fl_over':'no', 'fl_trim':'no', 'fl_bias':'no', 
        'fl_extract':'yes', 'fl_gsappwave':'yes', 'fl_wavtran':'no', 
        'fl_gscrrej':'no', 'fl_skysub':'no', 'fl_fluxcal':'no', 'fl_vardq':'yes', 
        'fl_inter':'no', 'rawpath':'./', 'logfile':'gfreduceLog.txt'
        }
    prefix='bprg'
    for epoch in ['S1']:
        for f in flat[epoch]:
            print "  -Extracting & normalizing file: %s" % (prefix+f)
            # outFile = prefix+f + '_flat'
            outFile = 'e'+prefix+f + '_flat'
            gmos.gfreduce (prefix+f, **flatFlags)
            # gmos.gfresponse (prefix+f, outimage=outFile, **responseFlags)
            gmos.gfresponse ('e'+prefix+f, outimage=outFile, **responseFlags)
    # Clean up
    iraf.imdele ("gS2012*.fits,rgS2012*.fits")
    
    print ("=== Processing Arc Files ===")
    gmos.gswavelength.unlearn()
    waveFlags = {
        'fwidth':8, 'minsep':2.5, 'fl_inter':'no', 'logfile':'gswaveLog.txt'
        }
    arcFlags = {
        'fl_addmdf':'yes', 'fl_over':'yes', 'fl_trim':'yes', 'fl_bias':'no', 
        'fl_extract':'yes', 'fl_gsappwave':'yes', 'fl_wavtran':'no', 
        'fl_gscrrej':'no', 'fl_skysub':'no', 'fl_fluxcal':'no', 'fl_vardq':'no', 
        'fl_fixgaps':'no', 'reference':'eprgS20120827S0069_trace', 
        'fl_inter':'no', 'rawpath':'./raw/', 'logfile':'gfreduceLog.txt'
        }
    arcs = fs.fileListQuery(dbFile, fs.createQuery('arc', qds1), qds1)
    for f in arcs:
        print "  -Processing file: %s" % (f)
        gmos.gfreduce (f, **arcFlags)
        gmos.gswavelength ('erg'+f, **waveFlags)

    print ("=== Processing Std Star Exposures ===")

    # Record some associative calibration information.
    std1 = ReductParam ('std', qds1, 'raw/obsLog.sqlite3', bias='MCbiasS1', 
                        trace='eprgS20120829S0062_trace', flat='ebprgS20120829S0062_flat', 
                        gaps='eprgS20120829S0062_gaps', wvtran='ergS20120829S0199')
    # std2 = ReductParam ('sciSpec', qds2, 'raw/obsLog.sqlite3', bias='MCbiasS2', 
    #                     trace='eprgS20121031S0025_trace', flat='eprgS20121031S0025_flat', 
    #                     gaps='eprgS20121031S0025_gaps', wvtran='ergS20120828S0199')

    print (" -- Basic Processing --")
    basicFlags = {
        'fl_addmdf':'yes', 'fl_over':'yes', 'fl_trim':'yes', 'fl_bias':'yes', 
        'fl_extract':'no', 'fl_gsappwave':'no', 'fl_wavtran':'no', 'fl_fixgaps':'no', 
        'fl_gscrrej':'no', 'fl_skysub':'no', 'fl_fluxcal':'no', 'fl_vardq':'yes', 
        'fl_inter':'no', 'rawpath':'./raw/', 'logfile':'gfreduceLog.txt'
        }
    prefix = 'rg'
    for s in [std1]:
        for f in s.fileList:
            print "  -Processing file: %s" % (prefix+f)
            gmos.gfreduce (f, bias=s.MCbias, **basicFlags)
            # Insert the BPM as the initial DQ extensions.
            rgf = prefix + f
            for i in ['1','2','3']:
                dq = '[DQ,'+i+']'
                dqow = '[DQ,'+i+',overwrite+]'
                iraf.imcopy(bpm+dq, rgf+'.fits'+dqow)

    # Clean up
    iraf.imdel('gS2012*.fits')
    prefix = 'rg'
    print (" -- Artifact Rejection --")
    gemtools.gemcrspec.unlearn()
    gemtools.gemcrspec.fl_vardq='yes'
    gemtools.gemcrspec.verbose='yes'
    gemtools.gemcrspec.logfile='gemcrspecLog.txt'

    for f in std1.fileList:
        print "  -CR rejecting file: %s" % (prefix+f)
        xFile = 'x' + prefix + f
        gemtools.gemcrspec (prefix+f, xFile, sigfrac=0.32, niter=4)
        gemtools.gemfix (xFile, 'p'+xFile, method='fixpix')

    print (" -- Scattered Light Correction --")
    prefix = 'pxrg'
    for s in [std1]:
        for f in s.fileList:
            fName = prefix + f
            print "  -Subtracting scattered light from file: %s" % (fName)
            gmos.gfscatsub (fName, s.MCgaps)

    # Apply QE correction here for GMOS-S CCDs.
    print (" -- Spectral Processing --")
    stdFlags = {
        'fl_addmdf':'no', 'fl_over':'no', 'fl_trim':'no', 'fl_bias':'no', 
        'fl_extract':'yes', 'trace':'no', 'recenter':'no', 'fl_gscrrej':'no', 
        'fl_fixgaps':'yes', 'fl_gsappwave':'yes', 'fl_wavtran':'yes', 
        'w1':5618., 'w2':'INDEF', 'dw':0.4622, 'nw':2822, 
        'fl_skysub':'yes', 'sepslits':'yes', 'fl_fluxcal':'no',  
        'fl_vardq':'yes', 'fl_inter':'no', 'rawpath':'./'
        }
    prefix = 'bpxrg'

    for s in [std1]:
        for f in s.fileList:
            print "  -Processing file: %s" % (prefix+f)
            gmos.gfreduce (prefix+f, reference=s.MCtrace, response=s.MCresponse, 
                           wavtraname=s.MCwavtran, **stdFlags)

    print (" -- Aperture Summation --")
    gmos.gfapsum.unlearn()
    gmos.gfapsum.fl_inter="no"
    gmos.gfapsum.logfile="gfapsumLog.txt"
    prefix = 'stebpxrg'

    # EG131:
    for f in std1.fileList:
        print "  -Summing file: %s" % (prefix+f)
        gmos.gfapsum (prefix+f, outimages='EG131sum_B6-625', lthreshold=0.)
    # LTT9239:
    # for f in std2.fileList:
    #     print "  -Summing file: %s" % (prefix+f)
    #     gmos.gfapsum (prefix+f, lthreshold=0.)

    # # Combine exposures
    # gemtools.gemcombine.unlearn()
    # gemtools.gemcombine.logfile="gemcombineLog.txt"
    # prefix = 'astepxrg'
    # gemtools.gemcombine (','.join(prefix+str(x) for x in std2.fileList), 
    #                      'LTT9239sum_B6-625', reject='none', scale='exposure')

    print (" -- Sensitivity Derivation --")

    gmos.gsstandard.unlearn()
    gsstdFlags = {
        'caldir':'./', 'extinction':'onedstds$ctioextinct.dat', 'order':11,
        'observatory':'Gemini-South', 'logfile':'gsstandardLog.txt', 
        'fl_inter':'yes'
        }
    gmos.gsstandard ('EG131sum_B6-625', 'std_B6-625', 
                     'Sens_B6-625', starname='eg131', **gsstdFlags)

    # Clean up
    iraf.imdel('gS2012*.fits,rgS2012*.fits')

    print ("=== Reducing Science Exposures ===")
    print (" -- Basic Processing --")
                             
    sci = ReductParam ('sciSpec', qds1, 'raw/obsLog.sqlite3', bias='MCbiasS1', 
                        trace='eprgS20120827S0069_trace', flat='ebprgS20120827S0069_flat', 
                        wvtran='ergS20120829S0199', sensfunc='Sens_B6-625')
                          
    for f in sci.fileList:
        print "  -Processing file: %s" % (prefix+f)
        gmos.gfreduce (f, bias=sci.MCbias, **basicFlags)
        # Insert the BPM as the initial DQ extensions.
        rgf = prefix + f
        for i in ['1','2','3']:
            dq = '[DQ,'+i+']'
            dqow = '[DQ,'+i+',overwrite+]'
            iraf.imcopy(bpm+dq, rgf+'.fits'+dqow)

    print (" -- Artifact Rejection --")

    # This step would be folded into the basic processing loop, execpt that the 
    # required CPU is of order an hour per exposure. 
    for f in sci.fileList:
        print "  -CR rejecting file: %s" % (prefix+f)
        xFile = 'x' + prefix + f
        gemtools.gemcrspec (prefix+f, xFile, sigfrac=0.32, niter=4)
        gemtools.gemfix (xFile, 'p'+xFile, method='fixpix')

    # Perform QE correction here. 
    print (" -- Spectral Processing --")
    sciFlags = {
        'fl_addmdf':'no', 'fl_over':'no', 'fl_trim':'no', 'fl_bias':'no', 
        'fl_extract':'yes', 'trace':'no', 'recenter':'no', 'fl_gscrrej':'no', 
        'fl_fixgaps':'yes', 'fl_gsappwave':'yes', 'fl_wavtran':'yes', 
        'w1':5618., 'w2':'INDEF', 'dw':0.4622, 'nw':2822, 
        'fl_skysub':'yes', 'sepslits':'yes', 'fl_fluxcal':'yes',  
        'fl_vardq':'yes', 'fl_inter':'no', 'rawpath':'./'
        }
    prefix = 'pxrg'
    for f in sci.fileList:
        print "  -Processing file: %s" % (prefix+f)
        gmos.gfreduce (prefix+f, reference=sci.MCtrace, response=sci.MCresponse, 
                        wavtraname=sci.MCwavtran, sfunction=sci.MCsens, **sciFlags)

    # Clean up
    iraf.imdel('gN2012*.fits')

    print (" -- Wavelength Zeropoint Correction -- ")

    rv.rvidlines.unlearn()
    rv.observatory='gemini-south'
    rvFlags = {
        'coordlist':'skylines.txt', 'ftype':'emission', 'nsum':1, 'threshold':7., 
        'maxfeatures':10, 'fwidth':10., 'cradius':10., 'minsep':5., 
        'logfile':'rvLog.txt'
        }
    rv.rvidlines ('stepxrgS20120827S0066.fits[SKY]', **rvFlags)
    # ...and so on for the other target exposures. 

    # Results of the RV analysis borrowed from JT measurements: 
    waveCorr = {
        'cstepxrgS20120827S0066':5618.164,
        'cstepxrgS20120827S0067':5618.268,
        'cstepxrgS20120827S0068':5618.387,
        'cstepxrgS20120827S0070':5618.560,
        'cstepxrgS20120827S0071':5618.667,
        'cstepxrgS20120827S0072':5618.787,
        }
    iraf.hedit.unlearn()
    iraf.hedit.update = 'yes'
    iraf.hedit.verify = 'no'

    for file,w0 in waveCorr.iteritems():
        iraf.hedit (file+'.fits[sci]', 'CRVAL1', w0)

    print ("=== Finished Calibration Processing ===")
    print (" ")
    print ("=== Create Advanced Products ===")

    print (" -- Create (x,y,lambda) datacube --")
    gmos.gfcube.unlearn()
    cubeFlags = {
        'fl_atmdisp':'yes', 'fl_flux':'yes', 'fl_var':'yes', 'fl_dq':'yes',
        'logfile':'gfcubeLog.txt'
    }
    prefix = 'cstepxrg'
    for f in sci.fileList:
        gmos.gfcube (prefix+f, **cubeFlags)

    print (" -- Align the datacubes spatially --")
    # The offsets were derived by James Turner using his pyfalign software. 
    iraf.hedit ('dcstepxrgS20120827S0070.fits[sci]', 'CRVAL1', 65.48042)
    iraf.hedit ('dcstepxrgS20120827S0071.fits[sci]', 'CRVAL1', 65.47042)
    iraf.hedit ('dcstepxrgS20120827S0071.fits[sci]', 'CRVAL2', 0.045000)
    iraf.hedit ('dcstepxrgS20120827S0072.fits[sci]', 'CRVAL1', 65.50042)
    iraf.hedit ('dcstepxrgS20120827S0072.fits[sci]', 'CRVAL2', 0.035000)

    print (" -- Combine the datacubes --")
    iraf.imcombine.unlearn()
    iraf.imcombine ('dcstepxrgS20120827S00*.fits[SCI]', 'j2240-097_cube.fits')
    
    print ("=== End of Advanced Processing ===")

if __name__ == "__main__":
    gmos_ifu_proc()
