
#!/usr/bin/env python

import os, shutil, sys, numpy as np, cPickle as pickle, shelve
import utils
import time
import reasonable as cg

CURRDIR = os.getcwd()
RoomTemp = 300.

PdbName = sys.argv[1]
BBType = sys.argv[2]
OutDir = os.path.abspath(sys.argv[3])
NReplica = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    
# parse paths etc
OutDir = os.path.join(OutDir, PdbName)
if not os.path.isdir(OutDir): os.system('mkdir -p %s' % OutDir)
NativePdb = utils.parseNative(PdbName, MasterDir = os.path.expanduser('~/protein_model/native_struct/mapped'))
try:
    AATopClustPdb = utils.parseNative(PdbName, MasterDir = os.path.expanduser('~/protein_model/native_struct/ff96_igb5_glghs_topclust_mapped'))
except IOError:
    print 'Utils Error: Requested Top clust pdb does not exist'
    AATopClustPdb = None
FF_File, FFMetadata = utils.parseBBFF(BBType)
Prefix = 'prot_' + PdbName

# parse forcefield parameters
MinBondOrd = FFMetadata['MinBondOrd']
NKnot = FFMetadata['NKnot']
SPCut = FFMetadata['Cut']
# Go parameters
NativeSigma = FFMetadata['NativeSigma']
NativeEpsilon = FFMetadata['NativeEpsilon']
NativeCut = FFMetadata['NativeCut']
NonNativeSigma = FFMetadata['NonNativeSigma']
NonNativeEpsilon = FFMetadata['NonNativeEpsilon']
NonNativeCut = FFMetadata['NonNativeCut']
    
# temp schedule
TLow = 270.
THigh = 500.
Temps = np.logspace(np.log10(TLow), np.log10(THigh), NReplica)
TempInd = np.argmin(abs(Temps - RoomTemp))
TempSet = Temps[TempInd]
    
# time-step
TimeStep = FFMetadata['TimeStep'] # femto-seconds
    
# MD iterations
NStepsMin = 10000                   # 10 ps
NStepsEquil = 50000000              # 50 ns
NStepsProd  = 20000000              # 20 ns
NStepsSwap = 1000                   # 1 ps
StepFreq = int(NStepsProd / 10000)  # need 10000 frames, 2 ps
    
# REMD script template
mdstr = '''
#!/usr/bin/env python
import os, sys, numpy as np, time
import reasonable as cg, utils

# prefix
Prefix = "%(PREFIX)s"

# check if only analysis needs to be done
FullPrefix = os.path.join(os.getcwd(), Prefix)
isDone = utils.checkGoREMD(FullPrefix, [%(TEMPSET)3.2f])
if isDone: exit()
        
# set up config object
cfg = cg.config.Config()

# backbone settings
cfg.MinBondOrd = %(MINBONDORD)d
cfg.NKnot = %(NKNOT)d
cfg.SPCut = %(SPCUT)g
cfg.NCOSType = 0

# native contacts
cfg.NativeType = 0
cfg.NativeEpsilon = %(NATIVEEPSILON)g
cfg.NativeSigma = None # enforce auto-sigma
# also ensure default value for Native Cut

# nonnative contacts
cfg.NonNativeType = 0
cfg.NonNativeSigma = %(NONNATIVESIGMA)g
cfg.NonNativeEpsilon = %(NONNATIVEEPSILON)g
cfg.NonNativeCut = None # enfore WCA

# timestep
cfg.TimeStep = %(TIMESTEP)g

# temp schedule
Temps = np.logspace(np.log10(%(TLOW)3.2f), np.log10(%(THIGH)3.2f), %(NREPLICA)d)
       
# compile Sys object
p, Sys = cg.cgmodel.makeGoSys(NativePdb = "%(NATIVEPDB)s", cfg = cfg, Prefix = Prefix, TempSet = %(TEMPSET)3.2f)

# load backbone forcefield file
cg.cgmodel.loadParam(Sys, "%(FF_FILE)s")

# compile REMD object
md = cg.mdsim.REMD(p, Sys, Prefix = Prefix, Temps = Temps,
                   NStepsMin = %(NSTEPSMIN)d, NStepsEquil = %(NSTEPSEQUIL)d, NStepsProd = %(NSTEPSPROD)d,
                   NStepsSwap = %(NSTEPSSWAP)d, StepFreq = %(STEPFREQ)d)

# run REMD
t1 = time.time()
md.runREMD()
t2 = time.time()
        
# reorder by temperature only at room temp
md.reorderTraj(ReorderTemps = [%(TEMPSET)3.2f] )
t3 = time.time()
        
# print stats
print "REMD time: ", (t2-t1), " seconds"
print "Reordering time: ", (t3-t2), " seconds"
'''
    
# job script template
jobstr = '''
#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N %(JOBNAME)s
        
export PYTHONPATH=$PYTHONPATH:~/protein_model
date
python remd.py
mkdir -p ./NativeAnalysis
mkdir -p ./AATopClustAnalysis
python ~/Go/analyze_go.py %(NATIVEPDB)s %(PREFIX)s ./ ./NativeAnalysis
python ~/Go/analyze_go.py %(AATOPCLUSTPDB)s %(PREFIX)s ./ ./AATopClustAnalysis
date
'''
# dict for filling md script template
d1 = {
      'PREFIX': Prefix,
      'NATIVEPDB': NativePdb,
        
      'FF_FILE': FF_File,
      'MINBONDORD': MinBondOrd,
      'NKNOT' : NKnot,
      'SPCUT' : SPCut,
      
      'NATIVEEPSILON': NativeEpsilon,
      'NATIVESIGMA': NativeSigma,
      'NATIVECUT': NativeCut,

      'NONNATIVEEPSILON': NonNativeEpsilon,
      'NONNATIVESIGMA': NonNativeSigma,
      'NONNATIVECUT' : NonNativeCut,
        
      'TLOW'    : TLow,
      'THIGH'   : THigh,
      'NREPLICA': NReplica,
      'TEMPSET' : TempSet,
        
      'TIMESTEP' : TimeStep,
      'NSTEPSMIN': NStepsMin,
      'NSTEPSEQUIL': NStepsEquil,
      'NSTEPSPROD' : NStepsProd,
      'NSTEPSSWAP' : NStepsSwap,
      'STEPFREQ': StepFreq
    }

# dict for filling job script template
d2 = {'JOBNAME': Prefix, 'NATIVEPDB': NativePdb, 'AATOPCLUSTPDB': AATopClustPdb, 'PREFIX': Prefix}

# run REMD job
mdscript = os.path.join(OutDir, 'remd.py')
jobscript = os.path.join(OutDir, 'remdjob.sh')
if not os.path.isfile(mdscript): file(mdscript, 'w').write(mdstr % d1)
file(jobscript, 'w').write(jobstr % d2)
os.chdir(OutDir)
os.system('qsub remdjob.sh')
os.chdir(CURRDIR)

