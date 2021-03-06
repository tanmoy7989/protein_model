
#!/usr/bin/env python

import os, shutil, sys, copy, numpy as np, cPickle as pickle, shelve
import utils
import time
import reasonable as cg

CURRDIR = os.getcwd()
RoomTemp = 300.

hasPseudoGLY = True

PdbName = sys.argv[1]
FFType = sys.argv[2]
OutDir = os.path.abspath(sys.argv[3])
DelFrac = float(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] else 0.0
NReplica = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    
# parse paths etc
OutDir = os.path.join(OutDir, '%1.2f' % DelFrac)
if not os.path.isdir(OutDir): os.system('mkdir -p %s' % OutDir)
if hasPseudoGLY:
    MasterDir_Native = os.path.expanduser('~/protein_model/native_struct/mapped_pseudoGLY')
    MasterDir_AATopClust = os.path.expanduser('~/protein_model/native_struct/ff96_igb5_glghs_topclust_mapped_pseudoGLY')
else:
    MasterDir_Native = os.path.expanduser('~/protein_model/native_struct/mapped')
    MasterDir_AATopClust = os.path.expanduser('~/protein_model/native_struct/ff96_igb5_glghs_topclust_mapped')

NativePdb = utils.parseNative(PdbName, MasterDir = MasterDir_Native)
try:
    AATopClustPdb = utils.parseNative(PdbName, MasterDir = MasterDir_AATopClust)
except IOError:
    print 'Utils Error: Requested Top clust pdb does not exist'
    AATopClustPdb = None
FF_File, FFMetadata = utils.parseGoFF(FFType)
Prefix = 'prot_' + PdbName

# backbone forcefield parameters
# Go parameters are loaded directly from the forcefield file
MinBondOrd = FFMetadata['MinBondOrd']
NKnot = FFMetadata['NKnot']
SPCut = FFMetadata['Cut']
NativeCut = FFMetadata['NativeCut']
Bonded_NCOSType = FFMetadata['Bonded_NCOSType']
NCOSType = FFMetadata['NCOSType']

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
NStepsEquil = 80000000              # 50 ns
NStepsProd  = 40000000              # 20 ns
NStepsSwap = 2000                   # 1 ps
StepFreq = int(NStepsProd / 10000)  # need 10000 frames, 2 ps
    
# REMD script template
mdstr = '''
#!/usr/bin/env python
import os, sys, numpy as np, time
import utils
from reasonable import const, config, cgmodel, md

# prefix
Prefix = "%(PREFIX)s"

# check if only analysis needs to be done
FullPrefix = os.path.join(os.getcwd(), Prefix)
isDone = utils.checkGoREMD(FullPrefix, [%(TEMPSET)3.2f])
if isDone: exit()
        
# set up config object
cfg = config.Config()

# pseudo glycine side chain
hasPseudoGLY = %(HASPSEUDOGLY)d
if hasPseudoGLY: cfg.AtomS['GLY'] = const.AtomS_GLY

# backbone settings (assumes single repulsive nonbonded BB-S potential by default)
cfg.MinBondOrd = %(MINBONDORD)d
cfg.NKnot = %(NKNOT)d
cfg.SPCut = %(SPCUT)g
cfg.Bonded_NCOSType = %(BONDED_NCOSTYPE)d
cfg.NCOSType = %(NCOSTYPE)d

# native contacts (spline)
cfg.NativeType = 1
cfg.NativeCut = %(NATIVECUT)g

# non-native contacts (WCA)
cfg.NonNativeType = 0

# timestep
cfg.TimeStep = %(TIMESTEP)g

# temp schedule
Temps = np.logspace(np.log10(%(TLOW)3.2f), np.log10(%(THIGH)3.2f), %(NREPLICA)d)
       
# compile Sys object
p, Sys = cgmodel.makeSplineGoSys(NativePdb = "%(NATIVEPDB)s", cfg = cfg, Prefix = Prefix, TempSet = %(TEMPSET)3.2f, DelFrac = %(DELFRAC)1.2f)

# load backbone forcefield file
cgmodel.loadParam(Sys, "%(FF_FILE)s")

# perform post-load Sys updates
cgmodel.UpdatePostLoad(Sys = Sys, cfg = cfg)

# show the forcefield
print Sys.ForceField.ParamString()

# compile REMD object
remd = md.REMD(p, Sys, Prefix = Prefix, Temps = Temps,
               NStepsMin = %(NSTEPSMIN)d, NStepsEquil = %(NSTEPSEQUIL)d, NStepsProd = %(NSTEPSPROD)d,
               NStepsSwap = %(NSTEPSSWAP)d, StepFreq = %(STEPFREQ)d)

# run REMD
t1 = time.time()
remd.runREMD()
t2 = time.time()
        
# reorder by temperature only at room temp
md.ReorderAll(ReorderTemps = [%(TEMPSET)3.2f], Prefix = "%(PREFIX)s", TempFile = "temps.txt", 
              NStepsEquil = %(NSTEPSEQUIL)d, NStepsProd = %(NSTEPSPROD)d,
              NStepsSwap = %(NSTEPSSWAP)d, StepFreq = %(STEPFREQ)d)
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
python ~/protein_model/analyze_go.py %(NATIVEPDB)s %(PREFIX)s ./ ./NativeAnalysis %(HASPSEUDOGLY)d
date
'''
# dict for filling md script template
d1 = {
      'PREFIX'          : Prefix,
      'NATIVEPDB'       : NativePdb,
        
      'FF_FILE'         : FF_File,
      'MINBONDORD'      : MinBondOrd,
      'NKNOT'           : NKnot,
      'SPCUT'           : SPCut,
      'BONDED_NCOSTYPE' : Bonded_NCOSType,
      'NCOSTYPE'        : NCOSType,
      
      'NATIVECUT'       : NativeCut,
      'DELFRAC'         : DelFrac,
      
      'HASPSEUDOGLY'    : hasPseudoGLY,

      'TLOW'            : TLow,
      'THIGH'           : THigh,
      'NREPLICA'        : NReplica,
      'TEMPSET'         : TempSet,
        
      'TIMESTEP'        : TimeStep,
      'NSTEPSMIN'       : NStepsMin,
      'NSTEPSEQUIL'     : NStepsEquil,
      'NSTEPSPROD'      : NStepsProd,
      'NSTEPSSWAP'      : NStepsSwap,
      'STEPFREQ'        : StepFreq
    }

# dict for filling job script template
d2 = {'JOBNAME': Prefix, 'NATIVEPDB': NativePdb, 'AATOPCLUSTPDB': AATopClustPdb, 'PREFIX': Prefix, 'HASPSEUDOGLY': hasPseudoGLY}

# run REMD job
mdscript = os.path.join(OutDir, 'remd.py')
jobscript = os.path.join(OutDir, 'remdjob.sh')
if not os.path.isfile(mdscript): file(mdscript, 'w').write(mdstr % d1)
file(jobscript, 'w').write(jobstr % d2)
os.chdir(OutDir)
os.system('qsub remdjob.sh')
os.chdir(CURRDIR)

