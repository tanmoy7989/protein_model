#!/usr/bin/env python

#### NOTE: PLEASE RUN FROM INSIDE THE INTENDED TARGET DIRECTORY ####

import os, sys, numpy as np, cPickle as pickle
from reasonable import config, cgmodel, md
import cgprotein as cglib
import utils
import argparse

# job template
jobscript = '''
#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -N %(JOBNAME)s
%(PARENV)s

date
%(CMD)s
date
'''

# Lammps reordering script
ReorderScript = '/home/cask0/home/tsanyal/mysoftware/mysim/reorderLammpsREMD.py'

# polypeptide analysis script
PolyAnalyzeScript = '/home/cask0/home/tsanyal/protein_model/analyze_poly.py'

# constants 
RoomTemp = 280.

# for histogramming etc
MinRg = 5.0
MaxRg = 15.0
MinRMSD = 5.0
MaxRMSD = 15.0

# MD iterations
NStepsMin = 10000 #10 ps
NStepsEquil = 400000000 #400 ns
NStepsProd = 100000000 #100 ns
NStepsSwap = 10000 #10 ps
NFrames = 20000
WriteFreq = int(NStepsProd / NFrames) # 5 ps

# get input
parser = argparse.ArgumentParser(description = 'CG Polypeptide simulations')
parser.add_argument('-r', '--resname', help = 'amino acid name')
parser.add_argument('-n', '--nmon', type = int, help = 'number of monomers')
parser.add_argument('-tf', '--tempfile', help = 'file with temp schedule')
parser.add_argument('-ff', '--forcefield', help = 'backbone forcefield type')
parser.add_argument('-o', '--outprefix', help = 'out prefix for all data')
parser.add_argument('-npdbs', '--nativepdbs', nargs = '+', help = 'native pdbs to compare against')
parser.add_argument('-op', '--operation', help = 'what to do?')

# parse input
args = parser.parse_args()
ResName = args.resname
NMon = args.nmon
TempFile = os.path.abspath(args.tempfile)
BBType = args.forcefield
Prefix = args.outprefix
Op = args.operation

# global data
OutDir = os.getcwd()
Temps = np.loadtxt(TempFile)
NProc = len(Temps)
InFile = Prefix + 'lammps.in'
LogFile = Prefix + 'lammps.log'
ScreenFile = Prefix + 'screen'
NativePdbs = []
if args.nativepdbs:
    NativePdbs = [os.path.abspath(i) for i in args.nativepdbs]


def RunSim():
    # set configuration object
    cfg = config.Config()
    # parse the forcefield
    ff_file, mdata = utils.parseBBFF(BBType = BBType)
    # backbone config
    cfg.MinBondOrd = mdata['MinBondOrd']
    cfg.NKnot = mdata['NKnot']
    cfg.SPCut = mdata['Cut']
    # backbone-sidechain config (single repulsive backbone-sidechain nonbonded potentials)
    cfg.Bonded_NCOSType = mdata['Bonded_NCOSType']
    # sidechain config
    cfg.NCOSType = mdata['NCOSType']                    
    # build the Sys object
    Seq = [ResName] * NMon
    TempSet = Temps[np.argmin(abs(Temps - RoomTemp))]
    p, Sys = cgmodel.makePolymerSys(Seq = Seq, cfg = cfg, TempSet = TempSet, Prefix = Prefix)   
    # load previously optimized bb interactions
    loadedPotentials = cgmodel.loadParam(Sys, ff_file)  
    # print the forcefield
    print '\n'
    print Sys.ForceField.ParamString()
    print '\n'
    # construct the REMD object
    remd = md.REMD(p = p, Sys = Sys, cfg = cfg, Prefix = Prefix, Temps = Temps,
                   NStepsMin = NStepsMin, NStepsEquil = NStepsEquil, NStepsProd = NStepsProd,
                   NStepsSwap = NStepsSwap, StepFreq = WriteFreq)
    # write lammps files
    print 'Writing Lammps Files...'
    remd.makeREMD()
    # write out jobscript
    parenv = '#$ -pe mpich %d' % NProc
    cmd = 'mpirun -np %d %s -partition %dx1 -in %s -log %s -screen %s' % (NProc, os.environ['LAMMPSEXEC'], NProc, InFile, LogFile, ScreenFile)
    d = {'NPROC': NProc, 'JOBNAME': Prefix + '_remd', 'PARENV': parenv, 'CMD': cmd}
    with open(os.path.join(OutDir, 'remd.sh'), 'w') as of: of.write(jobscript % d)
    return


def Reorder():
    parenv = '#$ -pe mpich %d' % NProc
    cmd = 'mpirun -np %d python %s %s %d %d -p %d -t %3.2f -v' % (NProc, ReorderScript, Prefix, NStepsSwap, WriteFreq, NStepsProd, RoomTemp)
    # write job script
    d = {'NPROC': NProc, 'JOBNAME': Prefix + '_reorder', 'PARENV': parenv, 'CMD': cmd}
    with open(os.path.join(OutDir, 'reorder.sh'), 'w') as of: of.write(jobscript % d)
    return


def Analyze():
    for i, x in enumerate(NativePdbs):
        this_OutDir = os.path.join(OutDir, 'Analysis%d' % (i+1))
        this_NativePdb = os.path.join(this_OutDir, x.split('/')[-1])
        this_T = RoomTemp
        os.system('mkdir -p %s' % this_OutDir)
        os.system('cp %s %s' % (x, this_NativePdb))
        # write analysis jobs
        parenv = ''
        cmd = 'python %s %s %s ../ ./ %3.2f %d %d %d' % (PolyAnalyzeScript, this_NativePdb.split('/')[-1], Prefix, 
                                                        this_T, NStepsProd, NStepsSwap, WriteFreq) 
        d = {'JOBNAME': Prefix + '_analyze', 'PARENV': parenv, 'CMD': cmd}
        with open(os.path.join(this_OutDir, 'analyze.sh'), 'w') as of: of.write(jobscript % d)
    return


#### MAIN ####
if Op == 'runsim':
    RunSim()
elif Op == 'reorder':
    Reorder()
elif Op == 'analyze':
    if NativePdbs: Analyze()
else:
    print 'ERROR: Uknown operation'
    exit()
