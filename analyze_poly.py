#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle
import cgprotein as lib, measure
from utils import FMT
import argparse

# min and max ranges for Rg and RMSD
# based on prior observations
MINRg = 5. #5.0
MAXRg = 15.0 #15.0
MINRMSD = 0. #0.0
MAXRMSD = 15. #15.0 # 15.0 for val

def NearestTemp(T, TList):
    ind = np.argmin(abs(TList - T))
    return TList[ind]

# user input
NativePdb = os.path.abspath(sys.argv[1])
Prefix = sys.argv[2]
TrajDir = os.path.abspath(sys.argv[3])
OutDir = os.path.abspath(sys.argv[4])
Temp = float(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] else None
NStepsProd = int(sys.argv[6]) if len(sys.argv) > 6 else None
NStepsSwap = int(sys.argv[7]) if len(sys.argv) > 7 else None
WriteFreq = int(sys.argv[8]) if len(sys.argv) > 8 else None

# trajectory locations
TrajPrefix = os.path.join(TrajDir, Prefix)
TempFile = os.path.join(TrajDir, 'temps.txt') 
Temps = np.loadtxt(TempFile)
if not Temp is None:
    TempSet = NearestTemp(Temp, Temps)
    TrajFn = FMT['TRAJ'] % (TrajPrefix, TempSet)
else:
    TempSet = None
    TrajFn = None

# output locations
OutPrefix = os.path.join(OutDir, Prefix)
if not os.path.isdir(OutDir): os.system('mkdir -p %s' % OutDir)

print '\n'
print 'ANALYZING POLYPEPTIDE'
print '====================='

# compute object
print '\nCreating Compute object'
calc = lib.Compute(NativePdb = NativePdb, TrajFn = TrajFn, Temp = TempSet, Prefix = OutPrefix)

# replica object
print '\nCreating Replica object'
rep = lib.Replica(NativePdb = NativePdb, TrajPrefix = TrajPrefix, Prefix = OutPrefix,
                  TempSet = TempSet, OrderParams = ['U', 'Rg', 'RMSD'],
                  NStepsProd = NStepsProd, NStepsSwap = NStepsSwap, WriteFreq = WriteFreq)

def Cluster():
    # top cluster
    print '\nClustering trajectory at Temp. %3.2f' % TempSet
    calc.Cluster()

def RamaChandran():
    # ramachandran plots from native
    print '\nGenerating Ramachandran plots from native'
    p = lib.ProteinNCOS(NativePdb)
    Phi, Psi = p.GetPhiPsi()
    picklename = OutPrefix + '_native.rama.pickle'
    with open(picklename, 'w') as of: pickle.dump((Phi, Psi), of)
    # ramachandran plots from entire traj
    print '\nRamachandran plots at Temp = %3.2fK' % TempSet
    calc.RamaChandran()

def FoldCurve():
    # folding curves
    print '\nCalculating folding curve from RMSD'
    rep.FoldCurve(MIN = MINRMSD, MAX = MAXRMSD)

def FieldMap():
    # free energy surface
    print '\nCalculating 2D PMF (Rg, RMSD) at %3.2fK' % TempSet
    rep.PMF2D(O1 = 'Rg', O2 = 'RMSD',
              MIN1 = MINRg, MAX1 = MAXRg, MIN2 = MINRMSD, MAX2 = MAXRMSD)
    

#### MAIN ####
FoldCurve()
FieldMap()
if not Temp is None:
    RamaChandran()
    #Cluster()
