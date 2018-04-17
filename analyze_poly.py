#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle
import cgprotein as lib, measure
from utils import FMT

def NearestTemp(T, TList):
    ind = np.argmin(abs(TList - T))
    return TList[ind]

# user input
NativePdb = os.path.abspath(sys.argv[1])
Prefix = sys.argv[2]
Temp = float(sys.argv[3])
TrajDir = os.path.abspath(sys.argv[4]) if sys.argv > 4 else os.getcwd()
OutDir = os.path.abspath(sys.argv[5]) if sys.argv > 5 else os.getcwd()

# temps
TempFile = os.path.join(TrajDir, 'temps.txt') 
Temps = np.loadtxt(TempFile)
TempSet = NearestTemp(Temp, Temps)

# parse paths, etc
TrajPrefix = os.path.join(TrajDir, Prefix)
TrajFn = FMT['TRAJ'] % (TrajPrefix, TempSet)
if not os.path.isdir(OutDir): os.system('mkdir -p %s' % OutDir)
OutPrefix = os.path.join(OutDir, Prefix)

print '\n'
print 'ANALYSING POLYMER'
print '-----------------\n'

# compute object
print 'Creating Compute object'
calc = lib.Compute(NativePdb = NativePdb, TrajFn = TrajFn, Temp = TempSet, Prefix = OutPrefix)

# top cluster
print 'Clustering trajectory at Temp. %3.2f' % TempSet
calc.Cluster()

# ramachandran plots from native
print 'Generating Ramachandran plot from native'
p = lib.ProteinNCOS(NativePdb)
Phi, Psi = p.GetPhiPsi()
picklename = OutPrefix + '_native.rama.pickle'
with open(picklename, 'w') as of: pickle.dump((Phi, Psi), of)

# ramachandran plots from entire traj
print 'Generating Ramachandran plot from trajectory'
calc.RamaChandran()

# replica object
print 'Creating Replica object'
rep = lib.Replica(NativePdb = NativePdb, TrajPrefix = TrajPrefix, Prefix = OutPrefix,
                  TempSet = TempSet, OrderParams = ['U', 'Rg', 'RMSD'])

# folding curves
rep.FoldCurve()

# free energy surface
rep.PMF2D(O1 = 'Rg', O2 = 'RMSD')
