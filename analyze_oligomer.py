#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle, shelve
import cgprotein as cg
from utils import FMT
import sim, pickleTraj

# room temp.
RoomTemp = 300.0

# user input
NativePdb = os.path.abspath(sys.argv[1])
Prefix = sys.argv[2]
TrajDir = os.path.abspath(sys.argv[3])
OutDir = os.path.abspath(sys.argv[4])
TempSet = float(sys.argv[5])
hasPseudoGLY = int(sys.argv[6]) if len(sys.argv) > 6 else False

# get Traj
if not os.path.isdir(OutDir): os.system('mkdir -p %s' % OutDir)
OutPrefix = os.path.join(OutDir, Prefix)
TrajPrefix = os.path.join(TrajDir, Prefix)
TrajFn = FMT['TRAJ'] % (TrajPrefix, TempSet)

# fancy intro
print 'ANALYSING MULTI-PROTEIN FOR %s at %3.2f K' % (Prefix, TempSet)
print '----------------------------------------------------------------'

# set up Compute object
print 'Creating Compute object'
calc = cg.Compute(NativePdb = NativePdb, TrajFn = TrajFn, Temp = TempSet, Prefix = OutPrefix, hasPseudoGLY = hasPseudoGLY)

# cluster the room temp traj
print 'Clustering trajectory'
calc.Cluster()

# top-cluster contact analysis
Trj = pickleTraj(TrajFn)
BoxL = Trj.FrameData['BoxL']
ClustPdb = FMT['CLUSTPDB'] % (calc.Prefix, calc.Temp)
p = cg.ProteinNCOS(ClustPdb, Model = 1)

# top cluster contact map
print 'Calculating contact map for top cluster'
cmap, cdist = p.GetResContacts(BoxL = BoxL)
picklename = FMT['RESCONTACTS'] % (calc.Prefix + '_topclust', calc.Temp)
with open(picklename, 'w') as of: pickle.dump(cmap, of)

# top cluster dihedrals
print 'Calculating dihedral angles for top cluster'
Phi, Psi = p.GetPhiPsi(BoxL = BoxL)
picklename = FMT['RAMA'] % (calc.Prefix + '_topclust', calc.Temp)
with open(picklename, 'w') as of: pickle.dump((Phi, Psi), of)

# traj contact map
print 'Calculating contact map from entire trajectory'
calc.ResContacts_frame()

# ramachandran plot
print 'Calculating Ramachandran plot from entire trajectory'
calc.RamaChandran()
