#!/usr/bin/env python

import os, sys, numpy as np

MasterDir = os.path.abspath(sys.argv[1])
Prefix = os.path.abspath(sys.argv[2])
hasPseudoGLY = int(sys.argv[3]) if len(sys.argv) > 3 else 0
FragLen = 17 if not len(sys.argv) > 4 else int(sys.argv[4])

TrajDir = os.path.join(MasterDir, 'r1-%da/data' % FragLen)
TempFile = os.path.join(MasterDir, 'r1-%da/data/temps.txt' % FragLen)
InPdb = os.path.join(MasterDir, 'r1-%da/data/0.current.pdb' % FragLen)

Temps = np.loadtxt(TempFile)
for i, T in enumerate(Temps):
    print 'Converting traj at temp %3.2f...\n' % T
    CGPrefix = '%s.%3.2f' % (Prefix, T)
    AATraj = os.path.join(TrajDir, '%d.mdtrj.crd.gz' % i)
    PrmTop = os.path.join(TrajDir, '%d.prmtop.parm7' % i)
    AmberEne = os.path.join(TrajDir, '%d.mdene.txt.gz' % i)
    LastNFrames = 20000
    cmd = 'python ~/protein_model/reasonable/mapNCOS.py map %s %s %d %s %s %s %d' % (InPdb, CGPrefix, hasPseudoGLY, AATraj, PrmTop, AmberEne, LastNFrames)
    os.system(cmd)

    
