#!/usr/bin/env python

import os, sys, numpy as np

ExpName = sys.argv[1]

# harmonic restraint experiment
if ExpName == 'harmonic_restraint_experiment':
    pset_type = sys.argv[2]
    FFType = sys.argv[3]
    OutDir = os.path.abspath(sys.argv[4])
    MetaDataFile = os.path.expanduser('~/protein_model/native_struct/native_metadata.txt')
    data = eval(file(MetaDataFile).read())
    PdbNames = data[pset_type]
    Script = os.path.expanduser('~/protein_model/%s.py' % ExpName)
    for i, p in enumerate(PdbNames):
        cmd = 'python %s %s %s %s' % (Script, p, FFType, OutDir)
        os.system(cmd)

# native-spline experiment
if ExpName == 'spline_experiment':
    pset_type = sys.argv[2]
    FFType = sys.argv[3]
    OutDir = os.path.abspath(sys.argv[4])
    NReplica = int(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] else 8
    MetaDataFile = os.path.expanduser('~/protein_model/native_struct/native_metadata.txt')
    data = eval(file(MetaDataFile).read())
    PdbNames = data[pset_type]
    Script = os.path.expanduser('~/protein_model/%s.py' % ExpName)
    for i, p in enumerate(PdbNames):
        cmd = 'python %s %s %s %s %d' % (Script, p, FFType, OutDir, NReplica)
        os.system(cmd)

# erodenative(-spline) experiment
if ExpName == 'erodenative_experiment':
    PdbName = sys.argv[2]
    FFType = sys.argv[3]
    OutDir = os.path.abspath(sys.argv[4])
    # del fractions
    DelFracList = []
    if len(sys.argv) >= 5:
        for i, frac in enumerate(sys.argv[5:]):
            DelFracList.append(float(frac))
    np.savetxt(os.path.join(OutDir, 'erodefrac.txt'), np.array(DelFracList), fmt = '%1.2f')
    if not DelFracList: DelFracList = [0.0]
    MetaDataFile = os.path.expanduser('~/protein_model/native_struct/native_metadata.txt')
    data = eval(file(MetaDataFile).read())
    Script = os.path.expanduser('~/protein_model/%s.py' % ExpName)
    for i, frac in enumerate(DelFracList):
        cmd = 'python %s %s %s %s % 1.2f' % (Script, PdbName, FFType, OutDir, frac)
        os.system(cmd)
