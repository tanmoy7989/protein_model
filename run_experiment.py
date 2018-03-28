#!/usr/bin/env python

import os, sys

ExpName = sys.argv[1]
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
