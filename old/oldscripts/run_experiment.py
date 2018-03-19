#!/usr/bin/env python

import os, sys, utils

print 'Please make sure you are running this from inside the experiment master directory'

Ex = sys.argv[1]
Op = sys.argv[2]
Pset_type = sys.argv[3]
if Op == 'run': FFType = sys.argv[4]

MdataFile = os.path.expanduser('~/Go/native_struct/native_metadata.txt')
Mdata = eval(file(MdataFile).read())
pset = Mdata[Pset_type]

# run templates
remd_template = 'python ~/Go/%s.py run %s %s ./'
plot_template = 'python ~/Go/%s.py plot %s %s ./'

if Op == 'run':
    for i, p in enumerate(pset):
        cmd = remd_template % (Ex, p, FFType)
        os.system(cmd)

if Op == 'plot':
    cmd = plot_template % (Ex, Pset_type)
    os.system(cmd)
