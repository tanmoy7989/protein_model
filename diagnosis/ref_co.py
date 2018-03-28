#!/usr/bin/env python

import numpy as np
import utils, cgprotein

polymers = ['ala15', 'leu15', 'val15']

for p in polymers:
    dict_AA = utils.parseAA(PolyPrefix = p)
    nativepdb = dict_AA['Pdb']
    trajfn = dict_AA['Traj']
    Temp = dict_AA['Temp']
    calc = cgprotein.Compute(NativePdb = nativepdb, TrajFn = trajfn, Temp = Temp, Prefix = p)
    nativeco, cohist = calc.CompareCO()
    avgco = np.sum(cohist[0] * cohist[1]) / np.sum(cohist[1])
    print 'Polymer: %s, Native CO: %g, TrajAvg CO: %g' % (p, nativeco, avgco)

