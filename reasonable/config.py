
#!/usr/bin/env python

''' simple class template to supply configs to the model
    '''

import numpy as np
from const import *

nice = lambda x : 0.0 if x is None else x

# create a config object and populate with defaults
class Config(object):
    def __init__(self):
        for k,v in DEFAULTS.iteritems():
            self.__setattr__(k, v)

    def __repr__(self):
        s1 = HELPSTR
        s2 = 'Current config:\n'
        for k, v in DEFAULTS.iteritems(): s2 += '%s: %s\n' % (k, self.__dict__[k])
        return s1 + '\n' + s2

    def Status(self):
        Backbone_status = 'Using MinBondOrd = %d, %d spline knots, Cutoff = %g A for backbone potentials, hasSpecialBBTorsions (for GLY and PRO) = %d' % (self.MinBondOrd, self.NKnot, self.SPCut, self.hasSpecialBBTorsions)
        Bonded_NCOS_status = {-1 : 'Ignoring backbone-sidechain bonded interactions (this can be fatal!!)',
                               0 : 'Using 21-alphabet backbone-sidechain bonded interactions',
                               1 : 'Using 1-alphabet backbone-sidechain bonded interactions',
                             }
        NCOS_status = {-1 : 'Ignoring backbone-sidechain nonbonded spline interactions',
                        0 : 'Using 21-alphabet backbone-sidechain nonbonded spline interactions',
                        1 : 'Using 1-alphabet backbone-sidechain nonbonded spline interactions',
                        2 : 'Using 1-alphabet constant soft core repulsive backbone-sidechain nonbonded interactions'
                      }

        Native_status =   {-1 : 'Ignoring native nonbonded interactions (not using a Go model)',
                            0 : 'Using single LJ native nonbonded interactions with Sigma = %2.2f A, Epsilon = %2.2f kT, Cutoff = %2.2f A' % (nice(self.NativeSigma), nice(self.NativeEpsilon), nice(self.NativeCut)),
                            1 : 'Using splined native nonbonded interactions with %d knots, Cutoff = %2.2f A' % (nice(self.NativeNKnot), nice(self.NativeCut)),
                            2 : 'Using harmonic restraint native nonbonded interactions with NativeFconst = %2.2f kT, HarmonicFluct = %2.2f A' % (nice(self.NativeFConst), nice(self.HarmonicFluct))
                          }

        NonNative_status = {-1 : 'Ignoring non-native nonbonded interactions',
                             0 : 'Using single WCA non-native nonbonded interactions with Sigma = %2.2f A, Epsilon = %2.2f kT, Cutoff = %2.2f A' % (nice(self.NonNativeSigma), nice(self.NonNativeEpsilon), nice(self.NonNativeCut)),
                             1 : 'Using splined non-native nonbonded interactions with %d knots, Cutoff = %2.2f A' % (nice(self.NonNativeNKnot), nice(self.NonNativeCut))
                            }
        s = '''
(R)elative (E)ntropy (A)ssisted, (S)tructure (O)ptimized (N)o (A)dded (B)informatic (LE)xicon
---------------------------------------------------------------------------------------------\n'''
        s += Backbone_status + '\n' + NCOS_status[self.NCOSType] + '\n' + Native_status[self.NativeType] + '\n' + NonNative_status[self.NonNativeType] + '\n'
        print s


