#!/usr/bin/env python

import sim
from const import *
import parsestruct as ps

''' Builds sidechain-sidechain nonbonded potentials.
    Functions are named numerically based on types in settings.
    Inputs are the NCOS cg protein and Sys objects
    LJ and WCA native and non-native energies should always be input in kT units'''

Verbose = True

class P_Sidechain(object):
    def __init__(self, p, Sys, cfg):
        ''' initialized a class for backbone potentials
            cfg is a configuration object'''
        # extract relevant settings
        self.MinBondOrd = cfg.MinBondOrd
        self.NKnot = cfg.NKnot
        self.SPCut = cfg.SPCut
        # also extract the complete config object (useful for Go models)
        self.cfg = cfg
        # extract cg protein and Sys objects
        self.p = p
        self.Sys = Sys
        # unpack sidechains
        if cfg.SSRefType == 'name':
            self.AtomS = p.AtomSbyRes
        elif cfg.SSRefType == 'number':
            self.AtomS = p.AtomSbyNum
        else:
            print 'Error: Unknown sidechain reference type, or not implemented yet'
            exit()

    def SS_0(self):
        '''21 alphabet nonbonded sidechain-sidechain spline potentials
        assumes sidechain reference by residue name'''
        if Verbose: print 'Generating 21-alphabet sidechain-sidechain nonbonded potentials. For splines, MinBondOrd = %d, %d knots, SPCut = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        Filter_SS_p = {}
        for i in range(self.p.NResTypes - 1):
            for j in range(i+1, self.p.NResTypes):
                r1 = self.p.ResTypes[i]
                r2 = self.p.ResTypes[j]
                if (self.AtomS[r1] is None) or (self.AtomS[r2] is None): continue
                Filter_SS_p[(r1, r2)] = sim.atomselect.PolyFilter([self.AtomS[r1], self.AtomS[r2]], MinBondOrd = self.MinBondOrd)
        Pair_SS = {}
        for i in range(self.p.NResTypes - 1):
            for j in range(i+1, self.p.NResTypes):
                r1 = self.p.ResTypes[i]
                r2 = self.p.ResTypes[j]
                if (self.AtomS[r1] is None) or (self.AtomS[r2] is None): continue
                Pair_SS[(r1, r2)] = sim.potential.PairSpline(self.Sys, Filter = Filter_SS_p[(r1,r2)], Label = 'NonBondS_%s_S_%s' % (r1,r2), NKnot = self.NKnot, Cut = self.SPCut)
        # populate
        ff = Pair_SS.values()
        return ff

    def SS_1(self):
        '''1 alphabet nonbonded sidechain-sidechain spline potentials
        assumes sidechain reference by residue number'''
        if Verbose: print 'Generating 1-alphabet sidechain-sidechain nonbonded potentials. For splines, MinBondOrd = %d, %d knots, SPCut = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        FilterS = sim.atomselect.Filter([self.AtomS[i] for i in range(self.p.NRes) if not self.AtomS[i] is None])
        Filter_SS_p = sim.atomselect.PolyFilter([FilterS, FilterS], MinBondOrd = self.MinBondOrd)
        Pair_SS = sim.potential.PairSpline(self.Sys, Filter = Filter_SS_p, Label = 'NonBondSS', NKnot = self.NKnot, Cut = self.SPCut)
        # populate
        ff = [Pair_SS]
        return ff

    def Go_native_0(self, ContactDict):
        '''1 alphabet Go like nonbonded native sidechain-sidechain LJ potentials
            with single sigma and epsilon'''
        if Verbose: print 'Generating sidechain-sidechain nonbonded LJ Go potentials between native contacts.'
        NativePairs, NonNativePairs, Topo2AID = ps.makeMatrix(self.p, ContactDict, self.Sys)
        Filter_Native = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = NativePairs, MinBondOrd = self.MinBondOrd)
        if self.cfg.NativeSigma is None:
            self.cfg.NativeSigma = ContactDict['d_native'].min() * (2**(-1/6.))
        NativeSigma = self.cfg.NativeSigma
        NativeEpsilon = self.cfg.NativeEpsilon * (kB * RoomTemp)
        NativeCut = self.cfg.NativeCut
        P = sim.potential.LJ(self.Sys, Filter = Filter_Native, Cut = NativeCut, Sigma = NativeSigma, Epsilon = NativeEpsilon, Shift = True, Label = 'NonBondNative')
        # populate
        if Verbose: print ' Using Sigma = %2.2f A, Epsilon = %2.2f kT, Cutoff = %2.2f A' % (NativeSigma, NativeEpsilon / (kB * RoomTemp), NativeCut)
        ff = [P]
        return ff

    def Go_native_1(self, ContactDict):
        '''1 alphabet Go like nonbonded native sidechain-sidechain spline potentials'''
        if Verbose: print 'Generating sidechain-sidechain nonbonded spline Go potentials between native contacts.'
        NativePairs, NonNativePairs, Topo2AID = ps.makeMatrix(self.p, ContactDict, self.Sys)
        Filter_Native = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = NativePairs, MinBondOrd = self.MinBondOrd)
        NativeCut = self.cfg.NativeCut
        NativeNKnot = self.cfg.NativeNKnot
        P = sim.potential.PairSpline(self.Sys, Filter = Filter_Native, Cut = NativeCut, NKnot = NativeNKnot, Label = 'NonBondNative')
        # populate
        if Verbose: print ' Using %d knots, Cutoff = %2.2f A' % (NativeNKnot, NativeCut)
        ff = [P]
        return ff

    def Go_native_2(self, ContactDict):
        '''1 alphabet Go like nonbonded native sidechain-sidechain harmonic restraints'''
        if Verbose: print 'Generating sidechain-sidechain nonbonded harmonic restraints between native contacts'
        NativePairs, NonNativePairs, Topo2AID_Native, Topo2AID_NonNative = ps.makeMatrix(self.p, ContactDict, self.Sys)
        if self.cfg.NativeFConst is None:
            self.cfg.NativeFConst = 0.5 * 1.0 / (2. * self.cfg.HarmonicFluct**2) # since sim uses k and not k/2
        NativeFConst = self.cfg.NativeFConst * (kB * RoomTemp)
        if Verbose: print ' Using FConst = %2.2f kT, HarmonicFluct = %2.2f A' % (NativeFConst / (kB * RoomTemp), self.cfg.HarmonicFluct) 
        P = []
        for k, (i,j) in enumerate(ContactDict['c_native']):
            # ignore residue pairs (without sidechains) not included in Topo2AID, while making filters
            if not Topo2AID_Native.__contains__( (i,j) ): continue
            if Verbose: print ' Applying harmonic restraint to : (%3d, %3d), (%3s, %3s) d0 = %2.2f A ' % (i, j, self.p.Seq[i], self.p.Seq[j], ContactDict['d_ss_native'][k])
            m, n = Topo2AID_Native[(i,j)]
            NAID = self.Sys.World.NAID
            this_NativePair = np.zeros([NAID, NAID], int)
            this_NativePair[m,n] = 1
            this_NativePair[n,m] = 1
            this_d0 = ContactDict['d_ss_native'][k]
            this_Filter = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = this_NativePair, MinBondOrd = self.MinBondOrd)
            this_P = sim.potential.Bond(self.Sys, FConst = NativeFConst, Dist0 = this_d0, Filter = this_Filter, Label = 'NonBondNative_%d_%d' % (i,j))
            P.append(this_P)
        # populate
        ff = P
        return ff

    def Go_nonnative_0(self, ContactDict):
        '''1 alphabet Go like nonbonded non-native sidechain-sidechain WCA potentials
            with single sigma and epsilon'''
        if Verbose: print 'Generating sidechain-sidechain nonbonded WCA potentials between non-native contacts.'
        NativePairs, NonNativePairs, Topo2AID_Native, Topo2AID_NonNative = ps.makeMatrix(self.p, ContactDict, self.Sys)
        Filter_NonNative = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = NonNativePairs, MinBondOrd = self.MinBondOrd)
        if self.cfg.NonNativeSigma is None:
            self.cfg.NonNativeSigma = ContactDict['d_native'].min() ** (2**(-1/6.))
        NonNativeSigma = self.cfg.NonNativeSigma
        NonNativeEpsilon = self.cfg.NonNativeEpsilon * (kB * RoomTemp)
        if self.cfg.NonNativeCut is None:
            self.cfg.NonNativeCut = NonNativeSigma * (2**(1/6.))
        NonNativeCut = self.cfg.NonNativeCut
        P = sim.potential.LJ(self.Sys, Filter = Filter_NonNative, Cut = NonNativeCut, Sigma = NonNativeSigma, Epsilon = NonNativeEpsilon, Shift = True, Label = 'NonBondNonNative')
        # populate
        if Verbose: print ' Using Sigma = %2.2f A, Epsilon = %2.2f kT, Cutoff = %2.2f A' % (NonNativeSigma, NonNativeEpsilon / (kB * RoomTemp), NonNativeCut)
        ff = [P]
        return ff

    def Go_nonnative_1(self, ContactDict):
        '''1 alphabet Go like nonbonded non-native sidechain-sidechain spline potentials'''
        if Verbose: print 'Generating sidechain-sidechain nonbonded spline potentials between non-native contacts.'
        NativePairs, NonNativePairs, Topo2AID_Native, Topo2AID_NonNative = ps.makeMatrix(self.p, ContactDict, self.Sys)
        Filter_NonNative = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = NonNativePairs, MinBondOrd = self.MinBondOrd)
        if self.cfg.NonNativeCut is None:
            # cutoff commensurate with a WCA whose corresponding LJ
            # has a min at min native-contact distance
            self.cfg.NonNativeCut = ContactDict['d_native'].min()
        NonNativeCut = self.cfg.NonNativeCut
        NonNativeNKnot = self.cfg.NonNativeNKnot
        P = sim.potential.PairSpline(self.Sys, Filter = Filter_NonNative, Cut = NonNativeCut, NKnot = NonNativeNKnot, Label = 'NonBondNonNative')
        # populate
        if Verbose: print ' Using %d knots, Cutoff = %2.2f A' % (NonNativeKnot, NonNativeCut)
        ff = [P]
        return ff




