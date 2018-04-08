#!/usr/bin/env python

import copy
import sim
from const import *
import parsestruct as ps

''' Builds sidechain-sidechain nonbonded potentials.
    Functions are named numerically based on types in settings.
    Inputs are the NCOS cg protein and Sys objects
    LJ and WCA native and non-native energies should always be input in kT units'''

Verbose = True

class P_Sidechain(object):
    def __init__(self, p, Sys, cfg = None, ContactDict = None):
        ''' initialized a class for backbone potentials
            cfg is a configuration object'''
        if cfg is None: cfg = p.cfg
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
        self.SSRefType = cfg.SSRefType
        if self.SSRefType == 'name':
            self.AtomS = p.AtomSbyRes
        elif self.SSRefType == 'number':
            self.AtomS = p.AtomSbyNum
        else:
            print 'Error: Unknown sidechain reference type, or not implemented yet'
            exit()
        # make native contact-filters if contactdict is supplied
        self.ContactDict = ContactDict
        if not self.ContactDict is None:
            NativePairs, NonNativePairs, Topo2AID_Native, Topo2AID_NonNative = ps.makeMatrix(self.p, self.ContactDict, self.Sys)
            self.NativePairs = NativePairs
            self.NonNativePairs = NonNativePairs
            self.Topo2AID_Native = Topo2AID_Native
            self.Topo2AID_NonNative = Topo2AID_NonNative

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
        assumes sidechain reference by residue name or number'''
        if Verbose: print 'Generating 1-alphabet sidechain-sidechain nonbonded potentials. For splines, MinBondOrd = %d, %d knots, SPCut = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        if self.SSRefType == 'number':
            FilterS = sim.atomselect.Filter([self.AtomS[i] for i in range(self.p.NRes) if not self.AtomS[i] is None])
        elif self.SSRefType == 'name':
            FilterS = sim.atomselect.Filter([self.AtomS[r] for r in self.p.ResTypes if not self.AtomS[r] is None])
        Filter_SS_p = sim.atomselect.PolyFilter([FilterS, FilterS], MinBondOrd = self.MinBondOrd)
        Pair_SS = sim.potential.PairSpline(self.Sys, Filter = Filter_SS_p, Label = 'NonBondSS', NKnot = self.NKnot, Cut = self.SPCut)
        # populate
        ff = [Pair_SS]
        return ff

    def Go_native_0(self, Sigma = None, Epsilon = None, Cut = None):
        '''1 alphabet Go like nonbonded native sidechain-sidechain LJ potentials
            with single sigma and epsilon'''
        Filter_Native = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = self.NativePairs, MinBondOrd = self.MinBondOrd)
        if Verbose: print 'Generating sidechain-sidechain nonbonded LJ Go potentials between native contacts.'
        if Sigma is None:
            if Verbose: print ' Native sigma not provided. Reverting to Auto-Sigma initially'
            Sigma = self.ContactDict['d_native'].min() * (2**(-1/6.))
        if Cut is None: Cut = self.cfg.NativeCut
        if Epsilon is None: Epsilon = self.cfg.NativeEpsilon
        P = sim.potential.LJ(self.Sys, Filter = Filter_Native, Cut = Cut, Sigma = Sigma, Epsilon = Epsilon, Shift = True, Label = 'NonBondNative')
        # populate
        if Verbose: print ' Using Sigma = %2.2f A, Epsilon = %2.2f kT, Cutoff = %2.2f A' % (Sigma, Epsilon / (kB * RoomTemp), Cut)
        ff = [P]
        return ff

    def Go_native_1(self, Cut = None):
        '''1 alphabet Go like nonbonded native sidechain-sidechain spline potentials'''
        Filter_Native = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = self.NativePairs, MinBondOrd = self.MinBondOrd)
        if Verbose: print 'Generating sidechain-sidechain nonbonded spline Go potentials between native contacts.'
        if Cut is None: Cut = self.cfg.NativeCut
        P = sim.potential.PairSpline(self.Sys, Filter = Filter_Native, Cut = Cut, NKnot = self.NKnot, Label = 'NonBondNative')
        # populate
        if Verbose: print ' Using Cutoff = %2.2f A' % Cut
        ff = [P]
        return ff

    def Go_native_2(self, FConst = None, HarmonicFluct = None):
        '''1 alphabet Go like nonbonded native sidechain-sidechain harmonic restraints'''
        if Verbose: print 'Generating sidechain-sidechain bonded harmonic restraints between native contacts'
        # determine harmonic fluct first
        if HarmonicFluct is None: HarmonicFluct = self.cfg.NativeHarmonicFluct
        # now determine force constant
        if FConst is None: FConst = 0.5 * 1.0 / (2. * HarmonicFluct**2) # since sim uses k and not k/2
        FConst *= (kB * RoomTemp)
        if Verbose: print ' Using FConst = %2.2f kT, HarmonicFluct = %2.2f A' % (FConst / (kB * RoomTemp), HarmonicFluct) 
        P = []
        for k, (i,j) in enumerate(self.ContactDict['c_native']):
            # ignore residue pairs not included in Topo2AID, while making filters
            # this ignores glycines if they are not included as pseudo side chains
            if not self.Topo2AID_Native.__contains__( (i,j) ): continue
            if Verbose: print ' Applying harmonic restraint to : (%3d, %3d), (%3s, %3s) d0 = %2.2f A ' % (i, j, self.p.Seq[i], self.p.Seq[j], self.ContactDict['d_ss_native'][k])
            m, n = self.Topo2AID_Native[(i,j)]
            # create filters
            NAID = self.Sys.World.NAID
            this_NativePair = np.zeros([NAID, NAID], int)
            this_NativePair[m,n] = 1
            this_NativePair[n,m] = 1
            this_d0 = self.ContactDict['d_ss_native'][k]
            this_Filter = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = this_NativePair, Bonded = True)
            # create native bonded potential
            this_P = sim.potential.Bond(self.Sys, FConst = FConst, Dist0 = this_d0, Filter = this_Filter, Label = 'Restraint_%d_%d' % (i,j))
            P.append(this_P)
        # populate
        ff = P
        return ff

    def Go_nonnative_0(self, Sigma = None, Epsilon = None, Cut = None):
        '''1 alphabet Go like nonbonded non-native sidechain-sidechain WCA potentials
            with single sigma and epsilon'''
        Filter_NonNative = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = self.NonNativePairs, MinBondOrd = self.MinBondOrd)
        if Verbose: print 'Generating sidechain-sidechain nonbonded WCA potentials between non-native contacts.'
        if Sigma is None:
            if Verbose: print ' Non-native sigma not provided. Reverting to Auto-Sigma initially'
            Sigma = self.ContactDict['d_native'].min() ** (2**(-1/6.))
        if Cut is None: Cut = Sigma * (2**(1/6.)) # WCA
        if Epsilon is None: Epsilon = self.cfg.NonNativeEpsilon
        P = sim.potential.LJ(self.Sys, Filter = Filter_NonNative, Cut = Cut, Sigma = Sigma, Epsilon = Epsilon, Shift = True, Label = 'NonBondNonNative')
        # populate
        if Verbose: print ' Using Sigma = %2.2f A, Epsilon = %2.2f kT, Cutoff = %2.2f A' % (Sigma, Epsilon / (kB * RoomTemp), Cut)
        ff = [P]
        return ff

    def Go_nonnative_1(self, Cut = None):
        '''1 alphabet Go like nonbonded non-native sidechain-sidechain spline potentials'''
        Filter_NonNative = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = self.NonNativePairs, MinBondOrd = self.MinBondOrd)
        if Verbose: print 'Generating sidechain-sidechain nonbonded spline potentials between non-native contacts.'
        # cutoff commensurate with a WCA whose corresponding LJ has a min at min native-contact distance
        if Cut is None: Cut = seld.ContactDict['d_native'].min()
        P = sim.potential.PairSpline(self.Sys, Filter = Filter_NonNative, Cut = Cut, NKnot = self.NKnot, Label = 'NonBondNonNative')
        # populate
        if Verbose: print ' Using Cutoff = %2.2f A' % Cut
        ff = [P]
        return ff




