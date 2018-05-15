#!/usr/bin/env python

import sim
from const import *

''' Builds backbone-sidechain bonded and nonbonded potentials.
    Functions are named numerically based on types in settings.
    Inputs are the NCOS cg protein and Sys objects'''

Verbose = True

class P_Backbone_Sidechain(object):
    def __init__(self, p, Sys, cfg = None):
        ''' initialized a class for backbone potentials
            cfg is a configuration object'''
        if cfg is None: cfg = p.cfg
        # extract relevant settings
        self.MinBondOrd = cfg.MinBondOrd
        self.NKnot = cfg.NKnot
        self.SPCut = cfg.SPCut
        # extract cg protein and Sys objects
        self.p = p
        self.Sys = Sys
        # unpack sidechain atoms
        self.SSRefType = cfg.SSRefType
        if self.SSRefType == 'name':
            self.AtomS = p.AtomSbyRes
        elif self.SSRefType == 'number':
            self.AtomS = p.AtomSbyNum
        else:
            print 'Error: Unknown sidechain reference type, or not implemented yet'
            exit()

    def BB_S_Bonded_0(self):
        '''21 alphabet bonded backbone-sidechain potentials
        assumes sidechain reference by residue name'''
        # create bonded potentials
        if Verbose: print 'Generating 21-alphabet bonded backbone-sidechain potentials. For splines, MinBondOrd = %d, %d knots, Cutoff = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        FilterC = sim.atomselect.Filter([AtomC, AtomC_GLY, AtomC_PRO])
        Filter_CS = dict( (r, sim.atomselect.PolyFilter([FilterC, self.AtomS[r]], Bonded = True)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Bond_CS = dict( (r, sim.potential.Bond(self.Sys, Filter = Filter_CS[r], Label = 'BondCS_%s' % r, Dist0 = 4.0, FConst = 1.0)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        # create angle potentials
        Filter_NCS = dict( (r, sim.atomselect.PolyFilter([AtomN, FilterC, self.AtomS[r]], Bonded = True)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Filter_SCO = dict( (r, sim.atomselect.PolyFilter([self.AtomS[r], FilterC, AtomO], Bonded = True)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Angle_NCS = dict( (r, sim.potential.AngleSpline(self.Sys, Filter = Filter_NCS[r], Label = 'AngleNCS_%s' % r, NKnot = self.NKnot)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Angle_SCO = dict( (r, sim.potential.AngleSpline(self.Sys, Filter = Filter_SCO[r], Label = 'AngleSCO_%s' % r, NKnot = self.NKnot)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        # create torsion potentials
        Filter_ONCS = dict( (r, sim.atomselect.PolyFilter([AtomO, AtomN, FilterC, self.AtomS[r]], Bonded = True)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Filter_SCON = dict( (r, sim.atomselect.PolyFilter([self.AtomS[r], FilterC, AtomO, AtomN], Bonded = True)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Torsion_ONCS = dict( (r, sim.potential.TorsionSpline(self.Sys, Filter = Filter_ONCS[r], Label = 'TorsionONCS_%s' % r, NKnot = self.NKnot)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Torsion_SCON = dict( (r, sim.potential.TorsionSpline(self.Sys, Filter = Filter_SCON[r], Label = 'TorsionSCON_%s' % r, NKnot = self.NKnot)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        # populate
        ff = Bond_CS.values() + \
            Angle_NCS.values() + Angle_SCO.values() + \
            Torsion_ONCS.values() + Torsion_SCON.values()
        return ff
    
    
    def BB_S_Bonded_1(self):
        '''1-alphabet bonded backbone-sidechain potentials
        sidechain reference can be residue name or number '''
        if Verbose: print 'Generating 1-alphabet bonded backbone-sidechain potentials. For splines, MinBondOrd = %d, %d knots, Cutoff = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        if self.SSRefType == 'number':
            FilterS = sim.atomselect.Filter([self.AtomS[i] for i in range(self.p.NRes) if not self.AtomS[i] is None])
        elif self.SSRefType == 'name':
            FilterS = sim.atomselect.Filter([self.AtomS[r] for r in self.p.ResTypes if not self.AtomS[r] is None])
        # create bonded potentials
        FilterC = sim.atomselect.Filter([AtomC, AtomC_GLY, AtomC_PRO])
        Filter_CS = sim.atomselect.PolyFilter([FilterC, FilterS], Bonded = True)
        Bond_CS = sim.potential.Bond(self.Sys, Filter = Filter_CS, Label = 'BondCS', Dist0 = 4.0, FConst = 1.0)
        # create angle potentials
        Filter_NCS = sim.atomselect.PolyFilter([AtomN, FilterC, FilterS], Bonded = True)
        Filter_SCO = sim.atomselect.PolyFilter([FilterS, FilterC, AtomO], Bonded = True)
        Angle_NCS = sim.potential.AngleSpline(self.Sys, Filter = Filter_NCS, Label = 'AngleNCS', NKnot = self.NKnot)
        Angle_SCO = sim.potential.AngleSpline(self.Sys, Filter = Filter_SCO, Label = 'AngleSCO', NKnot = self.NKnot)
        # create torsion potentials
        Filter_ONCS = sim.atomselect.PolyFilter([AtomO, AtomN, FilterC, FilterS], Bonded = True)
        Filter_SCON = sim.atomselect.PolyFilter([FilterS, FilterC, AtomO, AtomN], Bonded = True)
        Torsion_ONCS = sim.potential.TorsionSpline(self.Sys, Filter = Filter_ONCS, Label = 'TorsionONCS', NKnot = self.NKnot)
        Torsion_SCON = sim.potential.TorsionSpline(self.Sys, Filter = Filter_SCON, Label = 'TorsionSCON', NKnot = self.NKnot)
        # populate
        ff = [Bond_CS, \
              Angle_NCS, Angle_SCO, \
              Torsion_ONCS, Torsion_SCON]
        return ff
        
        
    def BB_S_Bonded_BondOnly(self):
        '''1-alphabet bonded backbone-sidechain potentials
        sidechain reference can be residue name or number '''
        if Verbose: print 'Generating backbone-sidechain bonds'
        if self.SSRefType == 'number':
            FilterS = sim.atomselect.Filter([self.AtomS[i] for i in range(self.p.NRes) if not self.AtomS[i] is None])
        elif self.SSRefType == 'name':
            FilterS = sim.atomselect.Filter([self.AtomS[r] for r in self.p.ResTypes if not self.AtomS[r] is None])
        # create bonded potentials
        FilterC = sim.atomselect.Filter([AtomC, AtomC_GLY, AtomC_PRO])
        Filter_CS = sim.atomselect.PolyFilter([FilterC, FilterS], Bonded = True)
        Bond_CS = sim.potential.Bond(self.Sys, Filter = Filter_CS, Label = 'BondCS', Dist0 = 4.0, FConst = 1.0)
        # populate
        ff = [Bond_CS]
        return ff


    def BB_S_0(self):
        '''21 alphabet nonbonded backbone-sidechain potentials
        assumes sidechain reference by residue name'''
        if Verbose: print 'Generating 21-alphabet backbone-sidechain nonbonded potentials. For splines, MinBondOrd = %d, %d knots, Cutoff = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        FilterC = sim.atomselect.Filter([AtomC, AtomC_GLY, AtomC_PRO])
        Filter_NS_p = dict( (r, sim.atomselect.PolyFilter([AtomN, self.AtomS[r]], MinBondOrd = self.MinBondOrd)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Filter_CS_p = dict( (r, sim.atomselect.PolyFilter([FilterC, self.AtomS[r]], MinBondOrd = self.MinBondOrd)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Filter_OS_p = dict( (r, sim.atomselect.PolyFilter([AtomO, self.AtomS[r]], MinBondOrd = self.MinBondOrd)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Pair_NS = dict( (r, sim.potential.PairSpline(self.Sys, Filter = Filter_NS_p[r], Label = 'NonBondNS_%s' % r, NKnot = self.NKnot, Cut = self.SPCut)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Pair_CS = dict( (r, sim.potential.PairSpline(self.Sys, Filter = Filter_CS_p[r], Label = 'NonBondCS_%s' % r, NKnot = self.NKnot, Cut = self.SPCut)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        Pair_OS = dict( (r, sim.potential.PairSpline(self.Sys, Filter = Filter_OS_p[r], Label = 'NonBondOS_%s' % r, NKnot = self.NKnot, Cut = self.SPCut)) for r in self.p.ResTypes if not self.AtomS[r] is None)
        # populate
        ff = Pair_NS.values() + Pair_CS.values() + Pair_OS.values()
        return ff
    
    
    def BB_S_1(self):
        '''1-alphabet nonbonded backbone-sidechain potentials
        sidechain reference can be residue name or number'''
        if Verbose: print 'Generating 1-alphabet backbone-sidechain nonbonded potentials. For splines, MinBondOrd = %d, %d knots, Cutoff = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        if self.SSRefType == 'number':
            FilterS = sim.atomselect.Filter([self.AtomS[i] for i in range(self.p.NRes) if not self.AtomS[i] is None])
        elif self.SSRefType == 'name':
            FilterS = sim.atomselect.Filter([self.AtomS[r] for r in self.p.ResTypes if not self.AtomS[r] is None])
        FilterC = sim.atomselect.Filter([AtomC, AtomC_GLY, AtomC_PRO])
        Filter_NS_p = sim.atomselect.PolyFilter([AtomN, FilterS], MinBondOrd = self.MinBondOrd)
        Filter_CS_p = sim.atomselect.PolyFilter([FilterC, FilterS], MinBondOrd = self.MinBondOrd)
        Filter_OS_p = sim.atomselect.PolyFilter([AtomO, FilterS], MinBondOrd = self.MinBondOrd)
        Pair_NS = sim.potential.PairSpline(self.Sys, Filter = Filter_NS_p, Label = 'NonBondNS', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_CS = sim.potential.PairSpline(self.Sys, Filter = Filter_CS_p, Label = 'NonBondCS', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_OS = sim.potential.PairSpline(self.Sys, Filter = Filter_OS_p, Label = 'NonBondOS', NKnot = self.NKnot, Cut = self.SPCut)
        # populate
        ff = [Pair_NS, Pair_CS, Pair_OS]
        return ff

    
    def BB_S_2(self):
        '''1-alphabet constant low repulsive spline potential
        sidechain reference can be residue name or number'''
        if Verbose: print 'Generating a single (repulsive) backbone-sidechain nonbonded potential with MinBondOrd = %d, %d knots, Cutoff = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        if self.SSRefType == 'number':
            FilterS = sim.atomselect.Filter([self.AtomS[i] for i in range(self.p.NRes) if not self.AtomS[i] is None])
        elif self.SSRefType == 'name':
            FilterS = sim.atomselect.Filter([self.AtomS[r] for r in self.p.ResTypes if not self.AtomS[r] is None])
        FilterC = sim.atomselect.Filter([AtomC, AtomC_GLY, AtomC_PRO])
        FilterBB = sim.atomselect.Filter([AtomN, FilterC, AtomO])
        Filter_BBS = sim.atomselect.PolyFilter([FilterBB, FilterS], MinBondOrd = self.MinBondOrd)
        Pair_BBS = sim.potential.PairSpline(self.Sys, Filter = Filter_BBS, Label = 'NonBondBBS', NKnot = self.NKnot, Cut = self.SPCut)
        # populate
        ff = [Pair_BBS]
        return ff
