#!/usr/bin/env python

import sim
from const import *

''' Builds intra-backbone bonded and nonbonded potentials.
    Functions are named numerically based on types in settings.
    Inputs are the NCOS cg protein and Sys objects'''

Verbose = True
class P_Backbone(object):
    def __init__(self, p, Sys, cfg):
        ''' initialized a class for backbone potentials
            cfg is a configuration object'''
        # extract relevant settings
        self.MinBondOrd = cfg.MinBondOrd
        self.NKnot = cfg.NKnot
        self.SPCut = cfg.SPCut
        self.hasSpecialBBTorsions = cfg.hasSpecialBBTorsions
        # extract cg protein and Sys objects
        self.p = p
        self.Sys = Sys
    
    def BB_0(self):
        ''' bonded and nonbonded potentials between N,C,O atoms in the backbone'''
        # create bonded potentials
        if Verbose: print 'Generating intra-backbone bonded and nonbonded potentials. For splines, MinBondOrd = %d, %d knots, Cutoff = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        Filter_NC = sim.atomselect.PolyFilter([AtomN, AtomC], Bonded = True)
        Filter_CO = sim.atomselect.PolyFilter([AtomC, AtomO], Bonded = True)
        Filter_ON = sim.atomselect.PolyFilter([AtomO, AtomN], Bonded = True)
        Bond_NC = sim.potential.Bond(self.Sys, Filter = Filter_NC, Label = 'BondNC', Dist0 = 4.0, FConst = 1.0)
        Bond_CO = sim.potential.Bond(self.Sys, Filter = Filter_CO, Label = 'BondCO', Dist0 = 4.0, FConst = 1.0)
        Bond_ON = sim.potential.Bond(self.Sys, Filter = Filter_ON, Label = 'BondON', Dist0 = 4.0, FConst = 1.0)
        # create angle potentials
        Filter_NCO = sim.atomselect.PolyFilter([AtomN, AtomC, AtomO], Bonded = True)
        Filter_CON = sim.atomselect.PolyFilter([AtomC, AtomO, AtomN], Bonded = True)
        Filter_ONC = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC], Bonded = True)
        Angle_NCO = sim.potential.AngleSpline(self.Sys, Filter = Filter_NCO, Label = 'AngleNCO', NKnot = self.NKnot)
        Angle_CON = sim.potential.AngleSpline(self.Sys, Filter = Filter_CON, Label = 'AngleCON', NKnot = self.NKnot)
        Angle_ONC = sim.potential.AngleSpline(self.Sys, Filter = Filter_ONC, Label = 'AngleONC', NKnot = self.NKnot)
        # create torsion potentials
        Filter_NCON = sim.atomselect.PolyFilter([AtomN, AtomC, AtomO, AtomN], Bonded = True)
        Filter_CONC = sim.atomselect.PolyFilter([AtomC, AtomO, AtomN, AtomC], Bonded = True)
        Filter_ONCO = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC, AtomO], Bonded = True)
        Torsion_NCON = sim.potential.TorsionSpline(self.Sys, Filter = Filter_NCON, Label = 'TorsionNCON', NKnot = self.NKnot)
        Torsion_CONC = sim.potential.TorsionSpline(self.Sys, Filter = Filter_CONC, Label = 'TorsionCONC', NKnot = self.NKnot)
        Torsion_ONCO = sim.potential.TorsionSpline(self.Sys, Filter = Filter_ONCO, Label = 'TorsionONCO', NKnot = self.NKnot)
        # created nonbonded potentials
        Filter_NN_p = sim.atomselect.PolyFilter([AtomN, AtomN], MinBondOrd = self.MinBondOrd)
        Filter_CC_p = sim.atomselect.PolyFilter([AtomC, AtomC], MinBondOrd = self.MinBondOrd)
        Filter_OO_p = sim.atomselect.PolyFilter([AtomO, AtomO], MinBondOrd = self.MinBondOrd)
        Filter_NC_p = sim.atomselect.PolyFilter([AtomN, AtomC], MinBondOrd = self.MinBondOrd)
        Filter_NO_p = sim.atomselect.PolyFilter([AtomN, AtomO], MinBondOrd = self.MinBondOrd)
        Filter_CO_p = sim.atomselect.PolyFilter([AtomC, AtomO], MinBondOrd = self.MinBondOrd)
        Pair_NN = sim.potential.PairSpline(self.Sys, Filter = Filter_NN_p, Label = 'NonBondNN', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_CC = sim.potential.PairSpline(self.Sys, Filter = Filter_CC_p, Label = 'NonBondCC', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_OO = sim.potential.PairSpline(self.Sys, Filter = Filter_OO_p, Label = 'NonBondOO', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_NC = sim.potential.PairSpline(self.Sys, Filter = Filter_NC_p, Label = 'NonBondNC', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_NO = sim.potential.PairSpline(self.Sys, Filter = Filter_NO_p, Label = 'NonBondNO', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_CO = sim.potential.PairSpline(self.Sys, Filter = Filter_CO_p, Label = 'NonBondCO', NKnot = self.NKnot, Cut = self.SPCut)
        # create special torsion potentials for gly and pro
        if self.hasSpecialBBTorsions:
            if Verbose: print 'Generating special torsional backbone potentials for GLY and PRO...'
            if self.p.Seq.__contains__('GLY'):
                Filter_NCON_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
                Filter_CONC_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
                Filter_ONCO_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
                Torsion_NCON_GLY = sim.potential.TorsionSpline(self.Sys, Filter = Filter_NCON_GLY, Label = 'TorsionNCON_GLY', NKnot = self.NKnot)
                Torsion_CONC_GLY = sim.potential.TorsionSpline(self.Sys, Filter = Filter_CONC_GLY, Label = 'TorsionCONC_GLY', NKnot = self.NKnot)
                Torsion_ONCO_GLY = sim.potential.TorsionSpline(self.Sys, Filter = Filter_ONCO_GLY, Label = 'TorsionONCO_GLY', NKnot = self.NKnot)
            if self.p.Seq.__contains__('PRO'):
                Filter_NCON_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
                Filter_CONC_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
                Filter_ONCO_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
                Torsion_NCON_PRO = sim.potential.TorsionSpline(self.Sys, Filter = Filter_NCON_PRO, Label = 'TorsionNCON_PRO', NKnot = self.NKnot)
                Torsion_CONC_PRO = sim.potential.TorsionSpline(self.Sys, Filter = Filter_CONC_PRO, Label = 'TorsionCONC_PRO', NKnot = self.NKnot)
                Torsion_ONCO_PRO = sim.potential.TorsionSpline(self.Sys, Filter = Filter_ONCO_PRO, Label = 'TorsionONCO_PRO', NKnot = self.NKnot)
        # populate
        ff = [Bond_NC, Bond_CO, Bond_ON, \
              Angle_NCO, Angle_CON, Angle_ONC,\
              Torsion_NCON, Torsion_CONC, Torsion_ONCO,
              Pair_NN, Pair_CC, Pair_OO, Pair_NC, Pair_NO, Pair_CO]
        if self.hasSpecialBBTorsions:
            if p.Seq.__contains__('GLY'): ff.extend([Torsion_NCON_GLY, Torsion_CONC_GLY, Torsion_ONCO_GLY])
            if p.Seq.__contains__('PRO'): ff.extend([Torsion_NCON_PRO, Torsion_CONC_PRO, Torsion_ONCO_PRO])
        return ff
