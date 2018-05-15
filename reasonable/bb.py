#!/usr/bin/env python

import sim
from const import *

''' Builds intra-backbone bonded and nonbonded potentials.
    Functions are named numerically based on types in settings.
    Inputs are the NCOS cg protein and Sys objects'''

Verbose = True
class P_Backbone(object):
    def __init__(self, p, Sys, cfg = None):
        ''' initialized a class for backbone potentials
            cfg is a configuration object'''
        if cfg is None: cfg = p.cfg
        # extract relevant settings
        self.MinBondOrd = cfg.MinBondOrd
        self.NKnot = cfg.NKnot
        self.SPCut = cfg.SPCut
        self.hasSpecialBBGLYAngles = cfg.hasSpecialBBGLYAngles
        self.hasSpecialBBGLYTorsions = cfg.hasSpecialBBGLYTorsions
        self.hasSpecialBBPROAngles = cfg.hasSpecialBBPROAngles
        self.hasSpecialBBPROTorsions = cfg.hasSpecialBBPROTorsions
        # extract cg protein and Sys objects
        self.p = p
        self.Sys = Sys
    
    def BB_0(self):
        ''' bonded and nonbonded potentials between N,C,O atoms in the backbone'''
        # create bonded potentials
        if Verbose: print 'Generating intra-backbone bonded and nonbonded potentials. For splines, MinBondOrd = %d, %d knots, Cutoff = %2.2f A' % (self.MinBondOrd, self.NKnot, self.SPCut)
        # create a common "C" atom to be used in all relevant potentials
        # except special angle and torsion potentials
        FilterC = sim.atomselect.Filter([AtomC, AtomC_GLY, AtomC_PRO])
        Filter_NC = sim.atomselect.PolyFilter([AtomN, FilterC], Bonded = True)
        Filter_CO = sim.atomselect.PolyFilter([FilterC, AtomO], Bonded = True)
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
        Filter_CONC = sim.atomselect.PolyFilter([FilterC, AtomO, AtomN, AtomC], Bonded = True)
        Filter_ONCO = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC, AtomO], Bonded = True)
        Torsion_NCON = sim.potential.TorsionSpline(self.Sys, Filter = Filter_NCON, Label = 'TorsionNCON', NKnot = self.NKnot)
        Torsion_CONC = sim.potential.TorsionSpline(self.Sys, Filter = Filter_CONC, Label = 'TorsionCONC', NKnot = self.NKnot)
        Torsion_ONCO = sim.potential.TorsionSpline(self.Sys, Filter = Filter_ONCO, Label = 'TorsionONCO', NKnot = self.NKnot)
        # created nonbonded potentials
        Filter_NN_p = sim.atomselect.PolyFilter([AtomN, AtomN], MinBondOrd = self.MinBondOrd)
        Filter_CC_p = sim.atomselect.PolyFilter([FilterC, FilterC], MinBondOrd = self.MinBondOrd)
        Filter_OO_p = sim.atomselect.PolyFilter([AtomO, AtomO], MinBondOrd = self.MinBondOrd)
        Filter_NC_p = sim.atomselect.PolyFilter([AtomN, FilterC], MinBondOrd = self.MinBondOrd)
        Filter_NO_p = sim.atomselect.PolyFilter([AtomN, AtomO], MinBondOrd = self.MinBondOrd)
        Filter_CO_p = sim.atomselect.PolyFilter([FilterC, AtomO], MinBondOrd = self.MinBondOrd)
        Pair_NN = sim.potential.PairSpline(self.Sys, Filter = Filter_NN_p, Label = 'NonBondNN', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_CC = sim.potential.PairSpline(self.Sys, Filter = Filter_CC_p, Label = 'NonBondCC', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_OO = sim.potential.PairSpline(self.Sys, Filter = Filter_OO_p, Label = 'NonBondOO', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_NC = sim.potential.PairSpline(self.Sys, Filter = Filter_NC_p, Label = 'NonBondNC', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_NO = sim.potential.PairSpline(self.Sys, Filter = Filter_NO_p, Label = 'NonBondNO', NKnot = self.NKnot, Cut = self.SPCut)
        Pair_CO = sim.potential.PairSpline(self.Sys, Filter = Filter_CO_p, Label = 'NonBondCO', NKnot = self.NKnot, Cut = self.SPCut)
        # create special angle and torsion potentials for gly and pro
        if self.hasSpecialBBGLYAngles:
            if 'GLY' in self.p.Seq:
                if Verbose: print 'Generating special backbone angle potentials for GLY...'
                Filter_NCO_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO], Bonded = True)
                Filter_CON_GLY = sim.atomselect.PolyFilter([AtomC_GLY, AtomO, AtomN], Bonded = True)
                Filter_ONC_GLY = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC_GLY], Bonded = True)
                Angle_NCO_GLY = sim.potential.AngleSpline(self.Sys, Filter = Filter_NCO_GLY, Label = 'AngleNCO_GLY', NKnot = self.NKnot)
                Angle_CON_GLY = sim.potential.AngleSpline(self.Sys, Filter = Filter_CON_GLY, Label = 'AngleCON_GLY', NKnot = self.NKnot)
                Angle_ONC_GLY = sim.potential.AngleSpline(self.Sys, Filter = Filter_ONC_GLY, Label = 'AngleONC_GLY', NKnot = self.NKnot)
        if self.hasSpecialBBGLYTorsions:
            if 'GLY' in self.p.Seq:
                if Verbose: print 'Generating special backbone torsional potentials for GLY...'
                Filter_NCON_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
                Filter_ONCO_GLY = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC_GLY, AtomO], Bonded = True)
                Torsion_NCON_GLY = sim.potential.TorsionSpline(self.Sys, Filter = Filter_NCON_GLY, Label = 'TorsionNCON_GLY', NKnot = self.NKnot)
                Torsion_ONCO_GLY = sim.potential.TorsionSpline(self.Sys, Filter = Filter_ONCO_GLY, Label = 'TorsionONCO_GLY', NKnot = self.NKnot)
        if self.hasSpecialBBPROAngles:
            if 'PRO' in self.p.Seq:
                if Verbose: print 'Generating special backbone angle potentials for PRO...'
                Filter_NCO_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO], Bonded = True)
                Filter_CON_PRO = sim.atomselect.PolyFilter([AtomC_PRO, AtomO, AtomN], Bonded = True)
                Filter_ONC_PRO = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC_PRO], Bonded = True)
                Angle_NCO_PRO = sim.potential.AngleSpline(self.Sys, Filter = Filter_NCO_PRO, Label = 'AngleNCO_PRO', NKnot = self.NKnot)
                Angle_CON_PRO = sim.potential.AngleSpline(self.Sys, Filter = Filter_CON_PRO, Label = 'AngleCON_PRO', NKnot = self.NKnot)
                Angle_ONC_PRO = sim.potential.AngleSpline(self.Sys, Filter = Filter_ONC_PRO, Label = 'AngleONC_PRO', NKnot = self.NKnot)
        if self.hasSpecialBBPROTorsions:
            if 'PRO' in self.p.Seq:
                if Verbose: print 'Generating special backbone torsion potentials for PRO...'
                Filter_NCON_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
                Filter_ONCO_PRO = sim.atomselect.PolyFilter([Atomo, AtomN, AtomC_PRO, AtomO], Bonded = True)
                Torsion_NCON_PRO = sim.potential.TorsionSpline(self.Sys, Filter = Filter_NCON_PRO, Label = 'TorsionNCON_PRO', NKnot = self.NKnot)
                Torsion_ONCO_PRO = sim.potential.TorsionSpline(self.Sys, Filter = Filter_ONCO_PRO, Label = 'TorsionONCO_PRO', NKnot = self.NKnot)
        # populate
        ff = [Bond_NC, Bond_CO, Bond_ON, \
              Angle_NCO, Angle_CON, Angle_ONC,\
              Torsion_NCON, Torsion_CONC, Torsion_ONCO,
              Pair_NN, Pair_CC, Pair_OO, Pair_NC, Pair_NO, Pair_CO]
        if self.hasSpecialBBGLYAngles:
            if 'GLY' in self.p.Seq: ff.extend([Angle_NCO_GLY, Angle_CON_GLY, Angle_ONC_GLY])
        if self.hasSpecialBBGLYTorsions:
            if 'GLY' in self.p.Seq: ff.extend([Torsion_NCON_GLY, Torsion_ONCO_GLY])
        if self.hasSpecialBBPROAngles:
            if 'PRO' in self.p.Seq: ff.extend([Angle_NCO_PRO, Angle_CON_PRO, Angle_ONC_PRO])
        if self.hasSpecialBBPROTorsions:
            if 'PRO' in self.p.Seq: ff.extend([Torsion_NCON_PRO, Torsion_ONCO_PRO])
        return ff
     
     
    def BB_BondOnly(self):
        ''' only bond between N,C,O atoms in the backbone'''
        # create bonded potentials
        if Verbose: print 'Generating intra-backbone bonds'
        FilterC = sim.atomselect.Filter([AtomC, AtomC_GLY, AtomC_PRO])
        Filter_NC = sim.atomselect.PolyFilter([AtomN, FilterC], Bonded = True)
        Filter_CO = sim.atomselect.PolyFilter([FilterC, AtomO], Bonded = True)
        Filter_ON = sim.atomselect.PolyFilter([AtomO, AtomN], Bonded = True)
        Bond_NC = sim.potential.Bond(self.Sys, Filter = Filter_NC, Label = 'BondNC', Dist0 = 4.0, FConst = 1.0)
        Bond_CO = sim.potential.Bond(self.Sys, Filter = Filter_CO, Label = 'BondCO', Dist0 = 4.0, FConst = 1.0)
        Bond_ON = sim.potential.Bond(self.Sys, Filter = Filter_ON, Label = 'BondON', Dist0 = 4.0, FConst = 1.0)
        # populate
        ff = [Bond_NC, Bond_CO, Bond_ON]
        return ff
