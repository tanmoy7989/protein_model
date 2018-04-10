
#!/usr/bin/env python

''' provides routines for relative entropy optimization in single and extended ensemble cases'''

import os, numpy as np
import sim, protein
import pickleTraj
from const import *

# settings
sim.srel.optimizetrajlammps.useREMD = True
sim.export.lammps.InnerCutoff = 0.02
sim.srel.optimizetraj.PlotFmt = 'svg'

class Srel(object):
    ''' runs relative entropy optimization for polymers, Go models in single and extended ensemble situations'''
    def __init__(self, p, Sys, AATraj, Temps, cfg = None, Prefix = None, OutDir = os.getcwd(), **kwargs):
        # read in protein, Sys and traj objects
        self.p = p
        self.Sys = Sys
        self.AATraj = AATraj
        # config object
        self.cfg = cfg if not cfg is None else self.p.cfg
        # prefix
        if Prefix is None: Prefix = p.Prefix
        self.Prefix = os.path.abspath(os.path.join(OutDir, Prefix))
        # read in the trajectory
        self.Trj = pickleTraj(self.AATraj)
        # initial conditions (choose to seed with native pos or not)
        # (native pos = RefPosInd-th frame of supplied AA Traj)
        self.useSeed = kwargs.get('useSeed', False)
        RefPosInd = kwargs.get('RefPosInd', None)
        if not self.useSeed: self.RefPosInd = None
        else: self.RefPosInd = RefPosInd
        sim.system.init.velocities.Canonical(self.Sys, Temp = self.Sys.TempSet)
        # 1-1 map (since supplied traj will be already mapped)
        Map = sim.atommap.PosMap()
        for (i, a) in enumerate(self.Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]
        self.Map = Map
        # temp schedule
        self.Temps = Temps
        # MD iterations
        self.NStepsMin = kwargs.get('NStepsMin', DEFAULTS['NStepsMin'])
        self.NStepsEquil = kwargs.get('NStepsEquil', DEFAULTS['NStepsEquil'])
        self.NStepsProd = kwargs.get('NStepsEquil', DEFAULTS['NStepsProd'])
        self.NStepsSwap = kwargs.get('NStepsProd', DEFAULTS['NStepsSwap'])
        self.StepFreq = kwargs.get('StepFreq', DEFAULTS['StepFreq'])
        self.TimeStep = kwargs.get('TimeStep', DEFAULTS['TimeStep'])
        # rescale the iterations according to the timestep
        for m in self.Sys.Int.Methods: m.TimeStep *= self.TimeStep
        self.NStepsMin = int(self.NStepsMin / self.TimeStep)
        self.NStepsEquil = int(self.NStepsEquil/ self.TimeStep)
        self.NStepsProd  = int(self.NStepsProd / self.TimeStep)
        self.NStepsSwap = int(self.NStepsSwap/ self.TimeStep)
        self.StepFreq = int(self.StepFreq  / self.TimeStep)
        # max srel iterations
        MaxIter = kwargs.get('MaxIter', dict(Bond = None, Angle = None, Torsion = None, NonBond = None, all = None) )
        self.MaxIter = MaxIter
        # Hessian
        self.useBondHessian = kwargs.get('useBondHessian', True)
        # permanently frozen potentials (constant potentials)
        self.PermaFrostList = kwargs.get('PermaFrostList', [])
        self.isPermaFrost = lambda P : self.PermaFrostList.__contains__(P.Name)
    
    
    def runPolymerSrel(self, OptStages = []):
        ''' optimize a NCOS backbone forcefield from a polypeptide traj'''
        # opt-stages
        if not OptStages: OptStages = ['Bond', 'Angle', 'Torsion', 'NonBond', 'all']
        # treat inner core (don't do this for frozen potentials)
        print '\nTreating inner core for nonbonded spline potentials...'
        for P in self.Sys.ForceField:
            # turn off constraints for PermaFrost potentials
            if self.isPermaFrost(P) and P.IsSpline: P.ConstrainSpline = False
            # apply constraints to all other potentials
            if P.Name.startswith("NonBond") and P.IsSpline and not self.isPermaFrost(P):
                P.EneInner = "20kT"
                P.EneSlopeInner = None
        # optimizer initialization (constraining all splines)
        print 'Starting Srel optimization for polypeptide...'
        Opt = sim.srel.OptimizeTrajClass(self.Sys, self.Map, Traj = self.Trj, SaveLoadArgData = True, FilePrefix = self.Prefix, Verbose = True, RefPosInd = self.RefPosInd)
        Opt = sim.srel.UseLammps(Opt)
        Opt.TempFileDir = os.getcwd()
        Opt.MinReweightFrames = None # need this to work with smaller mod traj
        # feed in all the settings to the Lammps export
        sim.export.lammps.LammpsExec = LAMMPSEXEC
        sim.export.lammps_REMD.NStepsSwap = self.NStepsSwap
        sim.export.lammps_REMD.TEMPS = self.Temps
        # stagewise optimization
        # bonds
        if OptStages.__contains__('Bond'):
            Opt.Iter = 0
            print 'Optimizing Bonds...'
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith("Bond") and not self.isPermaFrost(P):
                    P.Estimate() # initial guess based on Gaussian estimation
                    P.UnfreezeParam()
            if not self.useBondHessian: print 'Turning of Hessian Descent while optimizing bonds'
            Opt.UseHessian = self.useBondHessian # switch off hessian if bond distributions have multiple peaks
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['Bond'])
        # switch hessian back on regardless of if it was on or off during bond optimization
        Opt.UseHessian = True
        # angles
        if OptStages.__contains__('Angle'):
            Opt.Iter = 0
            print 'Optimizing Angles...'
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith("Angle") and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['Angle'])
        # torsions
        if OptStages.__contains__('Torsion'):
            Opt.Iter = 0
            print 'Optimizing Torsions...'
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith("Torsion") and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['Torsion'])
        # nonbond
        if OptStages.__contains__('NonBond'):
            Opt.Iter = 0
            print 'Optimizing NonBonds...'
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith("NonBond") and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['NonBond'])
        # simultaneous optimization
        if OptStages.__contains__('all'):
            Opt.Iter = 0
            print 'Optimizing all...'
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if not P.Name.startswith("Bond") and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['all'])

        del Opt


    def runGoSrel(self, OptNonNative = False):
        '''optimize Go potentials with a fixed backbone and a traj'''
        # turn off spline constraints for all backbone potentials
        OptSet = ['NonBondNative'] if not OptNonNative else ['NonBondNative', 'NonBondNonNative']
        print '\nTurning off constraints for backbone potentials...'
        for P in self.Sys.ForceField:
            if not P.IsSpline: continue
            if not OptSet.__contains__(P.Name):
                P.ConstrainSpline = False
            else:
                print 'Constraining potential: %s' % P.Name
        # treat inner core (don't do this for frozen potentials)
        print 'Treating inner core for spline Go potentials...'
        for P in self.Sys.ForceField:
            if ["NonBondNative", "NonBondNonNative"].__contains__(P.Name) and P.IsSpline and not self.isPermaFrost(P):
                P.EneInner = "20kT"
                P.EneSlopeInner = None
        # set up optimizer object
        print '\nStarting Srel optimization for Go model...'
        Opt = sim.srel.OptimizeTrajClass(self.Sys, self.Map, Traj = self.Trj, SaveLoadArgData = True, FilePrefix = self.Prefix, Verbose = True, RefPosInd = self.RefPosInd)
        Opt = sim.srel.UseLammps(Opt)
        Opt.TempFileDir = os.getcwd()
        Opt.MinReweightFrames = None # need this to work with smaller mod traj
        # feed in all the settings to the Lammps export
        sim.export.lammps.LammpsExec = LAMMPSEXEC
        sim.export.lammps_REMD.NStepsSwap = self.NStepsSwap
        sim.export.lammps_REMD.TEMPS = self.Temps
        # optimize
        Opt.Iter = 0
        if not OptNonNative:
            print 'Optimizing only native Go interactions...'
        else:
            print 'Optimizing ALL Go interactions...'
        for P in self.Sys.ForceField: P.FreezeParam()
        for P in self.Sys.ForceField:
            if OptSet.__contains__(P.Name) and not self.isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq)


    def runExtendedEnsemblePolymerSrel(self, OptStages = []):
        pass
