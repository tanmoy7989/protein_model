
#!/usr/bin/env python

''' provides routines for relative entropy optimization in single and extended ensemble cases'''

import os, numpy as np, copy
import sim, protein
import pickleTraj
from const import *

# common optimizer settings
sim.srel.optimizetrajlammps.useREMD = True
sim.export.lammps.InnerCutoff = 0.02
sim.srel.optimizetraj.PlotFmt = 'svg'

class Srel(object):
    ''' runs relative entropy optimization for polymers, Go models in single and extended ensemble situations'''
    def __init__(self, p, Sys, AATraj, Temps, cfg = None, Prefix = None, OutDir = os.getcwd(), **kwargs):
        # read in protein, Sys and traj objects
        print '\nRELATIVE ENTROPY OPTIMIZATION'
        print '=============================\n'
        self.p = p
        self.AATraj = AATraj
        # is this a single traj Srel optimization or extended ensemble?
        if isinstance(self.AATraj, list):
            self.__init_MultiSrel(Sys, AATraj)
        else:
            self.__init_Srel(Sys, AATraj)
        # config object
        self.cfg = cfg if not cfg is None else self.p.cfg
        # prefix
        if Prefix is None: Prefix = p.Prefix
        self.Prefix = os.path.abspath(os.path.join(OutDir, Prefix))
        # initial conditions (choose to seed with native pos or not)
        # (native pos = RefPosInd-th frame of supplied AA Traj)
        self.useSeed = kwargs.get('useSeed', False)
        RefPosInd = kwargs.get('RefPosInd', None)
        if not self.useSeed: self.RefPosInd = None
        else: self.RefPosInd = RefPosInd
        for Sys in self.AllSys:
            sim.system.init.velocities.Canonical(Sys, Temp = Sys.TempSet)
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
        for Sys in self.AllSys:
            for m in Sys.Int.Methods: m.TimeStep *= self.TimeStep
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
        # feed in all the settings to the Lammps export
        sim.export.lammps.LammpsExec = LAMMPSEXEC
        sim.export.lammps_REMD.NStepsSwap = self.NStepsSwap
        sim.export.lammps_REMD.TEMPS = self.Temps
    
    def __init_Srel(self, Sys, AATraj):
        ''' set Sys, Trj and Map for a single traj'''
        self.NTraj = 1
        self.Sys = Sys
        self.AllSys = [self.Sys]
        self.AATraj = AATraj
        # pickle supplied Traj
        self.Trj = pickleTraj(self.AATraj)
        # 1-1 map (since supplied traj will be already mapped)
        Map = sim.atommap.PosMap()
        for (i, a) in enumerate(self.Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]
        self.Map = Map
        
    def __init_MultiSrel(self, Sys, AATraj):
        '''set Sys, Trj and Map for extended ensemble optimization'''
        self.NTraj = len(AATraj)
        self.Sys = Sys if isinstance(Sys, list) else [Sys] * self.NTraj
        self.AllSys = self.Sys
        # list of pickled Trj objects
        self.Trj = [pickleTraj(Traj) for Traj in self.AATraj]
        # list of 1-1 maps (since all supplied trajs will be already mapped)
        Map = [sim.atommap.PosMap() for i in range(self.NTraj)] # need to create separate PosMap instances
        for i in range(self.NTraj):
            for (j, a) in enumerate(self.Sys[i].Atom): Map[i] += [sim.atommap.AtomMap(Atoms1 = j, Atom2 = a)]
        self.Map = Map
            
    def isPermaFrost(self, P):
        return all( [ (P in Sys.ForceField and P.Name in self.PermaFrostList) for Sys in self.AllSys ] )
        
    def ReportStage(self, Stage):
        ''' generate report strings at every stage, displaying frozen
        and optimizable potentials'''
        FrozenPotentials = []
        UnfrozenPotentials = []
        for P in self.AllSys[0].ForceField:
            Test_others = (not Stage == 'all') and P.Name.startswith(Stage) and not self.isPermaFrost(P)
            Test_all = (Stage == 'all') and not P.Name.startswith('Bond') and not self.isPermaFrost(P)
            # add other tests for new stages if necessary
            if (Test_others or Test_all): UnfrozenPotentials.append(P.Name)
            else: FrozenPotentials.append(P.Name)
        FrozenPotentials.sort()
        UnfrozenPotentials.sort()
        if not FrozenPotentials: FrozenPotentials = ['None']
        if not UnfrozenPotentials: UnfrozenPotentials = ['None']
        print '\n\nOPTIMIZING STAGE: %s' % Stage.upper()
        print '============================='
        print 'Freezing potentials: %s' % (', '.join(FrozenPotentials))
        print 'Optimzing potentials: %s' % (', '.join(UnfrozenPotentials))
            
    def SetConstraints(self):
        '''makes sure that frozen potentials are not constrained'''
        s = 'Removing constraints for potentials: '
        RemovedConstraintPotentials = []
        for i, Sys in enumerate(self.AllSys):
            for P in Sys.ForceField:
                Test = P.IsSpline and P.Name in self.PermaFrostList
                if Test:
                    P.ConstrainSpline = False
                    if i == 0: RemovedConstraintPotentials.append(P.Name) # construct reportstring only once
        RemovedConstraintPotentials.sort()
        if not RemovedConstraintPotentials: RemovedConstraintPotentials = ['None']
        print s + ', '.join(RemovedConstraintPotentials)
    
    def ConstrainInnerCore(self):
        '''constrains the inner core of nonbonded spline potentials
        for which ConstrainSpline is True. Uses a MaxEne of 20 kT and
        a variable slope: this strategy worked well for small tests'''
        # make sure that self.SetConstraints() has been used prior to this
        s = '\nConstraining inner core with EneMax = 20 kT and variable slope for nonbonded potentials: '
        ConstrainedPotentials = []
        for i, Sys in enumerate(self.AllSys): 
            for P in Sys.ForceField:
                Test1 = P.IsSpline and P.ConstrainSpline
                Test2 = P.Name.startswith('NonBond')
                Test = Test1 and Test2
                if Test:
                    if i == 0: ConstrainedPotentials.append(P.Name) # construct reportstring only once
                    P.EneInner = "20kT"
                    P.EneSlopeInner = None
        ConstrainedPotentials.sort()
        if not ConstrainedPotentials: ConstrainedPotentials = ['None']
        print s + ', '.join(ConstrainedPotentials)
    
    def PrepOpt(self, Opt):
        ''' encapsulates common settings for an optimizer object
        (post initializing the object)'''
        # Note: common optimizer settings defined towards the beginning
        # of this module need to be defined prior to creating the optimizer
        # so cannot be defined here
        Opt = sim.srel.UseLammps(Opt)
        Opt.TempFileDir = os.getcwd()
        Opt.MinReweightFrames = None # need this to work with smaller mod traj
    
            
    def runPolymerSrel(self, OptStages = []):
        ''' optimize a NCOS backbone forcefield from a polypeptide traj'''
        # opt-stages
        if not OptStages: OptStages = ['Bond', 'Angle', 'Torsion', 'NonBond', 'all']
        # set constraints
        self.SetConstraints()
        # treat inner core
        self.ConstrainInnerCore()
        # optimizer initialization (constraining all splines)
        Opt = sim.srel.OptimizeTrajClass(self.Sys, self.Map, Traj = self.Trj, SaveLoadArgData = True, FilePrefix = self.Prefix, Verbose = True, RefPosInd = self.RefPosInd)
        # prepare the optimizer object
        self.PrepOpt(Opt) 
        print '\nStarting stagewise optimization for polypeptide...\n'
        # stagewise optimization
        # bonds
        if 'Bond' in OptStages:
            Opt.Iter = 0
            self.ReportStage('Bond')
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith("Bond") and not self.isPermaFrost(P):
                    P.Estimate() # initial guess based on Gaussian estimation
                    P.UnfreezeParam()
            if not self.useBondHessian: print ' Turning of Hessian Descent while optimizing bonds'
            Opt.UseHessian = self.useBondHessian # switch off hessian if bond distributions have multiple peaks
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['Bond'])
        # switch hessian back on regardless of if it was on or off during bond optimization
        Opt.UseHessian = True
        # angles
        if 'Angle' in OptStages:
            Opt.Iter = 0
            self.ReportStage('Angle')
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith("Angle") and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['Angle'])
        # torsions
        if 'Torsion' in OptStages:
            Opt.Iter = 0
            self.ReportStage('Torsion')
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith("Torsion") and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['Torsion'])
        # nonbond
        if 'NonBond' in OptStages:
            Opt.Iter = 0
            self.ReportStage('NonBond')
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith("NonBond") and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['NonBond'])
        # simultaneous optimization
        if 'all' in OptStages:
            Opt.Iter = 0
            self.ReportStage('all')
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if not P.Name.startswith("Bond") and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq, MaxIter = self.MaxIter['all'])

        del Opt


    def runGoSrel(self, Constrain = True, OptNonNative = False):
        '''optimize Go potentials with a fixed backbone'''
        # turn off spline constraints for all backbone potentials
        OptStages = ['NonBondNative'] 
        if OptNonNative: OptStages += ['NonBondNonNative', 'Go_all']
        # set constraints
        for P in self.Sys.ForceField:
            Test1 = 'NonBondNative' in P.Name
            Test2 = 'NonBondNonNative' in P.Name and OptNonNative
            Test = Test1 or Test2
            if not Test: self.PermaFrostList.append(P.Name)
        self.SetConstraints()
        # treat inner core
        if Constrain: self.ConstrainInnerCore()
        # set up optimizer object
        print '\nStarting Go model optimization...\n'
        Opt = sim.srel.OptimizeTrajClass(self.Sys, self.Map, Traj = self.Trj, SaveLoadArgData = True, FilePrefix = self.Prefix, Verbose = True, RefPosInd = self.RefPosInd)
        # prepare the optimizer object
        self.PrepOpt(Opt)
        # stagewise optimization (usually just the native interactions)
        # native interactions
        if 'NonBondNative' in OptStages:
            Opt.Iter = 0
            self.ReportStage('NonBondNative')
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith('NonBondNative') and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq)
        # non-native interactions
        if 'NonBondNonNative' in OptStages:
            Opt.Iter = 0
            self.ReportStage('NonBondNonNative')
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if P.Name.startswith('NonBondNonNative') and not self.isPermaFrost(P):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq)
        # all Go interactions simultaneously
        if 'Go_all' in OptStages:
            Opt.Iter = 0
            self.ReportStage('Go_all')
            for P in self.Sys.ForceField: P.FreezeParam()
            for P in self.Sys.ForceField:
                if (not self.isPermaFrost(P)) and (P.Name.startswith('NonBondNative') or P.Name.startswith('NonBondNonNative')):
                    P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq)
        
        del Opt


    def runMutliPolymerSrel(self, OptStages = []):
        pass
    
    
    def runMultiGoSrel(self, Constrain = True, OptNonNative = False):
        '''optimize Go potentials with a fixed backbone from an ensemble of peptides'''
        # turn off spline constraints for all backbone potentials
        OptStages = ['NonBondNative'] 
        if OptNonNative: OptStages += ['NonBondNonNative', 'Go_all']
        # set constraints
        for P in self.Sys[0].ForceField:
            Test1 = 'NonBondNative' in P.Name
            Test2 = 'NonBondNonNative' in P.Name and OptNonNative
            Test = Test1 or Test2
            if not Test: self.PermaFrostList.append(P.Name)
        self.SetConstraints()
        # treat inner core
        if Constrain: self.ConstrainInnerCore()
        # set up optimizer object
        print '\nStarting extended ensemble Go model optimization...\n'
        print 'Creating ensemble of optimizers...'
        OptList = []
        for i in range(self.NTraj):
            print '\nOptimizer %d' % i
            # initialize the optimizer object
            thisOpt = sim.srel.OptimizeTrajClass(self.Sys[i], self.Map[i], Traj = self.Trj[i], SaveLoadArgData = True, 
                                             FilePrefix = self.Prefix + '_%d' % i, Verbose = True, RefPosInd = self.RefPosInd)
            # prepare the optimizer object
            self.PrepOpt(thisOpt)
            OptList.append(thisOpt)
        # initialize the multi-optimizer object
        Opt = sim.srel.OptimizeMultiTrajClass(OptList, FilePrefix = self.Prefix)
        # stagewise optimization (usually just the native interactions)
        # native interactions
        if 'NonBondNative' in OptStages:
            Opt.Iter = 0
            self.ReportStage('NonBondNative')
            for Sys in self.Sys:
                for P in Sys.ForceField: P.FreezeParam()
            for Sys in self.Sys:
                for P in Sys.ForceField:
                    if P.Name.startswith('NonBondNative') and not self.isPermaFrost(P):
                        P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq)
        # non-native interactions
        if 'NonBondNonNative' in OptStages:
            Opt.Iter = 0
            self.ReportStage('NonBondNonNative')
            for Sys in self.Sys:
                for P in Sys.ForceField: P.FreezeParam()
            for Sys in self.Sys:
                for P in Sys.ForceField:
                    if P.Name.startswith('NonBondNonNative') and not self.isPermaFrost(P):
                        P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq)
        # all Go interactions simultaneously
        if 'Go_all' in OptStages:
            Opt.Iter = 0
            self.ReportStage('Go_all')
            for Sys in self.Sys:
                for P in Sys.ForceField: P.FreezeParam()
            for Sys in self.Sys:
                for P in Sys.ForceField:
                    if (not self.isPermaFrost(P)) and (P.Name.startswith('NonBondNative') or P.Name.startswith('NonBondNonNative')):
                        P.UnfreezeParam()
            Opt.RunConjugateGradient(StepsEquil = self.NStepsEquil, StepsProd = self.NStepsProd, StepsStride = self.StepFreq)
        


