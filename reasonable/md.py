
#!/usr/bin/env python

''' provides routines for running MD / REMD simulations on the model and multiplexing trajectories based on temp
    using LAMMPS
    '''

import os, numpy as np
import sim, protein
import pickleTraj
from const import *

Verbose = True

# LAMMPS settings
sim.srel.optimizetrajlammps.useREMD = True
sim.export.lammps.InnerCutoff = 0.02

class REMD(object):
    ''' runs REMD and multiplexes trajectories according to temperature using LAMMPS
        input MD iterations and timestep are all measured in femtoseconds'''
    def __init__(self, p, Sys, Prefix = None, InitPdb = None, Temps = None, TempFile = None, OutDir = os.getcwd() , **kwargs):
        self.p = p
        self.Sys = Sys
        self.TempSet = self.Sys.TempSet
        self.OutDir = OutDir
        if Prefix is None: Prefix = self.p.Prefix
        self.Prefix = os.path.abspath(os.path.join(self.OutDir, Prefix))
        self.InitPdb = None
        # temp schedule
        self.Temps = Temps
        self.TempFile = TempFile
        if self.TempFile is None: self.TempFile = os.path.join(self.OutDir, 'temps.txt')
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

    def GenRandInitPos(self):
        '''generate random initial position based on /share/apps/scripts/template.pdb
            to run from a different seed structure just place a cg pdb called init.pdb with that structure'''
        tmpPdb = os.path.join(os.getcwd(), 'tmp.pdb')
        if self.InitPdb is None: self.InitPdb = 'init.pdb'
        if not os.path.isfile(self.InitPdb):
            if Verbose: print 'Generating fully extended initial structure...'
            # create an ALL-ATOM protein class object for the given sequence
            pobj = protein.ProteinClass(Seq = self.p.Seq)
            pobj.Dehydrogen()
            pobj = pobj.Decap()
            pobj.WritePdb(tmpPdb)
            # coarse grain all-atom proteinclass object
            os.system('python %s %s %s' % (MAPSCRIPT, tmpPdb, self.InitPdb.split('.')[0]))
            # remove all-atom pdb
            if os.path.isfile(tmpPdb): os.remove(tmpPdb)
            del pobj
        print 'Using init conf as generated in : %s' % self.InitPdb
        pobj = protein.ProteinClass(self.InitPdb)
        initpos = pobj.Pos
        return initpos

    def runREMD(self):
        ''' run REMD using LAMMPS'''
        # random initial structure
        sim.export.lammps.LammpsExec = LAMMPSEXEC
        self.Sys.Arrays.Pos = self.GenRandInitPos()
        sim.system.init.velocities.Canonical(self.Sys, Temp = self.Sys.TempSet)
        # record temp schedule
        if not os.path.isfile(self.TempFile):
            np.savetxt(self.TempFile, self.Temps)
        else:
            self.Temps = np.loadtxt(self.TempFile)
        # feed in all the settings to the Lammps export
        sim.export.lammps_REMD.NStepsSwap = self.NStepsSwap
        sim.export.lammps_REMD.TEMPS = self.Temps
        LammpsFiles, TrajFile = sim.export.lammps_REMD.MakeLammpsReplicaMD(self.Sys, Prefix = self.Prefix, TrajFile = '.lammpstrj',
                                                                           NStepsMin = self.NStepsMin, NStepsEquil = self.NStepsEquil,
                                                                           NStepsProd = self.NStepsProd, WriteFreq = self.StepFreq)
        # run Lammps REMD
        InFile, DataFile, TableFile, DihedralFile = LammpsFiles
        LogFile, ScreenFile, returncode = sim.export.lammps_REMD.RunLammpsReplica(InFile, Prefix = self.Prefix, Verbose = Verbose)
        return TrajFile, LogFile

    def reorderTraj(self, ReorderTemps = None):
        '''reordering trajectories by temp'''
        Temps = np.loadtxt(self.TempFile)
        if ReorderTemps is None: ReorderTemps = Temps
        TrajFile = self.Prefix + '.lammpstrj'
        LogFile = self.Prefix + 'lammps.log'
        for T in ReorderTemps:
            TempInd = [list(Temps).index(t) for t in list(Temps) if '%3.2f' % t  == '%3.2f' % T][0]
            MultiTrajFn = FMT['TRAJ'] % (self.Prefix, T)
            MultiEneFn = FMT['ENE'] % (self.Prefix, T)
            if os.path.isfile(MultiTrajFn) and os.path.isfile(MultiEneFn): continue
            if Verbose: print 'Collecting replica frames at temperature ', T
            RepIndsMaster = [np.where(x[1:] == TempInd)[0][0] for x in np.loadtxt(LogFile, skiprows = 3)]
            RepInds = RepIndsMaster[int(self.NStepsEquil / self.NStepsSwap) : ]
            RepInds = RepInds[:-1]
            this_Traj = {}; this_Ene = {}
            TrajList = [] ; EneList = []
        
            if self.StepFreq <= self.NStepsSwap: # assume mod(NStepsSwap, StepFreq) = 0
                for ii, i in enumerate(RepInds):
                    if not i in this_Traj.keys():
                        thisTrajFn = '%s.%d.gz' % (TrajFile, i)
                        thisLogFn = '%s.%d' % (LogFile, i)
                        this_Traj[i] = pickleTraj(thisTrajFn, LogFile = thisLogFn, LogFileToken = '#run production')
                        this_Ene[i] = this_Traj[i].ThermoDict['PEnergy']
                    this_Traj_slice = this_Traj[i][int(self.NStepsEquil / self.StepFreq) : ]
                    this_Ene_slice = this_Ene[i][int(self.NStepsEquil / self.StepFreq) : ]
                    start = ii * self.NStepsSwap / self.StepFreq ; stop = (ii+1) * self.NStepsSwap / self.StepFreq
                    TrajList.append(this_Traj_slice[start:stop])
                    EneList.extend(this_Ene_slice[start:stop])

            else:
                NSkip = self.StepFreq / self.NStepsSwap # assume mod(StepFreq, NStepsSwap) = 0
                for ii, i in enumerate(RepInds[0::NSkip]):
                    if not i in this_Traj.keys():
                        thisTrajFn = '%s.%d.gz' % (TrajFile, i)
                        thisLogFn = '%s.%d' % (LogFile, i)
                        this_Traj[i] = pickleTraj(thisTrajFn, LogFile = thisLogFn, LogFileToken = '#run production')
                        this_Ene[i] = this_Traj[i].ThermoDict['PEnergy']
                    this_Traj_slice = this_Traj[i][int(self.NStepsEquil / self.StepFreq) : ]
                    this_Ene_slice = this_Ene[i][int(self.NStepsEquil / self.StepFreq) : ]
                    TrajList.append(this_Traj_slice[ii:ii+1])
                    EneList.extend(this_Ene_slice[ii:ii+1])

        MultiTraj = sim.traj.Multi(TrajList)
        sim.traj.base.Convert(InTraj = MultiTraj, OutTrajClass = sim.traj.LammpsWrite, FileName = MultiTrajFn)
        np.savetxt(MultiEneFn, EneList, fmt = '%11.4e')

