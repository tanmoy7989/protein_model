
#!/usr/bin/env python

''' provides routines for running MD / REMD simulations on the model and multiplexing trajectories based on temp
    using LAMMPS
    '''

import os, numpy as np
import sim, protein
import pickleTraj

import topo, mapNCOS
from const import *

Verbose = True

# LAMMPS settings
sim.srel.optimizetrajlammps.useREMD = True
sim.export.lammps.InnerCutoff = 0.02
sim.export.lammps.MaxPairEnekBT = MaxLammpsPairEnekBT

class REMD(object):
    ''' runs REMD and multiplexes trajectories according to temperature using LAMMPS
        input MD iterations and timestep are all measured in femtoseconds'''
    def __init__(self, p, Sys, cfg = None, Prefix = None, InitPdb = None, Temps = None, TempFile = None, OutDir = os.getcwd() , **kwargs):
        # read in protein and Sys objects
        self.p = p
        self.Sys = Sys
        self.TempSet = self.Sys.TempSet
        # config object
        self.cfg = cfg if not cfg is None else self.p.cfg
        # Prefixes, Output locations etc 
        self.OutDir = OutDir
        if Prefix is None: Prefix = self.p.Prefix
        self.RunPrefix = Prefix
        self.Prefix = os.path.abspath(os.path.join(self.OutDir, Prefix))
        self.InitPdb = None
        # temp schedule
        self.Temps = Temps
        self.TempFile = TempFile
        if self.TempFile is None: self.TempFile = os.path.join(self.OutDir, 'temps.txt')
        # MD iterations
        self.NStepsMin = kwargs.get('NStepsMin', DEFAULTS['NStepsMin'])
        self.NStepsEquil = kwargs.get('NStepsEquil', DEFAULTS['NStepsEquil'])
        self.NStepsProd = kwargs.get('NStepsProd', DEFAULTS['NStepsProd'])
        self.NStepsSwap = kwargs.get('NStepsSwap', DEFAULTS['NStepsSwap'])
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
            if Verbose: print 'Generating fully extended initial AA structure...'
            # create an ALL-ATOM protein class object for the given sequence
            pobj = protein.ProteinClass(Seq = self.p.Seq)
            pobj.WritePdb(tmpPdb)
            # coarse grain all-atom proteinclass object
            mapNCOS.Map(InPdb = tmpPdb, CGPrefix = self.InitPdb.split('.pdb')[0], hasPseudoGLY = self.p.hasPseudoGLY)
            # map to polymer
            if self.cfg.Map2Polymer:
                mapNCOS.Map2Polymer(Pdb = self.InitPdb, PolyName = self.cfg.PolyName, AAPdb = tmpPdb, hasPseudoGLY = self.p.hasPseudoGLY, MappedPrefix = self.InitPdb.split('.pdb')[0])
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
        LammpsFiles, TrajFile = sim.export.lammps_REMD.MakeLammpsReplicaMD(self.Sys, Prefix = self.RunPrefix, TrajFile = '.lammpstrj',
                                                                           NStepsMin = self.NStepsMin, NStepsEquil = self.NStepsEquil,
                                                                           NStepsProd = self.NStepsProd, WriteFreq = self.StepFreq)
        # run Lammps REMD
        InFile, DataFile, TableFile, DihedralFile = LammpsFiles
        LogFile, ScreenFile, returncode = sim.export.lammps_REMD.RunLammpsReplica(InFile, Prefix = self.RunPrefix, Verbose = Verbose)
        return TrajFile, LogFile

    def reorderTraj(self, ReorderTemps = None):
        '''reordering trajectories by temp'''
        Temps = np.loadtxt(self.TempFile)
        if ReorderTemps is None: ReorderTemps = Temps
        # extract filenames
        TrajFile = self.Prefix + '.lammpstrj'
        LogFile = self.Prefix + 'lammps.log'
        # start reordering for each requested temp
        for T in ReorderTemps:
            # get array index of this temp in temp schedule
            TempInd = [list(Temps).index(t) for t in list(Temps) if '%3.2f' % t  == '%3.2f' % T][0]
            # set filenames for reordered traj and energy files at this temp
            MultiTrajFn = FMT['TRAJ'] % (self.Prefix, T)
            MultiEneFn = FMT['ENE'] % (self.Prefix, T)
            if os.path.isfile(MultiTrajFn) and os.path.isfile(MultiEneFn): continue
            if Verbose: print 'Collecting replica frames at temperature ', T
            # collect all replica indices at this temp
            RepIndsMaster = [np.where(x[1:] == TempInd)[0][0] for x in np.loadtxt(LogFile, skiprows = 3)]
            # keep only replica inds for the production run part
            RepInds = RepIndsMaster[int(self.NStepsEquil / self.NStepsSwap) : ]
            RepInds = RepInds[:-1]
            this_Traj = {}; this_Ene = {}
            TrajList = [] ; EneList = []
        
            # case 1
            if self.StepFreq <= self.NStepsSwap: # assume mod(NStepsSwap, StepFreq) = 0
                # extract sim-style trj and ene objects for all replicas if not already done
                for ii, i in enumerate(RepInds):
                    if not i in this_Traj.keys():
                        thisTrajFn = '%s.%d.gz' % (TrajFile, i)
                        thisLogFn = '%s.%d' % (LogFile, i)
                        this_Traj[i] = pickleTraj(thisTrajFn, LogFile = thisLogFn, LogFileToken = '#run production')
                        this_Ene[i] = this_Traj[i].ThermoDict['PEnergy']
                    # keep only the parts for the production run
                    this_Traj_slice = this_Traj[i][int(self.NStepsEquil / self.StepFreq) : ]
                    this_Ene_slice = this_Ene[i][int(self.NStepsEquil / self.StepFreq) : ]
                    # find start and stop indices for this block 
                    start = ii * self.NStepsSwap / self.StepFreq ; stop = (ii+1) * self.NStepsSwap / self.StepFreq
                    TrajList.append(this_Traj_slice[start:stop])
                    EneList.extend(this_Ene_slice[start:stop])

            # case 2
            else:
                NSkip = self.StepFreq / self.NStepsSwap # assume mod(StepFreq, NStepsSwap) = 0
                # etract sim-style trj and ene objects for all replicas if not already done
                for ii, i in enumerate(RepInds[0::NSkip]):
                    if not i in this_Traj.keys():
                        thisTrajFn = '%s.%d.gz' % (TrajFile, i)
                        thisLogFn = '%s.%d' % (LogFile, i)
                        this_Traj[i] = pickleTraj(thisTrajFn, LogFile = thisLogFn, LogFileToken = '#run production')
                        this_Ene[i] = this_Traj[i].ThermoDict['PEnergy']
                    # keep only the parts for the production run
                    this_Traj_slice = this_Traj[i][int(self.NStepsEquil / self.StepFreq) : ]
                    this_Ene_slice = this_Ene[i][int(self.NStepsEquil / self.StepFreq) : ]
                    # find start and stop indices for this block
                    TrajList.append(this_Traj_slice[ii:ii+1])
                    EneList.extend(this_Ene_slice[ii:ii+1])

        MultiTraj = sim.traj.Multi(TrajList)
        sim.traj.base.Convert(InTraj = MultiTraj, OutTrajClass = sim.traj.LammpsWrite, FileName = MultiTrajFn)
        np.savetxt(MultiEneFn, EneList, fmt = '%11.4e')

