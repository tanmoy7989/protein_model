
#!/usr/bin/env python

''' provides routines for running MD / REMD simulations on the model and multiplexing trajectories based on temp
    using LAMMPS
    '''

import os, numpy as np, time, cPickle as pickle
from multiprocessing import Pool
from functools import partial

import sim, protein
import pickleTraj

import topo, mapNCOS
from const import *

Verbose = True

# LAMMPS settings
sim.srel.optimizetrajlammps.useREMD = True
sim.export.lammps.InnerCutoff = 0.02
sim.export.lammps.MaxPairEnekBT = MaxLammpsPairEnekBT

# REMD base class
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
            if Verbose: print '\nGenerating fully extended initial AA structure...'
            # create an ALL-ATOM protein class object for the given sequence
            pobj = protein.ProteinClass(Seq = self.p.Seq)
            pobj.WritePdb(tmpPdb)
            # coarse grain all-atom proteinclass object
            mapNCOS.Map(InPdb = tmpPdb, CGPrefix = self.InitPdb.split('.pdb')[0], hasPseudoGLY = self.p.hasPseudoGLY)
            # remove all-atom pdb
            if os.path.isfile(tmpPdb): os.remove(tmpPdb)
            del pobj
        print '\nUsing init conf as generated in : %s' % self.InitPdb
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


# Re-ordering Lammps replica trajectories by temp
# need to be placed outside the base class to be
# parallelizable using multiprocessing
def Pickle(RepInd, Prefix):
    ''' pickles sim-style traj objects from each replica traj'''
    # filenames
    print 'Pickling Replica %d' % RepInd
    TrajFile = Prefix + '.lammpstrj'
    LogFile = Prefix + 'lammps.log'
    myTrajFn = TrajFile + '.%d.gz' % RepInd
    myLogFn = LogFile + '.%d' % RepInd
    pickleTraj(myTrajFn, LogFile = myLogFn, LogFileToken = '#run production')
    return

def Reorder(Temp, Prefix, TempFile, NStepsEquil, NStepsProd, NStepsSwap, StepFreq):
    ''' reorders trajectories at a particular given Temp.
    Temp needs to be the first input to enable use with partial() while spawning multiprocessing.Pool'''
    # temp schedule
    Temps = np.loadtxt(TempFile)
    TempInd = [list(Temps).index(t) for t in list(Temps) if '%3.2f' % t  == '%3.2f' % Temp][0]
   # Trj objects
    TrjList = []
    TrajFile = Prefix + '.lammpstrj'
    for i in range(len(Temps)):
        TrjPickle = TrajFile + '.%d.gz.pickle' % i
        with open(TrjPickle, 'r') as of: Trj = pickle.load(of)
        TrjList.append(Trj)
    # filenames
    MultiTrajFn = FMT['TRAJ'] % (Prefix, Temp)
    MultiEneFn = FMT['ENE'] % (Prefix, Temp)
    if os.path.isfile(MultiTrajFn) and os.path.isfile(MultiEneFn): return
    print 'Reordering replica frames at temperature %3.2f' %  Temp
    # collect all replica indices at this temp
    LogFile = Prefix + 'lammps.log'
    RepIndsMaster = [np.where(x[1:] == TempInd)[0][0] for x in np.loadtxt(LogFile, skiprows = 3)]
    # keep only replica inds for the production run part
    RepInds = RepIndsMaster[int(NStepsEquil / NStepsSwap) : ]
    RepInds = RepInds[:-1]
    # output lists
    OutTrjList = []
    OutEneList = []
    # case 1
    if StepFreq <= NStepsSwap: # assume mod(NStepsSwap, StepFreq) = 0
        for ii, i in enumerate(RepInds):
            # point to the Trj and Ene object for this replica
            this_Trj = TrjList[i]
            this_Ene = this_Trj.ThermoDict['PEnergy']
            # keep only the parts for the production run
            this_Trj_prod = this_Trj[int(NStepsEquil / StepFreq) : ]
            this_Ene_prod = this_Ene[int(NStepsEquil / StepFreq) : ]
            # find start and stop indices for this block 
            start = ii * NStepsSwap / StepFreq
            stop = (ii+1) * NStepsSwap / StepFreq
            # write into output lists
            OutTrjList.append(this_Trj_prod[start:stop])
            OutEneList.extend(this_Ene_prod[start:stop])
    # case 2
    else:
        NSkip = StepFreq / NStepsSwap # assume mod(StepFreq, NStepsSwap) = 0
        for ii, i in enumerate(RepInds[0::NSkip]):
            # point to the Trj and Ene object for this replica
            this_Trj = TrjList[i]
            this_Ene = this_Trj.ThermoDict['PEnergy'] 
            # keep only the parts for the production run
            this_Trj_prod = this_Trj[int(NStepsEquil / StepFreq) : ]
            this_Ene_prod = this_Ene[int(NStepsEquil / StepFreq) : ]
            # find start and stop indices for this block
            start = ii
            stop = ii+1
            OutTrjList.append(this_Trj_prod[start:stop])
            OutEneList.extend(this_Ene_prod[start:stop])
    MultiTrj = sim.traj.Multi(OutTrjList)
    sim.traj.base.Convert(InTraj = MultiTrj, OutTrajClass = sim.traj.LammpsWrite, FileName = MultiTrajFn)
    np.savetxt(MultiEneFn, OutEneList, fmt = '%11.4e')
    return 

def ReorderAll(Prefix, NStepsEquil, NStepsProd, NStepsSwap, StepFreq, TempFile = 'temps.txt', ReorderTemps = None):
    Temps = np.loadtxt(TempFile)
    if ReorderTemps is None: ReorderTemps = Temps
    # pickle all traj in parallel
    worker_1 = partial(Pickle, Prefix = Prefix)
    pool_1 = Pool(processes = len(Temps))
    pool_1.map(worker_1, range(len(Temps)))
    pool_1.close()
    # reorder all traj in parallel
    worker_2 = partial(Reorder, Prefix = Prefix, TempFile = TempFile,
                       NStepsEquil = NStepsEquil, NStepsProd = NStepsProd, 
                       NStepsSwap = NStepsSwap, StepFreq = StepFreq)
    pool_2 = Pool(processes = len(ReorderTemps))
    pool_2.map(worker_2, ReorderTemps)
    pool_2.close()
    return
