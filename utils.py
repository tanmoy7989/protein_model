#!/usr/bin/env python

''' This module mostly provides information
to set up paths, forcefields, etc throughout the project.
Basically it contains code snippets wrapped up in functions that need
to be used over and over again
'''

import os, sys, numpy as np, cPickle as pickle, shelve
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
    
# file formats throughout this project
FMT = {
       'GENERAL'        : '%s.%3.2f',
       'TRAJ'           : '%s.%3.2f.lammpstrj.gz',
       'PDB'            : '%s.%3.2f.pdb',
       'ENE'            : '%s.%3.2f.ene.dat.gz',
       'RMSD'           : '%s.%3.2f.RMSD.pickle',
       'RAMA'           : '%s.%3.2f.rama.pickle',
       'RESCONTACTS'    : '%s.%3.2f.rescontacts.pickle',
       'CLUSTPDB'       : '%s.%3.2f.clust.pdb',
       'CLUSTSUM'       : '%s.%3.2f.clustsum.txt',
       'DATASHELF'   :    '%s.ordparam.shelf',
       'FOLDCURVE'      : '%s_%s.foldcurve.pickle',
       'PMF1D'          : '%s.%3.2f.%s.pickle',
       'PMF2D'          : '%s.%3.2f.%s_%s.pickle'
      }

# forcefield types
BBTYPES = ['ff_ala_spc', 'ff_ala15', 'ff_leu15', 'ff_val15',
           'ff_leu15_bbonly', 'ff_val15_bbonly']
GoFFTYPES = ['ff_leu15_protg_lj', 'ff_leu15_protg_spline', 'ff_ala_spc_protg_lj', 'ff_ala_spc_protg_spline',
             'ff_leu15_bbonly_protg_lj', 'ff_leu15_bbonly_protg_spline',
             'ff_val15_bbonly_protg_lj', 'ff_val15_bbonly_protg_spline']

# important master paths and metadata files
FFDIR = os.path.expanduser('~/protein_model/cgff')
FFMETADATAFILE = os.path.expanduser('~/protein_model/cgff/ffs/ff_metadata.txt')
NATIVEDIR = os.path.expanduser('~/protein_model/native_struct/mapped')
NATIVEMETADATAFILE = os.path.expanduser('~/protein_model/native_struct/native_metadata.txt')


def parseAA(PolyPrefix, MasterDir = None, TempSet = None):
    if TempSet is None: TempSet = 300.0
    # extract the trajectory closest to TempSet
    if MasterDir is None: MasterDir = os.path.expanduser('~/protein_model/%s_AA' % PolyPrefix)
    TempFile = os.path.join(MasterDir, 'Lammps', 'temps.txt')
    Temps = np.loadtxt(TempFile)
    Ind = np.argmin(abs(Temps - TempSet))
    TrajTemp = Temps[Ind]
    TrajFn = os.path.join(MasterDir, 'Lammps', FMT['TRAJ'] % (PolyPrefix, TrajTemp))
    # extract energy data, pdb and clustered pdbs if present
    EneFn = os.path.join(MasterDir, 'Lammps', FMT['ENE'] % (PolyPrefix, TrajTemp))
    Pdb = os.path.join(MasterDir, 'Lammps', FMT['PDB'] % (PolyPrefix, TrajTemp))
    Pdb_generic = os.path.join(MasterDir, 'Lammps', PolyPrefix + '.pdb')
    # check if these exist
    if not os.path.isfile(TrajFn): TrajFn = None
    if not os.path.isfile(EneFn): EneFn = None
    if not os.path.isfile(Pdb):
        if os.path.isfile(Pdb_generic): Pdb = Pdb_generic
        else: Pdb = None
    # check if enough replica information is present
    hasReplica = True
    for T in Temps:
        this_TrajFn = os.path.join(MasterDir, 'Lammps', FMT['TRAJ'] % (PolyPrefix, T))
        this_EneFn = os.path.join(MasterDir, 'Lammps', FMT['ENE'] % (PolyPrefix, T)) 
        if not (os.path.isfile(this_TrajFn) and os.path.isfile(this_EneFn)):
            hasReplica = False
            print this_TrajFn, this_EneFn, ' does not exist'
            break
    
    ret = {'MasterDir': MasterDir, 'TempFile': TempFile, 'Traj': TrajFn, 
           'Ene': EneFn, 'Pdb': Pdb, 'Temp': TrajTemp, 'TempInd': Ind,
           'hasReplica': hasReplica}
    
    return ret

def parseCG(PolyPrefix, MasterDir = None, TempSet = None):
    if TempSet is None: TempSet = 300.0
    # extract the trajectory closest to TempSet
    if MasterDir is None: MasterDir = os.path.join(FFDIR, '%s_modtraj' % PolyPrefix)
    # CG simulations absent
    if not os.path.isdir(MasterDir):
        ret = {'hasReplica': False}
        return ret
    TempFile = os.path.join(MasterDir, 'temps.txt')
    Temps = np.loadtxt(TempFile)
    Ind = np.argmin(abs(Temps - TempSet))
    TrajTemp = Temps[Ind]
    TrajFn = os.path.join(MasterDir, FMT['TRAJ'] % (PolyPrefix, TrajTemp))
    # extract energy data, pdb and clustered pdbs if present
    EneFn = os.path.join(MasterDir, FMT['ENE'] % (PolyPrefix, TrajTemp))
    # check if these exist
    if not os.path.isfile(TrajFn): TrajFn = None
    if not os.path.isfile(EneFn): EneFn = None
    # check if enough replica information is present
    hasReplica = True
    for T in Temps:
        this_TrajFn = os.path.join(MasterDir, FMT['TRAJ'] % (PolyPrefix, T))
        this_EneFn = os.path.join(MasterDir,  FMT['ENE'] % (PolyPrefix, T)) 
        if not (os.path.isfile(this_TrajFn) and os.path.isfile(this_EneFn)):
            hasReplica = False
            break
            
    ret = {'MasterDir': MasterDir, 'TempFile': TempFile, 'Traj': TrajFn, 
           'Ene': EneFn, 'Temp': TrajTemp, 'TempInd': Ind, 'hasReplica': hasReplica}
    return ret

def parseBBFF(BBType, MasterDir = None):
    if not BBTYPES.__contains__(BBType):
        raise IOError('Requested backbone forcefield not present in repository')
    if MasterDir is None: MasterDir = os.path.join(FFDIR, 'ffs')
    FF_File = os.path.join(MasterDir, BBType + '_BB.dat')
    mdata = eval(file(FFMETADATAFILE).read())
    ret = FF_File, mdata[BBType]
    return ret

def parseGoFF(FFType, MasterDir = None):
    if not GoFFTYPES.__contains__(FFType):
        raise IOError('Requested Go forcefield not present in repository')
    if MasterDir is None: MasterDir = os.path.join(FFDIR, 'ffs')
    FF_File = os.path.join(MasterDir, FFType + '.dat')
    mdata = eval(file(FFMETADATAFILE).read())
    ret = FF_File, mdata[FFType]
    return ret

def parseNative(PdbName, MasterDir = None):
    if MasterDir is None: MasterDir = NATIVEDIR
    Pdb = os.path.join(MasterDir, PdbName+'.pdb')
    if not os.path.isfile(Pdb):
        raise IOError('Requested native pdb not present in repository')
    return Pdb

def hasCGData(Prefix, Key, OutDir = None):
    if OutDir is None: OutDir = os.getcwd()
    OutShelf = os.path.join(OutDir, Prefix + '.shelf')
    hasData = True
    if os.path.isfile(OutShelf):
        oshelf = shelve.open(OutShelf)
        if not oshelf.has_key(Key): hasData = False
        oshelf.close()
    else: hasData = False
    if hasData == True:
        oshelf = shelve.open(OutShelf)
        Data = oshelf[Key]
    else: Data = None
    return hasData, Data

def getPanelLayout(pset_type):
    native_metadata = eval(file(NATIVEMETADATAFILE).read())
    pset_small = native_metadata['small_set']
    pset_large = native_metadata['large_set']
    layout = {'small_set': {'NRows': 2, 'NCols': 5, 'pset': pset_small},
              'large_set': {'NRows': 3, 'NCols': 3, 'pset': pset_large}
              }
    return layout[pset_type]

def checkGoREMD(Prefix, Temps = None):
    Dir = os.path.dirname(Prefix)
    TempFile = os.path.join(Dir, 'temps.txt')
    if Temps is None: Temps = np.loadtxt(TempFile)
    isDone = True
    for T in Temps:
        thisTrajFn = FMT['TRAJ'] % (Prefix, T)
        thisEneFn = FMT['ENE'] % (Prefix, T)
        if not (os.path.isfile(thisTrajFn) and os.path.isfile(thisEneFn)):
            isDone = False
            break
    return isDone

def isConverged(TrajDir, PdbName, TempSet = 300, OutFile = None):
    if OutFile is None: OutFile = os.path.abspath('./rmsd_convergence.png')
    import cgprotein as lib, measure as m
    TempFile = os.path.join(TrajDir, 'temps.txt')
    Temps = np.loadtxt(TempFile)
    Ind = np.argmin(abs(Temps - TempSet))
    TrajTemp = Temps[Ind]
    Traj = os.path.join(TrajDir, 'prot_%s.%3.2f.lammpstrj.gz' % (PdbName, TrajTemp))
    NativePdb = parseNative(PdbName)
    LammpsREMDLog = os.path.join(TrajDir, 'prot_%slammps.log' % PdbName)
    
    calc = lib.Compute(TrajFn = Traj, NativePdb = NativePdb, Temp = TrajTemp)
    rmsd = calc.RMSD_frame()
    avgrmsd = rmsd.mean()
    stdrmsd = np.std(rmsd, ddof = 1)
    avgerr = stdrmsd / avgrmsd
    avgerr *= 100.
    m.NBins = 50 ; m.NBlocks = 4 ; m.Normalize = True
    rmsdhist = m.makeHist(rmsd)
    
    Walk = np.loadtxt(LammpsREMDLog, skiprows = 3)[:, 1:]
    
    fig = plt.figure(facecolor = 'w', edgecolor = 'w', figsize = (8, 4))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax1.plot(rmsd, lw = 1, color = 'blue', label = 'fluct from mean = %3.2f %%' % avgerr)
    ax1.axhline(avgrmsd, ls = '--', color = 'black', lw = 2)
    ax1.legend(loc = 'best')
    ax2.errorbar(rmsdhist[0], rmsdhist[1], yerr = rmsdhist[2], color = 'black', lw = 2, marker = 'o', markersize = 5)
    ax2.set_xlabel(r'$RMSD (\AA)$') ; ax2.set_ylabel('distribution')
    for i, T in enumerate(Temps):
        x = Walk[:,i]
        y = np.array([Temps[int(k)] for k in x])
        m.NBins = len(Temps) ; m.NBlocks = 1; m.Normalize = False
        walkhist = m.makeHist(y)
        ax3.plot(walkhist[0], walkhist[1], lw = 2, label = 'Replica: %d' % i)
        ax3.legend(loc = 'best', fontsize = 5)
    ax3.set_xlabel('Temp (K)') ; ax3.set_ylabel('distribution')
    fig.tight_layout()
    plt.savefig(OutFile)
    return
    

#### COMMAND LINE USAGE ####
if __name__ == '__main__':
    Op = sys.argv[1]

    # check rmsd convergence
    if Op == 'isconverged':
        TrajDir = os.path.abspath(sys.argv[2])
        PdbName = sys.argv[3]
        OutFile = os.path.abspath(sys.argv[4]) if len(sys.argv) > 4 else None
        isConverged(TrajDir = TrajDir, PdbName = PdbName, OutFile = OutFile)
