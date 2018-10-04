'''/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <tanmoy.7989@gmail.com> wrote this file.  As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return. Tanmoy Sanyal
 * ----------------------------------------------------------------------------
 */
'''

#!/usr/bin/env python

import os, sys, cPickle as pickle, shelve, numpy as np, copy
import sim, protein, measure, pickleTraj, utils
import cgproteinlib as lib
import pymbar, whamlib

# file formats
FMT = utils.FMT

# ascii magic for formatted output
UP_ONE_LINE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

# CG O-N "peptide" bond length
# (all forcefields roughly predict this value)
PeptideBondLen = 1.75 #A

# Boltzmann constant
kB = 0.001987

# Cutoff for metrics
OCut = {'RMSD': 2.0, 'Rg': 2.0, 'REE': 2.0}
MaxCluster = 10

# Histogramming
NBlocks = 4
NBins = 50

# LogFile reading style for getting energies
LogFileStyle = 'mystyle' # 'simstyle' for using sim.traj.Lammps



# CG protein object that defines basic topology and computation routines for order params
class ProteinNCOS(object):
    ResRadius = 8.0
    MinCO = 3
    
    def __init__(self, Pdb, Model = None, hasPseudoGLY = False):
        # protected (internal) proteinclass object
        self.Pdb = Pdb
        self.__p = protein.ProteinClass(Pdb, Model = Model)
        
        # copy over some attributes from the original proteinclass
        self.AtomNames = self.__p.AtomNames()
        self.AtomRes = self.__p.AtomRes()
        self.AtomResNum = self.__p.AtomResNum
        self.ChainResNums = self.__p.ChainResNums
        self.ChainAtomNums = self.__p.ChainAtomNums
        self.Chains = self.__p.Chains
        
        # sequence manipulation
        self.UpdateSeq()
        
        # set atom start indices
        self.hasPseudoGLY = hasPseudoGLY
        self.StartInds = []
        self.GetStartInds()
        
        # bacbone indices
        self.BBInds = self.GetBBInds()
        
        # Co-ordinates manipulation
        self.Pos = None
        self.UpdatePos()
        
        # chains
        self.NChains = len(self.ChainResNums) - 1
        if self.NChains == 1:
            self.AtomChain = [' '] * len(self.AtomNames)
        else:
            self.AtomChain = []
            for i, chain in enumerate(self.Chains):
                n = len(range(self.ChainAtomNums[i], self.ChainAtomNums[i+1]))
                self.AtomChain.extend( [chain] * n)
    
        
    def UpdateSeq(self):
        # updates sequence and residue type counts
        self.Seq = self.__p.Seq
        self.NRes = len(self.Seq)
        self.ResTypes = list(set(self.Seq))
        self.ResCount = dict( (x, self.Seq.count(x)) for x in set(self.Seq) )
    
    def WritePdb(self, Filename):
        self.__p.WritePdb(Filename)

    def UpdatePos(self, Pos = None):
        # updates coordinates and backbone positions
        # can be also called by resetting self.Pos directly
        if Pos is None:
            Pos = self.__p.Pos if self.Pos is None else self.Pos
        self.__p.Pos = Pos
        self.ResPos = self.__p.ResPos()
        self.Pos = Pos
        self.BBPos = self.GetBBPosSlice()
   
    def GetStartInds(self):
        # don't let this get calculated twice
        if len(self.StartInds) == self.NRes: return
        x = 0
        for i, r in enumerate(self.Seq):
            self.StartInds.append(x)
            if r == 'GLY' and not self.hasPseudoGLY: x+= 3
            else: x += 4
        return
    
    def GetBBInds(self, ResNums = None):
        if ResNums is None: ResNums = range(0, self.NRes)
        inds = []
        for i in ResNums:
            x = self.StartInds[i]
            inds.extend([x, x+1, x+2])
        return inds

    def GetPosSlice(self, ResNums = None):
        # returns coordinates of given residue-atoms
        if ResNums is None: ResNums = range(0, self.NRes)
        return self.__p[ResNums[0]: ResNums[-1]+1].Pos
    
    def GetBBPosSlice(self, ResNums = None): 
        # returns backbone co-ordinates of given residue atoms
        inds = self.GetBBInds(ResNums)
        return self.Pos[inds]
   
    def GetResContacts(self, BoxL = np.zeros(3)):
        # get residue contact map and distances
        ContactDist = np.zeros([self.NRes, self.NRes], np.float64)
        ContactMap = np.zeros([self.NRes, self.NRes], np.float64)
        ContactMap, ContactDist = lib.respairs_frame(respos = self.ResPos, resradius = self.ResRadius, minco = self.MinCO, boxl = BoxL)
        return ContactMap, ContactDist
    
    def GetCO(self, ContactMap = None, BoxL = np.zeros(3)):
        # not sure if defined for an oligomer
        if ContactMap is None:
            ContactMap, ContactDist = self.GetResContacts(BoxL = BoxL) 
        L = self.NRes # chain length
        N = np.sum(ContactMap) # number of contacts
        inds = np.nonzero(ContactMap)
        inds = zip(inds[0], inds[1])
        S = np.sum( np.array([abs(i-j) for (i,j) in inds]) )
        CO = S / (N * L) # Plaxco, Simmons, Baker, J. Mol. Biol. (1998) 277, 985-994 
        return CO
            
    def Rg(self, BoxL = np.zeros(3), BB = False):
        # calculates overall Rg (of the entire oligomer when appropriate)
        Pos = self.Pos if not BB else self.BBPos
        return lib.rg_frame(pos = Pos, boxl = BoxL)
    
    def REE(self, BoxL = np.zeros(3), BB = True):
        # calculates end-to-end distance along backbone (of the entire oligomer when appropriate)
        # setting BB = False would be bizarre
        Pos = self.Pos if not BB else self.BBPos
        return lib.ree_frame(pos = self.BBPos, boxl = BoxL) 
    
    def QuickRMSD(self, other):
        # returns backbone RMSD between 2 CG proteins
        Pos1 = other.BBPos
        Pos2 = self.BBPos
        return sim.geom.RMSD(Pos1, Pos2)
    
    def RMSD_Res(self, other, ResNums = None):
        # align complete structures but return RMSD
        # contribution only from requested residues
        if ResNums is None: ResNums = range(self.NRes)
        Pos1 = other.BBPos
        Pos2 = self.BBPos
        PosVec1, PosVec2, RotMat, Residuals = sim.geom.AlignmentRMSD(Pos1 = Pos1, Pos2 = Pos2)
        Pos1_ = Pos1 + PosVec1
        Pos2_ = np.dot(Pos2 + PosVec2, RotMat)
        # get per-atom RMSD for backbone atoms
        rmsd_atom = np.zeros(len(self.BBInds))
        for ii, i in enumerate(self.BBInds):
            rmsd_atom[ii] = sim.geom.dRMSD(Pos1_[i], Pos2_[i])
        # accumulate per-atom RMSDs over residues
        rmsd_res = 0.0 ; natoms = 0
        for i, r in enumerate(self.Seq):
            if not ResNums.__contains__(i): continue
            rmsd += rmsd_atom[i*3 : (i+1)*3]
            natoms += 1
            rmsd_res /= float(natoms)
        return rmsd_res
    
    def RMSD_ResType(self, other, ResTypes = None):
        # returns backbone RMSD averaged over ResTypes
        if ResTypes is None: ResTypes = self.ResTypes
        # align the 2 structures
        Pos1 = other.BBPos
        Pos2 = self.BBPos
        PosVec1, PosVec2, RotMat, Residuals = sim.geom.AlignmentRMSD(Pos1 = Pos1, Pos2 = Pos2)
        Pos1_ = Pos1 + PosVec1
        Pos2_ = np.dot(Pos2 + PosVec2, RotMat)
        # get per-atom RMSD for backbone atoms
        rmsd_atom = np.zeros(len(self.BBInds))
        for i in range(len(self.BBInds)):
            rmsd_atom[i] = sim.geom.dRMSD(Pos1_[i], Pos2_[i])
        # accumulate per-atom RMSDs over residues
        rmsd_restype = dict((r, 0.0) for r in ResTypes)
        for i, r in enumerate(self.Seq):
            if not ResTypes.__contains__(r): continue
            rmsd_restype[r] += np.mean(rmsd_atom[i*3 : (i+1)*3])
        # average over number of each residue
        for r in ResTypes:
            rmsd_restype[r] /= float(self.ResCount[r])
        return rmsd_restype
    
    def GetPhiPsi(self, ResNums = None, BoxL = np.zeros(3)):
        if ResNums is None: ResNums = range(self.NRes)
        # exclude residues at chain termini
        for i in (self.ChainResNums + [0, self.NRes-1]):
            if ResNums.__contains__(i): ResNums.remove(i)
        # get the entire backbone pos coordinates
        BBPos = self.GetBBPosSlice()
        # minimage for periodic box
        if all(BoxL): BBPos = lib.reimagechain(pos = BBPos, boxl = BoxL)
        Phi = np.zeros(len(ResNums))
        Psi = np.zeros(len(ResNums))
        for i, r in enumerate(ResNums):
            Pos = BBPos[r*3 : (r+1)*3]
            PosN = Pos[0] ; PosC = Pos[1] ; PosO = Pos[2]
            PosO_prev = BBPos[(r-1)*3 + 2]    
            PosN_next = BBPos[(r+1)*3]
            Phi[i] = sim.geom.Dihedral(PosO_prev, PosN, PosC, PosO)
            Psi[i] = sim.geom.Dihedral(PosN, PosC, PosO, PosN_next)
            Phi[i] = sim.geom.NormRad(Phi[i])
            Psi[i] = sim.geom.NormRad(Psi[i])
        return Phi, Psi
        
    def GetPhiPsiDiff(self, other, ResNums = None):
        if ResNums is None: ResNums = range(self.NRes)
        Phi1, Psi1 = other.GetPhiPsi(ResNums = ResNums)
        Phi, Psi = self.GetPhiPsi(ResNums = ResNums)
        PhiDiff = np.zeros(len(ResNums))
        PsiDiff = np.zeros(len(ResNums))
        for i in ResNums:
            PhiDiff[i] = sim.geom.NearestAngle(Phi[i] - Phi1[i], 0.0)
            PsiDiff[i] = sim.geom.NearestAngle(Psi[i] - Psi1[i], 0.0)
        Diff = np.sqrt(PhiDiff**2 + PsiDiff**2)
        return PhiDiff, PsiDiff, Diff
            
      
         
# class that defines computes over trajectories at a fixed temperature
# functions prefixed Quick, calculate globally over the entire backbone
# functions suffixed _frame only calculates per frame data and no histograms
class Compute(object):
    Recompute = False
    def __init__(self, NativePdb, TrajFn = None, EneFn = None, Temp = None, Prefix = 'compute', hasPseudoGLY = False):
        # parse out-prefix
        self.Prefix = os.path.abspath(Prefix)
        self.OutPrefix = self.Prefix.split('/')[-1]
        self.OutDir = os.path.dirname(self.Prefix)
        
        # placeholders
        self.Trj = None
        self.Ene = None
        self.NFrames = None
        self.FrameRange = None
        self.BoxL = 0.
        
        # parse traj
        if not TrajFn is None:
            self.TrajFn = os.path.abspath(TrajFn)
            self.Trj = pickleTraj(self.TrajFn)
            self.FrameRange = range(0, len(self.Trj))
            self.NFrames = len(self.FrameRange)
            self.BoxL = self.Trj.FrameData.get('BoxL', np.zeros(3))
        
        # energies
        if not EneFn is None:
            self.EneFn = EneFn
            self.Ene = np.loadtxt(self.EneFn)
            if self.FrameRange is None and self.NFrames is None:
                self.FrameRange = range(0, len(self.Ene))
                self.NFrames = len(self.FrameRange)
        
        # record temp
        self.Temp = Temp
        
        # cg protein objects for native and predicted
        self.pNative = ProteinNCOS(NativePdb, hasPseudoGLY = hasPseudoGLY)
        self.p = ProteinNCOS(NativePdb, hasPseudoGLY = hasPseudoGLY)
    
    def ProgressMonitor(self, OrderParamText):
        if not self.Temp is None:
            s = 'Calculating %s for %s at %3.2f K...' % (OrderParamText, self.OutPrefix, self.Temp)
        else:
            s = 'Calculating %s for %s...' % (OrderParamText, self.OutPrefix)
        return s
   
    def Update(self, TrajFn = None, Temp = None):
        # updates the compute object if traj or temp is changed
        # if native pdb is changed then better create another instance
        if not Temp is None: self.Temp = Temp
        if not TrajFn is None:
            self.TrajFn = os.path.abspath(TrajFn)
            self.Trj = pickleTraj(self.TrajFn)
            self.FrameRange = range(0, len(self.Trj))
            self.NFrames = len(self.FrameRange)
            self.BoxL = self.Trj.FrameData.get('BoxL', np.zeros(3))
    
    def Rg_frame(self):
        # per frame Rg calculator
        rg_frame = np.zeros(self.NFrames)
        # loop over frames
        pb = sim.utility.ProgressBar(Text = self.ProgressMonitor('Rg'), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            rg_frame[i] = self.p.Rg(BoxL = self.BoxL)
            pb.Update(i)
        return rg_frame
    
    def REE_frame(self):
        # per frame REE calculator
        ree_frame = np.zeros(self.NFrames)
        # loop over frames
        pb = sim.utility.ProgressBar(Text = self.ProgressMonitor('REE'), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            ree_frame[i] = self.p.REE(Boxl = self.BoxL)
            pb.Update(i)
        return ree_frame
            
    def RMSD_frame(self):
        # per frame rmsd calculator
        rmsd_frame = np.zeros(self.NFrames)
        # loop over frames
        pb = sim.utility.ProgressBar(Text = self.ProgressMonitor('RMSD'), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            rmsd_frame[i] = self.p.QuickRMSD(self.pNative)
            pb.Update(i)
        return rmsd_frame
    
    def PhiPsi_frame(self):
        # per frame dihedral angle calculator
        phi_frame = np.zeros(self.NFrames)
        psi_frame = np.zeros(self.NFrames)
        NExcludedAtoms = len(self.pNative.ChainResNums)
        Phi = np.zeros([self.NFrames, self.p.NRes-NExcludedAtoms], np.float64)
        Psi = np.zeros([self.NFrames, self.p.NRes-NExcludedAtoms], np.float64)
        pb = sim.utility.ProgressBar(Text = self.ProgressMonitor('Phi and Psi angles'), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            Phi[i, :], Psi[i, :] = self.p.GetPhiPsi(BoxL = self.BoxL)
            pb.Update(i)
        Phi = Phi.flatten()
        Psi = Psi.flatten()
        return Phi, Psi
    
    def ResContacts_frame(self):
        # makes sense to store the entire per frame data since
        # multiple computations can be done on this entire data set
        # and those calculations need not be saved to a pickle
        picklename = FMT['RESCONTACTS'] % (self.Prefix, self.Temp)
        if os.path.isfile(picklename) and not self.Recompute: return
        # loop over frames to get average CG contact distances   
        ContactMap = np.zeros([self.NFrames, self.p.NRes, self.p.NRes], np.float64)
        ContactDist = np.zeros([self.NFrames, self.p.NRes, self.p.NRes], np.float64)
        pb = sim.utility.ProgressBar(Text = self.ProgressMonitor('contacts'), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            x, y = self.p.GetResContacts(BoxL = self.BoxL)
            ContactMap[i,:,:] = x
            ContactDist[i,:,:] = y
            pb.Update(i)
        ret = (ContactMap, ContactDist)
        with open(picklename, 'w') as of: pickle.dump(ret, of)
        return
        
    def QuickRMSD(self):
        # get overall RMSD distribution, but don't save to file
        # should be very fast and callable from replica and other external routines
        rmsd_frame = np.zeros(self.NFrames)
        rmsd_hist = None
        # loop over frames
        pb = sim.utility.ProgressBar(Text = self.ProgressMonitor('RMSD'), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            rmsd_frame[i] = self.p.QuickRMSD(self.pNative)
            pb.Update(i)
        # histogram for overall RMSD
        measure.NBins = NBins
        measure.NBlocks = NBlocks
        measure.NFrames = self.NFrames
        rmsd_hist = measure.makeHist(rmsd_frame)
        return rmsd_hist
        
    def QuickRMSD_Res(self, ResNums = None):
        if ResNums == []: return None
        # get overall RMSD distribution only for selected residues
        rmsd_res_frame = np.zeros(self.NFrames)
        # loop over frames
        pb = sim.utility.ProgressBar(Text = self.ProgressMonitor('RMSD for select residues'), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            rmsd_res_frame[i] = self.p.RMSD_Res(self.pNative, ResNums = ResNums)
            pb.Update(i)
        return rmsd_res_frame
        
    def RMSD(self, ResTypes = None):
        # get overall RMSD distribution, and per-residue-type RMSD
        print '\n\n'
        if ResTypes is None: ResTypes = self.pNative.ResTypes
        picklename = FMT['RMSD'] % (self.Prefix, self.Temp) 
        if os.path.isfile(picklename) and not self.Recompute: return
        rmsd_hist = self.QuickRMSD()
        rmsd_restype_frame = dict((r, np.zeros(self.NFrames)) for r in ResTypes)
        rmsd_restype_block = dict((r, np.zeros(NBlocks)) for r in ResTypes)
        rmsd_restype = dict((r, None) for r in ResTypes)
        rmsd_restype_err = dict((r, 0.0) for r in ResTypes)
        # loop over frames
        pb = sim.utility.ProgressBar(Text = self.ProgressMonitor('detailed RMSD by residue number and type'), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            x = self.p.RMSD_ResType(self.pNative, ResTypes = ResTypes)
            for k in rmsd_restype_frame.keys(): 
                rmsd_restype_frame[k][i] = x[k]
            pb.Update(i)
        
        # block average per-residue-type RMSD
        BlockSize = int(self.NFrames / NBlocks)
        for b in range(NBlocks):
            start = b * BlockSize
            stop = (b+1) * BlockSize if not b == NBlocks - 1 else self.NFrames
            for r in ResTypes:
                rmsd_restype_block[r][b] = np.mean(rmsd_restype_frame[r][start:stop])
        for r in ResTypes: rmsd_restype[r] = np.mean(rmsd_restype_block[r])
        if NBlocks > 1:
            for r in ResTypes: rmsd_restype_err[r] = np.std(rmsd_restype_block[r], ddof = 1)
        
        # write to pickle
        ret = rmsd_hist, (rmsd_restype, rmsd_restype_err)
        with open(picklename, 'w') as of: pickle.dump(ret, of)
        return
        
    def Cluster(self):
        clustfilename = FMT['CLUSTPDB'] % (self.Prefix, self.Temp)
        sumfilename = FMT['CLUSTSUM'] % (self.Prefix, self.Temp)
        if os.path.isfile(clustfilename) and os.path.isfile(sumfilename) and not self.Recompute: return
        clustret = sim.cluster.ClusterMSS(Trj = self.Trj, Cutoff = OCut['RMSD'], MaxCluster = MaxCluster)
        # prep headvars
        HeadVars = {'NAtom': len(self.pNative.AtomNames), 
                    'AtomNames': self.pNative.AtomNames,
                    'AtomIsProtein': [True]*len(self.pNative.AtomNames),
                    'AtomRes': self.pNative.AtomRes,
                    'AtomResNum': self.pNative.AtomResNum,
                    'AtomChain': self.pNative.AtomChain}
        # write out clust results
        sim.cluster.WriteClustResults(self.Trj, clustret, sumfilename, clustfilename,
                                      HeadVars = HeadVars)
        # append connectivity to the end of the file
        s = file(self.pNative.Pdb, 'r').readlines()
        start = [s.index(i) for i in s if i.startswith('CONECT')][0]
        stop = len(s)
        connectstr = ''.join(s[start:stop])
        clustpdbstr = file(clustfilename).read()
        file(clustfilename, 'w').write(clustpdbstr + '\n' + connectstr)
        return

    def CompareContactMap(self, CompType = 'traj'):
        # get AA contact map
        NativeContactMap, NativeContactDist = self.pNative.GetResContacts()
        # get CG contact map
        if CompType == 'traj':
            self.ResContacts_frame()
            ContactPickle = FMT['RESCONTACTS'] % (self.Prefix, self.Temp)
            with open(ContactPickle, 'r') as of: data = pickle.load(of)
            ContactMap, ContactDist = data
            ClustContactMap = np.mean(ContactMap, axis = 0)
        if CompType == 'topclust':
            self.Cluster()
            ClustPdb = FMT['CLUSTPDB'] % (self.Prefix, self.Temp)
            pClust = ProteinNCOS(ClustPdb, Model = 1)
            ClustContactMap, ClustContactDist = pClust.GetResContacts()
        return NativeContactMap, ClustContactMap

    def GetFracNativeContacts(self):
        # get AA contact map
        NativeContactMap, NativeContactDist = self.pNative.GetResContacts()
        ind = (NativeContactMap == 1)
        # get CG contact map
        self.ResContacts_frame()
        # get per frame fraction of native contacts
        ContactPickle = FMT['RESCONTACTS'] % (self.Prefix, self.Temp)
        with open(ContactPickle, 'r') as of: data = pickle.load(of)
        ContactMap, ContactDist = data
        NFrames = ContactMap.shape[0]
        frac = np.zeros(NFrames)
        for i in range(NFrames):
            frac[i] = np.sum(ContactMap[i,:,:][ind]) / np.sum(NativeContactMap)
        measure.NBins = NBins
        measure.NBlocks = NBlocks
        measure.Normalize = True
        hist = measure.makeHist(frac)
        return hist
    
    def GetCO(self):
        self.ResContacts_frame()
        ContactPickle = FMT['RESCONTACTS'] % (self.Prefix, self.Temp)
        with open(ContactPickle, 'r') as of: data = pickle.load(of)
        ContactMap, ContactDist = data
        NFrames = ContactMap.shape[0]
        CO = np.zeros(NFrames)
        for i in range(NFrames):
            CO[i] = self.p.GetCO(ContactMap = ContactMap[i,:,:])
        # remove nans and infs, if any
        CO = CO[~np.isnan(CO)]
        CO = CO[~np.isinf(CO)]
        measure.NBins = NBins
        measure.NBlocks = NBlocks
        measure.Normalize = True
        hist = measure.makeHist(CO)
        return hist
    
    def CompareCO(self):
        # get native contact order
        NativeContactMap, NativeContactDist = self.pNative.GetResContacts()
        NativeCO = self.pNative.GetCO(ContactMap = NativeContactMap)
        ClustHist = self.GetCO()
        return NativeCO, ClustHist

    def ContactDistCorr(self):
        # get AA contact distances
        pass
    
    def RamaChandran(self, doRamaProb = True):
        measure.NBlocks = 1
        picklename = FMT['RAMA'] % (self.Prefix, self.Temp)
        if os.path.isfile(picklename) and not self.Recompute: return
        NExcludedAtoms = len(self.pNative.ChainResNums)
        Phi = np.zeros([self.NFrames, self.p.NRes-NExcludedAtoms], np.float64)
        Psi = np.zeros([self.NFrames, self.p.NRes-NExcludedAtoms], np.float64)
        pb = sim.utility.ProgressBar(Text = self.ProgressMonitor('Phi and Psi angles'), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            Phi[i, :], Psi[i, :] = self.p.GetPhiPsi(BoxL = self.BoxL)
            pb.Update(i)    
        if doRamaProb:
            measure.Normalize = True
            measure.NBins = NBins
            measure.NBlocks = NBlocks
            hist = measure.makeHist2D(x = Phi.flatten(), y = Psi.flatten())
            ret = (Phi, Psi, hist)
        else: ret = (Phi, Psi)
        with open(picklename, 'w') as of: pickle.dump(ret, of)
        return
    
    def SpecificHeat(self):
        Cv_block = (1. / (kB * self.Temp**2.)) * np.ones(NBlocks)
        BlockSize = int(self.NFrames / NBlocks)
        Ene = np.zeros(self.NFrames)
        for i, frame in enumerate(self.FrameRange): Ene[i] = self.Ene[frame]
        for b in range(NBlocks):
            start = b * BlockSize
            stop = (b+1) * BlockSize if b < NBlocks-1 else self.NFrames
            E = Ene[start:stop]
            DeltaE = E - np.mean(E)
            Cv_block[b] *= np.mean( (DeltaE)**2. )
        err = 0.0
        if NBlocks > 1: err = np.std(Cv_block, ddof = 1)
        Cv = np.mean(Cv_block)
        return Cv, err
        

# class that calculates order params over replicas 
class Replica(object):
    ReInitWeights = False
    ReCompute = False
    
    def __init__(self, NativePdb, TrajPrefix, Prefix, TempSet = None, OrderParams = None, hasPseudoGLY = False, 
                 NStepsProd = None, NStepsSwap = None, WriteFreq = None):
        # order param calc functions
        self.MasterOrderParamDict = {'U'        : self.U,
                                     'Rg'       : self.Rg,
                                     'REE'      : self.REE,
                                     'RMSD'     : self.RMSD}
        # native (reference pdb)
        self.NativePdb = NativePdb
        
        # parse input and output locations (prefixes should contain dir names)
        self.TrajPrefix = os.path.abspath(TrajPrefix)
        self.Prefix = os.path.abspath(Prefix)
        self.TrajDir = os.path.dirname(self.TrajPrefix)
        self.OutDir = os.path.dirname(self.Prefix)
        
        # temp schedule
        self.TempFile = os.path.join(self.TrajDir, 'temps.txt')
        self.Temps = np.loadtxt(self.TempFile)
        if TempSet is None: TempSet = 300.0
        self.UpdateTemp(TempSet)
        
        # MD iterations required to re-order order params
        self.NStepsProd = NStepsProd
        self.NStepsSwap = NStepsSwap
        self.WriteFreq = WriteFreq
        
        # parse replica, traj and ene fn list
        # replica are continuous (unordered in temp. space) trajectories
        # whereas "traj" are discontinuous trajectories reordered by temperature
        self.LogFile = os.path.join(self.TrajDir, self.TrajPrefix + 'lammps.log')
        self.ReplicaFnList = [os.path.join(FMT['REPLICA'] % (self.TrajPrefix, i)) for i in range(len(self.Temps))]
        self.TrajFnList = [os.path.join(FMT['TRAJ'] % (self.TrajPrefix, t)) for t in self.Temps]
        self.EneFnList = [os.path.join(FMT['ENE'] % (self.TrajPrefix, t)) for t in self.Temps]
       
        # detect reordering
        self.DetectReorderMode()
       
        # order param storage location
        self.OrderParams = OrderParams if not OrderParams is None else self.MasterOrderParamDict.keys()
        if not 'U' in self.OrderParams: self.OrderParams += ['U']
        self.DataShelf = os.path.join(FMT['DATASHELF'] % self.Prefix)
        self.RawDataShelf = os.path.join(FMT['DATASHELF'] % (self.Prefix + '_raw'))
        
        # link a compute object
        self.Calc = Compute(NativePdb = self.NativePdb, Prefix = self.Prefix, hasPseudoGLY = hasPseudoGLY)
        
        # get all orderparams from all replicas
        if self.isReordered: self.GetAllData_Ordered()
        else: self.GetAllData_Unordered()
        
        # get all config weights
        self.GetConfigWeights()
            
    def UpdateTemp(self, TempSet = None):
        if TempSet is None: return
        self.TempSet = self.Temps[np.argmin(abs(self.Temps - TempSet))]
        return        
                
    def genShelfKey(self, paramname, temp):
        # generate key for data shelf
        return '%s_%3.2f' % (paramname, temp)
    
    def genRawShelfKey(self, paramname, replicaind):
        # return key for raw data shelf with data for each replica
        return '%s_%d' % (paramname, replicaind)
    
    def U(self, TrajFn, Temp = None, ReplicaInd = None):
        # when energy files are present
        if self.isReordered:
            # ene file parser (TrajFn is essentially an EneFn)
            # notational abuse for consistency
            Ene = np.loadtxt(TrajFn)
        else:
            # when ene files absent, extract from traj if simstyle log files
            if LogFileStyle == 'simstyle':
                TrajName = TrajFn
                LogFile = self.LogFile + '.%d' % ReplicaInd
                Token = '#run production'
                Trj = pickleTraj(TrajName, LogFile = LogFile, LogFileToken = Token) 
                Ene = Trj.ThermoDict['PEnergy']
            
            # when ene files are absent, parse directly from log-files if mystyle log files
            if LogFileStyle == 'mystyle':
                LogFile = self.LogFile + '.%d' % ReplicaInd
                with open(LogFile, 'r') as ologf: lines = ologf.readlines()
                Token = 0
                for line in lines:
                    Token += 1
                    if 'Minimization stats' in line: break
                lines = lines[Token:]
                start = 1 + [lines.index(line) for line in lines if 'PotEng' in line][0]
                stop = [lines.index(line) for line in lines if 'Loop time' in line][0]
                lines = lines[start : stop]
                Ene = [float(line.split()[2]) for line in lines]
                Ene = np.array(Ene)
        return Ene
        
    def RMSD(self, TrajFn, Temp):
        # RMSD calculator
        self.Calc.Update(TrajFn = TrajFn, Temp = Temp)
        return self.Calc.RMSD_frame()
   
    def Rg(self, TrajFn, Temp):
        # Rg calculator
        self.Calc.Update(TrajFn = TrajFn, Temp = Temp)
        return self.Calc.Rg_frame()
    
    def REE(self, TrajFn, Temp):
        # REE calculator
        self.Calc.Update(TrajFn = TrajFn, Temp = Temp)
        return self.Calc.REE_frame()
    
    def DetectReorderMode(self):
        # detects if temp. ordered trajectories are present
        x = all([os.path.isfile(i) for i in self.TrajFnList])
        y = all([os.path.isfile(i) for i in self.EneFnList])
        z = all([os.path.isfile(i) for i in self.ReplicaFnList])
        Test1 = x and y
        Test2 = z
        if (not Test1) and (not Test2):
            print 'ERROR: No trajectories here!'
            exit()
        elif Test1:
            print 'Detected discontinous trajectories, but ordered by temperature'
            self.isReordered = True
            self.NFrames = len(np.loadtxt(self.EneFnList[0]))
        else:
            print 'Detected continuous trajectories, but not ordered by temperature'
            self.isReordered = False
            self.NFrames = int(self.NStepsProd / self.WriteFreq)
        return
        
    def Reorder(self, OrderParam, Temp):
        # reorders a given order parameter calculated from different replica
        # at a particular temp (present in the supplied temp. schedule)
        # assumes that all replicas have already been pickled
        if self.NStepsProd is None or self.NStepsSwap is None or self.WriteFreq is None:
            print 'ERROR: Reordering order parameters requires md iteration info'
            return
        TempInd = list(self.Temps).index(Temp)
        RepIndsMaster = [np.where(x[1:] == TempInd)[0][0] for x in np.loadtxt(self.LogFile, skiprows = 3)]
        RepInds = RepIndsMaster[ - int(self.NStepsProd / self.NStepsSwap) : ]
        # data type check
        if isinstance(OrderParam, list): OrderParam = np.array(OrderParam)
        # keep only production frames for this order parameter
        O1 = OrderParam[:, -int(self.NStepsProd / self.WriteFreq) : ]
        O2 = []
        # case 1
        if self.WriteFreq <= self.NStepsSwap:
            for ii, i in enumerate(RepInds):
                # point to the row for this replica
                thisO = O1[i, :]
                # find start and stop indices
                start = ii * self.NStepsSwap / self.WriteFreq
                stop = (ii+1) * self.NStepsSwap / self.WriteFreq
                # write to output
                O2.extend(thisO[start:stop])
        # case 2
        else:
            NSkip = self.WriteFreq / self.NStepsSwap
            for ii, i in enumerate(RepInds[0::NSkip]):
                # point to the row for this replica
                thisO = O1[i, :]
                # find start and stop indices
                start = ii
                stop = ii + 1
                O2.extend(thisO[start:stop])
        # return array
        O2 = np.array(O2)
        return O2
    
    
    def GetAllData_Ordered(self):
        # when trajectories are discontinuous i.e. ordered by temp.
        d = shelve.open(self.DataShelf)
        print '\nAssembling order params at each temperature...'
        for o in self.OrderParams:
            print '\nOrderParam = %s' % o
            for i, t in enumerate(self.Temps):
                Fn = self.EneFnList[i] if o == 'U' else self.TrajFnList[i]
                if not os.path.isfile(Fn):
                    print 'Skipping temp %3.2f K: Traj not found' % t
                    continue
                # skip if already computed
                key = self.genShelfKey(o, t)
                if d.has_key(key): continue
                print '  %s, %3.2f K\n' % (o, t)
                # if not computed, compute and add to shelf
                func = self.MasterOrderParamDict[o]
                d[key] = func(TrajFn = Fn, Temp = t)
                sys.stdout.write(ERASE_LINE)
                sys.stdout.write(UP_ONE_LINE)
            print '--------------------------------'
        d.close()
        return
        
    def GetAllData_Unordered(self):
        # when trajectories are continuous i.e. unordered in temp. space
        d = shelve.open(self.DataShelf)
        dr = shelve.open(self.RawDataShelf)
        print '\nAssembling order params from each replica...'
        for o in self.OrderParams:
            print '\nOrderParam = %s' % o
            # calculate the order param from each replica and accumulate
            thisO_raw = []
            for i in range(len(self.Temps)):
                Fn = self.ReplicaFnList[i]
                if not os.path.isfile(Fn):
                    print 'Skipping replica %d: Traj not found' % i
                    continue
                key = self.genRawShelfKey(o, i)
                # compute only if not present in shelf
                if not dr.has_key(key):
                    print '  %s, Replica %d\n' % (o, i)
                    func = self.MasterOrderParamDict[o]
                    if o == 'U': ret = func(TrajFn = Fn, Temp = None, ReplicaInd = i)
                    else: ret = func(TrajFn = Fn, Temp = None)
                    dr[key] = ret
                    sys.stdout.write(ERASE_LINE)
                    sys.stdout.write(UP_ONE_LINE)
                # if present access from shelf
                else:
                    ret = dr[key]
                # add to list (append not extend)
                thisO_raw.append(list(ret))
            # convert list of lists to 2D array
            thisO_raw = np.array(thisO_raw)
            # reorder by temperature
            print 'Reordering by temperature...'
            for i, t in enumerate(self.Temps):
                key = self.genShelfKey(o, t)
                if d.has_key(key): continue
                print '  %s, %3.2f K\n' % (o, t)
                thisO = self.Reorder(thisO_raw, t)
                sys.stdout.write(ERASE_LINE)
                sys.stdout.write(UP_ONE_LINE)
                d[key] = thisO
            print '--------------------------------'
        dr.close()
        d.close()
        return
    
    
    def GetConfigWeights(self):
        # calculate config weights
        # must be called only after the data shelf has been populated completely
        d = shelve.open(self.DataShelf)
        if not self.ReInitWeights and d.has_key('w_kn') and d.has_key('logw_kn'):
            d.close()
            return
        
        beta_k = 1.0 / (kB * self.Temps)
        K = len(self.Temps)
        N = self.NFrames
        
        # extract all energies
        U_kn = np.zeros([K, N], np.float64)
        for k, t in enumerate(self.Temps):
            key = self.genShelfKey('U', t)
            U_kn[k,:] = d[key]
        # loop over blocks
        w_kn = {}
        logw_kn = {}
        BlockSize = int(self.NFrames / NBlocks)
        print 'Calculating config weights...'
        for b in range(NBlocks):
            if NBlocks > 1: print 'Block: ', b
            start = b * BlockSize
            stop = (b+1) * BlockSize if not b == NBlocks - 1 else N
            # get reduced cross energies at all temps
            this_U_kn = U_kn[:, start:stop]
            this_N = this_U_kn.shape[1]
            this_N_k = this_N * np.ones(K, np.uint8) # data-type must be int to avoid error from mbar
            this_u_kln = np.zeros([K,K,this_N], np.float64)
            for k in range(K):
                for l in range(K):
                    this_u_kln[k, l, 0:this_N_k[k]] = beta_k[l] * this_U_kn[k, 0:this_N_k[k]]
            # initialize mbar object
            mbar = pymbar.mbar.MBAR(this_u_kln, this_N_k, verbose = False)
            # get free energies and shift by value at lowest temp
            this_f_k = mbar.f_k
            #this_f_k = whamlib.free_energy(ekn = this_U_kn, betak = beta_k, nbin = 100, niterbin = 100, niterall = 10) 
            # get log weights at all temps for this block        
            for k, t in enumerate(self.Temps):
                log_w_kn = whamlib.log_weight(ekn = this_U_kn, betak = beta_k, targetbeta = beta_k[k], fk = this_f_k)
                logw_kn[ (k,b) ] = log_w_kn
                w_kn[ (k,b) ] = np.exp(log_w_kn)
        # dump to shelf
        d['w_kn'] = w_kn
        d['logw_kn'] = logw_kn
        d.close()
        return 
        
    
    def FoldCurve(self, O = 'RMSD', MAX = None, MIN = None):
        # calculate folding curve w.r.t. chosen order param (usually rmsd)
        picklename = FMT['FOLDCURVE'] % (self.Prefix, O)
        if os.path.isfile(picklename): return
        # see if config weights need to be recalculated
        if self.ReInitWeights: self.GetConfigWeights()
        # retrieve all data for orderparam
        d = shelve.open(self.DataShelf)
        x_kn = np.zeros([len(self.Temps), self.NFrames], np.float64)
        for k, t in enumerate(self.Temps):
            key = self.genShelfKey(O, t)
            x_kn[k, :] = d[key]
        # get overall max and min and bin w.r.t them
        x_min = MIN
        x_max = MAX
        if x_min is None:
            x_min = x_kn.min() * (1. + measure.HistPadding)
        if x_max is None:
            x_max = x_kn.max() * (1. - measure.HistPadding)
        dx = (x_max - x_min) / float(NBins)
        x_centers =  x_min + dx * (0.5 + np.arange(NBins))
        cut_inds = (x_centers <= OCut[O])
        # computing folding fraction block by block
        BlockSize = int(self.NFrames / NBlocks)
        foldfrac_block = np.zeros([len(self.Temps), NBlocks])
        for k, t in enumerate(self.Temps):
            print '\nTarget Temp = %3.2f K' % t
            for b in range(NBlocks):
                if NBlocks > 1: print ' Block: ', b
                start = b * BlockSize
                stop = (b+1) * BlockSize if not b == NBlocks - 1 else self.NFrames
                # get config weights
                weights = d['w_kn'][ (k,b) ].flatten()
                # calculate 1D histograms (don't normalize)
                x = x_kn[:, start:stop].flatten()
                this_hist = np.zeros(NBins)
                for i, this_x in enumerate(x):
                    idx = int( (this_x - x_min) / dx)
                    if idx >= NBins: idx = NBins - 1
                    if idx < 0: idx = 0
                    this_hist[idx] += weights[i]
                #measure.NBins = NBins
                #measure.NBlocks = 1
                #measure.Normalize = False
                #this_bin_centers, this_hist, this_err = measure.makeHist(x, weights = weights, bintuple = (x_min, x_max, dx))
                
                # computing folding fraction at (temp, block) as cumulative histogram within a cutoff
                foldfrac_block[k, b] = this_hist[cut_inds].sum() / float(this_hist.sum())

        foldfrac = np.mean(foldfrac_block, axis = 1)
        if NBlocks > 1: err = np.std(foldfrac_block, axis = 1, ddof = 1)
        else: err = np.zeros(len(self.Temps))
        # write to pickle
        ret = (self.Temps, foldfrac, err)
        with open(picklename, 'w') as of: pickle.dump(ret, of)
        d.close()
        return
        
    
    def PMF(self, O):
        # calculates pmf w.r.t order param O
        picklename = FMT['PMF1D'] % (self.Prefix, self.TempSet, O)
        if os.path.isfile(picklename): return
        # see if config weights need to be recalculated
        if self.ReInitWeights: self.GetConfigWeights()
        # extract order parames and energies
        d = shelve.open(self.DataShelf)
        x_kn = np.zeros([len(self.Temps), self.NFrames], np.float64)
        for k, t in enumerate(self.Temps):
            key = self.genShelfKey(O, t)
            x_kn[k, :] = d[key]
        # get overall max and min and bin w.r.t them
        x_min = x_kn.min() * (1.0 - measure.HistPadding)
        x_max = x_kn.max() * (1.0 + measure.HistPadding)
        dx = (x_max - x_min) / float(NBins)
    
        # compute PMF block by block
        BlockSize = int(self.NFrames / NBlocks)
        pmf_block = np.zeros([NBlocks, NBins])
        TempInd = np.argmin(abs(self.Temps - self.TempSet))
        for b in range(NBlocks):
            if NBlocks > 1: print 'Block: ', b
            start = b * BlockSize
            stop = (b+1) * BlockSize if not b == NBlocks - 1 else self.NFrames
            x = x_kn[:, start:stop].flatten()
            weights = d['w_kn'][ (TempInd, b) ].flatten()
            measure.NBins = NBins
            measure.NBlocks = 1
            bintuple = (x_min, x_max, dx)
            bin_centers, this_hist, this_err = measure.makeHist(x, weights = weights, bintuple = bintuple)
            pmf_block[b, :] = - (kB * self.TempSet) * np.log(this_hist)
        
        # trim the pmf
        pmf = np.mean(pmf_block, axis = 0)
        if NBlocks > 1: err = np.std(pmf_block, axis = 0, ddof = 1)
        else: err = np.zeros(NBins)
        # write to pickle
        ret = (bin_centers, pmf, err)
        with open(picklename, 'w') as of: pickle.dump(ret, of)
        d.close()
        return

    
    def PMF2D(self, O1, O2, MIN1 = None, MIN2 = None, MAX1 = None, MAX2 = None):
        # calculate 2D pmf w.r.t order params O1 and O2
        picklename = FMT['PMF2D'] % (self.Prefix, self.TempSet, O1, O2)
        if os.path.isfile(picklename): return
        # see if config weights need to be recalculated
        if self.ReInitWeights: self.GetConfigWeights()
        # extract order parames and energies
        d = shelve.open(self.DataShelf)
        x_kn = np.zeros([len(self.Temps), self.NFrames], np.float64)
        y_kn = np.zeros([len(self.Temps), self.NFrames], np.float64)
        for k, t in enumerate(self.Temps):
            xkey = self.genShelfKey(O1, t)
            ykey = self.genShelfKey(O2, t)
            x_kn[k, :] = d[xkey]
            y_kn[k, :] = d[ykey]
        # get overall max and min and bin w.r.t them
        x_min = MIN1
        x_max = MAX1
        y_min = MIN2
        y_max = MAX2
        if x_min is None:
            x_min = x_kn.min() * (1.0 - measure.HistPadding)
        if x_max is None:
            x_max = x_kn.max() * (1.0 + measure.HistPadding)
        dx = (x_max - x_min) / float(NBins)
        if y_min is None:
            y_min = y_kn.min() * (1.0 - measure.HistPadding)
        if y_max is None:
            y_max = y_kn.max() * (1.0 + measure.HistPadding)
        dy = (y_max - y_min) / float(NBins)
        # bin centers
        x_centers =  x_min + dx * (0.5 + np.arange(NBins))
        y_centers =  y_min + dy * (0.5 + np.arange(NBins))
    
        # compute PMF block by block
        BlockSize = int(self.NFrames / NBlocks)
        pmf_block = np.zeros([NBlocks, NBins, NBins])
        TempInd = np.argmin(abs(self.Temps - self.TempSet))
        for b in range(NBlocks):
            if NBlocks > 1: print 'Block: ', b
            start = b * BlockSize
            stop = (b+1) * BlockSize if not b == NBlocks - 1 else self.NFrames
            x = x_kn[:, start:stop].flatten()
            y = y_kn[:, start:stop].flatten()
            weights = d['w_kn'][ (TempInd, b) ].flatten()
            # calculate histogram
            this_hist = np.zeros([NBins, NBins])
            for i in range(len(x)):
                idx = int( (x[i] - x_min) / dx)
                idy = int( (y[i] - y_min) / dy)
                if idx >= NBins: idx = NBins - 1
                if idx < 0: idx = 0
                if idy >= NBins: idy = NBins - 1
                if idy < 0: idy = 0
                this_hist[idx, idy] += weights[i]    
            # normalize histogram
            this_rownorm = dy * np.sum(this_hist, axis = 1)
            this_norm = dx * np.sum(this_rownorm)
            this_hist /= this_norm
            pmf_block[b, :, :] = -np.log(this_hist)
            #measure.NBins = NBins
            #measure.NBlocks = 1
            #bintuple = ( (x_min, y_min), (x_max, y_max), (dx, dy) )
            #bin_centers, this_hist, this_err = measure.makeHist2D(x, y, weights = weights, bintuple = bintuple)
            #pmf_block[b, :, :] =  - (kB * self.TempSet) * np.log(this_hist)
             
             
        # calculate errors
        pmf = np.mean(pmf_block, axis = 0)
        if NBlocks > 1: err = np.std(pmf_block, axis = 0, ddof = 1)
        else: err = np.zeros([NBins, NBins])
        # write to pickle
        ret = ( (x_centers, y_centers), pmf, err)
        with open(picklename, 'w') as of: pickle.dump(ret, of)
        d.close()
        return

       
