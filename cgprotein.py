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
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# file formats
FMT = utils.FMT

# CG O-N "peptide" bond length
# (all forcefields roughly predict this value)
PeptideBondLen = 1.75 #A

kB = 0.001987
Emax = 80* 0.6 # 20 kT at T = 300 K 

OCut = {'RMSD': 3.0, 'Rg': 3.0, 'REE': 3.0}
MaxCluster = 10

StepFreq = 1
NBlocks = 4
NBins = 50

# Misc utility functions
def TrimPMF(pmf, LastNBins = 10, Dim = 2):
    # trims a pmf to remove noise
    global Emax
    LastNBins = 10
    # hard cut for all finite but very high and potentially noisy) points on the free energy surface
    if Dim == 1:
        for i in range(len(pmf)):
            if np.isfinite(pmf[i]) and pmf[i] > Emax:
                pmf[i] = Emax
    if Dim == 2:
        for i in range(pmf.shape[0]):
            for j in range(pmf.shape[1]):
                if np.isfinite(pmf[i,j]) and pmf[i,j] > Emax:
                    pmf[i,j] = Emax     
    # shift the pmf so that the average of the smallest 
    # LastNBins is zero
    #pmf_ = pmf[np.isfinite(pmf)]
    #if Dim == 2: pmf_ = pmf_.flatten()
    #offset = np.mean(np.sort(pmf_)[0:LastNBins])
    #pmf -= offset
    return pmf


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
        if all(BoxL) :BBPos = lib.reimagechain(pos = self.Pos, boxl = BoxL)
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
    def __init__(self, NativePdb, TrajFn = None, Temp = None, Prefix = 'compute', hasPseudoGLY = False):
        # parse out-prefix
        self.Prefix = os.path.abspath(Prefix)
        self.OutPrefix = self.Prefix.split('/')[-1]
        self.OutDir = os.path.dirname(self.Prefix)
        # parse traj
        if not TrajFn is None:
            self.TrajFn = os.path.abspath(TrajFn)
            self.Trj = pickleTraj(self.TrajFn)
            self.FrameRange = range(0, len(self.Trj), StepFreq)
            self.NFrames = len(self.FrameRange)
            self.BoxL = self.Trj.FrameData.get('BoxL', np.zeros(3))
        # record temp
        if not Temp is None:
            self.Temp = Temp
        # cg protein objects for native and predicted
        self.pNative = ProteinNCOS(NativePdb, hasPseudoGLY = hasPseudoGLY)
        self.p = ProteinNCOS(NativePdb, hasPseudoGLY = hasPseudoGLY)
   
    def Update(self, TrajFn = None, Temp = None):
        # updates the compute object if traj or temp is changed
        # if native pdb is changed then better create another instance
        if not Temp is None: self.Temp = Temp
        if not TrajFn is None:
            self.TrajFn = os.path.abspath(TrajFn)
            self.Trj = pickleTraj(self.TrajFn)
            self.FrameRange = range(0, len(self.Trj), StepFreq)
            self.NFrames = len(self.FrameRange)
            self.BoxL = self.Trj.FrameData.get('BoxL', np.zeros(3))
    
    def Rg_frame(self):
        # per frame Rg calculator
        rg_frame = np.zeros(self.NFrames)
        # loop over frames
        pb = sim.utility.ProgressBar('Calculating overall Rg for %s at %3.2f K...' % (self.OutPrefix, self.Temp), Steps = self.NFrames)
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
        pb = sim.utility.ProgressBar('Calculating overall REE for %s at %3.2f K...' % (self.OutPrefix, self.Temp), Steps = self.NFrames)
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
        pb = sim.utility.ProgressBar('Calculating overall RMSD for %s at %3.2f K...' % (self.OutPrefix, self.Temp), Steps = self.NFrames)
        for i, frame in enumerate(self.FrameRange):
            Pos = self.Trj[frame]
            self.p.UpdatePos(Pos)
            rmsd_frame[i] = self.p.QuickRMSD(self.pNative)
            pb.Update(i)
        return rmsd_frame
    
    def ResContacts_frame(self):
        # makes sense to store the entire per frame data since
        # multiple computations can be done on this entire data set
        # and those calculations need not be saved to a pickle
        picklename = FMT['RESCONTACTS'] % (self.Prefix, self.Temp)
        if os.path.isfile(picklename) and not self.Recompute: return
        # loop over frames to get average CG contact distances   
        ContactMap = np.zeros([self.NFrames, self.p.NRes, self.p.NRes], np.float64)
        ContactDist = np.zeros([self.NFrames, self.p.NRes, self.p.NRes], np.float64)
        pb = sim.utility.ProgressBar(Text = 'Enumerating contacts at %3.2fK...' % self.Temp, Steps = self.NFrames)
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
        pb = sim.utility.ProgressBar('Calculating overall RMSD for %s at %3.2f K...' % (self.OutPrefix, self.Temp), Steps = self.NFrames)
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
        pb = sim.utility.ProgressBar('Calculating selected residue RMSD for %s at %3.2f K...' % (self.OutPrefix, self.Temp), Steps = self.NFrames)
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
        pb = sim.utility.ProgressBar('Calculating detailed RMSD for %s at %3.2f K...' % (self.OutPrefix, self.Temp), Steps = self.NFrames)
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
        pb = sim.utility.ProgressBar('Calculating Phi, Psi at %3.2f K...' % self.Temp, Steps = self.NFrames)
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
        
        

# class that calculates order params over replicas 
class Replica(object):
    ReInitWeights = False
    ReCompute = False
    
    def __init__(self, NativePdb, TrajPrefix, Prefix, TempSet = None, OrderParams = None, hasPseudoGLY = False):
        # order param calc functions
        self.MasterOrderParamDict = {'U': self.U,
                                     'Rg': self.Rg,
                                     'REE': self.REE,
                                     'RMSD': self.RMSD}
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
        T = TempSet if not TempSet is None else 300.0
        self.TempSet = self.Temps[np.argmin(abs(self.Temps - T))]
        # parse traj and ene fn list
        self.TrajFnList = [os.path.join(FMT['TRAJ'] % (self.TrajPrefix, t)) for t in self.Temps]
        self.EneFnList = [os.path.join(FMT['ENE'] % (self.TrajPrefix, t)) for t in self.Temps]
        # number of frames in each replica (assuming all replica have same number of frames)
        self.NFrames = len(np.loadtxt(self.EneFnList[0])[0::StepFreq])
        # order param storage location
        self.OrderParams = OrderParams if not OrderParams is None else self.MasterOrderParamDict.keys()
        if not self.OrderParams.__contains__('U'): self.OrderParams += ['U']
        self.DataShelf = os.path.join(FMT['DATASHELF'] % self.Prefix)
        # link a compute object
        self.Calc = Compute(NativePdb = self.NativePdb, Prefix = self.Prefix, hasPseudoGLY = hasPseudoGLY)
        # get all orderparams from all replicas
        self.GetAllData()
        # get all config weights
        self.GetConfigWeights()
       
    def genShelfKey(self, paramname, temp):
        # generate key for data shelf
        return '%s_%3.2f' % (paramname, temp)
    
    def U(self, TrajFn, Temp = None):
        # ene file parser (TrajFn is essentially an EneFn)
        # notational abuse for consistency
        Ene = np.loadtxt(TrajFn)
        return Ene[0::StepFreq]
        
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
    
    def GetAllData(self):
        # calculate all order params at all temperatures
        print 'Assembling order params for each Replica...'
        d = shelve.open(self.DataShelf)
        for o in self.OrderParams:
            for i, t in enumerate(self.Temps):
                Fn = self.EneFnList[i] if o == 'U' else self.TrajFnList[i]
                if not os.path.isfile(Fn):
                    print 'Replica at temp %3.2f K not found' % t
                    continue
                key = self.genShelfKey(o, t)
                if d.has_key(key): continue
                print '\n%s, %3.2f K' % (o, t)
                func = self.MasterOrderParamDict[o]
                d[key] = func(TrajFn = Fn, Temp = t)
        d.close()
        return 
    
    
    def GetConfigWeights(self):
        # calculate config weights
        # must be called only after the data shelf has been populated completely
        d = shelve.open(self.DataShelf)
        if not self.ReInitWeights and d.has_key('w_kn'):
            d.close()
            return
        
        import pymbar, whamlib
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
                    this_u_kln[k, l, :] = beta_k[l] * this_U_kn[k, :]
            # initialize mbar object
            mbar = pymbar.mbar.MBAR(this_u_kln, this_N_k, verbose = False)
            # get free energies and shift by value at lowest temp
            this_f_k = mbar.f_k
            this_f_k -= this_f_k.min()
            # get log weights at all temps for this block        
            for k, t in enumerate(self.Temps):
                log_w_kn = whamlib.log_weight(ekn = this_U_kn, betak = beta_k, targetbeta = beta_k[k], fk = this_f_k)
                w_kn[ (k, b) ] = np.exp(log_w_kn)
        # dump to shelf
        d['w_kn'] = w_kn
        d.close()
        return 
        
    
    def FoldCurve(self, O = 'RMSD'):
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
        x_max = x_kn.max() * (1. - measure.HistPadding)
        x_min = x_kn.min() * (1. + measure.HistPadding)
        dx = (x_max - x_min) / float(NBins)
        x_centers =  x_min + dx * (0.5 + np.arange(NBins))
        cut_inds = (x_centers <= OCut[O])
        # computing folding fraction block by block
        BlockSize = int(self.NFrames / NBlocks)
        foldfrac_block = np.zeros([len(self.Temps), NBlocks])
        for k, t in enumerate(self.Temps):
            print 'Target Temp = %3.2f K' % t
    	    for b in range(NBlocks):
    	        if NBlocks > 1: print ' Block: ', b
    	        start = b * BlockSize
    	        stop = (b+1) * BlockSize if not b == NBlocks - 1 else self.NFrames
    	        # get config weights
    	        weights = d['w_kn'][ (k,b) ].flatten()
    	        print 'Calculating histograms...'
    	        x = x_kn[:, start:stop].flatten()
    	        measure.NBins = NBins
    	        measure.NBlocks = 1
    	        this_bin_centers, this_hist, this_err = measure.makeHist(x, weights = weights, bintuple = (x_min, x_max, dx))
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
        for b in range(NBlocks):
            if NBlocks > 1: print 'Block: ', b
            start = b * BlockSize
            stop = (b+1) * BlockSize if not b == NBlocks - 1 else self.NFrames
            x = x_kn[:, start:stop].flatten()
            weights = d['w_kn'][ (k, b) ].flatten()
            measure.NBins = NBins
            measure.NBlocks = 1
            bintuple = (x_min, x_max, dx)
            bin_centers, this_hist, this_err = measure.makeHist(x, weights = weights, bintuple = bintuple)
            pmf_block[b, :] = - (kB * self.TempSet) * np.log(this_hist)
        
        # trim the pmf
        pmf = np.mean(pmf_block, axis = 0)
        pmf = TrimPMF(pmf, Dim = 1)
        if NBlocks > 1: err = np.std(pmf_block, axis = 0, ddof = 1)
        else: err = np.zeros(NBins)
        # write to pickle
        ret = (bin_centers, pmf, err)
        with open(picklename, 'w') as of: pickle.dump(ret, of)
        d.close()
        return

    
    def PMF2D(self, O1, O2):
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
        x_min = x_kn.min() * (1.0 - measure.HistPadding)
        x_max = x_kn.max() * (1.0 + measure.HistPadding)
        dx = (x_max - x_min) / float(NBins)
        y_min = y_kn.min() * (1.0 - measure.HistPadding)
        y_max = y_kn.max() * (1.0 + measure.HistPadding)
        dy = (y_max - y_min) / float(NBins)
    
        # compute PMF block by block
        BlockSize = int(self.NFrames / NBlocks)
        pmf_block = np.zeros([NBlocks, NBins, NBins])
        for b in range(NBlocks):
            if NBlocks > 1: print 'Block: ', b
            start = b * BlockSize
            stop = (b+1) * BlockSize if not b == NBlocks - 1 else self.NFrames
            x = x_kn[:, start:stop].flatten()
            y = y_kn[:, start:stop].flatten()
            weights = d['w_kn'][ (k, b) ].flatten()
            measure.NBins = NBins
            measure.NBlocks = 1
            bintuple = ( (x_min, y_min), (x_max, y_max), (dx, dy) )
            bin_centers, this_hist, this_err = measure.makeHist2D(x, y, weights = weights, bintuple = bintuple)
            pmf_block[b, :, :] = - (kB * self.TempSet) * np.log(this_hist)
            
        # trim the pmf
        x_centers, y_centers = bin_centers
        pmf = np.mean(pmf_block, axis = 0)
        pmf = TrimPMF(pmf, Dim = 2)
        if NBlocks > 1: err = np.std(pmf_block, axis = 0, ddof = 1)
        else: err = np.zeros([NBins, NBins])
        # write to pickle
        ret = ( (x_centers, y_centers), pmf, err)
        with open(picklename, 'w') as of: pickle.dump(ret, of)
        d.close()
        return

       
