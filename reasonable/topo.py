#!/usr/bin/env python

''' Builds the toplogy of the sim-style Sys object '''

import os, numpy as np
from const import *
import sim
import protein, mdsim # /share/apps/scripts/mdsim

Verbose = True
protein.MinCODflt = MinCO
protein.ResRadiusDflt = ResRadius
mdsim.RunDflts['MODEL'] = 'ff96glghs'

class ProteinNCOS(object):
    ''' A rudimentary object to create and store the NCOS toplogy 
        The atom ordering within a residue is N C O S
        a global config object must be supplied that contains
        sim-style AtomType objects
    '''    
    def __init__(self, cfg, Pdb = None, Seq = None, Prefix = 'ncos'):
        # set Prefix
        self.Prefix = Prefix
        # unpack sidechain atomtypes
        self.AtomS = cfg.AtomS
        # has special torsion potentials
        self.hasSpecialBBGLYTorsions = cfg.hasSpecialBBGLYTorsions
        self.hasSpecialBBPROTorsions = cfg.hasSpecialBBPROTorsions
        # extract the entire config object just in case
        self.cfg = cfg
        # sequence book-keeping
        if Seq is None:
            # internal proteinclass object
            self.p0 = protein.ProteinClass(Pdb)
            self.Seq = self.p0.Seq
            self.Pos = self.p0.Pos
            # get PdbName
            self.PdbName = Pdb.split('/')[-1].split('.pdb')[0]
        else:
            self.Seq = Seq
            self.PdbName = None
            self.Pos = None
        self.ResTypes = list(set(self.Seq))
        self.NRes = len(self.Seq)
        self.NResTypes = len(self.ResTypes)
        self.StartAtomInds = []
        self.BondPairs = []
        self.__SetStartAtomInds()
        self.__SetBonds()
   
    def __repr__(self):
        s = ' '.join(self.Seq)
        return s
        
    def __SetStartAtomInds(self):
        ''' set the index of first atoms for each residue'''
        x = 0
        for i, r in enumerate(self.Seq):
            self.StartAtomInds.append(x)
            if self.AtomS[r] is None: x += 3
            else: x += 4
        
    def __SetBonds(self):
        ''' create a simple bond pair list for sim'''
        for i, r in enumerate(self.Seq):
            AtomNInd = self.StartAtomInds[i]
            AtomCInd = AtomNInd + 1
            AtomOInd = AtomNInd + 2
            # intra-residue backbone bonds
            BondPairs = [(AtomNInd, AtomCInd), (AtomCInd, AtomOInd)]
            # backbone-sidechain bond
            if not self.AtomS[r] is None:
                # sidechain bonds only for residues whose sidechains are not None
                AtomSInd = AtomNInd + 3
                BondPairs.append( (AtomCInd, AtomSInd) )
            # peptide bond
            if i < self.NRes - 1:
                NextAtomNInd = self.StartAtomInds[i+1]
                BondPairs.append( (AtomOInd, NextAtomNInd) )
            # add to master list
            self.BondPairs.extend(BondPairs)
    
    def GetBBInds(self):
	''' extract the backbone indices of the structure'''
	BBInds = []
	for i, r in enumerate(self.Seq):
		k = self.StartAtomInds[i]
		BBInds.extend( [k, k+1, k+2] )
	return BBInds
   
    def GetSInds(self):
        ''' extract the sidechain indices of the structure
            returns alpha carbon index for residues with None sidechain'''
        SInds = []
        for i, r in enumerate(self.Seq):
            idx = 1 if self.AtomS[r] is None else 3
            SInds.append(self.StartAtomInds[i] + idx)
        return SInds

    def GetResPos(self):
        ''' get the COM of the residues'''
        ResPos = np.zeros([self.NRes, 3])
        for i in range(self.NRes):
            start = self.StartAtomInds[i]
            stop = self.StartAtomInds[i+1] if i < self.NRes-1 else len(self.Pos)
            ResPos[i, :] = np.mean(self.Pos[start:stop, :], axis = 0)
        return ResPos

    def GetResContactList(self):
        ''' find native contacts with ResRadius'''
        ResContactList = []
        ResPos = self.GetResPos()
        for i in range(self.NRes - 1):
            for j in range(i+1, self.NRes):
                # ignore adjacent residues
                if abs(i-j) < MinCO: continue
                # calculate distance between residues
                d_ij = ResPos[j] - ResPos[i]
                # check if contact
                dsq = np.sum(d_ij * d_ij)
                if dsq <= ResRadius * ResRadius: ResContactList.append( (i,j) )
        return ResContactList

    def Map2Polymer(self, PolyName, AAPdb = None, EneMin = False, DelTempPdb = True):
        ''' maps the given Pdb to a polymer of equivalent length 
        and returns a coarse grained protein obj for the mapped polymer
        Energy minimization not yet implemented'''
        if Verbose: print 'Mapping native structure to %s' % PolyName
        # read in unmapped Pdb
        if AAPdb is None: AAPdb = os.path.join(NATIVEPATH['Unmapped'], self.PdbName + '.pdb')
        # check if unmapped Pdb has multiple frames
        with open(AAPdb, 'r') as of: head = of.readlines(10)
        head = ''.join(head)
        if head.__contains__('Model'):  p_AA = protein.ProteinClass(AAPdb, Model = 1)
        else: p_AA = protein.ProteinClass(AAPdb)
        # map the polymer to this sequence
        p_AA = p_AA.Decap()
        NewSeq = [PolyName] * len(p_AA.Seq)
        p_AA = p_AA.MutateSeq(NewSeq)
        p_AA = p_AA.Cap() # add ACE and NME caps, else Amber starts yelling
        NewAAPdb = 'tmpAA.pdb'
        p_AA.WritePdb(NewAAPdb)
        # energy minimize this pdb
        RunPath = os.path.join(os.getcwd(), 'AmberEneMin')
        if EneMin:
            if Verbose: print 'Energy minimizing the mapped AA polymer...'
            mdsim.RunDflts['MODEL']= 'ff96glghs'
            x = mdsim.SimClass(RunPath)
            x.SysInitPdb(NewAAPdb)
            x.SysBuild()
            x.RunMin(500, 500) # 1000 total min steps
            os.remove(NewAAPdb)
            NewAAPdb = os.path.join(RunPath, 'current.pdb')
        # coarse grain this pdb
        NewCGPdb = 'tmpCG.pdb'
        cmd = 'python %s %s %s' % (MAPSCRIPT, NewAAPdb, NewCGPdb.split('.pdb')[0])
        os.system(cmd)
        # read in the cg ProteinNCOS object
        p_New = self.__class__(Pdb = NewCGPdb, cfg = self.cfg)
        # delete temp files
        for i in [NewAAPdb, NewCGPdb]: os.remove(i)
        if os.path.isdir(RunPath): os.system('rm -r ' + RunPath)
        return p_New

    

def MakeSys(p, cfg = None, NMols = 1):
    '''generate the topology for the sim Sys object
       currently implemented for single chains
       p is ProteinNCOS type object
       cfg is the global config object '''
    if Verbose == True:
        print 'Generating System topology...'
    # unpack sidechain atoms
    AtomS = p.AtomS
    # generate the molecule
    AtomList = []
    for i,r in enumerate(p.Seq):
        # backbone atoms
        if p.hasSpecialBBGLYTorsions:
            if r == 'GLY': res = [AtomN, AtomC_GLY, AtomO]
        elif p.hasSpecialBBPROTorsions:
            if r == 'PRO': res = [AtomN, AtomC_PRO, AtomO]
        else: res = [AtomN, AtomC, AtomO]
        # sidechain
        if not AtomS[r] is None: res.append(AtomS[r])
        AtomList.extend(res)
    Mol = sim.chem.MolType(p.Prefix, [i for i in AtomList])
    # generate bonds
    for (i,j) in p.BondPairs: Mol.Bond(i,j)
    # create System
    World = sim.chem.World([Mol], Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(World, Name = p.Prefix)
    #for i in range(NChains): Sys += Mol.New()
    Sys += Mol.New()
    return Sys



