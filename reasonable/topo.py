#!/usr/bin/env python

''' Builds the toplogy of the sim-style Sys object '''

import os, numpy as np, copy, string, cPickle as pickle
from const import *
import sim
import protein

Verbose = True
protein.MinCODflt = MinCO
protein.ResRadiusDflt = ResRadius

class ProteinNCOS(object):
    ''' A rudimentary object to create and store the NCOS toplogy 
        The atom ordering within a residue is N C O S
        a global config object must be supplied that contains
        sim-style AtomType objects
    '''    
    def __init__(self, cfg, Pdb = None, Seq = None, Model = None, Prefix = 'ncos'):
        if Verbose: print 'Initializing a NCOS protein object...'
        # set Prefix
        self.Prefix = Prefix
        # has special GLY Params
        self.hasSpecialBBGLYAngles = cfg.hasSpecialBBGLYAngles
        self.hasSpecialBBGLYTorsions = cfg.hasSpecialBBGLYTorsions
        self.hasSpecialGLYParams = cfg.hasSpecialGLYParams()
        # has special PRO params
        self.hasSpecialBBPROAngles = cfg.hasSpecialBBPROAngles
        self.hasSpecialBBPROTorsions = cfg.hasSpecialBBPROTorsions
        self.hasSpecialPROParams = cfg.hasSpecialPROParams()
        # has special sidechains for glycine
        self.hasPseudoGLY = cfg.hasPseudoGLY()
        if self.hasSpecialGLYParams and self.hasPseudoGLY:
            print 'Error: Cannot have pseudo GLY side chain and special GLY BB torsion simultaneously'
            exit()
        # sidechain referencing
        self.SSRefType = cfg.SSRefType
        # extract the entire config object just in case
        self.cfg = cfg
        # if a CG Pdb is provided instead of a Seq
        if Seq is None:
            # internal proteinclass object
            self.p0 = protein.ProteinClass(Pdb, Model = Model)
            self.Seq = self.p0.Seq
            self.Pos = self.p0.Pos
            self.ResChainInds = [self.p0.ResChain(i) for i, r in enumerate(self.Seq)]
            self.Chains = self.p0.Chains
            self.NChains = len(self.Chains)
        # if a Seq is provided
        else:
            self.Seq = None
            self.Pos = None
            self.ResChainInds = []
            self.Chains = []
            self.NChains = None
            self.__SetChains(Seq)
        print 'Sequence: %s' % (' '.join(self.Seq))
        # sequence book-keeping
        self.ResTypes = list(set(self.Seq))
        self.NRes = len(self.Seq)
        self.NResTypes = len(self.ResTypes)
        # unpack sidechain atomtypes and assign to this object
        self.AtomSbyNum = []
        self.AtomSbyRes = {}
        self.__SetSideChains()
        # build startatoms for quick backbone referencing
        self.StartAtomInds = []
        self.RelativeStartAtomInds = []
        self.__SetStartAtomInds()
        self.__SetRelativeStartAtomInds()
        # atomnames
        self.AtomNames = []
        self.__SetAtomNames()
        # generate bonds
        self.BondPairs = None
        self.__SetBonds()
   
    def __repr__(self):
        s = ' '.join(self.Seq)
        return s
        
    def __SetSideChains(self):
        ''' sets sidechain atoms for this protein object
        to maintain compatibility with both NCOS and Go type of systems
        i.e. both named according to residue number (as a list) 
        and residue name (as a dict)
        this method can be extended for other kinds of reduced alphabet schemes
        '''
        AtomS = self.cfg.AtomS
        for i, r in enumerate(self.Seq):
            if AtomS[r] is None:
                this_AtomSbyRes = None
                this_AtomSbyNum = None
            else:
                this_AtomSbyRes = AtomS[r]
                this_AtomSbyNum = copy.copy(self.cfg.AtomS[r])
                this_AtomSbyNum.Name = 'S_%d' % i
            self.AtomSbyRes[r] = this_AtomSbyRes
            self.AtomSbyNum.append(this_AtomSbyNum)
        
    def __SetStartAtomInds(self):
        ''' set the index of first atoms for each residue'''
        x = 0
        for i, r in enumerate(self.Seq):
            self.StartAtomInds.append(x)
            if self.AtomSbyNum[i] is None: x += 3
            else: x += 4
    
    def __SetRelativeStartAtomInds(self):
        ''' set the index of first atoms for each residue relative to the chain'''
        x = 0
        CurrentChain = 0
        for i, r in enumerate(self.Seq):
            thisChain = self.ResChainInds[i]
            # check if transitioning chains and reset the counter
            if thisChain > CurrentChain:
                x = 0
                CurrentChain = thisChain
            self.RelativeStartAtomInds.append(x)
            if self.AtomSbyNum[i] is None: x += 3
            else: x += 4

    def __SetChains(self, Seq):
        ''' figure out multiple chain related stuff'''
        self.Seq = []
        Test = [isinstance(x, list) for x in Seq]
        if not all(Test): Seq = [Seq]
        for ii, chain in enumerate(Seq):
            if len(Seq) > 1: self.Chains.append(string.ascii_uppercase[ii])
            else: self.Chains.append('')
            for i, res in enumerate(chain):
                self.Seq.append(res)
                self.ResChainInds.append(ii)
            self.NChains = len(self.Chains)
    
    def __SetAtomNames(self):
        '''sets atom names. Assumes that sidechains have been already named / numbered'''
        AtomS = self.cfg.AtomS
        for i, r in enumerate(self.Seq):
            self.AtomNames.extend(['N', 'C', 'O'])
            if not AtomS[r] is None:
                if self.SSRefType == 'name': self.AtomNames.append(self.AtomSbyRes[r].Name)
                if self.SSRefType == 'number': self.AtomNames.append(self.AtomSbyNum[i].Name)
        
    def __SetBonds(self):
        ''' create a simple bond pair list for sim
        assumes that Chain information has already been created'''
        self.BondPairs = [ [] for i in range(self.NChains) ]
        for i, r in enumerate(self.Seq):
            AtomNInd = self.RelativeStartAtomInds[i]
            AtomCInd = AtomNInd + 1
            AtomOInd = AtomNInd + 2
            # intra-residue backbone bonds
            BondPairs = [(AtomNInd, AtomCInd), (AtomCInd, AtomOInd)]
            # backbone-sidechain bond
            if not self.AtomSbyNum[i] is None:
                # sidechain bonds only for residues whose sidechains are not None
                AtomSInd = AtomNInd + 3
                BondPairs.append( (AtomCInd, AtomSInd) )
            # peptide bond 
            if i < self.NRes - 1:
                NextAtomNInd = self.RelativeStartAtomInds[i+1]
                # prevent inter-chain bonding
                if NextAtomNInd > AtomNInd: BondPairs.append( (AtomOInd, NextAtomNInd) )
            # add to master list
            ChainInd = self.ResChainInds[i]
            self.BondPairs[ChainInd].extend(BondPairs)
    
    def BondNativeContacts(self, ContactDict):
        if Verbose: print 'Bonding sidechains of native contacts for applying harmonic restraints'
        SInds = self.GetSInds()
        for i,j in ContactDict['c_native']:
            self.BondPairs.append( (SInds[i], SInds[j]) ) 
    
    def GetBBInds(self, ResNums = None):
        ''' extract the backbone indices of the structure'''
        BBInds = []
        if ResNums is None: ResNums = range(self.NRes)
        for i in ResNums:
            k = self.StartAtomInds[i]
            BBInds.extend( [k, k+1, k+2] )
        return BBInds
   
    def GetSInds(self, ResNums = None):
        ''' extract the sidechain indices of the structure
            returns alpha carbon index for residues with None sidechain'''
        SInds = []
        if ResNums is None: ResNums = range(self.NRes)
        for i in ResNums:
            idx = 1 if self.AtomSbyNum[i] is None else 3
            SInds.append(self.StartAtomInds[i] + idx)
        return SInds

    def GetResPos(self, ResNums = None):
        ''' get the COM of the residues'''
        if ResNums is None: ResNums = range(self.NRes)
        ResPos = np.zeros([len(ResNums), 3])
        for ii, i in enumerate(ResNums):
            start = self.StartAtomInds[i]
            stop = self.StartAtomInds[i+1] if i < self.NRes-1 else len(self.Pos)
            ResPos[ii, :] = np.mean(self.Pos[start:stop, :], axis = 0)
        return ResPos

    def GetResContactList(self, ResNums = None):
        ''' find native contacts with ResRadius'''
        if ResNums is None: ResNums = range(self.NRes)
        ResContactList = []
        ResPos = self.GetResPos(ResNums = ResNums)
        for i in ResNums:
            for j in ResNums:
                # ignore same residue
                if i == j: continue
                # ignore adjacent residues
                if abs(i-j) < MinCO: continue
                # calculate distance between residues
                d_ij = ResPos[j] - ResPos[i]
                # check if contact
                dsq = np.sum(d_ij * d_ij)
                if dsq <= ResRadius * ResRadius: ResContactList.append( (i,j) )
        return ResContactList

    def GetResSSList(self, ResNums = None):
        ''' find all residue pairs that are MinCO residues away'''
        if ResNums is None: ResNums = range(self.NRes)
        ResSSList = []
        for i in ResNums:
            for j in ResNums:
                # ignore same residue
                if i == j: continue
                # ignore adjacent residues
                if abs(i-j) < MinCO: continue
                ResSSList.append( (i,j) )
        return ResSSList
    

def MakeSys(p, cfg = None, NChains = 1):
    '''generate the topology for the sim Sys object
       currently implemented for single chains
       p is ProteinNCOS type object
       cfg is the global config object '''
    if Verbose == True:
        print 'Creating System...'
    # generate the molecule as a list of lists
    AtomList = [ [] for i in range(p.NChains)]
    s = 'Adding side chain atoms by %s: ' % p.SSRefType
    for i, r in enumerate(p.Seq):
        # backbone atoms
        res = [AtomN, AtomC, AtomO]
        if r == 'GLY':
            if p.hasSpecialGLYParams: res = [AtomN, AtomC_GLY, AtomO]
        if r == 'PRO':
            if p.hasSpecialPROParams: res = [AtomN, AtomC_PRO, AtomO]
        # sidechain atoms
        # ref by name
        if p.SSRefType == 'name':
            if not p.AtomSbyRes[r] is None:
                res.append(p.AtomSbyRes[r])
                s += '%s' % res[-1].Name
                if i < p.NRes-1: s += ', '
        # ref by number
        elif p.SSRefType == 'number':
            if not p.AtomSbyNum[i] is None:
                res.append(p.AtomSbyNum[i])
                s += '%s (%s)' % (res[-1].Name, r)
                if i < p.NRes-1: s+= ', '
        else:
            print 'Error: Unknown sidechain reference type, or not implemented yet'
            exit()
        # figure out correct chain
        ChainInd = p.ResChainInds[i]
        # add to AtomList maintaining correct chain
        AtomList[ChainInd].extend(res)
    if Verbose: print s
    # generate the sim-style MolTypes
    MolList = []
    for i in range(p.NChains):
        Mol =  sim.chem.MolType(p.Prefix, [atom for atom in AtomList[i]])
        # generate bonds for this molecule
        for (m, n) in p.BondPairs[i]: Mol.Bond(m, n)
        # populate MolList
        MolList.append(Mol)
    # create System
    World = sim.chem.World(MolList, Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(World, Name = p.Prefix)
    if p.NChains > 1: print 'Creating oligomer...'
    # add molecules to the system
    for i, Mol in enumerate(MolList):
        if Verbose and p.NChains > 1:
            print 'Added chain %d' % i
        Sys += Mol.New()
    return Sys

