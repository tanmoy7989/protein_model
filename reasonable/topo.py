#!/usr/bin/env python

''' Builds the toplogy of the sim-style Sys object '''

import os, numpy as np, copy
from const import *
import sim
import protein, mdsim

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
        if Verbose: print 'Creating a NCOS protein object...'
        # set Prefix
        self.Prefix = Prefix
        # has special torsion potentials
        self.hasSpecialBBGLYTorsions = cfg.hasSpecialBBGLYTorsions
        self.hasSpecialBBPROTorsions = cfg.hasSpecialBBPROTorsions
        # has special sidechains for glycine
        self.hasPseudoGLY = cfg.hasPseudoGLY()
        if self.hasSpecialBBGLYTorsions and self.hasPseudoGLY:
            print 'Error: Cannot have pseudo GLY side chain and special GLY BB torsion simultaneously'
            exit()
        # sidechain referencing
        self.SSRefType = cfg.SSRefType
        # extract the entire config object just in case
        self.cfg = cfg
        # sequence book-keeping
        if Seq is None:
            # internal proteinclass object
            self.p0 = protein.ProteinClass(Pdb, Model = Model)
            self.Seq = self.p0.Seq
            self.Pos = self.p0.Pos
        else:
            self.Seq = Seq
            self.Pos = None
        self.ResTypes = list(set(self.Seq))
        self.NRes = len(self.Seq)
        self.NResTypes = len(self.ResTypes)
        # unpack sidechain atomtypes and assign to this object
        self.AtomSbyNum = []
        self.AtomSbyRes = {}
        self.__SetSideChains()
        # build startatoms for quick backbone referencing
        self.StartAtomInds = []
        self.__SetStartAtomInds()
        # generate bonds
        self.BondPairs = []
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
        
    def __SetBonds(self):
        ''' create a simple bond pair list for sim'''
        for i, r in enumerate(self.Seq):
            AtomNInd = self.StartAtomInds[i]
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
                NextAtomNInd = self.StartAtomInds[i+1]
                BondPairs.append( (AtomOInd, NextAtomNInd) )
            # add to master list
            self.BondPairs.extend(BondPairs)
    
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
    

def MakeSys(p, cfg = None, NMols = 1):
    '''generate the topology for the sim Sys object
       currently implemented for single chains
       p is ProteinNCOS type object
       cfg is the global config object '''
    if Verbose == True:
        print 'Generating backbone topology...'
    # generate the molecule
    AtomList = []
    s = ' Added side chain atoms by %s: ' % p.SSRefType
    for i,r in enumerate(p.Seq):
        # backbone atoms
        if p.hasSpecialBBGLYTorsions:
            if r == 'GLY': res = [AtomN, AtomC_GLY, AtomO]
        elif p.hasSpecialBBPROTorsions:
            if r == 'PRO': res = [AtomN, AtomC_PRO, AtomO]
        else: res = [AtomN, AtomC, AtomO]
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
        # add to AtomList 
        AtomList.extend(res)
    s += '\n'
    if Verbose: print s
    Mol = sim.chem.MolType(p.Prefix, [i for i in AtomList])
    # generate bonds
    for (i,j) in p.BondPairs: Mol.Bond(i,j)
    # create System
    World = sim.chem.World([Mol], Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(World, Name = p.Prefix)
    #for i in range(NChains): Sys += Mol.New()
    Sys += Mol.New()
    return Sys



