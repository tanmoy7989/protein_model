#!/usr/bin/env python

''' Builds the toplogy of the sim-style Sys object '''

import os, numpy as np
from const import *
import sim
import protein

Verbose = True
protein.MinCODflt = MinCO
protein.ResRadiusDflt = ResRadius

class ProteinNCOS(object):
    ''' A rudimentary object to create and store the NCOS toplogy 
        The atom ordering within a residue is N C O S
    '''    
    def __init__(self, Pdb = None, Seq = None, Prefix = 'ncos'):
        self.Prefix = Prefix
        if Seq is None:
            # internal proteinclass object
            self.p0 = protein.ProteinClass(Pdb)
            self.Seq = self.p0.Seq
        else:
            self.Seq = Seq
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
            if r == 'GLY': x += 3
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
            if not r == 'GLY':
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
    
    def MakeSys(self, cfg, NMols = 1):
        '''generate the topology for the sim Sys object
           currently implemented for single chains
           cfg is the global config object '''
        hasSpecialBBTorsions = cfg.hasSpecialBBTorsions
        if Verbose == True:
            print 'Generating System topology...'
        # generate the molecule
        AtomList = []
        for i,r in enumerate(self.Seq):
            # backbone atoms
            if hasSpecialBBTorsions:
                if r == 'GLY': res = [AtomN, AtomC_GLY, AtomO]
                elif r == 'PRO': res = [AtomN, AtomC_PRO, AtomO]
            else: res = [AtomN, AtomC, AtomO]
            # sidechain
            if not r == 'GLY': res.append(AtomS[r])
            AtomList.extend(res)
        Mol = sim.chem.MolType(self.Prefix, [i for i in AtomList])
        # generate bonds
        for (i,j) in self.BondPairs:
            Mol.Bond(i,j)
        # create System
        World = sim.chem.World([Mol], Dim = 3, Units = sim.units.AtomicUnits)
        Sys = sim.system.System(World, Name = self.Prefix)
        #for i in range(NChains): Sys += Mol.New()
        Sys += Mol.New()
        return Sys
