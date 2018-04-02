#!/usr/bin/env python

''' Defines all sim-style global objects for this model'''

import os, numpy as np
import sim

# energy constants
kB = 0.001987 # Boltzmann constant in kcal/mol
RoomTemp = 300.0 # K
Escale = lambda T: kB * T # all energies are expressed in these units

# Atomic masses
# Ref: Wikipedia
AtomMass = {'N' : 14.006, 
            'H' :  1.008, 
            'C' : 12.010, 
            'O' : 15.999,
            'S' : 32.065
}

# Residue Masses (dry weight, i.e. minus weight of H2O)
# Ref: http://www.its.caltech.edu/~ppmal/sample_prep/work3.html
ResMass = {'ALA' :  71.079,
           'ARG' : 156.188, 
           'ASN' : 114.104, 
           'ASP' : 115.089, 
           'CYS' : 103.145, 
           'GLN' : 128.131, 
           'GLU' : 129.116, 
           'GLY' :  57.052, 
           'HIS' : 137.141, 
           'ILE' : 113.160, 
           'LEU' : 113.160, 
           'LYS' : 128.170, 
           'MET' : 131.199, 
           'PHE' : 147.177, 
           'PRO' :  97.117, 
           'SER' :  87.078, 
           'THR' : 101.105, 
           'TRP' : 186.213, 
           'TYR' : 163.176, 
           'VAL' :  99.130,
           'BB'  :  56.041 # this is the mass of 'NH-CH-CO'
}

# sim-style atom objects for backbone atoms
# ignoring hydrogens and coarse graining only heavy atoms
# special atom objects for GLY and PRO
# Sidechains can be configurable and are shifted into the DEFAULTS dict
AtomN = sim.chem.AtomType('N', Mass = AtomMass['N'])
AtomC = sim.chem.AtomType('C', Mass = AtomMass['C'])
AtomO = sim.chem.AtomType('O', Mass = AtomMass['C'] + AtomMass['O'])
AtomC_GLY = sim.chem.AtomType('C', Mass = AtomMass['C'])
AtomC_PRO = sim.chem.AtomType('C', Mass = AtomMass['C'])
AtomS_GLY = sim.chem.AtomType('S_GLY', Mass = 1.008)
DfltAtomS = {}
for r, mass in ResMass.iteritems():
    if r == 'BB': continue
    if r == 'GLY': DfltAtomS[r] = None
    else: DfltAtomS[r] = sim.chem.AtomType('S_%s' % r, Mass = mass - ResMass['BB'])

# contact prediction
ResRadius = 8.0 #A
MinCO = 3

# LAMMPS binary
LAMMPSEXEC = os.environ['LAMMPSEXEC']

# model defaults
# put anything here that you want to be dynamically configurable
DEFAULTS = dict(

            # sidechain naming scheme
            SSRefType = 'name',
            
            # sim-style sidechain objects (none for glycine by default)
            AtomS = DfltAtomS,

            # backbone-sidechain options (full 21-alphabet by default)
            Bonded_NCOSType = 0,
            NCOSType = 0,
                
            # Go model types (not a Go model by default)
            NativeType = -1,
            NonNativeType = -1,
                
            # backbone
            MinBondOrd = 5,
            NKnot = 40,
            SPCut = 10.0,
            hasSpecialBBGLYTorsions = False,
            hasSpecialBBPROTorsions = False,

            # Go models
            # 1) LJ and splined native interactions
            NativeEpsilon = 4 * kB * RoomTemp,
            NativeCut = 1.2 * ResRadius,
            
            # 2) Harmonic restraints
            NativeHarmonicFluct = 1.0,
            NativeFConst = None,
            Map2Polymer = False,
            PolyName = None,

            # 3) LJ and splined nonnative interactions
            NonNativeEpsilon = 4 * kB * RoomTemp,
            NonNativeCut = None, # force a WCA-type cutoff
        
            # timestep
            TimeStep = 1.0, #fs
                
            # MD iterations
            NStepsMin = 10000,
            NStepsEquil = 10000000,
            NStepsProd = 10000000,
            NStepsSwap = 1000,
            StepFreq = 1000
)

# file formats
FMT = {'TRAJ'  : '%s.%3.2f.lammpstrj.gz',
       'PDB'   : '%s.%3.2f.pdb',
       'ENE'   : '%s.%3.2f.ene.dat.gz'
}

# native struct locations
NATIVEPATH = {'Unmapped'  : os.path.expanduser('~/protein_model/native_struct/unmapped'),
              'Mapped'    : os.path.expanduser('~/protein_model/native_struct/mapped')} 

# mapping script
MAPSCRIPT = os.path.expanduser('~/protein_model/reasonable/map.py')
HELPSTR = '''
CG Model Types:
SSRefType: controls whether sidechains are named according to name or number of the residue
(number required for Go models)
SSRefType = 'name' -> sidechains named as 'S_<resname>' # DEFAULT
SSRefType = 'number' -> sidechains names as 'S_<number>'

Bonded_NCOSType: controls backbone-sidechain bonded potentials'
Bonded_NCOSType = -1 -> no backbone-sidechain bonded potentials (why ?!)
Bonded_NCOSType = 0 -> 21-alphabet backbone-sidechain bonded potentials
Bonded_NCOSType = 1 -> 1-alphabet backbone-sidechain bonded potentials

NCOSType: controls backbone-sidechain nonbonded potentials
NCOSType = -1 -> no backbone-sidechain nonbonded potentials
NCOSType = 0 -> 21-alphabet backbone-sidechain nonbonded potentials
NCOSType = 1 -> 1-alphabet backbone-sidechain nonbonded potentials
NCOSType = 2 -> 1-alphabet constant soft core repulsive backbone-sidechain nonbonded potentials

NativeType and NonNativeType: controls sidechain-sidechain potentials in Go models
NativeType = -1 -> no interaction between native contacts
NativeType = 0 -> single LJ potential between native contacts
NativeType = 1 -> spline potential between native contacts
NativeType = 2 -> harmonic restraints between native contacts

NonNativeType = -1 -> no interaction between nonnative contacts
NonNativeType = 0 -> single WCA potential between native contacts
NonNativeType = 1 -> spline potential commensurate with a WCA between nonnative contacts
'''

