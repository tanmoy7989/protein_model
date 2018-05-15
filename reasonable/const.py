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
AtomC_GLY = sim.chem.AtomType('C_GLY', Mass = AtomMass['C'])
AtomC_PRO = sim.chem.AtomType('C_PRO', Mass = AtomMass['C'])
AtomS_GLY = sim.chem.AtomType('S_GLY', Mass = 1.008)
DfltAtomS = {}
for r, mass in ResMass.iteritems():
    if r == 'BB': continue
    if r == 'GLY': DfltAtomS[r] = None
    else: DfltAtomS[r] = sim.chem.AtomType('S_%s' % r, Mass = mass - ResMass['BB'])

# contact prediction
ResRadius = 8.0 #A
MinCO = 3

# Miyazawa-Jernigan interaction matrix (in RT units)
# citation: Miyazawa, Jernigan, Macromolecules 1985, 18, 434-552
MJResNames = ['CYS', 'MET', 'PHE', 'ILE', 'LEU', 
              'VAL', 'TRP', 'TYR', 'ALA', 'GLY', 
              'THR', 'SER', 'GLN', 'ASN', 'GLU', 
              'ASP', 'HIS', 'ARG', 'LYS', 'PRO']
MJDATA = (kB * RoomTemp) * np.array([ \
[-5.44, -5.05, -5.63, -5.03, -5.03, -4.46, -4.76, -3.89, -3.38, -3.16, -2.88, -2.86, -2.73, -2.59, -2.08, -2.66, -3.63, -2.70, -1.54, -2.92],
[ 0.00, -6.06, -6.68, -6.33, -6.01, -5.52, -6.37, -4.92, -3.99, -3.75, -3.73, -3.55, -3.17, -3.50, -3.19, -2.90, -3.31, -3.49, -3.11, -4.11],
[ 0.00,  0.00, -6.85, -6.39, -6.26, -5.75, -6.02, -4.95, -4.36, -3.72, -3.76, -3.56, -3.30, -3.55, -3.51, -3.31, -4.61, -3.54, -2.83, -3.73],
[ 0.00,  0.00,  0.00, -6.22, -6.17, -5.58, -5.64, -4.63, -4.41, -3.65, -3.74, -3.43, -3.22, -2.99, -3.23, -2.91, -3.76, -3.33, -2.70, -3.47],
[ 0.00,  0.00,  0.00,  0.00, -5.79, -5.38, -5.50, -4.26, -3.96, -3.43, -3.43, -3.16, -3.09, -2.99, -2.91, -2.59, -3.84, -3.15, -2.63, -3.06],
[ 0.00,  0.00,  0.00,  0.00,  0.00, -4.94, -5.05, -4.05, -3.62, -3.06, -2.95, -2.79, -2.67, -2.36, -2.56, -2.25, -3.38, -2.78, -1.95, -2.96],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -5.42, -4.44, -3.93, -3.37, -3.31, -2.95, -3.16, -3.11, -2.94, -2.91, -4.02, -3.56, -2.49, -3.66],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -3.55, -2.85, -2.50, -2.48, -2.30, -2.53, -2.47, -2.42, -2.25, -3.33, -2.75, -2.01, -2.80],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -2.51, -2.15, -2.15, -1.89, -1.70, -1.44, -1.51, -1.57, -2.09, -1.50, -1.10, -1.81],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -2.17, -2.03, -1.70, -1.54, -1.56, -1.22, -1.62, -1.94, -1.68, -0.84, -1.72],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -1.72, -1.59, -1.59, -1.51, -1.45, -1.66, -2.31, -1.97, -1.02, -1.66],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -1.48, -1.37, -1.31, -1.48, -1.46, -1.94, -1.22, -0.83, -1.35],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -0.89, -1.36, -1.33, -1.26, -1.85, -1.85, -1.02, -1.73],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -1.59, -1.43, -1.33, -2.01, -1.41, -0.91, -1.43],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -1.18, -1.23, -2.27, -2.07, -1.60, -1.40],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -0.96, -2.14, -1.98, -1.32, -1.19],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -2.78, -2.12, -1.09, -2.17],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -1.39, -0.06, -1.85],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.13, -0.67],
[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -1.18]
])

# Accessing the MJ matrix data (upper triangular elements)
MJMATRIX = {}
for i in range(len(MJResNames)):
    for j in range(i, len(MJResNames)):
        res_i = MJResNames[i]
        res_j = MJResNames[j]
        MJMATRIX[ (res_i, res_j) ] = -MJDATA[i,j]
        MJMATRIX[ (res_j, res_i) ] = -MJDATA[i,j]

# LAMMPS settings
LAMMPSEXEC = os.environ['LAMMPSEXEC']
MaxLammpsPairEnekBT = 20 # might need to increase to prevent sim crashes

# model defaults
# put anything here that you want to be dynamically configurable
DEFAULTS = dict(

            # sidechain naming scheme
            SSRefType = 'name',
            
            # sim-style sidechain objects (none for glycine by default)
            AtomS = DfltAtomS,

            # backbone-sidechain options (1-alphabet for bonded and constant repulsive for non-bonded by default)
            Bonded_NCOSType = 1,
            NCOSType = 2,
                
            # Go model types (not a Go model by default)
            NativeType = -1,
            NonNativeType = -1,
                
            # backbone
            MinBondOrd = 5,
            NKnot = 40,
            SPCut = 10.0,
            hasSpecialBBGLYAngles = False,
            hasSpecialBBGLYTorsions = False,
            hasSpecialBBPROAngles = False,
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
       
            # MJ model
            MJSigma = 4.0, # A
            MJSigmas = [],
            MJAutoSigma = False,

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
MAPSCRIPT = os.path.expanduser('~/protein_model/reasonable/mapNCOS.py')
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

