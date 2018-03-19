'''
/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <tanmoy.7989@gmail.com> wrote this file.  As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return. Tanmoy Sanyal
 * ----------------------------------------------------------------------------
 */
'''
 
#!/usr/bin/env python

import os, sys, numpy as np, time
from multiprocessing import Pool
import sim, protein, parse_potential as pp, pickleTraj, utils

#### GLOBALS ####

# output location
Prefix = 'Go'
OutputDir = os.getcwd()

# file format info
FMT = utils.FMT

# native struct info
InPdb = None
Seq = None
NRes = None
ResTypes = None
NResTypes = None
ProtObj = None # proteinclass object created through protein.py
ResContactList = None
ResMass = {'ALA': 71.079, 'ARG': 156.188, 'ASN': 114.104, 'ASP': 115.089, 'CYS': 103.145, 'GLN': 128.131, 'GLU': 129.116, 'GLY': 57.052, 'HIS': 137.141, 'ILE': 113.160, 'LEU': 113.160, 'LYS': 128.17, 'MET': 131.199, 'PHE': 147.177, 'PRO': 97.117, 'SER': 87.078, 'THR': 101.105, 'TRP': 186.213, 'TYR': 163.176, 'VAL': 99.13, 'BB': 56.0}

# contact enumeration settings
ResRadius = 8.0
MinCO = 3

# ambient conditions
TempSet = 300.00
BoxL = 0.0

# backbone optimization types
# NCOSType = 0 for retaining backbone-sidechain potentials (NS, OS, CS)
# NCOSType = 1 for ignoring these potentials
NCOSType = 0

# backbone forcefield parameters
NCOS_forcefield_file = None
BB_forcefield_file = None
MinBondOrd = 5
NKnot = 20
SPCut = 10.0

# Common Go params
NativeCut = None
NonNativeCut = None

# Go Types (Native Types)
# GoType = -1 for random coils with no side chain sites and interactions
# GoType = 0 for random coils with side chain sites but no interactions
# GoType = 1 for Native LJ with single sigma
# GoType = 2 for Native LJ with varying sigma for each residue pair
# GoType = 3 for Native Spline
# GoType = 4 for harmonic restraint with a single force constant
# IncludeGly = True includes native contact pairs one or both member/(s) of which is/are GLY 
GoType = 1
IncludeGLY = False

# Native LJ  params
NativeSigma = None
NativeEpsilon = None

# Native Spline params
NativeNKnot = 20
NativeKnots = None

# Native Harmonic Restraint params
NativeFConst = None
HarmonicFluct = 1.0 # A (thermal fluctuations of S-S distance under an ENM type potential at room temp.)

# Non-Native Types
# NonNativeType = -1 # no non-native potentials
# NonNative Type = 0 # choose non-native sigma and epsilon same as native
# NonNativeType = 1 # choose non-native sigma as longest non-native contact, epsilon same as native
# NonNativeType = 2 # make it a spline, oh yeah!
hasNonNative = True
NonNativeType = 0

# NonNative WCA params
NonNativeSigma = None
NonNativeEpsilon = None

# NonNative Spline params
NonNativeNKnot = 20
NonNativeKnots = None

# REMD settings
NStepsMin = 100000
NStepsEquil = 50000000
NStepsProd = 50000000
NStepsSwap = 1000
StepFreq = 10000
Temps = np.logspace(np.log10(280), np.log10(500), 24)

# Lammps settings
LammpsExec = 'lmp_mpich2'
sim.srel.optimizetrajlammps.useREMD = True
sim.export.lammps.InnerCutoff = 0.02

# Srel settings
useBondHessian = True
sim.srel.optimizetraj.PlotFmt = 'svg'
OptStages = ['Bond', 'Angle', 'Torsion', 'NonBond', 'all' , 'Go1' , 'Go2']
sim.srel.base.DiffEneFracTol = 0.1 # relaxed tolerance to prevent crash due to lammps-sim mismatch
LoadedPotentials = [] # potentials that were loaded from supplied forcefield file
PermaFrost = [] # potentials whose initial estimates are good enough, so no need to optimize
MaxIter = {'Bond': None, 'Angle': None, 'Torsion': None, 'NonBond': None, 'all': None}

# Initial Seed Pos manipulation flags
# if this is true, make sure an appropriate RefPosInd is used
# might not necessarily be frame 0
InitPdb = None
useSeed = False
RefPosInd = None

# Extended ensemble settings
# ExType = 0 for no extended ensemble
# ExType = 1 for NCOS extended ensemble from ensembles of homopolymers
# ExType = 2 for Go extended ensemble from ensemble of different peptides (keeping backbone frozen)
ExType = 0

######## UTILITY FUNCTIONS ########

# convert from residue to absolute atom indices
def res2atom(resnum, resatomind):
    if resatomind > 3: raise TypeError('Each residue has 4 CG beads numbered 0,1,2,3')
    startind = 0
    for i,r in enumerate(Seq[:resnum]):
        if r == 'GLY': startind += 3
        else: startind += 4
    return startind + resatomind
    
    
# parse native pdb and get contact details
def parseNative():
    global OutputDir
    global InPdb, Seq, NRes, ResTypes, NResTypes, ResContactList
    # set contact settings
    protein.MinCODflt = MinCO
    protein.ResRadiusDflt = ResRadius
    # read native structure and parse seq
    if Seq is None or (not InPdb is None and os.path.isfile(InPdb)):
        p = protein.ProteinClass(Pdb = InPdb)
        p.Dehydrogen() ; p = p.Decap()
        Seq = p.Seq
        # record contact distances and contact distance histogram
        ResContactList = p.ResContactList()
        p_cg = p.CoarseGrained()
        ResPos = p_cg.ResPos()
        x_native = []
        x_nonnative = []
        x_native_ss = []
        print 'Collecting COM and sidechain distances between native and non-native contacts'
        for i in range(len(Seq)-1):
            for j in range(i+1, len(Seq)):
                # flags
                status = ''
                UsedAlphaC = ''
                # record COM distances between residues
                d_ij = ResPos[j] - ResPos[i]
                d = np.sqrt(np.sum(d_ij * d_ij))  
                # record S-S distances between residues
                ind_m = 3
                ind_n = 3
                # use alpha-C positions if GLY
                if Seq[i] == 'GLY':
                    ind_m = 1
                    UsedAlphaC += '%s ' % i
                if Seq[j] == 'GLY':
                    ind_n = 1
                    UsedAlphaC += '%s ' % j
                m = res2atom(i, ind_m)
                n = res2atom(j, ind_n)
                d_ss_ij = p.Pos[n] - p.Pos[m]
                d_ss = np.sqrt(np.sum(d_ss_ij * d_ss_ij))
                # native contacts
                if ResContactList.__contains__( (i,j) ):
                    x_native.append(d)
                    x_native_ss.append(d_ss)
                    status = 'native'
                # non-native contacts
                else: 
                    if abs(i-j) >= MinCO:
                        x_nonnative.append(d) # all contacts should be 3 apart at least
                        status = 'non-native'
                    else:
                        status = 'adjacent'
                # print stats
                if not UsedAlphaC: UsedAlphaC = 'none' 
                print 'Contact: (%3d, %3d), Residues: (%3s, %3s), Status: %12s, Used alpha-C for residue: %3s' % (i,j, Seq[i],Seq[j], status, UsedAlphaC)
        np.savetxt(os.path.join(OutputDir ,'nativecontactdist.txt'), np.array(x_native))
        np.savetxt(os.path.join(OutputDir, 'nativessdist.txt'), np.array(x_native_ss))
        np.savetxt(os.path.join(OutputDir, 'nonnativecontactdist.txt'), np.array(x_nonnative))
        file(os.path.join(OutputDir, 'nativecontactlist.txt'), 'w').write(str(ResContactList))
    NRes = len(Seq)
    ResTypes = list(set(Seq))
    NResTypes = len(ResTypes)
    return

    
# get backbone positions for CG peptides
def getBBInds(resnums = None):
    global Seq
    if resnums is None: resnums = range(len(Seq))
    BBInds = []
    for r in resnums:
        inds = [res2atom(r,x) for x in [0,1,2]]
        BBInds.extend(inds)
    return BBInds


# generate random initial position based on /share/apps/scripts/template.pdb
def genRandInitPos():
    global Seq, InitPdb
    tmpPdb = 'tmp.pdb'
    if InitPdb is None: InitPdb = 'init.pdb'
    if not os.path.isfile(InitPdb):
        print 'Generating Initial structure...'
        p = protein.ProteinClass(Seq = Seq)
        p.Dehydrogen()
        p = p.Decap()
        p.WritePdb(tmpPdb)
        MAPSCRIPT = '/home/cask0/home/tsanyal/Go/map.py'
        os.system('python %s %s %s' % (MAPSCRIPT, tmpPdb, InitPdb.split('.')[0]))
        if os.path.isfile(tmpPdb): os.remove(tmpPdb)
    print 'Using init conf as generated in : %s' % InitPdb
    p = protein.ProteinClass(InitPdb)
    initpos = p.Pos
    return initpos



####### GENERATE GO POTENTIALS ########

def genGoPotentials(Sys):
    #TODO: write for multiple chains
    global OutputDir
    global NRes, ResContactList, ResRadius, MinBondOrd
    global GoType, IncludeGLY, NonNativeType
    global NativeCut, NonNativeCut
    global NativeSigma, NativeEpsilon, NativeNKnot, NativeKnots
    global NativeFConst, HarmonicFluct
    global NonNativeSigma, NonNativeEpsilon, NonNativeNKnot, NonNativeKnots
    
    # extract native and non-native distances
    d_native = np.loadtxt(os.path.join(OutputDir, 'nativecontactdist.txt'))
    d_nonnative = np.loadtxt(os.path.join(OutputDir, 'nonnativecontactdist.txt'))
    d_ss = np.loadtxt(os.path.join(OutputDir, 'nativessdist.txt'))
    
    # create side chain pair potential filters
    # do this only if has LJ type native interactions and/or non-native interactions
    Test1 = not (GoType == -1 or GoType == 0 or GoType == 4)
    Test2 = not NonNativeType == -1
    if (Test1 or Test2):
        print 'Creating side chain potential filters'
        NativePairs = np.zeros([Sys.World.NAID, Sys.World.NAID], int)
        NonNativePairs = np.zeros([Sys.World.NAID, Sys.World.NAID], int)
        dict_d_native = {}
        AID2Res = {}
        for i in range(NRes-1):
            for j in range(i+1, NRes):
                if Seq[i] == 'GLY' or Seq[j] == 'GLY':
                    if IncludeGLY:
                        # include GLY
                        print 'Native Contact Pair: (%d, %d) Retaining C-alpha for Glycines' % (i,j)
                    else:
                        # don't include glycines
                        print 'Native Contact Pair: (%d, %d) Ignoring pairs with Glycines' % (i,j)
                        continue
                if IncludeGLY:
                    idx_i = 1 if Seq[i] == 'GLY' else 3
                    idx_j = 1 if Seq[j] == 'GLY' else 3
                else:
                    idx_i = 3
                    idx_j = 3
                m = Sys.Mol[0][res2atom(i,idx_i)].AID
                n = Sys.Mol[0][res2atom(j,idx_j)].AID
                # contact residues (LJ)
                if ResContactList.__contains__((i,j)):
                    NativePairs[m, n] = 1
                    NativePairs[n, m] = 1
                    ind = ResContactList.index( (i,j) )
                    AID2Res[ (m,n) ] = (i,j)
                    dict_d_native[ (m,n)] = d_native[ind] * (2 ** (-1/6.))
                # non-contact residues (WCA)
                else:
                    if abs(i-j) >= MinCO:
                        NonNativePairs[m, n] = 1
                        NonNativePairs[n, m] = 1
        Filter_Native = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = NativePairs, MinBondOrd = MinBondOrd)
        Filter_NonNative = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = NonNativePairs, MinBondOrd = MinBondOrd)
    
    # No Native interactions
    if GoType == -1 or GoType == 0:
        P_n = None

    # Native LJ with all sidechains having same sigma
    # that gives the min inter-res contact distance
    if GoType == 1:
        if NativeSigma is None: NativeSigma = d_native.min() * (2 **(-1/6.))
        if NativeCut is None: NativeCut = 1.2 * ResRadius # since no native pair is at greater inter-res distance than ResRadius (then add 20% as safety)
        print 'Using LJ Native interactions with sigma = %g A, Epsilon = %g kT, Cutoff = %gA' % (NativeSigma, NativeEpsilon / 0.6, NativeCut)
        P_n = sim.potential.LJ(Sys, Filter = Filter_Native, Cut = NativeCut, 
                               Sigma = NativeSigma, Epsilon = NativeEpsilon, 
                               Shift = True, Label = 'NonBondNative')
    
    # Native LJ with sidechains having different sigma
    # based on their separate inter-res contact distances
    if GoType == 2:
        if NativeCut is None: NativeCut = 1.2 * ResRadius # since no native pair is at greater inter-res distance than ResRadius (then add 20% as safety)
        LJA = np.zeros([Sys.World.NAID, Sys.World.NAID], np.float64)
        LJB = np.zeros([Sys.World.NAID, Sys.World.NAID], np.float64)
        print 'Populating LJ Matrix'
        for i in range(Sys.World.NAID - 1):
            for j in range(i+1, Sys.World.NAID):
                if NativePairs[i,j] == 1:
                    print 'Native Contact Pair: ', AID2Res[ (i,j) ], 'LJ-Sigma = %g A' % dict_d_native[ (i,j) ]
                    LJA[i,j] = LJA[j,i] = 4 * NativeEpsilon * dict_d_native[ (i,j) ] ** 12.0
                    LJB[i,j] = LJB[j,i] = 4 * NativeEpsilon * dict_d_native[ (i,j) ] ** 6.0
        P_n = sim.potential.LJMatrix(Sys, Filter = Filter_Native, Cut = NativeCut, LJA = LJA, LJB = LJB,
                                     Shift = True, Label = 'NonBondNative')
    
    # Native Spline
    if GoType == 3:
        if NativeCut is None: NativeCut = 1.2 * ResRadius # since no native pair is at greater inter-res distance than ResRadius (then add 20% as safety)
        print 'Using Spline Native interactions with Cutoff = %g A' % NativeCut
        P_n = sim.potential.PairSpline(Sys, Filter = Filter_Native, NKnot = NativeNKnot, Cut = NativeCut, Label = 'NonBondNative')
        if not NativeKnots is None: P_n.SetParam(Knots = NativeKnots)
        
    # Native harmonic constraints
    if GoType == 4:
        if NativeFConst is None:
            # generate a weak force-constant
            # from s.d. of a gaussian with this bias energy,
            # Fconst = (k/2) = (kT / 2 * variance), variance = HarmonicFluct
            # increasing the variance decreases the strength of the harmonic bias
            NativeFConst = 0.6 / (2. * HarmonicFluct **2.) # kT = 0.6 at room temp.
        print 'Using Harmonic Native restraints with k = %g kcal/(mol A^2)' % NativeFConst
        # add separate bond potentials between side-chains of native residues
        # could have written using filters developed above, but faster as written here.
        # Also, ideally a "BondMatrix" type potential like the "LJMatrix" would be go easier on memory
        # but until that is implemented, we roll with this.
        P_n = []
        for r, (i,j) in enumerate(ResContactList):
            if Seq[i] == 'GLY' or Seq[j] == 'GLY':
                if IncludeGLY:
                    # include glycines
                    print 'Native Contact Pair: (%d, %d) , d0 = %g A, Retaining C-alpha for Glycines' % (i, j, d_ss[r])
                else:
                    # dont include glycines
                    print 'Native Contact Pair: (%d, %d), d0 = %g A, Ignoring pairs with Glycines' % (i, j, d_ss[r])
                    continue
            else:
                print 'Native Contact Pair: (%d, %d), d0 = %g A' % (i, j, d_ss[r])
            this_NativePair = np.zeros([Sys.World.NAID, Sys.World.NAID], int)
            if IncludeGLY:
                idx_i = 1 if Seq[i] == 'GLY' else 3
                idx_j = 1 if Seq[j] == 'GLY' else 3
            else:
                idx_i = 3
                idx_j = 3
            m = Sys.Mol[0][res2atom(i,idx_i)].AID
            n = Sys.Mol[0][res2atom(j,idx_j)].AID
            this_NativePair[m,n] = 1
            this_NativePair[n,m] = 1
            this_d0 = d_ss[r]
            this_Filter = sim.atomselect.PolyFilter([sim.atomselect.All, sim.atomselect.All], AIDPairs = this_NativePair, MinBondOrd = MinBondOrd) 
            this_P = sim.potential.Bond(Sys, FConst = NativeFConst, Dist0 = this_d0, Filter = this_Filter, Label = 'NonBondNative_%d_%d' % (i,j))
            P_n.append(this_P)

    # No NonNative interactions
    if NonNativeType == -1: P_nn = None

    # NonNative WCA with sigma = shortest non-native contact distance
    if NonNativeType == 0:
        if NonNativeSigma is None: NonNativeSigma = d_nonnative.min()
        if NonNativeCut is None: NonNativeCut = NonNativeSigma * (2 ** (1/6.))
        if NonNativeEpsilon is None:
            if not NativeEpsilon is None: NonNativeEpsilon = NativeEpsilonT
        print 'Using WCA Non-Native interactions with sigma = %gA, Epsilon = %g kT, Cutoff = %gA' % (NonNativeSigma, NonNativeEpsilon / 0.6, NonNativeCut)
        P_nn = sim.potential.LJ(Sys, Filter = Filter_NonNative, Cut = NonNativeCut,
                                Sigma = NonNativeSigma, Epsilon = NonNativeEpsilon,
                                Shift = True, Label = 'NonBondNonNative')

    # NonNative WCA with sigma = native sigma
    if NonNativeType == 1:
        if NonNativeSigma is None:
            if not NativeSigma is None: NonNativeSigma = NativeSigma
            else: NonNativeSigma = d_native.min() * (2**(-1/6.))
        if NonNativeCut is None: NonNativeCut = NonNativeSigma * (2**(1/6.))
        if NonNativeEpsilon is None:
            if not NativeEpsilon is None: NonNativeEpsilon = NativeEpsilon
        print 'Using WCA Non-Native interactions with sigma = %gA, Epsilon = %g kT, Cutoff = %gA' % (NonNativeSigma, NonNativeEpsilon / 0.6, NonNativeCut)
        P_nn = sim.potential.LJ(Sys, Filter = Filter_NonNative, Cut = NonNativeCut,
                                Sigma = NonNativeSigma, Epsilon = NonNativeEpsilon,
                                Shift = True, Label = 'NonBondNonNative')

    # NonNative spline with cutoff commensurate with a WCA corresponding to the native LJ potential
    # with sigma that attributes attractive interactions to the closest native-contact pair
    if NonNativeType == 2:
        if NonNativeCut is None:
            if not NativeSigma is None: NonNativeCut = NativeSigma * (2**(1/6.))
            else: NonNativeCut = d_native.min()
        print 'Using Spline Non-Native interactions with Cutoff = %gA' % NonNativeCut
        P_nn = sim.potential.PairSpline(Sys, Filter = Filter_NonNative, NKnot = NonNativeNKnot, Cut = NonNativeCut, Label = 'NonBondNonNative')
        if not NonNativeKnots is None: P_nn.SetParam(Knots = NonNativeKnots) 

    # add the potentials to the forcefield
    GoSet = []
    if not P_n is None:
        if not GoType == 4: GoSet.append(P_n)
        else: GoSet.extend(P_n)
    if not P_nn is None:
        GoSet.append(P_nn)
    return GoSet
    


####### SIM CG MODEL ########

# generate the Go system
def makeGoSys(hasSpecialBBParams = False, NChains = 1):
    global Prefix
    global Seq, NRes
    global BB_forcefield_file, Go_forcefield_file, GoSplineCut, MinBondOrd, LoadedPotentials
    global NKnot, SPCut
    global Temps, TempSet, BoxL
    global NCOSType
    
    # read Seq info
    parseNative()
        
    print 'System Temp: %g' % TempSet
    print 'REMD schedule: %gK - %gK, %d replicas' % (Temps[0], Temps[-1], len(Temps))
    
    ## create atomtypes
    AtomN = sim.chem.AtomType('N', Mass = 14.0)
    AtomC = sim.chem.AtomType('C', Mass = 12.0)
    if hasSpecialBBParams:
        if Seq.__contains__('GLY'): AtomC_GLY = sim.chem.AtomType('C', Mass = 12.0)
        if Seq.__contains__('PRO'): AtomC_PRO = sim.chem.AtomType('C', Mass = 12.0)
    AtomO = sim.chem.AtomType('O', Mass = 28.0)
    AtomS = []
    for i,r in enumerate(Seq):
        if not r == 'GLY':
            thisMass = ResMass[r] - ResMass['BB']
            thisName = 'S%d' % i
            thisAtomS = sim.chem.AtomType(thisName, Mass = thisMass)
            AtomS.append(thisAtomS)
        else: AtomS.append(None)

    # create molecule
    ResList = []
    AtomList = []
    for i,r in enumerate(Seq):
        if hasSpecialBBParams:
            if r == 'GLY': res = [AtomN, AtomC_GLY, AtomO, AtomS[i]]
            elif r == 'PRO': res = [AtomN, AtomC_PRO, AtomO, AtomS[i]]
            else: res = [AtomN, AtomC, AtomO, AtomS[i]]
        else:
            res = [AtomN, AtomC, AtomO, AtomS[i]]
        if r == 'GLY': res.remove(res[-1])
        ResList.append(res)
        AtomList.extend(res)
    Mol = sim.chem.MolType(Prefix, [x for x in AtomList])
    
    # create bonds
    for i, r in enumerate(Seq):
        # intra-res bonds
        resbondpairs = [(0,1), (1,2)] if r == 'GLY' else [(0,1), (1,2), (1,3)]
        for b in resbondpairs:
            m = res2atom(i, b[0])
            n = res2atom(i, b[1])
            Mol.Bond(m,n)
        # peptide bond
        if i < len(Seq) - 1:
            m = res2atom(i, 2)
            n = res2atom(i+1,0)
            Mol.Bond(m,n)

    # create System
    World = sim.chem.World([Mol], Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(World, Name = Prefix)
    for i in range(NChains): Sys += Mol.New()

    # create backbone potential filters
    FilterS = sim.atomselect.Filter([x for x in AtomS if not x is None]) 

    Filter_NC = sim.atomselect.PolyFilter([AtomN, AtomC], Bonded = True)
    Filter_CO = sim.atomselect.PolyFilter([AtomC, AtomO], Bonded = True)
    Filter_ON = sim.atomselect.PolyFilter([AtomO, AtomN], Bonded = True)
    Filter_CS = sim.atomselect.PolyFilter([AtomC, FilterS], Bonded = True)
    
    Filter_NCO = sim.atomselect.PolyFilter([AtomN, AtomC, AtomO], Bonded = True)
    Filter_CON = sim.atomselect.PolyFilter([AtomC, AtomO, AtomN], Bonded = True)
    Filter_ONC = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC], Bonded = True)
    Filter_NCS = sim.atomselect.PolyFilter([AtomN, AtomC, FilterS], Bonded = True)
    Filter_SCO = sim.atomselect.PolyFilter([FilterS, AtomC, AtomO], Bonded = True)
    
    Filter_NCON = sim.atomselect.PolyFilter([AtomN, AtomC, AtomO, AtomN], Bonded = True)
    Filter_CONC = sim.atomselect.PolyFilter([AtomC, AtomO, AtomN, AtomC], Bonded = True)
    Filter_ONCO = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC, AtomO], Bonded = True)
    Filter_ONCS = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC, FilterS], Bonded = True)
    Filter_SCON = sim.atomselect.PolyFilter([FilterS, AtomC, AtomO, AtomN], Bonded = True)
    if hasSpecialBBParams:
        if Seq.__contains__('GLY'):
            Filter_NCON_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
            Filter_CONC_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
            Filter_ONCO_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
        if Seq.__contains__('PRO'):
            Filter_NCON_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
            Filter_CONC_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
            Filter_ONCO_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
    
    Filter_NN_p = sim.atomselect.PolyFilter([AtomN, AtomN], MinBondOrd = MinBondOrd)
    Filter_CC_p = sim.atomselect.PolyFilter([AtomC, AtomC], MinBondOrd = MinBondOrd)
    Filter_OO_p = sim.atomselect.PolyFilter([AtomO, AtomO], MinBondOrd = MinBondOrd)
    Filter_NC_p = sim.atomselect.PolyFilter([AtomN, AtomC], MinBondOrd = MinBondOrd)
    Filter_NO_p = sim.atomselect.PolyFilter([AtomN, AtomO], MinBondOrd = MinBondOrd)
    Filter_CO_p = sim.atomselect.PolyFilter([AtomC, AtomO], MinBondOrd = MinBondOrd)
    Filter_NS_p = sim.atomselect.PolyFilter([AtomN, FilterS], MinBondOrd = MinBondOrd)
    Filter_CS_p = sim.atomselect.PolyFilter([AtomC, FilterS], MinBondOrd = MinBondOrd)
    Filter_OS_p = sim.atomselect.PolyFilter([AtomO, FilterS], MinBondOrd = MinBondOrd)
    
    # create backbone and backbone-sidechain potentials
    Bond_NC = sim.potential.Bond(Sys, Filter = Filter_NC, Label = 'BondNC', Dist0 = 4.0, FConst = 1.0)
    Bond_CO = sim.potential.Bond(Sys, Filter = Filter_CO, Label = 'BondCO', Dist0 = 4.0, FConst = 1.0)
    Bond_ON = sim.potential.Bond(Sys, Filter = Filter_ON, Label = 'BondON', Dist0 = 4.0, FConst = 1.0)
    Bond_CS = sim.potential.Bond(Sys, Filter = Filter_CS, Label = 'BondCS', Dist0 = 4.0, FConst = 1.0)
    
    Angle_NCO = sim.potential.AngleSpline(Sys, Filter = Filter_NCO, Label = 'AngleNCO', NKnot = NKnot)
    Angle_CON = sim.potential.AngleSpline(Sys, Filter = Filter_CON, Label = 'AngleCON', NKnot = NKnot)
    Angle_ONC = sim.potential.AngleSpline(Sys, Filter = Filter_ONC, Label = 'AngleONC', NKnot = NKnot)
    Angle_NCS = sim.potential.AngleSpline(Sys, Filter = Filter_NCS, Label = 'AngleNCS', NKnot = NKnot)
    Angle_SCO = sim.potential.AngleSpline(Sys, Filter = Filter_SCO, Label = 'AngleSCO', NKnot = NKnot)
    
    Torsion_NCON = sim.potential.TorsionSpline(Sys, Filter = Filter_NCON, Label = 'TorsionNCON', NKnot = NKnot)
    Torsion_CONC = sim.potential.TorsionSpline(Sys, Filter = Filter_CONC, Label = 'TorsionCONC', NKnot = NKnot)
    Torsion_ONCO = sim.potential.TorsionSpline(Sys, Filter = Filter_ONCO, Label = 'TorsionONCO', NKnot = NKnot)
    Torsion_ONCS = sim.potential.TorsionSpline(Sys, Filter = Filter_ONCS, Label = 'TorsionONCS', NKnot = NKnot)
    Torsion_SCON = sim.potential.TorsionSpline(Sys, Filter = Filter_SCON, Label = 'TorsionSCON', NKnot = NKnot)
    if hasSpecialBBParams:
        if Seq.__contains__('GLY'):
            Torsion_NCON_GLY = sim.potential.TorsionSpline(Sys, Filter = Filter_NCON_GLY, Label = 'TorsionNCON_GLY', NKnot = NKnot)
            Torsion_CONC_GLY = sim.potential.TorsionSpline(Sys, Filter = Filter_CONC_GLY, Label = 'TorsionCONC_GLY', NKnot = NKnot)
            Torsion_ONCO_GLY = sim.potential.TorsionSpline(Sys, Filter = Filter_ONCO_GLY, Label = 'TorsionONCO_GLY', NKnot = NKnot)
        if Seq.__contains__('PRO'):
            Torsion_NCON_PRO = sim.potential.TorsionSpline(Sys, Filter = Filter_NCON_PRO, Label = 'TorsionNCON_PRO', NKnot = NKnot)
            Torsion_CONC_PRO = sim.potential.TorsionSpline(Sys, Filter = Filter_CONC_PRO, Label = 'TorsionCONC_PRO', NKnot = NKnot)
            Torsion_ONCO_PRO = sim.potential.TorsionSpline(Sys, Filter = Filter_ONCO_PRO, Label = 'TorsionONCO_PRO', NKnot = NKnot)
    
    Pair_NN = sim.potential.PairSpline(Sys, Filter = Filter_NN_p, Label = 'NonBondNN', NKnot = NKnot, Cut = SPCut)
    Pair_CC = sim.potential.PairSpline(Sys, Filter = Filter_CC_p, Label = 'NonBondCC', NKnot = NKnot, Cut = SPCut)
    Pair_OO = sim.potential.PairSpline(Sys, Filter = Filter_OO_p, Label = 'NonBondOO', NKnot = NKnot, Cut = SPCut)
    Pair_NC = sim.potential.PairSpline(Sys, Filter = Filter_NC_p, Label = 'NonBondNC', NKnot = NKnot, Cut = SPCut)
    Pair_NO = sim.potential.PairSpline(Sys, Filter = Filter_NO_p, Label = 'NonBondNO', NKnot = NKnot, Cut = SPCut)
    Pair_CO = sim.potential.PairSpline(Sys, Filter = Filter_CO_p, Label = 'NonBondCO', NKnot = NKnot, Cut = SPCut)
    Pair_NS = sim.potential.PairSpline(Sys, Filter = Filter_NS_p, Label = 'NonBondNS', NKnot = NKnot, Cut = SPCut)
    Pair_CS = sim.potential.PairSpline(Sys, Filter = Filter_CS_p, Label = 'NonBondCS', NKnot = NKnot, Cut = SPCut)
    Pair_OS = sim.potential.PairSpline(Sys, Filter = Filter_OS_p, Label = 'NonBondOS', NKnot = NKnot, Cut = SPCut)
    
    # create side chain-side chain Go potentials
    GoSet = genGoPotentials(Sys)
    
    # populate backbone forcefield
    P_backbone = [Bond_NC, Bond_CO, Bond_ON,
                  Angle_NCO, Angle_CON, Angle_ONC,
                  Torsion_NCON, Torsion_CONC, Torsion_ONCO,
                  Pair_NN, Pair_CC, Pair_OO, Pair_NC, Pair_NO, Pair_CO]
    
    # special backbone parameters for GLY and PRO
    P_backbone_special = []              
    if hasSpecialBBParams:
        if Seq.__contains__('GLY'): P_backbone_special.extend([Torsion_NCON_GLY, Torsion_CONC_GLY, Torsion_ONCO_GLY])
        if Seq.__contains__('PRO'): P_backbone_special.extend([Torsion_NCON_PRO, Torsion_CONC_PRO, Torsion_ONCO_PRO])
    
    # choose how to populate backbone-sidechain bonded forcefield
    P_sidechain = []
    # GoType = -1 is a random coil with no side chain potentials at all
    if not GoType == -1: P_sidechain += [Bond_CS, Angle_NCS, Angle_SCO, Torsion_ONCS, Torsion_SCON] 
    
    # choose how to populate backbone-sidechain nonbonded forcefield
    if NCOSType == 0:
        # NCOSType = 0 includes nonbonded backbone-sidechain interactions
        P_sidechain += [Pair_NS, Pair_CS, Pair_OS]
    
    # populate sidechain-sidechain forcefield
    P_sidechain += GoSet
    
    # populate overall forcefield
    Sys.ForceField.extend(P_backbone)
    Sys.ForceField.extend(P_backbone_special)
    Sys.ForceField.extend(P_sidechain)
                                       
    # load backbone if available
    if not BB_forcefield_file is None: LoadedPotentials = pp.loadParam(Sys, BB_forcefield_file)
   
    print Sys.ForceField.ParamString()
    for P in Sys.ForceField:
        if P.Name.startswith('NonBond'):
            print P.Name, 'Cutoff = %2.2f' % P.Cut

    # set up histograms
    for P in Sys.ForceField:
        P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)
    # load the system
    Sys.Load()
    
    # set integrator properties
    Int = Sys.Int
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    
    # set ambient conditions
    Sys.Temp = TempSet
    Sys.TempSet = TempSet
    Sys.BoxL = BoxL
    
    return Sys
    

# generate the NCOS system
def makeNCOSSys(hasSpecialBBParams = False):
    global Prefix
    global Seq, NRes, ResTypes, NResTypes
    global NCOS_forcefield_file, MinBondOrd, LoadedPotentials
    global NKnot, SPCut
    global Temps, TempSet, BoxL
    global NCOSType
    global ExType
    
    print '\n\nSystem Temp: %g' % TempSet
    print 'REMD schedule: %gK - %gK, %d replicas' % (Temps[0], Temps[-1], len(Temps))
    if NCOSType == 0: print 'Retaining N-S, C-S, O-S nonbonded interactions'
    if NCOSType == 1: print 'Ignoring N-S, C-S, O-S nonbonded interactions'
    print 'Generating NCOS Model'
    
    # read Seq info
    parseNative()
    
    ## create atomtypes
    AtomN = sim.chem.AtomType('N', Mass = 14.0)
    AtomC = sim.chem.AtomType('C', Mass = 12.0)
    if hasSpecialBBParams:
        if Seq.__contains__('GLY'): AtomC_GLY = sim.chem.AtomType('C', Mass = 12.0)
        if Seq.__contains__('PRO'): AtomC_PRO = sim.chem.AtomType('C', Mass = 12.0)
    AtomO = sim.chem.AtomType('O', Mass = 28.0)
    AtomS = {}
    for r in Seq:
        if not r == 'GLY':
            thisMass = ResMass[r] - ResMass['BB']
            thisName = 'S_%s' % r
            thisAtomS = sim.chem.AtomType(thisName, Mass = thisMass)
            AtomS[r] = thisAtomS

    # create molecule
    ResList = []
    AtomList = []
    for r in Seq:
        if hasSpecialBBParams:
            if r == 'GLY': res = [AtomN, AtomC_GLY, AtomO, AtomS[r]]
            elif r == 'PRO': res = [AtomN, AtomC_PRO, AtomO, AtomS[r]]
            else: res = [AtomN, AtomC, AtomO, AtomS[r]]
        else:
            res = [AtomN, AtomC, AtomO, AtomS[r]]
        if r == 'GLY': res.remove(res[-1])
        ResList.append(res)
        AtomList.extend(res)
    Mol = sim.chem.MolType(Prefix, [x for x in AtomList])
    
    # create bonds
    for i, r in enumerate(Seq):
        # intra-res bonds
        resbondpairs = [(0,1), (1,2)] if r == 'GLY' else [(0,1), (1,2), (1,3)]
        for b in resbondpairs:
            m = res2atom(i, b[0])
            n = res2atom(i, b[1])
            Mol.Bond(m,n)
        # peptide bond
        if i < len(Seq) - 1:
            m = res2atom(i, 2)
            n = res2atom(i+1,0)
            Mol.Bond(m,n)

    # create System
    World = sim.chem.World([Mol], Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(World, Name = Prefix)
    Sys += Mol.New()

    # create backbone potential filters
    Filter_NC = sim.atomselect.PolyFilter([AtomN, AtomC], Bonded = True)
    Filter_CO = sim.atomselect.PolyFilter([AtomC, AtomO], Bonded = True)
    Filter_ON = sim.atomselect.PolyFilter([AtomO, AtomN], Bonded = True)
    Filter_CS = dict( (r, sim.atomselect.PolyFilter([AtomC, AtomS[r]], Bonded = True)) for r in ResTypes if not r == 'GLY')
    
    Filter_NCO = sim.atomselect.PolyFilter([AtomN, AtomC, AtomO], Bonded = True)
    Filter_CON = sim.atomselect.PolyFilter([AtomC, AtomO, AtomN], Bonded = True)
    Filter_ONC = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC], Bonded = True)
    Filter_NCS = dict( (r, sim.atomselect.PolyFilter([AtomN, AtomC, AtomS[r]], Bonded = True)) for r in ResTypes if not r == 'GLY')
    Filter_SCO = dict( (r, sim.atomselect.PolyFilter([AtomS[r], AtomC, AtomO], Bonded = True)) for r in ResTypes if not r == 'GLY')
    
    Filter_NCON = sim.atomselect.PolyFilter([AtomN, AtomC, AtomO, AtomN], Bonded = True)
    Filter_CONC = sim.atomselect.PolyFilter([AtomC, AtomO, AtomN, AtomC], Bonded = True)
    Filter_ONCO = sim.atomselect.PolyFilter([AtomO, AtomN, AtomC, AtomO], Bonded = True)
    Filter_ONCS = dict( (r, sim.atomselect.PolyFilter([AtomO, AtomN, AtomC, AtomS[r]], Bonded = True)) for r in ResTypes if not r == 'GLY')
    Filter_SCON = dict( (r, sim.atomselect.PolyFilter([AtomS[r], AtomC, AtomO, AtomN], Bonded = True)) for r in ResTypes if not r == 'GLY')
    if hasSpecialBBParams:
        if Seq.__contains__('GLY'):
            Filter_NCON_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
            Filter_CONC_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
            Filter_ONCO_GLY = sim.atomselect.PolyFilter([AtomN, AtomC_GLY, AtomO, AtomN], Bonded = True)
        if Seq.__contains__('PRO'):
            Filter_NCON_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
            Filter_CONC_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
            Filter_ONCO_PRO = sim.atomselect.PolyFilter([AtomN, AtomC_PRO, AtomO, AtomN], Bonded = True)
    
    Filter_NN_p = sim.atomselect.PolyFilter([AtomN, AtomN], MinBondOrd = MinBondOrd)
    Filter_CC_p = sim.atomselect.PolyFilter([AtomC, AtomC], MinBondOrd = MinBondOrd)
    Filter_OO_p = sim.atomselect.PolyFilter([AtomO, AtomO], MinBondOrd = MinBondOrd)
    Filter_NC_p = sim.atomselect.PolyFilter([AtomN, AtomC], MinBondOrd = MinBondOrd)
    Filter_NO_p = sim.atomselect.PolyFilter([AtomN, AtomO], MinBondOrd = MinBondOrd)
    Filter_CO_p = sim.atomselect.PolyFilter([AtomC, AtomO], MinBondOrd = MinBondOrd)
    Filter_NS_p = dict( (r, sim.atomselect.PolyFilter([AtomN, AtomS[r]], MinBondOrd = MinBondOrd)) for r in ResTypes if not r == 'GLY')
    Filter_CS_p = dict( (r, sim.atomselect.PolyFilter([AtomC, AtomS[r]], MinBondOrd = MinBondOrd)) for r in ResTypes if not r == 'GLY')
    Filter_OS_p = dict( (r, sim.atomselect.PolyFilter([AtomO, AtomS[r]], MinBondOrd = MinBondOrd)) for r in ResTypes if not r == 'GLY')
    
    # create side chain pair potential filters
    Filter_SS_p = {}
    for i in range(NResTypes):
        for j in range(NResTypes):
            if i > j: continue
            r1 = ResTypes[i]
            r2 = ResTypes[j]
            if (r1 == 'GLY') or (r2 == 'GLY'): continue
            Filter_SS_p[(r1, r2)] = sim.atomselect.PolyFilter([AtomS[r1], AtomS[r2]], MinBondOrd = MinBondOrd)
       
    # create backbone and backbone-sidechain potentials
    Bond_NC = sim.potential.Bond(Sys, Filter = Filter_NC, Label = 'BondNC', Dist0 = 4.0, FConst = 1.0)
    Bond_CO = sim.potential.Bond(Sys, Filter = Filter_CO, Label = 'BondCO', Dist0 = 4.0, FConst = 1.0)
    Bond_ON = sim.potential.Bond(Sys, Filter = Filter_ON, Label = 'BondON', Dist0 = 4.0, FConst = 1.0)
    Bond_CS = dict( (r, sim.potential.Bond(Sys, Filter = Filter_CS[r], Label = 'BondCS_%s' % r, Dist0 = 4.0, FConst = 1.0)) for r in ResTypes if not r == 'GLY')
    
    Angle_NCO = sim.potential.AngleSpline(Sys, Filter = Filter_NCO, Label = 'AngleNCO', NKnot = NKnot)
    Angle_CON = sim.potential.AngleSpline(Sys, Filter = Filter_CON, Label = 'AngleCON', NKnot = NKnot)
    Angle_ONC = sim.potential.AngleSpline(Sys, Filter = Filter_ONC, Label = 'AngleONC', NKnot = NKnot)
    Angle_NCS = dict( (r, sim.potential.AngleSpline(Sys, Filter = Filter_NCS[r], Label = 'AngleNCS_%s' % r, NKnot = NKnot)) for r in ResTypes if not r == 'GLY')
    Angle_SCO = dict( (r, sim.potential.AngleSpline(Sys, Filter = Filter_SCO[r], Label = 'AngleSCO_%s' % r, NKnot = NKnot)) for r in ResTypes if not r == 'GLY')
    
    Torsion_NCON = sim.potential.TorsionSpline(Sys, Filter = Filter_NCON, Label = 'TorsionNCON', NKnot = NKnot)
    Torsion_CONC = sim.potential.TorsionSpline(Sys, Filter = Filter_CONC, Label = 'TorsionCONC', NKnot = NKnot)
    Torsion_ONCO = sim.potential.TorsionSpline(Sys, Filter = Filter_ONCO, Label = 'TorsionONCO', NKnot = NKnot)
    Torsion_ONCS = dict( (r, sim.potential.TorsionSpline(Sys, Filter = Filter_ONCS[r], Label = 'TorsionONCS_%s' % r, NKnot = NKnot)) for r in ResTypes if not r == 'GLY')
    Torsion_SCON = dict( (r, sim.potential.TorsionSpline(Sys, Filter = Filter_SCON[r], Label = 'TorsionSCON_%s' % r, NKnot = NKnot)) for r in ResTypes if not r == 'GLY')
    if hasSpecialBBParams:
        if Seq.__contains__('GLY'):
            Torsion_NCON_GLY = sim.potential.TorsionSpline(Sys, Filter = Filter_NCON_GLY, Label = 'TorsionNCON_GLY', NKnot = NKnot)
            Torsion_CONC_GLY = sim.potential.TorsionSpline(Sys, Filter = Filter_CONC_GLY, Label = 'TorsionCONC_GLY', NKnot = NKnot)
            Torsion_ONCO_GLY = sim.potential.TorsionSpline(Sys, Filter = Filter_ONCO_GLY, Label = 'TorsionONCO_GLY', NKnot = NKnot)
        if Seq.__contains__('PRO'):
            Torsion_NCON_PRO = sim.potential.TorsionSpline(Sys, Filter = Filter_NCON_PRO, Label = 'TorsionNCON_PRO', NKnot = NKnot)
            Torsion_CONC_PRO = sim.potential.TorsionSpline(Sys, Filter = Filter_CONC_PRO, Label = 'TorsionCONC_PRO', NKnot = NKnot)
            Torsion_ONCO_PRO = sim.potential.TorsionSpline(Sys, Filter = Filter_ONCO_PRO, Label = 'TorsionONCO_PRO', NKnot = NKnot)
    
    Pair_NN = sim.potential.PairSpline(Sys, Filter = Filter_NN_p, Label = 'NonBondNN', NKnot = NKnot, Cut = SPCut)
    Pair_CC = sim.potential.PairSpline(Sys, Filter = Filter_CC_p, Label = 'NonBondCC', NKnot = NKnot, Cut = SPCut)
    Pair_OO = sim.potential.PairSpline(Sys, Filter = Filter_OO_p, Label = 'NonBondOO', NKnot = NKnot, Cut = SPCut)
    Pair_NC = sim.potential.PairSpline(Sys, Filter = Filter_NC_p, Label = 'NonBondNC', NKnot = NKnot, Cut = SPCut)
    Pair_NO = sim.potential.PairSpline(Sys, Filter = Filter_NO_p, Label = 'NonBondNO', NKnot = NKnot, Cut = SPCut)
    Pair_CO = sim.potential.PairSpline(Sys, Filter = Filter_CO_p, Label = 'NonBondCO', NKnot = NKnot, Cut = SPCut)
    Pair_NS = dict( (r, sim.potential.PairSpline(Sys, Filter = Filter_NS_p[r], Label = 'NonBondNS_%s' % r, NKnot = NKnot, Cut = SPCut)) for r in ResTypes if not r == 'GLY')
    Pair_CS = dict( (r, sim.potential.PairSpline(Sys, Filter = Filter_CS_p[r], Label = 'NonBondCS_%s' % r, NKnot = NKnot, Cut = SPCut)) for r in ResTypes if not r == 'GLY')
    Pair_OS = dict( (r, sim.potential.PairSpline(Sys, Filter = Filter_OS_p[r], Label = 'NonBondOS_%s' % r, NKnot = NKnot, Cut = SPCut)) for r in ResTypes if not r == 'GLY')
    
    # create sidechain-sidechain nonbonded potentials
    Pair_SS = {}
    for i in range(NResTypes):
        for j in range(NResTypes):
            if i > j: continue
            r1 = ResTypes[i]
            r2 = ResTypes[j]
            if (r1 == 'GLY') or (r2 == 'GLY'): continue
            Pair_SS[(r1, r2)] = sim.potential.PairSpline(Sys, Filter = Filter_SS_p[(r1,r2)], Label = 'NonBondS_%s_S_%s' % (r1,r2), Cut = SPCut)
    
    # populate backbone forcefield
    P_backbone = [Bond_NC, Bond_CO, Bond_ON,
                  Angle_NCO, Angle_CON, Angle_ONC,
                  Torsion_NCON, Torsion_CONC, Torsion_ONCO,
                  Pair_NN, Pair_CC, Pair_OO, Pair_NC, Pair_NO, Pair_CO]
    
    # special backbone parameters for GLY and PRO
    P_backbone_special = []              
    if hasSpecialBBParams:
        if Seq.__contains__('GLY'): P_backbone_special.extend([Torsion_NCON_GLY, Torsion_CONC_GLY, Torsion_ONCO_GLY])
        if Seq.__contains__('PRO'): P_backbone_special.extend([Torsion_NCON_PRO, Torsion_CONC_PRO, Torsion_ONCO_PRO])
    
    # populate backbone-sidechain bonded forcefield
    P_sidechain = Bond_CS.values() + \
                  Angle_NCS.values() + Angle_SCO.values() + \
                  Torsion_ONCS.values() + Torsion_SCON.values()
                  
    # choose how to populate backbone-sidechain nonbonded forcefield
    if NCOSType == 0:
        # NCOSType = 0 includes nonbonded backbone-sidechain interactions
        P_sidechain += Pair_NS.values() + Pair_CS.values() + Pair_OS.values()
    
    # populate sidechain-sidechain nonbonded forcefield
    P_sidechain += Pair_SS.values()
            
    # populate overall forcefield
    Sys.ForceField.extend(P_backbone)
    Sys.ForceField.extend(P_backbone_special)
    Sys.ForceField.extend(P_sidechain)
                                                 
    # load supplied forcefield 
    if not NCOS_forcefield_file is None:
        LoadedPotentials = pp.loadParam(Sys, NCOS_forcefield_file)
    
    # set up histograms
    for P in Sys.ForceField:
        P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)
    
    # load the system
    Sys.Load()
    
    # set integrator properties
    Int = Sys.Int
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    Int.Method.LangevinGamma = 0.01
    
    # set ambient conditions
    Sys.Temp = TempSet
    Sys.TempSet = TempSet
    Sys.BoxL = 0.0
    
    return Sys



######## REMD ########

# wrapper on lammps_REMD        
def runREMD(Sys):
    global OutputDir, Prefix
    global NStepsMin, NStepsEquil, NStepsProd, NStepsSwap, StepFreq, Temps, TempSet

    # random initial structure
    sim.export.lammps.LammpsExec = LammpsExec
    Sys.Arrays.Pos = genRandInitPos()
    sim.system.init.velocities.Canonical(Sys, Temp = Sys.TempSet)
 
    # record temp schedule
    np.savetxt(os.path.join(OutputDir, 'temps.txt'), Temps)

    # feed in all the settings to the Lammps export
    sim.export.lammps_REMD.NStepsSwap = NStepsSwap
    sim.export.lammps_REMD.TEMPS = Temps
    LammpsFiles, TrajFile = sim.export.lammps_REMD.MakeLammpsReplicaMD(Sys, Prefix = Prefix, TrajFile = '.lammpstrj',
                                                                       NStepsMin = NStepsMin, NStepsEquil = NStepsEquil, 
                                                                       NStepsProd = NStepsProd, WriteFreq = StepFreq)
    # run Lammps REMD
    InFile, DataFile, TableFile, DihedralFile = LammpsFiles
    LogFile, ScreenFile, returncode = sim.export.lammps_REMD.RunLammpsReplica(InFile, Prefix = Prefix, Verbose = True)
    return TrajFile, LogFile


# reordering trajectories by temp
def reorderTraj(ReorderTemps = None):
    global Prefix, OutputDir
    global NStepsEquil, NStepsProd, NStepsSwap, StepFreq, Temps
    Temps = np.loadtxt(os.path.join(OutputDir, 'temps.txt'))
    if ReorderTemps is None: ReorderTemps = Temps
    TrajFile = os.path.join(OutputDir, Prefix + '.lammpstrj')
    LogFile = os.path.join(OutputDir, Prefix + 'lammps.log')
    
    for T in ReorderTemps:
        TempInd = [list(Temps).index(t) for t in list(Temps) if '%3.2f' % t  == '%3.2f' % T][0]
        MultiTrajFn = FMT['TRAJ'] % (Prefix, T)
        MultiEneFn = FMT['ENE'] % (Prefix, T)
        if os.path.isfile(MultiTrajFn) and os.path.isfile(MultiEneFn): continue
        print 'Collecting replica frames at temperature ', T
        RepIndsMaster = [np.where(x[1:] == TempInd)[0][0] for x in np.loadtxt(LogFile, skiprows = 3)]
        RepInds = RepIndsMaster[int(NStepsEquil/NStepsSwap) : ]
        RepInds = RepInds[:-1]
        this_Traj = {}; this_Ene = {}
        TrajList = [] ; EneList = []
       
        if StepFreq <= NStepsSwap: # assume mod(NStepsSwap, StepFreq) = 0
            for ii, i in enumerate(RepInds):
                if not i in this_Traj.keys():
                    thisTrajFn = '%s.%d.gz' % (TrajFile, i)
                    thisLogFn = '%s.%d' % (LogFile, i)
                    this_Traj[i] = pickleTraj(thisTrajFn, LogFile = thisLogFn, LogFileToken = '#run production') 
                    this_Ene[i] = this_Traj[i].ThermoDict['PEnergy']
                this_Traj_slice = this_Traj[i][int(NStepsEquil/StepFreq) : ]
                this_Ene_slice = this_Ene[i][int(NStepsEquil/StepFreq) : ]    
                start = ii * NStepsSwap / StepFreq ; stop = (ii+1) * NStepsSwap / StepFreq
                TrajList.append(this_Traj_slice[start:stop])
                EneList.extend(this_Ene_slice[start:stop])
        
        else:
            NSkip = StepFreq / NStepsSwap # assume mod(StepFreq, NStepsSwap) = 0
            for ii, i in enumerate(RepInds[0::NSkip]):
                if not i in this_Traj.keys():
                    thisTrajFn = '%s.%d.gz' % (TrajFile, i)
                    thisLogFn = '%s.%d' % (LogFile, i)
                    this_Traj[i] = pickleTraj(thisTrajFn, LogFile = thisLogFn, LogFileToken = '#run production') 
                    this_Ene[i] = this_Traj[i].ThermoDict['PEnergy']
                this_Traj_slice = this_Traj[i][int(NStepsEquil/StepFreq) : ]
                this_Ene_slice = this_Ene[i][int(NStepsEquil/StepFreq) : ]    
                TrajList.append(this_Traj_slice[ii:ii+1])
                EneList.extend(this_Ene_slice[ii:ii+1])
    
        MultiTraj = sim.traj.Multi(TrajList)
        sim.traj.base.Convert(InTraj = MultiTraj, OutTrajClass = sim.traj.LammpsWrite, FileName = MultiTrajFn)
        np.savetxt(MultiEneFn, EneList, fmt = '%11.4e')
    

    
######## SREL OPTIMIZATION ########

def runSrel(Sys, AATraj):
    global OutputDir, Prefix
    global NStepsMin, NStepsEquil, NStepsProd, NStepsSwap, StepFreq, Temps, TempSet
    global OptStages, LoadedPotentials, PermaFrost, MaxIter
    global useSeed, RefPosInd
    
    # read in AA Traj
    Trj = pickleTraj(os.path.abspath(AATraj))
    
    # initial conditions (choose to seed with native pos or not)
    # (native pos = RefPosInd-th frame of supplied AA Traj)
    if useSeed is None: RefPosInd = None
    sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
    
    # treat inner core (don't do this for Go optimizations)
    if not OptStages == ['Go']:
        print '\nTreating inner core for nonbonded spline potentials...'
        for P in Sys.ForceField:
            if P.Name.startswith("NonBond") and P.Names.__contains__('spline'):
                P.EneInner = "20kT"
                P.EneSlopeInner = None
    
    # 1-1 map
    Map = sim.atommap.PosMap()
    for (i, a) in enumerate(Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]
    
    # optimizer initialization
    print '\nStarting Srel optimization...'
    Opt = sim.srel.OptimizeTrajClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix, Verbose = True, RefPosInd = RefPosInd)
    Opt = sim.srel.UseLammps(Opt)
    Opt.TempFileDir = os.getcwd()
    Opt.MinReweightFrames = None # need this to work with smaller mod traj
    
    # feed in all the settings to the Lammps export
    sim.export.lammps.LammpsExec = LammpsExec
    sim.export.lammps_REMD.NStepsSwap = NStepsSwap
    sim.export.lammps_REMD.TEMPS = Temps
        
    # check if any potentials are permanently frozen
    isPermaFrost = lambda P : PermaFrost.__contains__(P.Name)
    PermaFrostList = [P.Name for P in Sys.ForceField if isPermaFrost(P)]
    if PermaFrostList:
        print '\nPermanently frozen potentials:-'
        for x in PermaFrostList: print x
        
    # stagewise optimization
    if OptStages.__contains__('Bond'):
        Opt.Iter = 0
        print '\nOptimizing Bonds...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name.startswith("Bond"):
                if not isPermaFrost(P): P.UnfreezeParam()
                if not P.Name in LoadedPotentials:
                    print 'Estimating bonded potential ', P.Name
                    P.Estimate() # estimate only bonded potentials that are not supplied
        if not useBondHessian: print 'Turning of Hessian Descent while optimizing bonds'
        Opt.UseHessian = useBondHessian # switch off hessian if bond distributions have multiple peaks
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['Bond'])
    
    # switch hessian back on
    Opt.UseHessian = True
                  
    if OptStages.__contains__('Angle'):
        Opt.Iter = 0
        print '\n\nOptimizing Angles...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name.startswith("Angle") and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['Angle'])
    
    if OptStages.__contains__('Torsion'):
        Opt.Iter = 0
        print '\n\nOptimizing Torsions...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name.startswith("Torsion") and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['Torsion'])    
    
    if OptStages.__contains__('NonBond'):
        Opt.Iter = 0
        print '\n\nOptimizing NonBonds...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name.startswith("NonBond") and not isPermaFrost(P) and not P.Name == 'NonBondNonNative':
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['NonBond'])
    
    # simultaneous optimization
    if OptStages.__contains__('all'):
        Opt.Iter = 0
        print '\n\nOptimizing all...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if not P.Name.startswith("Bond") and not isPermaFrost(P) and not P.Name == 'NonBondNonNative':
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['all'])
    
    # optimize only the Native Go potentials
    if OptStages.__contains__('Go1'):
        Opt.Iter = 0
        print '\n\nOptimizing Native Go...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name.startswith('NonBondNative') and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq)

    # optimize all Go potentials
    if OptStages.__contains__('Go2'):
        Opt.Iter = 0
        print '\n\nOptimizing Native and Non-Native Go...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if ( P.Name.startswith('NonBondNative') or P.Name.startswith('NonBondNonNative') ) and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq)


######## EXTENDED ENSEMBLE SREL OPTIMIZATION ########

def runSrelExtendedEnsemble_NCOS(Sys, AATrajList):
    global OutputDir, Prefix
    global NStepsMin, NStepsEquil, NStepsProd, NStepsSwap, StepFreq, Temps, TempSet
    global OptStages, LoadedPotentials, PermaFrost, MaxIter
    global useSeed, RefPosInd
    
    # read in AA Traj
    TrjList = []
    for i, AATraj in enumerate(AATrajList):
        TrjList.append(pickleTraj(os.path.abspath(AATraj)))
        
    # initial conditions (choose to seed with native pos or not)
    # (native pos = RefPosInd-th frame of supplied AA Traj)
    if useSeed is None: RefPosInd = None
    sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
    
    # treat inner core (don't do this for Go optimizations)
    if not OptStages == ['Go']:
        print '\nTreating inner core for nonbonded spline potentials...'
        for P in Sys.ForceField:
            if P.Name.startswith("NonBond") and P.Names.__contains__('spline'):
                P.EneInner = "20kT"
                P.EneSlopeInner = None
    
    # 1-1 map
    Map = sim.atommap.PosMap()
    for (i, a) in enumerate(Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]
    
    # optimizer ensemble initialization
    print '\nCreating ensemble of optimizer objects...'
    OptList = []
    for i, Trj in enumerate(TrjList):
        print '\n%d)  AA Traj: %s' % (i, AATrajList[i].split('/')[-1])
        this_Opt = sim.srel.OptimizeTrajClass(Sys, Map, Traj = Trj, RefPosInd = RefPosInd,
                                              SaveLoadArgData = True, FilePrefix = Prefix + '_%d' % i, 
                                              Verbose = True)
        
        this_Opt = sim.srel.UseLammps(this_Opt)
        this_Opt.TempFileDir = os.getcwd()
        this_Opt.MinReweightFrames = None # need this to work with smaller mod traj
        OptList.append(this_Opt)
    
    # multi-optimizer initialization
    print '\nStarting relative entropy optimization...'
    sim.srel.optimizemultitraj.doParallel = True
    Opt = sim.srel.OptimizeMultiTrajClass(OptList, FilePrefix = Prefix)
    
    # feed in all the settings to the Lammps export
    sim.export.lammps.LammpsExec = LammpsExec
    sim.export.lammps_REMD.NStepsSwap = NStepsSwap
    sim.export.lammps_REMD.TEMPS = Temps
        
    # check if any potentials are permanently frozen
    isPermaFrost = lambda P : PermaFrost.__contains__(P.Name)
    if PermaFrost:
        print '\nPermanently frozen potentials:-'
        for x in PermaFrost: print x
        
    # stagewise optimization
    if OptStages.__contains__('Bond'):
        Opt.Iter = 0
        print '\nOptimizing Bonds...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name.startswith("Bond"):
                if not isPermaFrost(P): P.UnfreezeParam()
                if not P.Name in LoadedPotentials:
                    print 'Estimating bonded potential ', P.Name
                    P.Estimate() # estimate only bonded potentials that are not supplied
        if not useBondHessian: print 'Turning of Hessian Descent while optimizing bonds'
        Opt.UseHessian = useBondHessian # switch off hessian if bond distributions have multiple peaks
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['Bond'])
    
    # switch hessian back on
    Opt.UseHessian = True
                  
    if OptStages.__contains__('Angle'):
        Opt.Iter = 0
        print '\n\nOptimizing Angles...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name.startswith("Angle") and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['Angle'])
    
    if OptStages.__contains__('Torsion'):
        Opt.Iter = 0
        print '\n\nOptimizing Torsions...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name.startswith("Torsion") and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['Torsion'])    
    
    if OptStages.__contains__('NonBond'):
        Opt.Iter = 0
        print '\n\nOptimizing NonBonds...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name.startswith("NonBond") and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['NonBond'])
    
    # simultaneous optimization
    if OptStages.__contains__('all'):
        Opt.Iter = 0
        print '\n\nOptimizing all...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if not P.Name.startswith("Bond") and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq, MaxIter = MaxIter['all'])
    
    # optimize only the Native Go potentials
    if OptStages.__contains__('Go1'):
        print '\n\nOptimizing Native Go...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if P.Name == 'NonBondNative' and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq)

    # optimize all Go potentials
    if OptStages.__contains__('Go2'):
        print '\n\nOptimizing Native and Non-Native Go...'
        for P in Sys.ForceField: P.FreezeParam()
        for P in Sys.ForceField:
            if ['NonBondNative', 'NonBondNonNative'].__contains__(P.Name) and not isPermaFrost(P):
                P.UnfreezeParam()
        Opt.RunConjugateGradient(StepsEquil = NStepsEquil, StepsProd = NStepsProd, StepsStride = StepFreq)

