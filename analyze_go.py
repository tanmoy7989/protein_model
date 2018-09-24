#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle, shelve
import cgprotein as cg, utils
import sim

# constants and flags, etc
FMT = utils.FMT
RoomTemp = 300.0
Recompute = True

# user input
NativePdb = os.path.abspath(sys.argv[1])
Prefix = sys.argv[2]
TrajDir = os.path.abspath(sys.argv[3])
OutDir = os.path.abspath(sys.argv[4])

hasPseudoGLY = int(sys.argv[5]) if len(sys.argv) > 5 else False

calcFoldTemp = bool(sys.argv[6]) if len(sys.argv) > 6 and sys.argv[6] else False
NStepsProd = int(sys.argv[7]) if len(sys.argv) > 7 else None
NStepsSwap = int(sys.argv[8]) if len(sys.argv) > 8 else None
WriteFreq = int(sys.argv[9]) if len(sys.argv) > 9 else None

# get tempset
TempFile = os.path.join(TrajDir, 'temps.txt')
Temps = np.loadtxt(TempFile)
TempSet = Temps[np.argmin(abs(Temps - RoomTemp))]

# get Traj
if not os.path.isdir(OutDir): os.system('mkdir -p %s' % OutDir)
OutPrefix = os.path.join(OutDir, Prefix)
TrajPrefix = os.path.join(TrajDir, Prefix)
TrajFn = FMT['TRAJ'] % (TrajPrefix, TempSet)

# output shelf
OutShelf = os.path.join(OutPrefix + '.shelf')
Keys = ['rmsd', 
        'ramaprob',
        'ramaerr', 
        'clustpdb', 
        'rescontacts', 
        'contactmap_traj', 'contactmap_topclust',
        'fracnativecontacts',
        'co',
        'contactdistcorr',
        'foldcurve'
       ] 

# if any of the following does not exist / has not been created, exit
Required = [NativePdb, TrajFn]
if not all([os.path.isfile(i) for i in Required]):
    print 'analysis IO Error: Insufficient starting information. Traj or native struct missing'
    exit()

# fancy intro
print '\nANALYSING GO MODEL RESULTS FOR %s at %3.2f K' % (Prefix, TempSet)
print '----------------------------------------------------------------'

# set up Compute object
print '\nCreating Compute object'
calc = cg.Compute(NativePdb = NativePdb, TrajFn = TrajFn, Temp = TempSet, Prefix = OutPrefix, hasPseudoGLY = hasPseudoGLY)

# calculate overall rmsd
def RMSD(CompType = 'traj'):
    oshelf = shelve.open(OutShelf)
    if CompType == 'traj':
        print '\nComputing overall RMSD distribution'
        rmsdhist = calc.QuickRMSD()
        oshelf['rmsd'] = rmsdhist
    
    if CompType == 'topclust':
        print '\nComputing RMSD from top cluster'
        pNative = cg.ProteinNCOS(NativePdb)
        pClust = cg.ProteinNCOS(oshelf['clustpdb'])
        rmsd = pClust.QuickRMSD(pNative)
        oshelf['rmsd_topclust'] = rmsd
    oshelf.close()
    return

# cluster the room temp traj
def Cluster():
    print '\nClustering trajectory'
    calc.Cluster()
    oshelf = shelve.open(OutShelf)
    oshelf['clustpdb'] = FMT['CLUSTPDB'] % (OutPrefix, TempSet)
    oshelf.close()
    return

# calculate phi psi errors
def PhiPsiErr(ErrType = 'traj'):
    print '\nCalculating Ramachandran plot errors with native struct'
    # get native phi, psi
    pNative = cg.ProteinNCOS(NativePdb)
    Phi_Native, Psi_Native = pNative.GetPhiPsi()
    Err = np.zeros(len(Phi_Native), float)
    
    # errors from top cluster
    if ErrType == 'topclust':
        ClustPdb = FMT['CLUSTPDB'] % (OutPrefix, TempSet)
        if not os.path.isfile(ClustPdb):
            raise IOError('Top cluster for %s missing' % Prefix)
        p = cg.ProteinNCOS(ClustPdb, Model = 1)
        PhiErr, PsiErr, Err = pNative.GetPhiPsiDiff(p)
    
    # errors from entire traj
    if ErrType == 'traj':
        calc.RamaChandran()
        RamaPickle = FMT['RAMA'] % (OutPrefix, TempSet)
        with open(RamaPickle, 'r') as of: data = pickle.load(of)
        Phi, Psi, hist = data
        for i in range(len(Phi_Native)):
            # per frame phi-psi errors
            this_Phidiff = sim.geom.NearestAngle(Phi_Native[i] - Phi[:, i], 0.0)
            this_Psidiff = sim.geom.NearestAngle(Psi_Native[i] - Psi[:, i], 0.0)
            Err[i] = np.mean( np.sqrt(this_Phidiff**2 + this_Psidiff**2) ) # need to do a np.mean over frames
    
    oshelf = shelve.open(OutShelf)
    oshelf['ramaprob'] = RamaPickle
    oshelf['ramaerr'] = (Phi_Native, Psi_Native, Err)
    oshelf.close()
    return   

def ContactMap(CompType = 'traj'):
    print '\nComparing Native vs. Predicted contact maps based on: ', CompType
    NativeContactMap, ClustContactMap = calc.CompareContactMap(CompType = CompType)
    oshelf = shelve.open(OutShelf)
    oshelf['rescontacts'] = FMT['RESCONTACTS'] % (OutPrefix, TempSet)
    oshelf['contactmap_%s' % CompType] = (NativeContactMap, ClustContactMap)
    oshelf.close()

def FracNativeContacts():
    print '\nCalculating fraction of native contacts'
    frachist = calc.GetFracNativeContacts()
    oshelf = shelve.open(OutShelf)
    oshelf['fracnativecontacts'] = frachist
    oshelf.close()

def ContactOrder():
    print '\nComparing native vs. Predicted  contact orders'
    nativeCO, clustCOhist = calc.CompareCO()
    oshelf = shelve.open(OutShelf)
    oshelf['co'] = (nativeCO, clustCOhist)
    oshelf.close()

def FoldCurve():
    print '\nComputing folding curves'
    print 'Creating replica object...'
    rep = cg.Replica(NativePdb = NativePdb, TrajPrefix = TrajPrefix, TempSet = 300.0,
                     OrderParams = ['U', 'RMSD'], Prefix = OutPrefix,
                     NStepsProd = NStepsProd, NStepsSwap = NStepsSwap,
                     WriteFreq = WriteFreq)
    rep.FoldCurve()
    oshelf = shelve.open(OutShelf)
    oshelf['foldcurve'] = FMT['FOLDCURVE'] % (OutPrefix, 'RMSD')
    oshelf.close()



#### MAIN ####
RMSD('traj')
Cluster()
PhiPsiErr('traj')
RMSD('topclust')
ContactMap('traj')
ContactMap('topclust')
FracNativeContacts()
ContactOrder()
if calcFoldTemp: FoldCurve()
