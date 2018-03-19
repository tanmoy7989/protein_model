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

# get tempset
TempFile = os.path.join(TrajDir, 'temps.txt')
Temps = np.loadtxt(TempFile)
TempSet = Temps[np.argmin(abs(Temps - RoomTemp))]

# get Traj
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
        'co'
        'contactdistcorr'
       ] 

# fancy intro
print 'ANALYSING GO MODEL RESULTS FOR %s at %3.2f K' % (Prefix, TempSet)
print '----------------------------------------------------------------'

# set up Compute object
print 'Creating Compute object'
calc = cg.Compute(NativePdb = NativePdb, TrajFn = TrajFn, Temp = TempSet, Prefix = OutPrefix)

# calculate overall rmsd
def RMSD():
    print 'Computing overall RMSD distribution'
    rmsdhist = calc.QuickRMSD()
    oshelf = shelve.open(OutShelf)
    oshelf['rmsd'] = rmsdhist
    oshelf.close()
    return

# cluster the room temp traj
def Cluster():
    print 'Clustering trajectory'
    calc.Cluster()
    oshelf = shelve.open(OutShelf)
    oshelf['clustpdb'] = FMT['CLUSTPDB'] % (OutPrefix, TempSet)
    oshelf.close()
    return

# calculate phi psi errors
def PhiPsiErr(ErrType = 'traj'):
    print 'Calculating Ramachandran plot errors with native struct'
    # get native phi, psi
    pNative = cg.ProteinNCOS(NativePdb)
    Phi_Native, Psi_Native = pNative.GetPhiPsi(RamaType = 'Generic')
    Phi_Native, Psi_Native = cg.TrimDihedrals(Phi_Native, Psi_Native, RamaType = 'Generic', ResNums = range(pNative.NRes) )
    Err = np.zeros([pNative.NRes - 1], float)
    
    # errors from top cluster
    if ErrType == 'topclust':
        ClustPdb = FMT['CLUSTPDB'] % (OutPrefix, TempSet)
        if not os.path.isfile(ClustPdb):
            raise IOError('Top cluster for %s missing' % Prefix)
        p = cg.ProteinNCOS(ClustPdb, Model = 1)
        PhiErr, PsiErr, Err = pNative.GetPhiPsiDiff(p, RamaType = 'Generic')
        # remove undefined dihedrals and recompute Err
        PhiErr = PhiErr[1:] ; PsiErr = PsiErr[:-1]
        Err = np.sqrt(PhiErr**2. + PsiErr**2.)
    
    # errors from entire traj
    if ErrType == 'traj':
        calc.RamaChandran()
        RamaPickle = FMT['RAMA'] % (OutPrefix, TempSet)
        with open(RamaPickle, 'r') as of: data = pickle.load(of)
        Phi, Psi, hist = data['Generic']
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
    print 'Comparing Native vs. Predicted contact maps based on: ', CompType
    NativeContactMap, ClustContactMap = calc.CompareContactMap(CompType = CompType)
    oshelf = shelve.open(OutShelf)
    oshelf['rescontacts'] = FMT['RESCONTACTS'] % (OutPrefix, TempSet)
    oshelf['contactmap_%s' % CompType] = (NativeContactMap, ClustContactMap)
    oshelf.close()

def FracNativeContacts():
    print 'Calculating fraction of native contacts'
    frachist = calc.GetFracNativeContacts()
    oshelf = shelve.open(OutShelf)
    oshelf['fracnativecontacts'] = frachist
    oshelf.close()

def ContactOrder():
    print 'Comparing native vs. Predicted  contact orders'
    nativeCO, clustCOhist = calc.CompareCO()
    oshelf = shelve.open(OutShelf)
    oshelf['co'] = (nativeCO, clustCOhist)
    oshelf.close()
    

#### MAIN ####
RMSD()
PhiPsiErr('traj')
Cluster()
ContactMap('traj')
ContactMap('topclust')
FracNativeContacts()
ContactOrder()

