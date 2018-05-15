#!/usr/bin/env python

''' main routine that creates different types of models
    '''

import numpy as np, copy, random
import protein

from const import *
import topo, bb, bb_s, ss, mapNCOS
import parsestruct as ps

Verbose = True

def Preamble(s = None):
    if s is None:
        s = '''
\t\t\t\t================================================================================================
\t\t\t\t(R)elative (E)ntropy (A)ssisted, (S)tructure (O)ptimized, (N)o (A)dded (B)ioinformatic (LE)xicon
\t\t\t\t================================================================================================
'''
    return s


def CheckSys(p, Sys, cfg):
    ''' checks the system to ensure filters have been applied correctly''' 
    #TODO:
    pass

def PrepSys(Sys, TempSet = RoomTemp):
    # set up histogram
    for P in Sys.ForceField:
        P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)
    # set integrator properties
    Int = Sys.Int
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    Int.Method.LangevinGamma = 0.01
    # set ambient conditions
    Sys.Temp = TempSet
    Sys.TempSet = TempSet
    Sys.BoxL = 0.0

def loadParam(Sys, FF_file):
    # loads only those interactions in the current system
    # if they are also present in the supplied ff file
    print '\nLoading parameters from %s' % FF_file.split('/')[-1]
    FileParamDict = {}
    lines = file(FF_file).readlines()
    inds = [lines.index(line) for line in lines if line.startswith('>>> POTENTIAL')]
    for ii, i in enumerate(inds):
        k = lines[i].split('>>> POTENTIAL ')[-1].strip()
        if ii < len(inds)-1:
            paramlines = lines[(i+1):inds[ii+1]]
        else:
            paramlines = lines[(i+1):]
        s = ''.join(x for x in paramlines)
        v = eval(s)
        FileParamDict[k] = v
    # see if Sys has these potentials
    hasPotentials = []
    for P in Sys.ForceField:
        if P.Name in FileParamDict.keys():
            P.SetParam(**(FileParamDict[P.Name]))
            hasPotentials.append(P.Name)
    s = ' Over-writing potentials %s' % (', '.join(hasPotentials))
    print s
    return hasPotentials

class UpdatePostLoad(object):
    ''' updates Sys object after loading params from forcefield files'''
    def __init__(self, Sys, cfg):
        self.Sys = Sys
        self.cfg = cfg
        self.UpdateNonNativeWCA()
        self.UpdateNonNativeWCA_MJ()
        # other update methods can be added as required
        return None

    def UpdateNonNativeWCA(self):
        print 'Updating Non-native WCA potentials for Go model'
        for P in self.Sys.ForceField:
            if P.Name == 'NonBondNonNative' and self.cfg.NonNativeType == 0:
                P.Cut = P.Sigma[0] * 2**(1/6.)
        self.Sys.ForceField.Update()
        return
   
    def UpdateNonNativeWCA_MJ(self):
        # check if all native sigmas are the same
        # that means a single sigma was used for all native contacts
        # set the non-native sigma to this value
        if not self.cfg.MJSigmas: return
        Test = self.cfg.MJSigmas.count(self.cfg.MJSigmas[0]) == len(self.cfg.MJSigmas)
        if Test:
            NativeSigma = self.cfg.MJSigmas[0]
            print 'Updating Non-Native WCA for MJ Go model to use the common Native sigma = %2.2f A' % NativeSigma
            for P in self.Sys.ForceField:
                if P.Name == 'NonBondNonNative' and self.cfg.NonNativeType == 0:
                    P.SetParam(Sigma = NativeSigma)
                    P.Cut = P.Sigma[0] * 2**(1/6.)
            self.Sys.ForceField.Update()
        return


def ErodeNativeContacts(ContactDict, DelFrac):
    ''' removes DelFrac fraction of native contacts
    from the given residue contact '''
    ResContactList = ContactDict['c_native']
    print 'ResContactList: ', ResContactList
    print 'Deleting %2.1f %% random native contacts: ' % (100. * DelFrac), 
    n = int(round(DelFrac * len(ResContactList)))
    # take the floor and not the nearest int, if it deletes all contacts when not requested
    if n == len(ResContactList):
        if not DelFrac == 1.0: n = int(DelFrac * len(ResContactList))
    NewResContactList = copy.copy(ResContactList)
    Deleted = []
    for i in range(n):
        ind = random.randrange(len(NewResContactList))
        Deleted.append(NewResContactList[ind])
        NewResContactList.pop(ind)
    print Deleted
    del ContactDict
    # recreate contactdict with these contacts
    print '\n'
    print 'Re-creating contact information with reduced contact list'
    return NewResContactList


def makePolymerSys(Seq, cfg, Prefix = None, TempSet = RoomTemp, NChains = 1):
    print Preamble()
    # create system topology
    print '\n'
    p = topo.ProteinNCOS(Seq = Seq, cfg = cfg, Prefix = Prefix)
    # ensure that siechains are referenced according to residue name
    cfg.SSRefType = 'name'
    print '\n'
    Sys = topo.MakeSys(p = p, cfg = cfg, NChains = NChains)
    ff = []
    # create backbone potentials
    BB = bb.P_Backbone(p, Sys, cfg)
    print '\n'
    ff.extend(BB.BB_0())
    # create backbone-sidechain potentials
    BB_S = bb_s.P_Backbone_Sidechain(p, Sys, cfg)
    # 1-alphabet bonded potentials
    print '\n'
    if cfg.Bonded_NCOSType == 0:
        print 'Error: Multi-alphabet models not implemented yet'
        exit()
    if cfg.Bonded_NCOSType == 1: ff.extend(BB_S.BB_S_Bonded_1())
    # 1-alphabet or constant repulsive nonbonded potentials
    print '\n'
    if cfg.NCOSType == 0:
        print 'Error: Multi-alphabet models not implemented yet'
        exit()
    if cfg.NCOSType == 1: ff.extend(BB_S.BB_S_1())
    if cfg.NCOSType == 2: ff.extend(BB_S.BB_S_2())
    # create sidechain-sidechain potentials (1-alphabet)
    SS = ss.P_Sidechain(p, Sys, cfg)
    print '\n'
    ff.extend(SS.SS_1())
    # populate forcefield
    Sys.ForceField.extend(ff)
    # set up other system properties
    PrepSys(Sys, TempSet = TempSet)
    # compile
    if Verbose: print '\nCompiling model...'
    Sys.Load()
    return p, Sys


def makeHarmonicGoSys(NativePdb, cfg, Prefix = None, TempSet = RoomTemp):
    print Preamble()
    # map NativePdb to polymer
    if cfg.Map2Polymer:
        print '\n'
        PdbName = NativePdb.split('/')[-1].split('.pdb')[0]
        AAPdb = os.path.join(NATIVEPATH['Unmapped'], PdbName + '.pdb')
        MappedPdb = mapNCOS.Map2Polymer(Pdb = NativePdb, PolyName = cfg.PolyName, AAPdb = AAPdb, hasPseudoGLY = cfg.hasPseudoGLY()) 
    # ensure that sidechains are referenced according to residue number
    cfg.SSRefType = 'number'
    # create system topology
    print '\n'
    p = topo.ProteinNCOS(Pdb = NativePdb, cfg = cfg, Prefix = Prefix)
    # parse native struct for native contacts in given pdb
    print '\n'
    ContactDict = ps.ParsePdb(p)
    # using these native contacts get distances between contact sidechains in the polymer-mapped Pdb
    if cfg.Map2Polymer:
        print '\n'
        print 'Repopulating native contact database with intra-sidechain distances from %s mapped sequence' % cfg.PolyName.upper()
        p_Mapped = topo.ProteinNCOS(Pdb = MappedPdb, cfg = cfg, Prefix = Prefix + '_map2%s' % cfg.PolyName.lower())
        MappedContactDict = ps.ParsePdb(p_Mapped, ResContactList = ContactDict['ResContactList'])
        ContactDict = MappedContactDict
    # bond sidechains of native contacts for harmonic restraints
    # must be done prior to creating the Sys object
    print '\n'
    p.BondNativeContacts(ContactDict)
    # create Sys object 
    print '\n'
    Sys = topo.MakeSys(p = p, cfg = cfg)
    ff = []
    # create backbone potentials
    BB = bb.P_Backbone(p, Sys, cfg = cfg)
    print '\n'
    ff.extend(BB.BB_0())
    # create backbone-sidechain potentials
    BB_S = bb_s.P_Backbone_Sidechain(p, Sys, cfg = cfg)
    # 1-alphabet bonded potentials
    print '\n'
    if cfg.Bonded_NCOSType == 0:
        print 'Error: Multi-alphabet models not implemented yet'
        exit()
    if cfg.Bonded_NCOSType == 1: ff.extend(BB_S.BB_S_Bonded_1())
    # constant repulsive nonbonded potentials
    print '\n'
    if cfg.NCOSType == 0:
        print 'Error: Multi-alphabet models not implemented yet'
        exit()
    if cfg.NCOSType == 1: ff.extend(BB_S.BB_S_1())
    if cfg.NCOSType == 2: ff.extend(BB_S.BB_S_2())
    # create sidechain-sidechain potentials
    print '\n'
    SS = ss.P_Sidechain(p, Sys, cfg = cfg, ContactDict = ContactDict)
    # harmonic restraints between native contacts
    print '\n'
    cfg.NativeType = 2
    ff.extend(SS.Go_native_2(FConst = cfg.NativeFConst, HarmonicFluct = cfg.NativeHarmonicFluct))
    # populate forcefield
    Sys.ForceField.extend(ff)
    # set up other system properties
    PrepSys(Sys, TempSet = TempSet)
    # compile
    if Verbose: print '\nCompiling model...'
    Sys.Load()
    return p, Sys
    
    
def makeSplineGoSys(NativePdb, cfg, Prefix = None, TempSet = RoomTemp, DelFrac = None):
    print Preamble()
    # ensure that sidechains are referenced according to residue number
    cfg.SSRefType = 'number'
    # create system topology
    p = topo.ProteinNCOS(Pdb = NativePdb, cfg = cfg, Prefix = Prefix)
    # parse native struct for native contacts in given pdb
    print '\n'
    if DelFrac: ps.Verbose = False
    ContactDict = ps.ParsePdb(p)
    if DelFrac:
        NewResContactList = ErodeNativeContacts(ContactDict, DelFrac)
        ps.Verbose = True
        # recreate the contact information
        del ContactDict
        ContactDict = ps.ParsePdb(p, ResContactList = NewResContactList)
    # create Sys object 
    print '\n'
    Sys = topo.MakeSys(p = p, cfg = cfg)
    ff = []
    # create backbone potentials
    BB = bb.P_Backbone(p, Sys, cfg = cfg)
    print '\n'
    ff.extend(BB.BB_0())
    # create backbone-sidechain potentials
    BB_S = bb_s.P_Backbone_Sidechain(p, Sys, cfg = cfg)
    # 1-alphabet bonded potentials
    print '\n'
    if cfg.Bonded_NCOSType == 0:
        print 'Multi-alphabet models not implemented yet'
        exit()
    if cfg.Bonded_NCOSType == 1:ff.extend(BB_S.BB_S_Bonded_1())
    # 1-alphabet or constant repulsive nonbonded potentials
    print '\n'
    if cfg.NCOSType == 0:
        print 'Multi-alphabet models not implemented yet'
        exit()
    if cfg.NCOSType == 1: ff.extend(BB_S.BB_S_1())
    if cfg.NCOSType == 2: ff.extend(BB_S.BB_S_2())
    # create sidechain-sidechain potentials
    print '\n'
    SS = ss.P_Sidechain(p, Sys, cfg = cfg, ContactDict = ContactDict)
    # native contacts
    print '\n'
    cfg.NativeType = 1
    ff.extend(SS.Go_native_1(Cut = cfg.NativeCut))
    # non-native contacts
    # Note: the non-native cutoff needs to be supplied carefully to be compatible
    # as a WCA with the supplied sigma
    print '\n'
    if not cfg.NonNativeType == -1:
        cfg.NonNativeType = 0
        ff.extend(SS.Go_nonnative_0())
    # populate forcefield
    Sys.ForceField.extend(ff)
    # set up other system properties
    PrepSys(Sys, TempSet = TempSet)
    # compile
    if Verbose: print '\nCompiling model...'
    Sys.Load()
    return p, Sys


def makeSplineGoMultiSys(NativePdb, cfg, NSys = None, MasterPrefix = None, TempSet = None):
    print Preamble()
    # parse all lists
    if NSys is None: NSys = len(NativePdb)
    else:
        if not isinstance(NativePdb, list):
            NativePdb = [NativePdb] * NSys 
    if not isinstance(cfg, list): cfg = [cfg] * NSys
    if not isinstance(TempSet, list):
        if TempSet is None: TempSet = [RoomTemp] * NSys
        else: TempSet = [TempSet] * NSys
    # ensure that sidechains are referenced according to residue number
    for i in range(NSys): cfg[i].SSRefType = 'number'
    print 'Creating %d-system Go model ensemble...' % NSys 
    # create the systems
    pList = []
    SysList = []
    for i in range(NSys):
        print '\nENSEMBLE SYSTEM: %d' % i
        print '==================='
        # set Prefix
        Prefix = MasterPrefix + '_%d' % i
        # create system topology
        p = topo.ProteinNCOS(Pdb = NativePdb[i], cfg = cfg[i], Prefix = Prefix)
        pList.append(p)
        # parse native struct for native contacts in given pdb
        ps.Verbose = False
        ContactDict = ps.ParsePdb(p)
        # create Sys object 
        Sys = topo.MakeSys(p = p, cfg = cfg[i])
        ff = []
        # create backbone potentials
        BB = bb.P_Backbone(p, Sys, cfg = cfg[i])
        ff.extend(BB.BB_0())
        # create backbone-sidechain potentials
        BB_S = bb_s.P_Backbone_Sidechain(p, Sys, cfg = cfg[i])
        # 1-alphabet bonded potentials
        if cfg[i].Bonded_NCOSType == 0:
            print 'ERROR: Multi-alphabet models not implemented yet'
            exit()
        if cfg[i].Bonded_NCOSType == 1: ff.extend(BB_S.BB_S_Bonded_1())
        # 1-alphabet or constant repulsive nonbonded potentials
        if cfg[i].NCOSType == 0:
            print 'ERROR: Multi-alphabet models not implemented yet'
            exit()
        if cfg[i].NCOSType == 1: ff.extend(BB_S.BB_S_1())
        if cfg[i].NCOSType == 2: ff.extend(BB_S.BB_S_2())
        # create sidechain-sidechain potentials
        SS = ss.P_Sidechain(p, Sys, cfg = cfg[i], ContactDict = ContactDict)
        # native contacts
        cfg[i].NativeType = 1
        ff.extend(SS.Go_native_1(Cut = cfg[i].NativeCut))
        # non-native contacts
        # Note: the non-native cutoff needs to be supplied carefully to be compatible
        # as a WCA with the supplied sigma
        if not cfg[i].NonNativeType == -1:
            cfg[i].NonNativeType = 0
            ff.extend(SS.Go_nonnative_0())
        # populate forcefield
        Sys.ForceField.extend(ff)
        # set up other system properties
        PrepSys(Sys, TempSet = TempSet[i])
        SysList.append(Sys)
    # compile
    if Verbose: print '\nCompiling extended ensemble model...'
    for i, Sys in enumerate(SysList):
        print '\nSystem: %d' % i
        Sys.Load()
    return pList, SysList


def makeMJGoSys(NativePdb, cfg, Prefix = None, TempSet = RoomTemp, Sigma = None):
    print Preamble()
    # ensure that sidechains are referenced according to residue number
    cfg.SSRefType = 'number'
    # create system topology
    p = topo.ProteinNCOS(Pdb = NativePdb, cfg = cfg, Prefix = Prefix)
    # parse native struct for native contacts in given pdb
    print '\n'
    ContactDict = ps.ParsePdb(p)
    # create Sys object 
    print '\n'
    Sys = topo.MakeSys(p = p, cfg = cfg)
    ff = []
    # create backbone potentials
    BB = bb.P_Backbone(p, Sys, cfg = cfg)
    print '\n'
    ff.extend(BB.BB_0())
    # create backbone-sidechain potentials
    BB_S = bb_s.P_Backbone_Sidechain(p, Sys, cfg = cfg)
    # 1-alphabet bonded potentials
    print '\n'
    if cfg.Bonded_NCOSType == 0:
        print 'Multi-alphabet models not implemented yet'
        exit()
    if cfg.Bonded_NCOSType == 1:ff.extend(BB_S.BB_S_Bonded_1())
    # 1-alphabet or constant repulsive nonbonded potentials
    print '\n'
    if cfg.NCOSType == 0:
        print 'Multi-alphabet models not implemented yet'
        exit()
    if cfg.NCOSType == 1: ff.extend(BB_S.BB_S_1())
    if cfg.NCOSType == 2: ff.extend(BB_S.BB_S_2())
    # create sidechain-sidechain potentials
    print '\n'
    SS = ss.P_Sidechain(p, Sys, cfg = cfg, ContactDict = ContactDict)
    # native contacts
    print '\n'
    ff.extend(SS.Go_native_MJ(Cut = cfg.NativeCut, Sigma = Sigma))
    # non-native contacts
    # Note: the non-native cutoff needs to be supplied carefully to be compatible
    # as a WCA with the supplied sigma
    print '\n'
    if not cfg.NonNativeType == -1:
        cfg.NonNativeType = 0
        ff.extend(SS.Go_nonnative_0())
    # populate forcefield
    Sys.ForceField.extend(ff)
    # set up other system properties
    PrepSys(Sys, TempSet = TempSet)
    # compile
    if Verbose: print '\nCompiling model...'
    Sys.Load()
    return p, Sys


def makeNonBondOnlySplineGoSys(NativePdb, cfg, Prefix = None, TempSet = RoomTemp):
    print Preamble()
    # ensure that sidechains are referenced according to residue number
    cfg.SSRefType = 'number'
    # create system topology
    p = topo.ProteinNCOS(Pdb = NativePdb, cfg = cfg, Prefix = Prefix)
    # parse native struct for native contacts in given pdb
    print '\n'
    ContactDict = ps.ParsePdb(p)
    print '\n'
    # create Sys object 
    print '\n'
    Sys = topo.MakeSys(p = p, cfg = cfg)
    ff = []
    # create only bonded backbone potentials
    BB = bb.P_Backbone(p, Sys, cfg = cfg)
    print '\n'
    ff.extend(BB.BB_BondOnly())
    # create backbone-sidechain potentials
    BB_S = bb_s.P_Backbone_Sidechain(p, Sys, cfg = cfg)
    # 1-alphabet bonded potentials
    print '\n'
    ff.extend(BB_S.BB_S_Bonded_BondOnly())
    # create sidechain-sidechain potentials
    print '\n'
    SS = ss.P_Sidechain(p, Sys, cfg = cfg, ContactDict = ContactDict)
    # native contacts
    print '\n'
    cfg.NativeType = 1
    ff.extend(SS.Go_native_1(Cut = cfg.NativeCut))
    # non-native contacts
    # Note: the non-native cutoff needs to be supplied carefully to be compatible
    # as a WCA with the supplied sigma
    print '\n'
    if not cfg.NonNativeType == -1:
        cfg.NonNativeType = 0
        ff.extend(SS.Go_nonnative_0())
    # populate forcefield
    Sys.ForceField.extend(ff)
    # set up other system properties
    PrepSys(Sys, TempSet = TempSet)
    # compile
    if Verbose: print '\nCompiling model...'
    Sys.Load()
    return p, Sys

