#!/usr/bin/env python

''' main routine that creates different types of models
    '''

import numpy as np
import protein

from const import *
import topo, bb, bb_s, ss, mapNCOS
import parsestruct as ps

Verbose = True

def Preamble(s = None):
    if s is None:
        s = '''
================================================================================================
(R)elative (E)ntropy (A)ssisted, (S)tructure (O)ptimized, (N)o (A)dded (B)ioinformatic (LE)xicon
================================================================================================
'''
    return s

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

class UpdatePostLoad(object):
    ''' updates Sys object after loading params from forcefield files'''
    def __init__(self, Sys, cfg):
        self.Sys = Sys
        self.cfg = cfg
        self.UpdateNonNativeWCA()
        # other update methods can be added as required
        return None

    def UpdateNonNativeWCA(self):
        print 'Updating Non-native WCA potentials for Go model'
        for P in self.Sys.ForceField:
            if P.Name == 'NonBondNonNative' and self.cfg.NonNativeType == 0:
                P.Cut = P.Sigma[0] * 2**(1/6.)
        self.Sys.ForceField.Update()
        return

def loadParam(Sys, FF_file):
    # loads only those interactions in the current system
    # if they are also present in the supplied ff file
    print 'Loading parameters from %s' % FF_file.split('/')[-1]
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
        if FileParamDict.keys().__contains__(P.Name):
            P.SetParam(**(FileParamDict[P.Name]))
            hasPotentials.append(P.Name)
    s = ' Over-writing potentials %s' % (', '.join(hasPotentials))
    print s
    return hasPotentials

def CheckSys(p, Sys, cfg):
    ''' checks the system to ensure filters have been applied correctly''' 
    #TODO:
    pass

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
    cfg.Bonded_NCOSType = 1
    print '\n'
    ff.extend(BB_S.BB_S_Bonded_1())
    # 1-alphabet or constant repulsive nonbonded potentials
    print '\n'
    cfg.NCOSType == 2
    ff.extend(BB_S.BB_S_2())
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
    cfg.Bonded_NCOSType = 1
    ff.extend(BB_S.BB_S_Bonded_1())
    # constant repulsive nonbonded potentials
    print '\n'
    cfg.NCOSType = 2
    ff.extend(BB_S.BB_S_2())
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
    
    
def makeSplineGoSys(NativePdb, cfg, Prefix = None, TempSet = RoomTemp):
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
    cfg.Bonded_NCOSType = 1
    print '\n'
    ff.extend(BB_S.BB_S_Bonded_1())
    # 1-alphabet or constant repulsive nonbonded potentials
    print '\n'
    cfg.NCOSType == 2
    ff.extend(BB_S.BB_S_2())
    # create sidechain-sidechain potentials
    print '\n'
    SS = ss.P_Sidechain(p, Sys, cfg = cfg, ContactDict = ContactDict)
    # native contacts
    print '\n'
    cfg.NativeType == 1
    ff.extend(SS.Go_native_1(Cut = cfg.NativeCut))
    # non-native contacts
    # Note: the non-native cutoff needs to be supplied carefully to be compatible
    # as a WCA with the supplied sigma
    print '\n'
    cfg.NonNativeType == 0
    ff.extend(SS.Go_nonnative_0())
    # populate forcefield
    Sys.ForceField.extend(ff)
    # set up other system properties
    PrepSys(Sys, TempSet = TempSet)
    # compile
    if Verbose: print '\nCompiling model...'
    Sys.Load()
    return p, Sys


