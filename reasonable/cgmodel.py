#!/usr/bin/env python

''' main routine that creates different types of models
    '''

import numpy as np
import protein

from const import *
import topo, bb, bb_s, ss
import parsestruct as ps

Verbose = True

def Preamble(s = None):
    if s is None:
        s = '''
================================================================================================
(R)elative (E)ntropy (A)ssisted, (S)tructure (O)ptimized, (N)o (A)dded (B)ioinformatic (LE)xicon
================================================================================================ \n
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
    return hasPotentials
 
def makePolymerSys(Seq, cfg, Prefix = None, TempSet = RoomTemp):
    print Preamble()
    # create system topology
    p = topo.ProteinNCOS(Seq = Seq, cfg = cfg, Prefix = Prefix)
    Sys = topo.MakeSys(p = p, cfg = cfg)
    ff = []
    # create backbone potentials
    BB = bb.P_Backbone(p, Sys, cfg)
    ff.extend(BB.BB_0())
    # create backbone-sidechain potentials
    BB_S = bb_s.P_Backbone_Sidechain(p, Sys, cfg)
    # bonded potentials
    if cfg.Bonded_NCOSType == 0: ff.extend(BB_S.BB_S_Bonded_0())
    if cfg.Bonded_NCOSType == 1: ff.extend(BB_S.BB_S_Bonded_1())
    # nonbonded potentials
    if cfg.NCOSType == 0: ff.extend(BB_S.BB_S_0())
    if cfg.NCOSType == 1: ff.extend(BB_S.BB_S_1())
    if cfg.NCOSType == 2: ff.extend(BB_S.BB_S_2())
    # create sidechain-sidechain potentials (1-alphabet)
    SS = ss.P_Sidechain(p, Sys, cfg)
    ff.extend(SS.SS_1())
    # populate forcefield
    Sys.ForceField.extend(ff)
    # set up other system properties
    PrepSys(Sys, TempSet = TempSet)
    # compile
    if Verbose: print 'Compiling model...'
    Sys.Load()
    return p, Sys

def makeGoSys(NativePdb, cfg, Prefix = None, TempSet = RoomTemp):
    print Preamble()
    # create system topology
    p = topo.ProteinNCOS(Pdb = NativePdb, cfg = cfg, Prefix = Prefix)
    Sys = topo.MakeSys(p = p, cfg = cfg)
    ff = []
    # create backbone potentials
    BB = bb.P_Backbone(p, Sys, cfg)
    ff.extend(BB.BB_0())
    # create 1-alphabet backbone-sidechain potentials
    BB_S = bb_s.P_Backbone_Sidechain(p, Sys, cfg)
    # bonded potentials
    cfg.Bonded_NCOSType = 1
    ff.extend(BB_S.BB_S_Bonded_1())
    # nonbonded potentials
    if cfg.NCOSType == 0: ff.extend(BB_S.BB_S_0())
    if cfg.NCOSType == 1: ff.extend(BB_S.BB_S_1())
    if cfg.NCOSType == 2: ff.extend(BB_S.BB_S_2())
    # create sidechain-sidechain potentials
    # parse native struct
    ContactDict = ps.ParsePdb(p)
    # map to polymer
    if cfg.Map2Polymer:
        MappedContactDict = ps.Map2Polymer(p = p, PolyName = cfg.PolyName, ContactDict = ContactDict)
        ContactDict = MappedContactDict
    SS = ss.P_Sidechain(p, Sys, cfg)
    # native contacts
    if cfg.NativeType == 0: ff.extend(SS.Go_native_0(ContactDict))
    if cfg.NativeType == 1: ff.extend(SS.Go_native_1(ContactDict))
    if cfg.NativeType == 2: ff.extend(SS.Go_native_2(ContactDict))
    # non-native contacts
    if cfg.NonNativeType == 0: ff.extend(SS.Go_nonnative_0(ContactDict))
    if cfg.NonNativeType == 1: ff.extend(SS.Go_nonnative_1(ContactDict))
    # populate forcefield
    Sys.ForceField.extend(ff)
    # set up other system properties
    PrepSys(Sys, TempSet = TempSet)
    # compile
    if Verbose: print 'Compiling model...'
    Sys.Load()
    return p, Sys


