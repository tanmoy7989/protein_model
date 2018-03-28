#!/usr/bin/env python

import os, numpy as np
import reasonable as cg

cfg = cg.config.Config()
cfg.BondedNCOSType = 1
cfg.NCOSType = 2
cfg.AtomS['GLY'] = cg.const.AtomS_GLY
cfg.SSRefType = 'number'
cfg.NativeType = 2
cfg.NonNativeType = -1

pdb = os.path.expanduser('~/protein_model/native_struct/mapped_pseudoGLY/1l2y.pdb')
p = cg.topo.ProteinNCOS(Pdb = pdb, cfg = cfg)

Sys = cg.topo.MakeSys(p)

ff = []
p_bb = cg.bb.P_Backbone(p = p, Sys = Sys)
ff.extend(p_bb.BB_0())

p_bbs = cg.bb_s.P_Backbone_Sidechain(p = p, Sys = Sys)
ff.extend(p_bbs.BB_S_Bonded_1())
ff.extend(p_bbs.BB_S_2())

ContactDict = cg.parsestruct.ParsePdb(p)
p_ss = cg.ss.P_Sidechain(p = p, Sys = Sys)
ff.extend(p_ss.Go_native_2(ContactDict))

Sys.ForceField.extend(ff)
cg.cgmodel.PrepSys(Sys)

Sys.Load()
