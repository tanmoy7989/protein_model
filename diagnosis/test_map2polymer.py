#!/usr/bin/env python

import os, sys, numpy as np
import reasonable as cg

PdbName = sys.argv[1]

cfg = cg.config.Config()
cfg.BondedNCOSType = 1
cfg.NCOSType = 2
cfg.AtomS['GLY'] = cg.const.AtomS_GLY
cfg.SSRefType = 'number'

pdb = os.path.expanduser('~/protein_model/native_struct/mapped_pseudoGLY/%s.pdb' % PdbName)
p0 = cg.topo.ProteinNCOS(Pdb = pdb, cfg = cfg)
p1 = p0.Map2Polymer(PolyName = 'LEU', DelTempPdb = False)

