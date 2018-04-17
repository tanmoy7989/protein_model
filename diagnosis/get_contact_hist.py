#!/usr/bin/env python

import os, numpy as np
import matplotlib #; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import measure
from reasonable import topo, parsestruct as ps, config, const

Map2Polymer = True
PolyName = 'VAL'

cfg = config.Config()
cfg.AtomS['GLY'] = const.AtomS_GLY

measure.NBins = 25
measure.Normalize = True

def getData(pdbname, ax, of):
    print '\nPdb: %s\n-----------' % pdbname
    pdbfile = os.path.abspath('../native_struct/mapped_pseudoGLY/%s.pdb' % pdbname)
    p0 = topo.ProteinNCOS(Pdb = pdbfile, cfg = cfg, Prefix = pdbfile)
    p = p0.Map2Polymer(PolyName = PolyName) if Map2Polymer else p0
    cdict = ps.ParsePdb(p)
    d_native = cdict['d_native']
    d_ss = cdict['d_ss_native']
    h_native = measure.makeHist(d_native)
    h_ss = measure.makeHist(d_ss)
    of.write('%7s %14.3f %14.3f %14.3f %14.3f\n' % (pdbname, d_native.min(), d_native.max(), d_ss.min(), d_ss.max()))
    ax.plot(h_native[0], h_native[1], 'b-', lw = 2, label = 'COM-COM')
    ax.plot(h_ss[0], h_ss[1], 'g-', lw = 2, label = 'S-S')
    ax.legend()
    ax.set_title(pdbname)
    delflist = [pdbname+'_nativecontact.txt', pdbname+'_nonnativecontact.txt']
    for i in delflist:
        f = os.path.join(os.getcwd(), i)
        if os.path.isfile(f): os.remove(f)
    return
    
    
#### MAIN ####
native_mdatafile = os.path.abspath('../native_struct/native_metadata.txt')
mdata = eval(file(native_mdatafile).read())
pset = mdata['small_set'] + mdata['large_set']

outfile = os.path.abspath('./nativecontactdist')
if Map2Polymer: outfile += '_map2polymer_%s' % PolyName
outfile += '.dat'
of = open(outfile, 'w')
of.write('%7s % 14s %14s %14s %14s\n' % ('PdbName', 'min_com', 'max_com', 'min_ss', 'max_ss'))

fig = plt.figure(figsize = (30, 16), facecolor = 'w', edgecolor = 'w')
for i, pdbname in enumerate(pset):
    ax = fig.add_subplot(4,5,i+1)
    getData(pdbname, ax, of)
fig.tight_layout()

figname = 'nativecontact'
if Map2Polymer: figname += '_map2polymer_%s' % PolyName
figname += '.png'
fig.savefig(figname, bbox_inches = 'tight')
    






