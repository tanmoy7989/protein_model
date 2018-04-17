#!/usr/bin/env python

import os, sys, numpy as np
from reasonable import topo, parsestruct as ps, mapNCOS, config, const

PolyName = sys.argv[1]
TrajDir = os.path.abspath(sys.argv[2])

# settings
topo.Verbose = False
ps.Verbose = False
Temp = 294.84

# get the config object ready
cfg = config.Config()
cfg.AtomS['GLY'] = const.AtomS_GLY

# paths, etc
NativeDir_mapped = os.path.abspath('../native_struct/mapped_pseudoGLY')
NativeDir_unmapped = os.path.abspath('../native_struct/unmapped')
pset = ['1cb3', '1l2y', '2i9m', 'c_pep', 'ek_pep', '15beta', '1e0q', '1gb1', '1hrx', '1j4m']

def getData(structype):
    c = dict( (p, None) for p in pset)
    d_com = dict( (p, None) for p in pset)
    d_ss = dict( (p, None) for p in pset)
    for p in pset:
        if structype == 'native':
            # map 2 polymer
            AAPdb = os.path.join(NativeDir_unmapped, p + '.pdb')
            CGPdb = os.path.join(NativeDir_mapped, p + '.pdb')
            MappedPrefix = os.path.join(os.getcwd(), p + '_mapped2%s' % PolyName.lower())
            MappedPdb = mapNCOS.Map2Polymer(Pdb = CGPdb, AAPdb = AAPdb, PolyName = PolyName.upper(), hasPseudoGLY = True, MappedPrefix = MappedPrefix)
            # obtain the native contacts
            p0 = topo.ProteinNCOS(Pdb = CGPdb, cfg = cfg)
            c0 = ps.ParsePdb(p0)
            # obtain the distances between these native contacts from the mapped Pdb
            p1 = topo.ProteinNCOS(Pdb = MappedPdb, cfg = cfg)
            c1 = ps.ParsePdb(p0, ResContactList = c0['c_native'])
            # populate the dicts
            c[p] = c0['c_native']
            d_com[p] = c1['d_native']
            d_ss[p] = c1['d_ss_native']
            del p0, p1, c0, c1
            for f in [MappedPdb, 'ncos_nativecontact.txt', 'ncos_nonnativecontact.txt']:
                if os.path.isfile(f): os.remove(f)
        
        elif structype == 'pred':
            CGPdb = os.path.join(TrajDir, p, 'NativeAnalysis', 'prot_%s.%3.2f.clust.pdb' % (p, Temp))
            p0 = topo.ProteinNCOS(Pdb = CGPdb, cfg = cfg)
            c0 = ps.ParsePdb(p0)
            # populate the dicts
            c[p] = c0['c_native']
            d_com[p] = c0['d_native']
            d_ss[p] = c0['d_ss_native']
            del p0, c0
            for f in ['ncos_nativecontact.txt', 'ncos_nonnativecontact.txt']:
                if os.path.isfile(f): os.remove(f)

    return c, d_com, d_ss

print '\nGetting Native data'
c_native, d_com_native, d_ss_native = getData('native')
print '\n Getting predicted data...'
c_pred, d_com_pred, d_ss_pred = getData('pred')

os.system('clear')
print 'CONTACTS PREDICTED FOR %s' % TrajDir.split('/')[-1]
print '===============================================================================================\n'
print '%10s %20s %15s %15s %15s %15s %10s' % ('Peptide', 'contact', 'd_com_expt', 'd_com_pred', 'd_ss_expt', 'd_ss_pred', 'status')
for p in pset:
    NativePdb = os.path.join(NativeDir_mapped, p + '.pdb')
    pobj = topo.ProteinNCOS(Pdb = NativePdb, cfg = cfg)
    Seq = pobj.Seq
    for k, c in enumerate(c_native[p]):
        i, j = c
        status = 'retained' if c_pred[p].__contains__(c) else 'missed'
        item1 = p if k == 0 else ' '
        item2 = '(%d, %d), (%s, %s)' % (i, j, Seq[i], Seq[j])
        item3 = '%5.2f' % d_com_native[p][k]
        item5 = '%5.2f' % d_ss_native[p][k]
        item4 = ' ' 
        item6 = ' '
        if c_pred[p].__contains__(c):
            ind = c_pred[p].index(c)
            item4 = '%5.2f' % d_com_pred[p][ind]
            item6 = '%5.2f' % d_ss_pred[p][ind]
        item7 = status
        s = '%10s %20s %15s %15s %15s %15s %10s' % (item1, item2, item3, item4, item5, item6, item7)
        if status == 'missed':
            CRED = '\033[41m'
            CEND = '\033[0m'
            s = CRED + s + CEND
        print s
    raw_input()
    print '\n'










