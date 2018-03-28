#!/usr/bin/env python

import os, sys, numpy as np
from scipy.stats import linregress
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import protein
import reasonable as cg
import cgprotein as lib

PolyName = sys.argv[1]

def correlation(PdbName, ax1, ax2):
    Pdb = os.path.abspath('../native_struct/mapped/%s.pdb' % PdbName)
    AAPdb = os.path.abspath('../native_struct/unmapped/%s.pdb' % PdbName)
    # form protein struct by including glycines
    cfg = cg.config.Config()
    p = cg.topo.ProteinNCOS(cfg = cfg, Pdb = Pdb)
    # calculate contacts in original structure
    cdict0 = cg.parsestruct.ParsePdb(p)
    # remap to polymer
    cdict1 = cg.parsestruct.Map2Polymer(p = p, PolyName = PolyName, AAPdb = AAPdb, ContactDict = cdict0, EneMin = False, ReCalcContacts = False)
    #plot
    x_com = cdict0['d_native']
    y_com = cdict1['d_native']
    x_ss = cdict0['d_ss_native']
    y_ss = cdict1['d_ss_native']
    ax1.scatter(x_com, y_com, c = 'red', s = 50)
    ax1.plot(x_com, x_com, 'k-', lw = 2, Label = 'COM')
    ax1.legend(loc = 'best')
    
    ax2.scatter(x_ss, y_ss, c = 'red', s = 50)
    ax2.plot(x_ss, x_ss, 'k-', lw = 2, label = 'SS')
    ax2.legend(loc = 'best')
    
    for ax in [ax1, ax2]:
        ax.set_xlabel('actual pdb (cg mapped)')
        ax.set_ylabel('mapped to %s' % PolyName)
    return

def fracnative(PdbName, ax):
    Pdb = os.path.abspath('../native_struct/mapped/%s.pdb' % PdbName)
    AAPdb = os.path.abspath('../native_struct/unmapped/%s.pdb' % PdbName)
    # form protein struct by including glycines
    cfg = cg.config.Config()
    p = cg.topo.ProteinNCOS(cfg = cfg, Pdb = Pdb)
    # calculate contacts in original structure
    cdict0 = cg.parsestruct.ParsePdb(p)
    # remap to polymer but recalculate native contacts
    cdict1 = cg.parsestruct.Map2Polymer(p = p, PolyName = PolyName, AAPdb = AAPdb, ContactDict = cdict0, EneMin = False, ReCalcContacts = True)
    # contact map
    cmap0 = np.zeros([len(p.Seq), len(p.Seq)], float)
    cmap1 = np.zeros([len(p.Seq), len(p.Seq)], float)
    for k, (i,j) in enumerate(cdict0['c_native']): cmap0[i,j] = 1
    for k, (i,j) in enumerate(cdict1['c_native']): cmap1[i,j] = 1
    # get fraction of native contacts
    ind = (cmap0 == 1)
    frac = np.sum(cmap1[ind]) / np.sum(cmap0)
    hasMoreContacts = ( len(cdict1['c_native']) > len(cdict0['c_native']) )
    # plot
    ind0 = np.nonzero(cmap0)
    ind1 = np.nonzero(cmap1)
    ax.scatter(ind0[0], ind0[1], s = 100, c = 'blue', edgecolor = 'black', lw = 3)
    ax.scatter(ind1[1], ind1[0], s = 100, c = 'red', edgecolor = 'black', lw = 3)
    ax.plot(range(p.NRes+1), 'k--', lw = 2)
    ax.set_xticks(range(p.NRes))
    ax.set_yticks(range(p.NRes))
    ax.set_xlim([0, p.NRes-1])
    ax.set_ylim([0, p.NRes-1])
    ax.grid('on')
    return frac, hasMoreContacts
    
# create correlations for all sequences in all test sets
small_set = ['1cb3', '1l2y', '2i9m', 'c_pep', 'ek_pep', '15beta', '1e0q', '1gb1', '1hrx', '1j4m']
fig1 = plt.figure(figsize = (25, 12), facecolor = 'w', edgecolor = 'w')
fig2 = plt.figure(figsize = (25, 12), facecolor = 'w', edgecolor = 'w')
fig3 = plt.figure(figsize = (25, 12), facecolor = 'w', edgecolor = 'w')
fig4 = plt.figure(figsize = (10, 5), facecolor = 'w', edgecolor = 'w')
frac = np.zeros(len(small_set))
for i, pdbname in enumerate(small_set):
    print pdbname
    ax1 = fig1.add_subplot(2,5,i+1)
    ax2 = fig2.add_subplot(2,5,i+1)
    ax3 = fig3.add_subplot(2,5,i+1)
    for ax in [ax1, ax2, ax3]: ax.set_title(pdbname)
    correlation(pdbname, ax1, ax2)
    f, c = fracnative(pdbname, ax3)
    frac[i] = f

ax3 = fig4.add_subplot(1,2,1)
ax4 = fig4.add_subplot(1,2,2)
l = np.zeros(len(small_set))
for i, pdbname in enumerate(small_set):
    pdb = os.path.abspath('../native_struct/mapped/%s.pdb' % pdbname)
    p = lib.ProteinNCOS(pdb)
    l[i] = len(p.Seq)
ind0 = np.argsort(l[0:5])
ind1 = np.argsort(l[5:10])
print ind0, ind1
ax3.plot(l[0:5][ind0], frac[0:5][ind0], 'k-', lw = 3, marker = 'o', markersize = 6)
ax4.plot(l[5:10][ind1], frac[5:10][ind1], 'k-', lw = 3, marker = 'o', markersize = 6)
for ax in [ax3, ax4]:
    ax.set_xlabel('sequence length')
    ax.set_ylabel('frac. native contact. retained')
ax3.set_title('helical peptides')
ax4.set_title('hairpin peptides')

for fig in [fig1, fig2, fig3, fig4]: fig.tight_layout()
fig1.savefig('correlation_%s_COM.png' % PolyName, bbox_inches = 'tight')
fig2.savefig('correlation_%s_SS.png' % PolyName, bbox_inches = 'tight')
fig3.savefig('contactmap_%s.png' % PolyName, bbox_inches = 'tight')
fig4.savefig('fracnativecontact_%s.png' % PolyName, bbox_inches = 'tight')

    
