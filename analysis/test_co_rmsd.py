#!/usr/bin/env python

import os, sys, numpy as np, shelve
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import protein

Dir = os.path.abspath(sys.argv[1])
Clrs = {'ff_ala_spc': 'cyan', 'ff_ala15': 'red', 'ff_leu15': 'green', 'ff_val15': 'black'}

PSET = ['1cb3', '1l2y', '2i9m', 'c_pep', 'ek_pep', '15beta', '1e0q', '1gb1', '1hrx', '1j4m']


def f(ax1, ax2, pset, fftype, plotnative = True):
    N = len(pset)
    L = np.zeros(N)
    native_co = np.zeros(N)
    traj_co = np.zeros(N)
    traj_rmsd = np.zeros(N)

    for i, p in enumerate(pset):
        nativepdb = os.path.expanduser('~/Go/native_struct/mapped/%s.pdb' % p)
        prot =  protein.ProteinClass(nativepdb)
        L[i] = len(prot.Seq)
        shelf = os.path.join(Dir, fftype, p, 'prot_%s.shelf' % p)
        d = shelve.open(shelf)
        rmsdhist = d['rmsd']
        avgrmsd = np.sum(rmsdhist[0] * rmsdhist[1]) / np.sum(rmsdhist[1])
        nativeco, cohist = d['co']
        avgco = np.sum(cohist[0] * cohist[1]) / np.sum(cohist[1])
        native_co[i] = nativeco
        traj_co[i] = avgco
        traj_rmsd[i] = avgrmsd
        d.close()

    ind = np.argsort(L)
    s = ', '.join([pset[i] for i in ind])
    print s
    if plotnative:
        ax1.plot(L[ind], native_co[ind], color = 'blue', marker = 'o', label = 'Native')
    ax1.plot(L[ind], traj_co[ind], color = Clrs[fftype], marker = 'o', label = fftype)
    ax2.plot(L[ind], traj_rmsd[ind], color = Clrs[fftype], marker = 'o', markersize = 8, label = fftype)
    ax1.legend()
    ax2.legend()

fig1 = plt.figure()
fig2 = plt.figure()
ax11 = fig1.add_subplot(1,2,1)
ax12 = fig2.add_subplot(1,2,1)
#f(ax11, ax12, PSET[0:5], 'ff_ala_spc', plotnative = True)
f(ax11, ax12, PSET[0:5], 'ff_ala15',  plotnative = True)
f(ax11, ax12, PSET[0:5], 'ff_leu15', plotnative = False)
f(ax11, ax12, PSET[0:5], 'ff_val15', plotnative = False)
ax11.set_title('helical peptides')
ax12.set_title('helical peptides')

ax21 = fig1.add_subplot(1,2,2)
ax22 = fig2.add_subplot(1,2,2)
#f(ax21, ax22, PSET[5:10], 'ff_ala_spc', plotnative = True)
f(ax21, ax22, PSET[5:10], 'ff_ala15', plotnative = True)
f(ax21, ax22, PSET[5:10], 'ff_leu15', plotnative = False)
f(ax21, ax22, PSET[5:10], 'ff_val15', plotnative = False)
ax21.set_title('hairpin peptides')
ax22.set_title('hairpin peptides')

for ax in [ax11, ax21]:
    ax.set_xlabel('sequence length')
    ax.set_ylabel('relative contact order')

for ax in [ax12, ax22]:
    ax.set_xlabel('sequence length')
    ax.set_ylabel('avg. RMSD ' + r'$(\AA)$')

fig1.tight_layout()
fig2.tight_layout()

fig1.savefig('avgco.png', bbox_inches = 'tight')
fig2.savefig('avgrmsd.png', bbox_inches = 'tight')
