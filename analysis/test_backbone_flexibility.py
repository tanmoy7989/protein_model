#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle
import matplotlib ; matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmap
import matplotlib.pyplot as plt
import protein, cgprotein as cg, utils

fftypes = ['ala15', 'leu15', 'val15']
NRes = 15

# create ideal helix pdb
p = protein.ProteinClass(['ALA']*15)
for i in range(15): p.RotateToPhiPsi(ResNum = i, Phi = -60, Psi = -45)
p.WritePdb('ideal_helix_unmapped.pdb')
cmd = 'python ~/Go/map.py ideal_helix_unmapped.pdb ./ideal_helix'
os.system(cmd)

# get ideal beta-hairpin (take the topclust from val15_AA simulations)
hairpin_pdb = os.path.abspath('../val15_AA/topclust.pdb')

# Trajectories
Traj = {'ala_spc': os.path.abspath('../ala_spc_AA/Lammps/ala_spc.300.00.lammpstrj.gz'),
        'ala15': os.path.abspath('../ala15_AA/Lammps/ala15.299.00.lammpstrj.gz'),
        'leu15': os.path.abspath('../leu15_AA/Lammps/leu15_wham.407.00.lammpstrj.gz'),
        'val15': os.path.abspath('../val15_AA/Lammps/val15_wham.367.00.lammpstrj.gz')
       }

Temp = {'ala_spc': 300.0, 'ala15': 299.00, 'leu15': 407.00, 'val15': 367.00}
NativePdb = {'ala_spc': 'ideal_helix.pdb', 'ala15': 'ideal_helix.pdb', 'leu15': 'ideal_helix.pdb', 'val15': hairpin_pdb}

fig = plt.figure(figsize = (22, 15), facecolor = 'w', edgecolor = 'w')
ind = 1
gs = gridspec.GridSpec(3,4)
for i, fftype in enumerate(fftypes):
    print '-------- %s --------' % fftype
    # get native Ramachandran angles
    pNative = cg.ProteinNCOS(NativePdb[fftype])
    PhiNative, PsiNative = pNative.GetPhiPsi(RamaType = 'Generic')
    PhiNative, PsiNative = cg.TrimDihedrals(PhiNative, PsiNative, ResNums = range(NRes))
    PhiNative *= (180 / np.pi)
    PsiNative *= (180 / np.pi)
    
    # create compute object
    calc = cg.Compute(NativePdb = NativePdb[fftype], TrajFn = Traj[fftype], Temp = Temp[fftype], Prefix = fftype)
    
    # get cluster ramachandran angle probabilities
    calc.RamaChandran()
    RamaPickle = utils.FMT['RAMA'] % (fftype, Temp[fftype])
    with open(RamaPickle, 'r') as of: data = pickle.load(of)['Generic']
    (x,y), ramahist, err = data[-1]
    x *= (180 / np.pi)
    y *= (180 / np.pi)
    pmf = -np.log(ramahist)
    pmf = pmf.clip(min = -5, max = 5)
    
    # get cluster contact maps
    NativeCM, TrajCM = calc.CompareContactMap()
    
    # get fraction of native contacts
    frac, frachist, err = calc.GetFracNativeContacts()
    
    # plot ramachandran plot
    ax1 = fig.add_subplot(gs[0, i])
    im = ax1.imshow(np.transpose(pmf), origin = 'lower', aspect = 'auto', interpolation = 'gaussian', cmap = cmap.Reds,
                    extent = [x.min(), x.max(), y.min(), y.max()])
    ax1.scatter(PhiNative, PsiNative, s = 100, c = 'blue', marker = 'o', edgecolor = 'k', lw = 4)
    ax1.axhline(0., color = 'black', lw = 2)
    ax1.axvline(0., color = 'black', lw = 2)
    ax1.set_xlim([-190, 190]) ; ax1.set_ylim([-190, 190])
    ax1.set_xlabel(r'$\phi$')
    ax1.set_ylabel(r'$\psi$')
    ax1.set_title(fftype)
    
    # plot contactmaps
    ax2 = fig.add_subplot(gs[1,i])
    inds = np.nonzero(NativeCM)
    ax2.scatter(inds[0], inds[1], s = 100, c = 'blue', edgecolor = 'black', lw = 3)
    im = ax2.imshow(TrajCM, cmap = cmap.Reds, origin = 'lower', aspect = 'auto', interpolation = 'gaussian', extent = [0, NRes-1, 0, NRes-1], alpha = 0.8)
    ax2.plot(range(NRes+1), 'k--', lw = 2)
    ax2.set_xticks(range(NRes)) ; ax2.set_yticks(range(NRes))
    ax2.set_xlim([0, NRes]) ; ax2.set_ylim([0, NRes])
    ax2.grid('on')
    ax2.set_xlabel('residue')
    ax2.set_ylabel('residue')

    # plot frac of native contacts
    ax3 = fig.add_subplot(gs[2,i])
    avgfrac = np.sum(frac * frachist) / np.sum(frachist)
    ax3.errorbar(frac, frachist, yerr = err, lw = 2, marker = 'o', markersize = 4, label = 'Avg frac = %1.2f' % avgfrac)
    ax3.set_xlim([0,1])
    ax3.set_xlabel('fraction of native contacts')
    ax3.set_ylabel('distribution')

fig.tight_layout()
figname = 'backbone_flex.png'
plt.savefig(figname, bbox_inches = 'tight')
for x in ['ideal_helix_unmapped.pdb', 'ideal_helix.pdb']: os.remove(x)
