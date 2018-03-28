#!/usr/bin/env python
import os, sys, numpy as np, shelve
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import measure
measure.NBins = 25
measure.Normalize = True

pset = ['1cb3', '1l2y', '2i9m', 'c_pep', 'ek_pep', '15beta', '1e0q', '1gb1', '1hrx', '1j4m']
N = len(pset)
sigmas = np.zeros(N)
nativecos = np.zeros(N)
avgdist = np.zeros(N)
sigma0 = 3.78

for i, p in enumerate(pset):
    nativedistfile = os.path.abspath('../experiments/small_set_lj/ff_leu15_protg_lj_autosigma_nn1/%s/nativecontactdist.txt' % p)
    oshelf = os.path.abspath('../experiments/small_set_lj/ff_leu15_protg_lj_autosigma_nn1/%s/prot_%s.shelf' % (p, p))
    data = shelve.open(oshelf)
    nativeco, co = data['co']
    d = np.loadtxt(nativedistfile)
    sigmas[i] = d.min() * (2**(-1/6.))
    nativecos[i] = nativeco
    avgdist[i] = np.mean(d)

fig = plt.figure(facecolor = 'w', edgecolor = 'w')
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ind1 = np.argsort(avgdist[0:5])
ax1.plot(avgdist[0:5][ind1], sigmas[0:5][ind1], marker = 'o', markersize = 6, lw = 2)
ax1.axhline(sigma0, ls = '--', lw = 2, color = 'black')
ax1.set_xlabel('avg. native contact dist ' + r'$(\AA)$')
ax1.set_ylabel(r'$\sigma (\AA)$')
ax1.set_title('helical peptides')

ind2 = np.argsort(avgdist[5:10])
ax2.plot(avgdist[5:10][ind2], sigmas[5:10][ind2], marker = 'o', markersize = 6, lw = 2)
ax2.axhline(sigma0, ls = '--', lw = 2, color = 'black')
ax2.set_xlabel('avg. native contact dist ' + r'$(\AA)$')
ax2.set_ylabel(r'$\sigma (\AA)$')
ax2.set_title('hairpin peptides')

fig.tight_layout()
plt.savefig('sigmas.png', bbox_inches = 'tight')

