#!/usr/bin/env python

import os, numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cgprotein as lib

pdbname = '1l2y'
gly_dir = os.path.abspath('../native_struct/mapped_pseudoGLY')
nogly_dir = os.path.abspath('../native_struct/mapped')

p_gly = lib.ProteinNCOS(os.path.join(gly_dir, pdbname + '.pdb'), hasPseudoGLY = True)
phi_gly, psi_gly = p_gly.GetPhiPsi()

p_nogly = lib.ProteinNCOS(os.path.join(nogly_dir, pdbname + '.pdb'))
phi_nogly, psi_nogly = p_nogly.GetPhiPsi()

fig = plt.figure(figsize = (10, 5), facecolor = 'w', edgecolor = 'w')
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
ax1.scatter(phi_gly * (180./np.pi), psi_gly * (180./np.pi), c = 'red', s = 50)
ax2.scatter(phi_nogly * (180./np.pi), psi_nogly * (180./np.pi), c = 'red', s = 50)
ax1.set_title('peptide %s, pseudo gly' % pdbname)
ax2.set_title('peptide %s' % pdbname)
for ax in [ax1, ax2]:
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\psi$')
    ax.axhline(0.0, ls = '-', lw = 2, color = 'black')
    ax.axvline(0.0, ls = '-', lw = 2, color = 'black')
    ax.set_xlim([-180, 180])
    ax.set_ylim([-180, 180])

fig.tight_layout()
plt.show()
