#!/usr/bin/env python
import os, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import measure

sigma_lj = {'ala_spc': 3.7811, 'leu15': 3.7845}
epsilon_lj = {'ala_spc': 2.6623, 'leu15': 2.4763}

nativedistfile = os.path.abspath('../cgff/leu15_protg/nativelj/nativecontactdist.txt')
d = np.loadtxt(nativedistfile)
sigma0 = d.min() * (2**(-1/6.)) 

measure.NBins = 25
measure.NBlocks = 1
measure.Normalize = True
x, y, err = measure.makeHist(d)

plt.plot(x, y, 'b-', lw = 2)

l1 = 'ala_spc\n' + r'$\sigma = %2.2f \AA$' % sigma_lj['ala_spc'] + ' , ' + r'$\epsilon = %2.2f k_B T$' % (epsilon_lj['ala_spc'] / 0.6)
plt.axvline(sigma_lj['ala_spc'], ls = '--', lw = 2, color = 'red', label = l1)

l2 = 'leu15\n' + r'$\sigma = %2.2f \AA$' % sigma_lj['leu15'] + ' , ' + r'$\epsilon = %2.2f k_B T$' % (epsilon_lj['leu15'] / 0.6)
plt.axvline(sigma_lj['leu15'], ls = '--', lw = 2, color = 'black', label = l2)

l3 = '$\sigma_0 = %2.2f \AA$' % sigma0
plt.axvline(sigma0, ls = '--', lw = 2, color = 'blue', label = l3)

plt.legend()
plt.xlabel('res-res distance between native contacts ' + r'$ (\AA)$')
plt.ylabel('distribution')
plt.tight_layout()
plt.savefig('protg_sigma.png', bbox_inches='tight')
