#!/usr/bin/env python

import os, numpy as np, cPickle as pickle
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# retrieve data
ClustPrefix = os.path.abspath('./Analysis/leu15_bbs_topclust.294.84')
TrajPrefix = os.path.abspath('./Analysis/leu15_bbs.294.84')
with open(ClustPrefix + '.rescontacts.pickle', 'r') as of: ClustCM = pickle.load(of)
with open(ClustPrefix + '.rama.pickle', 'r') as of: PhiClust, PsiClust = pickle.load(of)
with open(TrajPrefix + '.rescontacts.pickle', 'r') as of: TrajCM, TrajDist = pickle.load(of)
with open(TrajPrefix + '.rama.pickle', 'r') as of: Phi, Psi, hist = pickle.load(of)

# find number of contacts over time
# and avg. TrajMap
NFrames = TrajCM.shape[0]
NContacts = np.zeros(NFrames)
AvgTrajCM = np.mean(TrajCM, axis = 0)
for i in range(NFrames):
    NContacts[i] = np.sum(TrajCM[i])

# make contact maps symmetric
AvgTrajCM += np.transpose(AvgTrajCM)
ClustCM += np.transpose(ClustCM)
NRes = ClustCM.shape[0]

# plot
fig = plt.figure(figsize = (30, 10), facecolor = 'w', edgecolor = 'w')

ax1 = fig.add_subplot(1,3,1)
im = ax1.imshow(AvgTrajCM, cmap = cm.Reds, origin = 'lower', aspect = 'auto', interpolation = 'gaussian', extent = [0, NRes-1, 0, NRes-1], alpha = 0.8)
inds = np.nonzero(ClustCM)
ax1.scatter(inds[0], inds[1], s = 100, c = 'blue', edgecolor = 'black', lw = 3)
ax1.plot(range(NRes+1), 'k--', lw = 2)
ax1.set_xticks(range(NRes))
ax1.set_yticks(range(NRes))
ax1.set_xlim([0, NRes-1])
ax1.set_ylim([0, NRes-1])
ax1.grid('on')

ax2 = fig.add_subplot(1,3,2)
ax2.plot(NContacts, lw = 0, marker = 'o', markersize = 10, color = 'black')
ax2.set_xlabel('iteration step')
ax2.set_ylabel('number of contacts')

ax3 = fig.add_subplot(1,3,3)
(PhiCenters, PsiCenters), h, err = hist
RamaPmf = - np.log(h)
# convert all angles to degrees
PhiClust *= (180. / np.pi)
PsiClust *= (180. / np.pi)
PhiCenters *= (180. / np.pi)
PsiCenters *= (180. / np.pi)
# trim pmf
#RamaPmf = RamaPmf.clip(min = -5, max = 5) # restrict within -5 kT and 5 kT
im = ax3.imshow(np.transpose(RamaPmf), origin = 'lower', aspect = 'auto', interpolation = 'gaussian', cmap = cm.Reds,
               extent = [PhiCenters.min(), PhiCenters.max(), PsiCenters.min(), PsiCenters.max()])
ax3.scatter(PhiClust, PsiClust, s = 100, c = 'blue', marker = 'o', edgecolor = 'k', lw = 4)
ax3.axhline(0., color = 'black', lw = 2)
ax3.axvline(0., color = 'black', lw = 2)
ax3.set_xlabel(r'$\phi$')
ax3.set_ylabel(r'$\psi$')
ax3.set_xlim([-180, 180])
ax3.set_ylim([-180, 180])

fig.tight_layout()
fig.savefig('contacts.png', bbox_inches = 'tight')
