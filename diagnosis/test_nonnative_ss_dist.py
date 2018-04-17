#!/usr/bin/env python
import os, numpy as np, pickle
import matplotlib # ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sim
import reasonable as cg, pickleTraj, measure, parse_potential as pp

measure.NBins = 100
measure.Normalize = True

pdbfile = os.path.abspath('../native_struct/mapped_pseudoGLY/1pgb.pdb')
trajfile = os.path.abspath('../protg_AA/Lammps_pseudoGLY/protg.301.59.lammpstrj.gz')

cfg = cg.config.Config()
cfg.AtomS['GLY'] = cg.const.AtomS_GLY
p = cg.topo.ProteinNCOS(Pdb = pdbfile, cfg = cfg, Prefix = 'protg_test')
ret = cg.parsestruct.ParsePdb(p)

c_nonnative = ret['c_nonnative']
x, y = zip(*c_nonnative)
resnums = list(x) + list(y)
resnums = list(set(resnums))
SInds_nonnative = p.GetSInds(ResNums = resnums)

d_nonnative_pdb = ret['d_nonnative']
d_ss_nonnative_pdb = ret['d_ss_nonnative']

trj = pickleTraj(trajfile)
NFrames = len(trj)
if not os.path.isfile('protg_traj_nonnative_ss.txt'):
    d_ss_nonnative_traj = np.zeros( [NFrames, len(c_nonnative)] )
    pb = sim.utility.ProgressBar(Text = 'Calculating s-s distances for non-native contacts...', Steps = NFrames * len(c_nonnative))
    count = 0
    for n in range(NFrames):
        Pos = trj[n]
        for k, (i, j) in enumerate(c_nonnative):
            si = SInds_nonnative[i]
            sj = SInds_nonnative[j]
            d_si_sj = Pos[sj] - Pos[si]
            d = np.sqrt(np.sum(d_si_sj * d_si_sj))
            d_ss_nonnative_traj[n, k] = d
            pb.Update(count)
            count += 1
    np.savetxt('protg_traj_nonnative_ss.txt', d_ss_nonnative_traj)
            
d_ss_nonnative_traj = np.loadtxt('protg_traj_nonnative_ss.txt')
# histograming
h1 = measure.makeHist(d_ss_nonnative_pdb)
h2 = measure.makeHist(d_ss_nonnative_traj)
ret = d_ss_nonnative_traj

# plot
fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
ax = fig.add_subplot(1,1,1)
ax.plot(h1[0], h1[1], 'b-', lw = 2, label = 'Native Struct')
ax.plot(h2[0], h2[1], 'r-', lw = 2, label = 'AA Traj')
ax.legend()
ax.set_xlabel('s-s distance between nonnative contacts')
ax.set_ylabel('distribution')
fig.tight_layout()
fig.savefig('test_protg.png', bbox_inches = 'tight')
plt.show()

flist = ['protg_test_nativecontact.txt', 'protg_nonnativecontact.txt']
for x in flist:
    if os.path.isfile(x): os.remove(x)

    

