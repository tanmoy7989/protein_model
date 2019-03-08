#!/usr/bin/env python

import os, sys, numpy as np, subprocess, shelve, pickle
from scipy.misc import logsumexp
from tqdm import trange
import sim, protein, whamlib

NBins = 50
NBlocks = 4
UseFrames = 200


TrajDir = os.path.abspath(sys.argv[1])

DataDir = os.path.join(TrajDir, 'data')

TempFile = os.path.join(DataDir, 'temps.txt')
Temps = np.loadtxt(TempFile)

SamplePdb = os.path.join(DataDir, '0.current.pdb')
p = protein.ProteinClass(SamplePdb)
NRes = len(p) - 2 # all AA simulations were capped

WhamPickle = os.path.join(TrajDir, 'anal', 'wham.pickled.dat')
AllEneFn = os.path.join(TrajDir, 'anal', 'wham.eikn.npy')
u_kn = np.load(AllEneFn)[0]
NFrames = u_kn.shape[1]
FrameRange = range(0, NFrames, int(NFrames / UseFrames)) 

ScoreShelf = os.path.join(TrajDir, 'ss.shelf')
of = shelve.open(ScoreShelf)


get_wsum = lambda x : np.exp(logsumexp(x))

def get_ss(pdb):
    cmd = 'stride ' + pdb
    proc = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    out, err = proc.communicate() ; retcode = proc.returncode
    HelixScore = 0.
    BetaScore = 0.
    CoilScore = 0.
    if not retcode:
        ss = [i for i in out.split('\n') if 'ASG' in i]
        for i in ss:
            struct = i.split()[5]
            if struct in ['H', 'G', 'I']: HelixScore += 1.0
            elif struct in ['E', 'B', 'b']: BetaScore += 1.0
            else: CoilScore += 1.0
        HelixScore /= float(NRes)
        BetaScore /= float(NRes)
        CoilScore /= float(NRes)
    return HelixScore, BetaScore, CoilScore


def get_ss_traj(TempInd):
    print 'Calculating STRIDE scores for secondary structures at %3.2f K' % Temps[TempInd]
    keys = ['HelixScore_%d' % TempInd, 'BetaScore_%d' % TempInd, 'CoilScore_%d' % TempInd]
    if all([k in of for k in keys]): return
    HelixScore = []
    BetaScore = []
    CoilScore = []
    TrajFn = os.path.join(DataDir, '%d.mdtrj.crd.gz' % TempInd)
    PrmTopFn = os.path.join(DataDir, '%d.prmtop.parm7' % TempInd)
    Trj = sim.traj.Amber(TrjFile = TrajFn, PrmtopFile = PrmTopFn)
    Trj[0]
    Trj = Trj[-NFrames : ]
    for ii in trange(UseFrames):
        i = FrameRange[ii]
        p.Pos = Trj[i]
        p1 = p.Decap()
        tmpPdb = 'this.pdb'
        p1.WritePdb(tmpPdb)
        a,b,c = get_ss(tmpPdb)
        HelixScore.append(a)
        BetaScore.append(b)
        CoilScore.append(c)
        del p1
        os.remove(tmpPdb)
    of[keys[0]] = np.array(HelixScore)
    of[keys[1]] = np.array(BetaScore)
    of[keys[2]] = np.array(CoilScore)
    return


def get_logw():
    u_kn_use = np.zeros([len(Temps), UseFrames])
    for ii, i in enumerate(FrameRange):
        u_kn_use[:, ii] = u_kn[:, i]
    print 'Calculating log-weights at all temperatures'
    with open(WhamPickle, 'r') as op: whamdata = pickle.load(op)
    fk = whamdata['FK']
    beta_k = whamdata['BetaK']
    for i, T in enumerate(Temps):
        key = 'logw_%d' % i
        print '%3.2f' % T 
        logw = whamlib.log_weight(ekn = u_kn_use, betak = beta_k, targetbeta = beta_k[i], fk = fk)
        of[key] = logw
    return


def get_foldcurve(SSType):
    print 'Calculating folding curves for ', SSType
    
    # extract scores as a 2-D array
    x_kn = []
    for i in range(len(Temps)):
        x_kn.append(of['%sScore_%d' % (SSType, i)])
    x_kn = np.array(x_kn)

    score_block = np.zeros([NBlocks, len(Temps)])
    blocksize = int(UseFrames / NBlocks)
    # for each block
    for b in range(NBlocks):
        # for all desired temps
        for k in range(len(Temps)):
            # for all temps and frames
            this_x = x_kn[:, b*blocksize : (b+1)*blocksize].flatten()
            this_logw = of['logw_%d' % k][:, b*blocksize : (b+1)*blocksize].flatten()
            this_score = np.sum(this_x * np.exp(this_logw)) / get_wsum(this_logw)
            # add to block array
            score_block[b,k] = this_score
    
    # errors
    score = np.mean(score_block, axis = 0)
    err = np.std(score_block, axis = 0, ddof = 1)
    of['foldfrac_%s' % SSType] = (Temps, score, err)
    return



######## MAIN ########
print 'Using last %d frames' % NFrames
print 'Using %d effective frames' % UseFrames

# get log-weights
get_logw()

# get scores for each traj
for i in range(len(Temps)): get_ss_traj(i)

# get folding curves
get_foldcurve('Helix')
get_foldcurve('Beta')
get_foldcurve('Coil')

of.close()






