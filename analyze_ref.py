#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cgprotein as lib, utils, measure

# declare globalls
global c_AA, c_CG, r_AA, r_CG

# file formats
FMT = utils.FMT

# user input
polymer = sys.argv[1]
OutPrefix = os.path.abspath(sys.argv[2])

# extract AA and CG paths
print 'Setting ouput locations'
AAPrefix = OutPrefix + '_AA'
CGPrefix = OutPrefix + '_CG'
OutDir = os.path.dirname(OutPrefix)
if not os.path.isdir(OutDir): os.mkdir(OutDir)

# nativepdb w.r.t all calculations are done
NativePdb = os.path.expanduser('~/Go/%s_AA/topclust.pdb' % polymer)
if not os.path.isfile(NativePdb):
    raise IOError('AA top cluster struct missing.')

# global compute and replica objects for AA and CG
print 'Creating global Compute objects for AA and CG'
c_AA = lib.Compute(NativePdb = NativePdb, Prefix = AAPrefix)
c_CG = lib.Compute(NativePdb = NativePdb, Prefix = CGPrefix)

# global replica objects for AA
print 'Creating global Replica objects for AA'
fAA = utils.parseAA(polymer)
if not fAA['hasReplica']:
    print 'Error: Not enough data to create AA Replica object'
    r_AA = None
else:
    r_AA = lib.Replica(NativePdb = NativePdb, Prefix = AAPrefix,
                       TrajPrefix = os.path.join(fAA['MasterDir'], 'Lammps', polymer))

# global replica objects for CG
print 'Creating global Replica objects for CG'
fCG = utils.parseCG(polymer)
if not fCG['hasReplica']:
    print 'Error: Not enough data to create CG Replica object'
    r_CG = None
else:
    r_CG = lib.Replica(NativePdb = NativePdb, Prefix = CGPrefix,
                       TrajPrefix = os.path.join(fCG['MasterDir'], polymer))

# check if any replica (AA or CG) present
hasReplica = fAA['hasReplica'] or fCG['hasReplica']
hasAAReplica = fAA['hasReplica']
hasCGReplica = fCG['hasReplica']


# folding temp. estimation
def getFoldTemp(t, f):
    ind = np.argmin(abs(f-0.5))
    return t[ind]

def makeFoldCurve():
    print 'FOLD CURVES'
    print '-----------'
    # get folding fractions
    global r_AA, r_CG
    if hasAAReplica: r_AA.FoldCurve()
    if hasCGReplica: r_CG.FoldCurve()
    # get folding temps and save to file
    AApickle = os.path.join(AAPrefix, FMT['FOLDCURVE'] % (AAPrefix, 'RMSD'))
    CGpickle = os.path.join(CGPrefix, FMT['FOLDCURVE'] % (CGPrefix, 'RMSD'))
    s = ''
    if hasAAReplica:
        with open(AApickle, 'r') as of: t_AA, f_AA, err_AA = pickle.load(of)
        tfold_AA = getFoldTemp(t_AA, f_AA)
        s += 'tfold_AA: %3.2f\n' % tfold_AA
    if hasCGReplica:
        with open(CGpickle, 'r') as of: t_CG, f_CG, err_CG = pickle.load(of)
        tfold_CG = getFoldTemp(t_CG, f_CG)
        s += 'tfold_CG: %3.2f' % tfold_CG
    datafile = OutPrefix + '.tfold.dat'
    with open(datafile, 'w') as of: of.write(s)
    # plot fold curves
    fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
    ax = fig.add_subplot(1,1,1)
    if hasAAReplica:
        ax.errorbar(t_AA, f_AA, yerr = err_AA, color = 'blue', lw = 3, markersize = 6, label = 'AA')
        ax.axvline(tfold_AA, ls = '--', color = 'blue', lw = 1)
    if hasCGReplica:
        ax.errorbar(t_CG, f_CG, yerr = err_CG, color = 'red', lw = 3, markersize = 6, label = 'CG')
        ax.axvline(tfold_CG, ls = '--', color = 'red', lw = 1)
    ax.legend()
    ax.set_xlabel('temp (K)', fontsize = 15)
    ax.set_ylabel('folding fraction (f)', fontsize = 15)
    fig.tight_layout()
    figname = OutPrefix + '.foldcurve.png'
    plt.savefig(figname, bbox_inches = 'tight')


def makePMF(TempAA, TempCG):
    print 'Rg REE PMFs'
    print '-----------'
    # get pmfs
    global r_AA, r_CG
    if hasAAReplica:
        fAA = utils.parseAA(polymer, TempSet = TempAA)
        TempAA = fAA['Temp']
        r_AA.TempSet = TempAA
        r_AA.PMF2D('Rg', 'RMSD')
        AApickle = FMT['PMF2D'] % (AAPrefix, TempAA, 'Rg', 'RMSD')
        with open(AApickle, 'r') as of: (x1, y1), pmf1, err1 = pickle.load(of)
    if hasCGReplica:
        fCG = utils.parseCG(polymer, TempSet = TempCG)
        TempCG = fCG['Temp']
        r_CG.TempSet = TempCG
        r_CG.PMF2D('Rg', 'RMSD')
        CGpickle = FMT['PMF2D'] % (CGPrefix, TempCG, 'Rg', 'RMSD')
        with open(CGpickle, 'r') as of: (x2, y2), pmf2, err2 = pickle.load(of)
    # plot pmfs
    if hasAAReplica and hasCGReplica:
        fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (10,7), facecolor = 'w', edgecolor = 'w')
        ax1 = axes.flat[0]
        ax2 = axes.flat[1]
        im1 = ax1.imshow(np.transpose(pmf1), aspect = 'equal', interpolation = 'gaussian', origin = 'lower',
                         cmap = cm.nipy_spectral, extent = [x1.min(), x1.max(), y1.min(), y1.max()])
        im2 = ax2.imshow(np.transpose(pmf2), aspect = 'equal', interpolation = 'gaussian', origin = 'lower',
                     cmap = cm.nipy_spectral, extent = [x2.min(), x2.max(), y2.min(), y2.max()])
        
        fig.colorbar(im1, ax=axes.ravel().tolist(), label = 'pmf (kcal/mol)', orientation = 'vertical')
        ax1.set_title('AA: %3.2f K' % TempAA)
        ax2.set_title('CG: %3.2f K' % TempCG)
        for ax in [ax1, ax2]: ax.set_xlabel(r'$R_g (\AA)$', fontsize = 12)
        ax1.set_ylabel(r'RMSD (\AA)$', fontsize = 12)
    
    elif hasAAReplica:
        fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
        ax = fig.add_subplot(1,1,1)
        im = ax.imshow(np.transpose(pmf1), aspect = 'equal', origin = 'lower',
                       cmap = cm.jet, extent = [x1.min(), x1.max(), y1.min(), y1.max()])
        fig.colorbar(im, label = 'pmf (kcal/mol)', orientation = 'vertical')
        ax.set_title('AA: %3.2f K' % TempAA)
        ax.set_xlabel(r'$R_g (\AA)$', fontsize = 12)
        ax.set_ylabel(r'$RMSD (\AA)$', fontsize = 12)

    elif hasCGReplica:
        fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
        ax = fig.add_subplot(1,1,1)
        im = ax.imshow(np.transpose(pmf2), aspect = 'equal', origin = 'lower',
                      cmap = cm.jet, extent = [x2.min(), x2.max(), y2.min(), y2.max()])
        fig.colorbar(im, label = 'pmf (kcal/mol)', orientation = 'vertical')
        ax.set_title('CG: %3.2f K' % TempCG)
        ax.set_xlabel(r'$R_g (\AA)$', fontsize = 12)
        ax.set_ylabel(r'$RMSD (\AA)$', fontsize = 12)

    else: return
    figname = OutPrefix + '.%3.2f.Rg_RMSD.png' % TempAA
    plt.savefig(figname)#, bbox_inches = 'tight')


def ClusterAA(Temp):
    print 'Clustering AA traj at %3.2f K' % Temp
    print '------------------------------'
    global c_AA
    # parse AA filenames and update Compute object
    fAA = utils.parseAA(polymer, TempSet = Temp)
    c_AA.Update(TrajFn = fAA['Traj'], Temp = fAA['Temp'])
    # run clustering algorithm
    c_AA.Cluster()
    

def ClusterCG(Temp):
    print 'Clustering CG traj at %3.2f K' % Temp
    print '------------------------------'
    global c_CG
    # parse CG filenames and update Compute object
    fCG = utils.parseCG(polymer, TempSet = Temp)
    c_CG.Update(TrajFn = fCG['Traj'], Temp = fCG['Temp'])
    # run clustering algorithm
    c_CG.Cluster()
    

def rad2deg(x):
    x *= (180. / np.pi)
    return x


def makeRamaPlot(TempAA, TempCG):
    print 'RAMACHANDRAN PLOTS'
    print '------------------'
    global c_AA, c_CG
    # parse AA and CG filenames
    fAA = utils.parseAA(polymer, TempSet = TempAA)
    fCG = utils.parseCG(polymer, TempSet = TempCG)
    TempAA = fAA['Temp']
    TempCG = fCG['Temp']
    # get top cluster ramachandran plots
    AAClustPdb = FMT['CLUSTPDB'] % (AAPrefix, fAA['Temp'])
    if not os.path.isfile(AAClustPdb):
        print 'Error: AA top cluster struct missing. Cannot compute top cluster dihderals'
        Phi_AA = None
        Psi_AA = None
    else:
        pAA = lib.ProteinNCOS(AAClustPdb)
        Phi_AA, Psi_AA = pAA.GetPhiPsi()
        Phi_AA, Psi_AA = lib.TrimDihedrals(Phi_AA, Psi_AA, range(pAA.NRes))
    CGClustPdb = FMT['CLUSTPDB'] % (CGPrefix, fCG['Temp'])
    if not os.path.isfile(CGClustPdb):
        print 'Error: CG top cluster struct missing. Cannot compute top cluster dihderals'
        Phi_CG = None
        Psi_CG = None
    else:
        pCG = lib.ProteinNCOS(CGClustPdb)
        Phi_CG, Psi_CG = pCG.GetPhiPsi()
        Phi_CG, Psi_CG = lib.TrimDihedrals(Phi_CG, Psi_CG, range(pCG.NRes))
    # get whole traj based ramachandran plots
    c_AA.Update(TrajFn = fAA['Traj'], Temp = TempAA)
    c_CG.Update(TrajFn = fCG['Traj'], Temp = TempCG)
    c_AA.RamaChandran()
    c_CG.RamaChandran()
    # plot
    AApickle = FMT['RAMA'] % (AAPrefix, TempAA)
    CGpickle = FMT['RAMA'] % (CGPrefix, TempCG)
    with open(AApickle, 'r') as of: phiaa, psiaa, haa = pickle.load(of)['Generic']
    with open(CGpickle, 'r') as of: phicg, psicg, hcg = pickle.load(of)['Generic']
    pmf_AA = -np.log(haa[1])
    pmf_CG = -np.log(hcg[1])
    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (10,7), facecolor = 'w', edgecolor = 'w')
    ax1 = axes.flat[0]
    ax2 = axes.flat[1]
    X1 = rad2deg(haa[0][0]) ; Y1 = rad2deg(haa[0][1])
    X2 = rad2deg(hcg[0][0]) ; Y2 = rad2deg(hcg[0][1])
    im1 = ax1.imshow(np.transpose(pmf_AA), aspect = 'equal', interpolation = 'none', origin = 'lower',
                     cmap = cm.jet, extent = [X1.min(), X1.max(), Y1.min(), Y1.max()])
    im2 = ax2.imshow(np.transpose(pmf_CG), aspect = 'equal', interpolation = 'none', origin = 'lower',
                     cmap = cm.jet, extent = [X2.min(), X2.max(), Y2.min(), Y2.max()])
    #TODO: scatter the top cluster dihedrals on this plot
    fig.colorbar(im1, ax=axes.ravel().tolist(), label = 'pmf (kT)', orientation = 'vertical')
    ax1.set_title('AA : %3.2f K' % TempAA)
    ax2.set_title('CG : %3.2f K' % TempCG)
    for ax in [ax1, ax2]:
        ax.set_xlabel(r'$\phi$')
        ax.set_ylabel(r'$\psi$')
        ax.set_xlim([-190, 190])
        ax.set_ylim([-190, 190])
        ax.axvline(0, ls = '-', lw = 1, color = 'black')
        ax.axhline(0, ls = '-', lw = 1, color = 'black')
    figname = OutPrefix + '.%3.2f.ramaplot.png' % TempAA
    plt.savefig(figname, bbox_inches = 'tight')
    

######## MAIN ########

measure.NBins = 60
measure.NBlocks = 4

# Replica processing
if hasReplica:
    # generate folding curve
    makeFoldCurve()
    # extract AA and CG folding temps
    tfold = file(OutPrefix + '.tfold.dat').readlines()
    tfold_AA = float(tfold[0].split()[-1])
    tfold_CG = float(tfold[1].split()[-1])
    # make 2D pmfs at room and folding temperatures
    makePMF(300., 300.)
    makePMF(tfold_AA, tfold_AA)
    makePMF(tfold_CG, tfold_CG)


# cluster AA and CG trajectories at room and folding temp
ClusterAA(300.)
#ClusterAA(tfold_AA)
#ClusterAA(tfold_CG)
#ClusterCG(300.0)
#ClusterCG(tfold_AA)
#ClusterCG(tfold_CG)

# compute Ramachandran plots at room and folding temp
makeRamaPlot(300., 300.)
if os.path.isfile(OutPrefix + '.tfold.dat'):
    makeRamaPlot(tfold_AA, tfold_AA)
    #makeRamaPlot(tfold_CG, tfold_CG)




    
    
