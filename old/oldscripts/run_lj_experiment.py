#!/usr/bin/env python

import os, shutil, sys, numpy as np, cPickle as pickle, shelve
import utils, parse_potential as pp
import time

CURRDIR = os.getcwd()
RoomTemp = 300. 
AutoSigma = True # Turning this off makes the model recorded native sigma value from forcefield file (mostly keep this turned on)

#### WORKHORSE THAT RUNS REMD AND ANALYSIS ####
def runExperiment():
    import Go
    PdbName = sys.argv[2]
    FFType = sys.argv[3]
    OutDir = os.path.abspath(sys.argv[4])
    NReplica = int(sys.argv[5]) if len(sys.argv) > 5 else 8

    # parse paths etc
    OutDir = os.path.join(OutDir, PdbName)
    if not os.path.isdir(OutDir): os.system('mkdir -p %s' % OutDir)
    NativePdb = utils.parseNative(PdbName, MasterDir = os.path.expanduser('~/Go/native_struct/mapped'))
    try:
        AATopClustPdb = utils.parseNative(PdbName, MasterDir = os.path.expanduser('~/Go/native_struct/ff96_igb5_glghs_topclust_mapped'))
    except IOError:
        print 'Utils Error: Requested Top clust pdb does not exist'
        AATopClustPdb = None
    Prefix = 'prot_' + PdbName

    # parse forcefield parameters
    GoFF_File, FFMetadata = utils.parseGoFF(FFType)
    
    # backbone 
    BBType = FFMetadata['BBType']
    BBFF_File = utils.parseBBFF(BBType)[0]
    MinBondOrd = FFMetadata['MinBondOrd']
    NKnot = FFMetadata['NKnot']
    Cut = FFMetadata['Cut']
    
    # sidechain
    NCOSType = FFMetadata['NCOSType']
    GoType = FFMetadata['GoType']
    NonNativeType = FFMetadata['NonNativeType']
    
    # get LJ / spline params
    # Note: NonNative params are recorded from the given
    # forcefield file but not used currently, as these vary
    # depending on the sequence
    if not os.path.isfile(GoFF_File):
        print 'I/O Error: forcefield file does not exist'
        exit()
    # single LJ/Sigma for all native interactions
    if GoType == 1:
        NativeSigma = FFMetadata['NativeSigma']
        NativeEpsilon = FFMetadata['NativeEpsilon']
        NativeCut = FFMetadata['NativeCut']
    # LJ matrix for native interactions
    if GoType == 2:
        print 'Error: Not implemented yet'
        exit()
    # splined potential for all native interactions
    if GoType == 3:
        NativeKnots = FFMetadata['NativeKnots']
        NativeNKnot = len(NativeKnots)
        NativeCut = FFMetadata['NativeCut']
    if NonNativeType == 0:
        print 'Error: Will generate too much repulsion'
        exit()
    if NonNativeType == 1:
        NonNativeSigma = FFMetadata['NonNativeSigma']
        NonNativeEpsilon = FFMetadata['NonNativeEpsilon']
        NonNativeCut = FFMetadata['NonNativeCut']
    if NonNativeType == 2:
        NonNativeKnots = FFMetadata['NonNativeKnots']
        NonNativeNKnot = len(NonNativeKnots)
        NonNativeCut = FFMetadata['NonNativeCut']
    
    # temp schedule
    TLow = 270.
    THigh = 500.
    Temps = np.logspace(np.log10(TLow), np.log10(THigh), NReplica)
    TempInd = np.argmin(abs(Temps - RoomTemp))
    TempSet = Temps[TempInd]

    # time-step
    Delta_T = FFMetadata['TimeStep'] # femto-seconds

    # MD iterations
    NStepsMin = int(10000 / Delta_T)            # 10 ps
    NStepsEquil = int(50000000 / Delta_T)       # 50 ns
    NStepsProd  = int(20000000 / Delta_T)       # 20 ns
    NStepsSwap = int(1000 / Delta_T)            # 1 ps
    StepFreq = int(NStepsProd / 10000)          # need 10000 frames, 2 ps

    # REMD setup template
    mdstr_setup = '''
#!/usr/bin/env python
import os, sys, numpy as np, time
import Go, utils

# check if only analysis needs to be done
FullPrefix = os.path.join(os.getcwd(), "%(PREFIX)s")
isDone = utils.checkGoREMD(FullPrefix, [%(TEMPSET)3.2f])
if isDone: exit()

#### REMD ####
# basic settings
Go.Prefix = "%(PREFIX)s"
Go.InPdb ="%(NATIVEPDB)s"
Go.Temps = np.logspace(np.log10(%(TLOW)3.2f), np.log10(%(THIGH)3.2f), %(NREPLICA)d)
Go.TempSet = %(TEMPSET)3.2f
    
# backbone settings
Go.MinBondOrd = %(MINBONDORD)d
Go.NKnot = %(NKNOT)d
Go.SPCut = %(CUT)g
Go.BB_forcefield_file = "%(BBFF_FILE)s"
Go.NCOSType = %(NCOSTYPE)d

# md iterations
Go.NStepsMin = %(NSTEPSMIN)d
Go.NStepsEquil = %(NSTEPSEQUIL)d
Go.NStepsProd = %(NSTEPSPROD)d
Go.NStepsSwap = %(NSTEPSSWAP)d
Go.StepFreq = %(STEPFREQ)d
    '''
    
    # side-chain interactions template
    mdstr_sidechain = {
    'native_LJ': '''
# native contacts
Go.GoType = 1
Go.NativeSigma = None if %(AUTOSIGMA)d else %(NATIVESIGMA)g
Go.NativeEpsilon = %(NATIVEEPSILON)g
Go.NativeCut = %(NATIVECUT)g
    ''',
    
    'native_Spline': '''
# native contacts
Go.GoType = 3
Go.NativeNKnot = %(NATIVENKNOT)d
Go.NativeKnots = %(NATIVEKNOTS)s
Go.NativeCut = %(NATIVECUT)g
    ''',

    'nonnative_LJ': '''
# non-native contacts
Go.NonNativeType = 1
Go.NonNativeSigma = %(NONNATIVESIGMA)g
Go.NonNativeEpsilon = %(NONNATIVEEPSILON)g
Go.NonNativeCut = %(NONNATIVECUT)g
    ''',

    'nonnative_Spline': '''
# non-native contacts
Go.NonNativeType = 2
Go.NonNativeNKnot = %(NONNATIVENKNOT)d
Go.NonNativeKnots = %(NONNATIVEKNOTS)s
Go.NonNativeCut = %(NONNATIVECUT)g
    '''}


    # run template
    mdstr_run = '''
# compile Sys object
Sys = Go.makeGoSys()
for m in Sys.Int.Methods: m.TimeStep *= %(TIMESTEP)g
    
# run REMD
t1 = time.time()
TrajFile, LogFile = Go.runREMD(Sys)
t2 = time.time()
    
# reorder by temperature only at room temp
Go.reorderTraj(ReorderTemps = [%(TEMPSET)3.2f] )
t3 = time.time()
    
# print stats
print "REMD time: ", (t2-t1), " seconds"
print "Reordering time: ", (t3-t2), " seconds"
'''

    # job script template
    jobstr = '''
#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N %(JOBNAME)s

export PYTHONPATH=$PYTHONPATH:~/Go
date
python remd.py
mkdir -p ./NativeAnalysis
mkdir -p ./AATopClustAnalysis
python ~/Go/analyze_go.py %(NATIVEPDB)s %(PREFIX)s ./ ./NativeAnalysis
python ~/Go/analyze_go.py %(AATOPCLUSTPDB)s %(PREFIX)s ./ ./AATopClustAnalysis
date
'''
    # dict for filling md script template
    d_setup_run = {
        'PREFIX'        : Prefix,
        'NATIVEPDB'     : NativePdb,
 
        'TLOW'          : TLow,
        'THIGH'         : THigh,
        'NREPLICA'      : NReplica,
        'TEMPSET'       : TempSet,
        'TEMPIND'       : TempInd,
     
        'BBFF_FILE'     : BBFF_File,
        'MINBONDORD'    : MinBondOrd,
        'NKNOT'         : NKnot,
        'CUT'           : Cut,
        'NCOSTYPE'      : NCOSType,

         'TIMESTEP'     : Delta_T,
        'NSTEPSMIN'     : NStepsMin,
        'NSTEPSEQUIL'   : NStepsEquil,
        'NSTEPSPROD'    : NStepsProd,
        'NSTEPSSWAP'    : NStepsSwap,
        'STEPFREQ'      : StepFreq
        }
    
    if GoType == 1:
        d_native = {
        'AUTOSIGMA'     : int(AutoSigma),
        'NATIVESIGMA'   : NativeSigma,
        'NATIVEEPSILON' : NativeEpsilon,
        'NATIVECUT'     : NativeCut, 
       }
    
    if GoType == 2:
        print 'Error: Not yet implemented'
        exit()

    if GoType == 3:
        d_native = {
        'NATIVENKNOT'   : NativeNKnot,
        'NATIVEKNOTS'   : str(NativeKnots),
        'NATIVECUT'     : NativeCut
        }

    if NonNativeType == 0:
        # verified this with simulations of protein-G
        print 'Error: Will generate too much non-native repulsion'
        exit()

    if NonNativeType == 1:
        d_nonnative = {
        'NONNATIVESIGMA'    : NativeSigma,
        'NONNATIVEEPSILON'  : NativeEpsilon,
        'NONNATIVECUT'      : NonNativeCut,
        }

    if NonNativeType == 2:
        d_nonnative = {
        'NONNATIVENKNOT'    : NonNativeNKnot,
        'NONNATIVEKNOTS'    : NonNativeKnots,
        'NONNATIVECUT'      : NonNativeCut
        }

    # extract complete dic
    d1 = {}
    for x in [d_setup_run, d_native, d_nonnative]:
        for k, v in x.iteritems(): d1[k] = v

    # dict for filling job script template
    d2 = {'JOBNAME': Prefix, 'NATIVEPDB': NativePdb, 'AATOPCLUSTPDB': AATopClustPdb, 'PREFIX': Prefix}

    # extract complete template
    s = mdstr_setup
    if GoType == 1:
        s += mdstr_sidechain['native_LJ']
    if GoType == 2:
        print 'Error: Not implemented yet'
        exit()
    if GoType == 3:
        s += mdstr_sidechain['native_Spline']
    if NonNativeType == 0:
        print 'Error: Will generate too much non-native repulsion'
        exit()
    if NonNativeType == 1:
        s += mdstr_sidechain['nonnative_LJ']
    if NonNativeType == 2:
        s += mdstr_sidechain['nonnative_Spline']
    s += mdstr_run

    # fill template and submit job
    print 'Using backbone type:' , BBType
    print 'Using GoType: ', GoType
    print 'Using Non-Native Type: ', NonNativeType
    mdscript = os.path.join(OutDir, 'remd.py')
    jobscript = os.path.join(OutDir, 'remdjob.sh')
    if not os.path.isfile(mdscript): file(mdscript, 'w').write(s % d1)
    file(jobscript, 'w').write(jobstr % d2)
    os.chdir(OutDir)
    os.system('qsub remdjob.sh')
    os.chdir(CURRDIR)
    return



#### PREDICTED STRUCTURES PANEL ####
def PlotPanel(DataDir, NativeDir = None, Prefix = None):
    print 'PLOTTING PANELS'
    print '---------------'
    import matplotlib ; matplotlib.use('Agg')
    import vis
    pset_type = sys.argv[2]
    OutDir = os.path.abspath(sys.argv[3])
    layout = utils.getPanelLayout(pset_type)
    pset = layout['pset']
    NativePdbs = []
    ClustPdbs = []
    Labels = []
    for i, p in enumerate(pset):
        print 'Aligning  ', p
        # extract native pdb
        nativepdb = utils.parseNative(p, MasterDir = NativeDir)
        NativePdbs.append(nativepdb)
        # see if cluster pdb is present
        hasclustpdb, clustpdb = utils.hasCGData(Prefix = 'prot_'+p, OutDir = os.path.join(OutDir, p, DataDir), Key = 'clustpdb')
        ClustPdbs.append(clustpdb)
        # extract rmsds for labels
        if not hasclustpdb: Labels.append(None)
        else:
            hasrmsdhist, rmsdhist = utils.hasCGData(Prefix = 'prot_'+p, OutDir = os.path.join(OutDir, p, DataDir), Key = 'rmsd')
            if hasrmsdhist:
                rmsdavg = np.sum(rmsdhist[0] * rmsdhist[1]) / np.sum(rmsdhist[1])
                label = p + ' ' + r'$RMSD = %1.1f \AA$' % rmsdavg
                Labels.append(label)
            else: Labels.append(None)
    # plot panel
    print 'Rendering...'
    vis.Panel(NativePdbs, ClustPdbs, layout['NRows'], layout['NCols'], Labels, OutDir = OutDir, PanelPrefix = Prefix)
    return


#### RAMACHANDRAN PLOT PANEL ####
def PlotRamaProb(DataDir, Prefix = None):
    print 'PLOTTING RAMACHANDRAN MAPS'
    print '--------------------------'
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    pset_type = sys.argv[2]
    OutDir = os.path.abspath(sys.argv[3])
    layout = utils.getPanelLayout(pset_type)
    pset = layout['pset']
    fig = plt.figure(figsize = (5*layout['NCols'], 4*layout['NRows']), facecolor = 'w', edgecolor = 'w')
    for i, p in enumerate(pset):
        ax = fig.add_subplot(layout['NRows'], layout['NCols'], i+1)
        ax.set_xlim([-190, 190])
        ax.set_ylim([-190, 190])
        ax.set_title(p)
        # extract data
        print 'Extracting Ramachandran maps for ', p
        hasRamaProb, RamaPickle = utils.hasCGData(Prefix = 'prot_'+p, OutDir = os.path.join(OutDir, p, DataDir), Key = 'ramaprob')
        hasRamaErr, RamaErr = utils.hasCGData(Prefix = 'prot_'+p, OutDir = os.path.join(OutDir, p, DataDir), Key = 'ramaerr')
        if not (hasRamaProb and hasRamaErr): 
            print 'Rama plot not found for %s not found' % p
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            continue
        PhiNative, PsiNative, Err = RamaErr
        with open(RamaPickle, 'r') as of: data = pickle.load(of)
        PhiClust, PsiClust, RamaHist = data['Generic']
        (PhiCenters, PsiCenters), h, err = RamaHist
        RamaPmf = - np.log(h)
        # convert all angles to degrees
        PhiNative *= (180. / np.pi)
        PsiNative *= (180. / np.pi)
        PhiCenters *= (180. / np.pi)
        PsiCenters *= (180. / np.pi)
        # trim pmf
        RamaPmf = RamaPmf.clip(min = -5, max = 5) # restrict within -5 kT and 5 kT
        # plot
        im = ax.imshow(np.transpose(RamaPmf), origin = 'lower', aspect = 'auto', interpolation = 'gaussian', cmap = cm.Blues,
                       extent = [PhiCenters.min(), PhiCenters.max(), PsiCenters.min(), PsiCenters.max()])
        ax.scatter(PhiNative, PsiNative, s = 100, c = 'red', marker = 'o', edgecolor = 'k', lw = 4)
        ax.axhline(0., color = 'black', lw = 2)
        ax.axvline(0., color = 'black', lw = 2)
    if Prefix is None: Prefix = 'ramaprob'
    figname = os.path.join(OutDir, Prefix+'.png')
    fig.tight_layout()
    plt.savefig(figname, bbox_inches='tight')
    return


#### RAMA ERRORS PANEL ####
def PlotPhiPsiErr(DataDir, Prefix = None):
    print 'PLOTTING RAMACHANDRAN ERRORS'
    print '----------------------------'
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pset_type = sys.argv[2]
    OutDir = os.path.abspath(sys.argv[3])
    layout = utils.getPanelLayout(pset_type)
    pset = layout['pset']
    fig = plt.figure(figsize = (5*layout['NCols'], 4*layout['NRows']), facecolor = 'w', edgecolor = 'w')
    for i, p in enumerate(pset):
        ax = fig.add_subplot(layout['NRows'], layout['NCols'], i+1)
        ax.set_xlim([-190, 190])
        ax.set_ylim([-190, 190])
        ax.set_title(p, fontsize = 8)
        # extract data
        print 'Extracting Ramachandran errors for ', p
        hasRamaErr, RamaErr = utils.hasCGData(Prefix = 'prot_'+p, OutDir = os.path.join(OutDir, p, DataDir), Key = 'ramaerr')
        if not hasRamaErr:
            print 'RamaErrors for %s not found' % p
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            continue
        PhiNative, PsiNative, Err = RamaErr
        # scale the errors
        PhiNative *= (180. / np.pi)
        PsiNative *= (180. / np.pi)
        Err *= (180. / np.pi)
        Err = (Err - Err.min()) / (Err.max() - Err.min())
        # plot
        plt.scatter(PhiNative, PsiNative, s = 500*Err, c = 'red', alpha=0.4, edgecolors="grey", linewidth=2)
        ax.axhline(0., color = 'black', lw = 2)
        ax.axvline(0., color = 'black', lw = 2)
    if Prefix is None: Prefix = 'ramaerr'
    figname = os.path.join(OutDir, Prefix+'.png')
    fig.tight_layout()
    plt.savefig(figname, bbox_inches='tight')
    return


#### CONTACT MAPS  ####
def PlotContactMap(DataDir, Prefix = None):
    print 'PLOTTING CONTACT MAPS'
    print '---------------------'
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cmap
    pset_type = sys.argv[2]
    OutDir = os.path.abspath(sys.argv[3])
    layout = utils.getPanelLayout(pset_type)
    pset = layout['pset']
    fig = plt.figure(figsize = (5*layout['NCols'], 4*layout['NRows']), facecolor = 'w', edgecolor = 'w')
    for i, p in enumerate(pset):
        ax = fig.add_subplot(layout['NRows'], layout['NCols'], i+1)
        ax.set_title(p, fontsize = 8)
        # extract data
        print 'Extracting contact maps for ', p
        hasCM_Traj, d1 = utils.hasCGData(Prefix = 'prot_'+p, OutDir = os.path.join(OutDir, p, DataDir), Key = 'contactmap_traj')
        hasCM_topclust, d2 = utils.hasCGData(Prefix = 'prot_'+p, OutDir = os.path.join(OutDir, p, DataDir), Key = 'contactmap_topclust')
        if not (hasCM_Traj and hasCM_topclust):
            print 'Contact maps for %s not found' % p
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            continue
        NativeCM_, TrajCM = d1
        NativeCM, TopClustCM = d2
        NRes = NativeCM.shape[0]
        inds_Native = np.nonzero(NativeCM)
        inds_TopClust = np.nonzero(TopClustCM)
        # plot
        im = ax.imshow(TrajCM, cmap = cmap.Blues, origin = 'lower', aspect = 'auto', interpolation = 'gaussian', extent = [0, NRes-1, 0, NRes-1], alpha = 0.8)
        ax.scatter(inds_Native[0], inds_Native[1], s = 100, c = 'blue', edgecolor = 'black', lw = 3)
        ax.scatter(inds_TopClust[1], inds_TopClust[0], s = 100, c = 'red', edgecolor = 'black', lw = 3)
        ax.plot(range(NRes+1), 'k--', lw = 2)
        ax.set_xticks(range(NRes))
        ax.set_yticks(range(NRes))
        ax.set_xlim([0, NRes-1])
        ax.set_ylim([0, NRes-1])
        ax.grid('on')
    if Prefix is None: Prefix = 'contactmap'
    figname = os.path.join(OutDir, Prefix+'.png')
    fig.tight_layout()
    plt.savefig(figname, bbox_inches='tight')
    return


#### AA-CG CONTACT DISTANCE CORRELATIONS ####
def PlotContactCorr(DataDir, Prefix = None):
    pass


#### FRACTION OF NATIVE CONTACT DISTRIBUTIONS ####
def PlotFracNativeContact(DataDir, Prefix = None):
    print 'PLOTTING FRACTION OF NATIVE CONTACTS'
    print '---------------------'
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cmap
    pset_type = sys.argv[2]
    OutDir = os.path.abspath(sys.argv[3])
    layout = utils.getPanelLayout(pset_type)
    pset = layout['pset']
    fig = plt.figure(figsize = (5*layout['NCols'], 5*layout['NRows']), facecolor = 'w', edgecolor = 'w')
    for i, p in enumerate(pset):
        ax = fig.add_subplot(layout['NRows'], layout['NCols'], i+1)
        ax.set_title(p, fontsize = 8)
        # extract data
        print 'Extracting native contact fraction for ', p
        hasFracNativeContact, FracNativeContact = utils.hasCGData(Prefix = 'prot_'+p, OutDir = os.path.join(OutDir, p, DataDir), Key = 'fracnativecontacts')
        if not hasFracNativeContact:
            print 'Contact maps for %s not found' % p
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            continue
        frac, hist, err = FracNativeContact
        # plot
        ax.errorbar(frac, hist, yerr = err, lw = 2, markersize = 4)
        ax.set_xlim([0, 1])
    if Prefix is None: Prefix = 'fracnativecontact'
    figname = os.path.join(OutDir, Prefix+'.png')
    fig.tight_layout()
    plt.savefig(figname, bbox_inches='tight')
    return


#### CONTACT ORDER ####
def PlotContactOrder(DataDir, Prefix = None):
    print 'PLOTTING CONTACT ORDERS'
    print '---------------------'
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cmap
    pset_type = sys.argv[2]
    OutDir = os.path.abspath(sys.argv[3])
    layout = utils.getPanelLayout(pset_type)
    pset = layout['pset']
    fig = plt.figure(figsize = (5*layout['NCols'], 4*layout['NRows']), facecolor = 'w', edgecolor = 'w')
    for i, p in enumerate(pset):
        ax = fig.add_subplot(layout['NRows'], layout['NCols'], i+1)
        ax.set_title(p, fontsize = 8)
        # extract data
        print 'Extracting native contact fraction for ', p
        hasCO, CO = utils.hasCGData(Prefix = 'prot_'+p, OutDir = os.path.join(OutDir, p, DataDir), Key = 'co')
        if not hasCO:
            print 'Contact orders for %s not found' % p
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            continue
        nativeCO, TrajCOHist = CO
        TrajCO, hist, err = TrajCOHist
        AvgCO = np.sum(TrajCO * hist) / np.sum(hist)
        # plot
        ax.errorbar(TrajCO, hist, yerr = err, lw = 2, markersize = 4)
        ax.axvline(nativeCO, ls = '-', lw = 2, color = 'red')
        ax.axvline(AvgCO, ls = '-', lw = 2, color = 'blue')
    figname = os.path.join(OutDir, 'co.png')
    fig.tight_layout()
    plt.savefig(figname, bbox_inches='tight')
    return



#### COMMAND LINE RUNNING ####
HelpStr = '''
    USAGE: python ~/Go/run_lj_experiment.py run PdbName BBType OutDir [NReplica]
    USAGE: python ~/Go/run_lj_experiment.py plot pset_type OutDir
    '''
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print HelpStr
        exit()
    
    elif sys.argv[1] == 'plot':
        PlotPanel(NativeDir = None, DataDir = 'NativeAnalysis', Prefix = 'vispanel_native')
        PlotPanel(NativeDir = os.path.expanduser('~/Go/native_struct/ff96_igb5_glghs_topclust_mapped'), DataDir = 'AATopClustAnalysis', Prefix = 'vispanel_topclust')

        PlotRamaProb(DataDir = 'NativeAnalysis', Prefix = 'ramaprob_native')
        PlotRamaProb(DataDir = 'AATopClustAnalysis', Prefix = 'ramaprob_topclust')

        PlotPhiPsiErr(DataDir = 'NativeAnalysis', Prefix = 'ramaerr_native')
        PlotPhiPsiErr(DataDir = 'AATopClustAnalysis', Prefix = 'ramaerr_topclust')

        PlotContactMap(DataDir = 'NativeAnalysis', Prefix = 'contactmap_native')
        PlotContactMap(DataDir = 'AATopClustAnalysis', Prefix = 'contactmap_topclust')

        PlotContactCorr(DataDir = 'NativeAnalysis', Prefix = 'contactcorr_native')
        PlotContactCorr(DataDir = 'AATopClustAnalysis', Prefix = 'contactcorr_topclust')

        PlotFracNativeContact(DataDir = 'NativeAnalysis', Prefix = 'fracnativecontact_native')
        PlotFracNativeContact(DataDir = 'AATopClustAnalysis', Prefix = 'fracnativecontact_topclust')

        PlotContactOrder(DataDir = 'NativeAnalysis', Prefix = 'co_native')
        PlotContactOrder(DataDir = 'AATopClustAnalysis', Prefix = 'co_topclust')
    
    elif sys.argv[1] == 'run':
        runExperiment()
    
    else:
        print HelpStr
        exit()

