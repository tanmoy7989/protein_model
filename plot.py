#!/usr/bin/env python

import os, shutil, sys, numpy as np, cPickle as pickle, shelve
import utils, vis

#### PREDICTED STRUCTURES PANEL ####
def PlotPanel(DataDir, NativeDir = None, Prefix = None):
    print 'PLOTTING PANELS'
    print '---------------'
    import matplotlib ; matplotlib.use('Agg')
    import vis
    pset_type = sys.argv[1]
    OutDir = os.path.abspath(sys.argv[2])
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
    pset_type = sys.argv[1]
    OutDir = os.path.abspath(sys.argv[2])
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
    pset_type = sys.argv[1]
    OutDir = os.path.abspath(sys.argv[2])
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
    pset_type = sys.argv[1]
    OutDir = os.path.abspath(sys.argv[2])
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
    pset_type = sys.argv[1]
    OutDir = os.path.abspath(sys.argv[2])
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
    pset_type = sys.argv[1]
    OutDir = os.path.abspath(sys.argv[2])
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


######## MAIN ########
HELPSTR = '''USAGE: python ~/protein_model/plot.py pset_type OutDir hasPseudoGLY (0 or 1)'''
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print HELPSTR
        exit()
    
    hasPseudoGLY = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    if hasPseudoGLY:
        NativeDir = os.path.expanduser('~/protein_model/native_struct/mapped_pseudoGLY')
        AATopClustDir = os.path.expanduser('~/protein_model/native_struct/ff96_igb5_glghs_topclust_mapped_pseudoGLY')
    else:
        NativeDir = os.path.expanduser('~/protein_model/native_struct/mapped')
        AATopClustDir = os.path.expanduser('~/protein_model/native_struct/ff96_igb5_glghs_topclust_mapped')


    PlotPanel(NativeDir = NativeDir, DataDir = 'NativeAnalysis', Prefix = 'vispanel_native')
    PlotPanel(NativeDir = AATopClustDir, DataDir = 'AATopClustAnalysis', Prefix = 'vispanel_topclust')
    
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

