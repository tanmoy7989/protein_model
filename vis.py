#!/usr/bin/env python
import os, sys, numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.image as mpimg, matplotlib.pyplot as plt
import sim, protein, cgprotein

sys.path.append('/home/cask0/home/tsanyal/protein_model')
from reasonable import mapNCOS

PYMOLEXEC = os.environ['PYMOL_PATH']
IMAGEMAGICKEXEC = 'convert' # imagemagick tool
doRotate = True

def ExtractConect(NativePdb, Pdb):
    # clustered pdbs come with all cluster structs
    # followed by CONECT records at the end. This parses 
    # CONECT records from the native pdb and places it after 
    # the top cluster struct
    p = protein.ProteinClass(Pdb)
    pdbstr = p.GetPdb()
    s = file(NativePdb, 'r').readlines()
    start = [s.index(i) for i in s if i.startswith('CONECT')][0]
    stop = len(s)
    conectstr = ''.join(s[start:stop])
    return pdbstr + '\n' + conectstr
    
def RotateProteinClass(p):
    p = p.Copy()
    #first rotate so smallest eigenvector is pointing to screen
    Pos = p.ResPos(ResAtom = "C")
    evals, evecs = np.linalg.eig(np.cov(Pos, rowvar = False, ddof = 0))
    v = evecs[:, np.argmin(abs(evals))]
    Vec, Ang = sim.geom.GetVecMapping(np.array([0.,0.,1.]), v)
    p.Rotate(Vec, Ang)
    #next rotate so larget eigenvector is pointing horozontal
    Pos = p.ResPos(ResAtom = "C")
    evals, evecs = np.linalg.eig(np.cov(Pos, rowvar = False, ddof = 0))
    v = evecs[:, np.argmax(abs(evals))]
    Vec, Ang = sim.geom.GetVecMapping(np.array([1.,0.,0.]), v)
    p.Rotate(Vec, Ang)
    #determine if we need to flip
    Pos = p.ResPos(ResAtom = "C")
    avgz = np.mean(Pos[:,2])
    frac = (Pos[:,2] > avgz).astype(float).sum() / len(Pos)
    if frac > 0.5: p.Rotate(np.array([1.,0.,0.]), 180.)
    return p

def AlignProtein(p, pNative, BBInds):
    #align using the Kabsch algorithm
    p = p.Copy()
    pNative = pNative.Copy()
    Vec1, Vec2, RotMat, Residuals = sim.geom.AlignmentRMSD(Pos1 = pNative.Pos[BBInds], Pos2 = p.Pos[BBInds])
    newNativeBBPos = pNative.Pos[BBInds] + Vec1
    newBBPos = np.dot(p.Pos[BBInds] + Vec2, RotMat)
    pNative.Pos[BBInds] = newNativeBBPos
    p.Pos[BBInds] = newBBPos
    return p, pNative

def Overlay(NativePdb, Pdb, OutPrefix = None, Label = '', SinglePlot = False, hasPseudoGLY = False):
    global doRotate
    if OutPrefix is None: OutPrefix = 'go'
    # parse BBInds
    p_cg = cgprotein.ProteinNCOS(NativePdb, hasPseudoGLY = hasPseudoGLY)
    BBInds = p_cg.GetBBInds()
    # read pdbs
    pNative = protein.ProteinClass(NativePdb)
    p = protein.ProteinClass(Pdb, Model = 1)
    # rotate the native pdb
    if doRotate: pNative = RotateProteinClass(pNative)
    # align with rotated native struct
    p, pNative = AlignProtein(p, pNative, BBInds)
    # write to first set of tmp pdb files
    tmpNativePdb = os.path.join(os.getcwd(), '%s_tmpnative.pdb' % OutPrefix)
    tmpPdb = os.path.join(os.getcwd(), '%s_tmp.pdb' % OutPrefix)
    pNative.WritePdb(tmpNativePdb)
    p.WritePdb(tmpPdb)
    # now reverse map these pdbs to generate an approximate carbonyl group
    # so that STRIDE can assign secondary structures
    mapNCOS.ReverseMap(CGPdb = tmpNativePdb, Prefix = tmpNativePdb.split('.pdb')[0])
    mapNCOS.ReverseMap(CGPdb = tmpPdb, Prefix = tmpPdb.split('.pdb')[0])
    # prepare pymol script
    s = '''
# set background and display style
bg_color white
set antialias, 1
set orthoscopic, on
set depth_cue, 0
set ray_trace_mode, 1
ray renderer=0
# load pdbs
load %(TMPNATIVEPDB)s, native
load %(TMPPDB)s, predicted
# color pdbs
color blue, native
color red, predicted
# display only backbone
select bb, name N+CA+C+O
hide all
show cartoon
#cartoon tube
# save to file
zoom complete=1
png %(FILENAME)s, width = 1200, height = 1200, dpi = 300, ray = 1
'''
    d = {'TMPNATIVEPDB': tmpNativePdb, 'TMPPDB': tmpPdb}
    d['FILENAME'] = OutPrefix + '.png' if not SinglePlot else OutPrefix + '_tmp0.png'
    tmpPml = OutPrefix + '.pml'
    file(OutPrefix + '.pml', 'w').write(s % d)
    cmdstr1 = '%s -Qc %s' % (PYMOLEXEC, tmpPml)
    os.system(cmdstr1)
    for x in [tmpNativePdb, tmpPdb, tmpPml]: os.remove(x)
    # if single plot with supplied labels is requested (mostly when used from the command line)
    if SinglePlot:
        if not Label:
            rmsd = sim.geom.RMSD(pNative.Pos[BBInds], p.Pos[BBInds]) ; print rmsd
            Label = r'$RMSD = %2.2f \AA$' % rmsd
        pic0 = OutPrefix + '_tmp0.png'
        pic1 = OutPrefix + '_tmp1.png'
        cmdstr2 = '%s %s -trim -bordercolor white -background white -border 50x50 -quality 100 %s' % (IMAGEMAGICKEXEC, pic0, pic1)
        os.system(cmdstr2)
        pic = mpimg.imread(pic1)
        fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
        ax = fig.add_subplot(1,1,1)
        ax.imshow(pic, aspect = 'auto')
        ax.set_title(Label, fontsize = 8)
        ax.set_xticks([]) ; ax.set_yticks([])
        ax.set_xticklabels([]) ; ax.set_yticklabels([])
        figname = OutPrefix + '.png'
        plt.savefig(figname, bbox_inches = 'tight')
        for x in [pic0, pic1]: os.remove(x)
    return

def Panel(NativePdbs, Pdbs, NRows, NCols, Labels = [], OutDir = None, PanelPrefix = None, hasPseudoGLY = False, DelOverlayPng = True, OverlayPngDir = None):
    if OutDir is None: OutDir = os.path.getcwd()
    if not Labels:
        Labels = [x.split('/')[-1].split('.pdb')[0] for x in NativePdbs]
    fig = plt.figure(facecolor = 'w', edgecolor = 'w', figsize = (4*NCols, 5*NRows))
    N = len(NativePdbs)
    for i in range(N):
        ax = fig.add_subplot(NRows, NCols, i+1)
        ax.set_xticks([]) ; ax.set_yticks([])
        ax.set_xticklabels([]) ; ax.set_yticklabels([])
        NativePdb = NativePdbs[i]
        Pdb = Pdbs[i]
        # check if native struct present
        if NativePdb is None: NativePdb = '<unknown>.pdb'
        if not os.path.isfile(NativePdb):
            print 'Native structure for %s not found' % NativePdb.split('/')[-1].split('.pdb')[0]
        # check if predicted struct present
        if Pdb is None or (not Pdb is None and not os.path.isfile(Pdb)):
            print 'Predicted structure for %s not found' % NativePdb.split('/')[-1].split('.pdb')[0]
            ax.set_title(NativePdb.split('/')[-1].split('.pdb')[0])
            continue
        # align and overlay in pymol
        OutPrefix0 = os.path.join(OutDir, './this0')
        OutPrefix1 = os.path.join(OutDir, './this1')
        Overlay(NativePdb, Pdb, OutPrefix =  OutPrefix0, hasPseudoGLY = hasPseudoGLY)
        # crop pic and make it nice using imagemagick
        pic0 = OutPrefix0 + '.png'
        pic1 = OutPrefix1 + '.png'
        cmdstr = '%s %s -trim -bordercolor white -background white -border 50x50 -quality 100 %s' % (IMAGEMAGICKEXEC, pic0, pic1)
        os.system(cmdstr)
        # plot panel
        pic = mpimg.imread(pic1)
        ax.imshow(pic, aspect = 'auto')
        ax.set_title(Labels[i])
        # delete intermediate pngs
        os.remove(pic0)
        if not DelOverlayPng:
            if OverlayPngDir is None:
                OverlayPngDir = os.path.join(OutDir, 'overlay_png')
            if not os.path.isdir(OverlayPngDir): os.system('mkdir -p %s' % OverlayPngDir)
            PdbName = NativePdb.split('/')[-1].split('.pdb')[0]
            newpic1 = os.path.join(OverlayPngDir, PdbName + '_overlay.png')
            os.system('mv %s %s' % (pic1, newpic1))
        else:
            os.remove(pic1)
    # save panel
    plt.subplots_adjust(wspace = 0)
    if PanelPrefix is None: PanelPrefix = 'vispanel'
    figname = os.path.join(OutDir, PanelPrefix + '.png')
    plt.savefig(figname, bbox_inches = 'tight')


#### COMMAND LINE USAGE ####
if __name__ == '__main__':
    # 1) help
    helpstr = 'USAGE: python ~/protein_model/vis.py NativePdb Pdb OutPrefix [hasPseudoGLY] [DelFinalPng]'
    if len(sys.argv) < 3: print helpstr
    
    #2) for single protein overlay
    if len(sys.argv) >= 4:
        NativePdb = os.path.abspath(sys.argv[1])
        Pdb = os.path.abspath(sys.argv[2])
        OutPrefix = os.path.abspath(sys.argv[3])
        hasPseudoGLY = int(sys.argv[4]) if len(sys.argv) > 4 else 0
        Overlay(NativePdb, Pdb, OutPrefix = OutPrefix, hasPseudoGLY = hasPseudoGLY, SinglePlot = True)
