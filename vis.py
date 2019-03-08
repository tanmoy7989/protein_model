#!/usr/bin/env python
import os, sys, numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.image as mpimg, matplotlib.pyplot as plt
import sim, protein, cgprotein
from reasonable import mapNCOS

Renderer='pymol'

IMAGEMAGICKEXEC = 'convert' # imagemagick tool
doRotate = True

# Renderers
# Pymol
PYMOLEXEC = os.environ['PYMOL_PATH']
s_pymol = '''
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

# VMD
VMDEXEC = 'vmd'
TACHYONEXEC = 'tachyon_LINUXAMD64'
s_vmd = '''
# set background and display style
color Display Background white
display projection Orthographic
axes location Off
display nearclip set 0.01
# load pdbs
mol addrep 0
mol new {%(TMPNATIVEPDB)s} type {pdb}
mol addrep 1
mol new {%(TMPPDB)s} type {pdb}
# align pdbs
set ref_sel [atomselect 0 "not name S"]
set sel [atomselect 1 "not name S"]
set rotmat [measure fit $sel $ref_sel]
set move_sel0 [atomselect 0 "all"]
set move_sel1 [atomselect 1 "all"]
set com [measure center $move_sel0 weight mass]
$move_sel1 move $rotmat
$move_sel0 moveby [vecscale -1.0 $com]
$move_sel1 moveby [vecscale -1.0 $com]
# cartoon reps
mol modstyle 0 0 NewCartoon 0.3 20.0 4.1 0
mol modcolor 0 0 ColorID 0
mol modmaterial 0 0 Diffuse
mol modstyle 0 1 NewCartoon 0.3 20.0 4.1 0
mol modcolor 0 1 ColorID 1
mol modmaterial 0 1 Diffuse
#render
render Tachyon %(FILEPREFIX)s "%(TACHYONEXEC)s" -aasamples 12 %%s -format TARGA -o %%s.tga
quit
'''

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
    # align with rotated native struct (produces weird results with vmd, so let vmd do its own )
    if not Renderer == 'vmd':
        p, pNative = AlignProtein(p, pNative, BBInds)
    # write to first set of tmp pdb files
    tmpNativePdb = os.path.join(os.getcwd(), '%s_tmpnative.pdb' % OutPrefix)
    tmpPdb = os.path.join(os.getcwd(), '%s_tmp.pdb' % OutPrefix)
    pNative.WritePdb(tmpNativePdb)
    p.WritePdb(tmpPdb)
    # now reverse map these pdbs to generate an approximate carbonyl group
    # so that STRIDE can assign secondary structures
    mapNCOS.ReverseMap(CGPdb = tmpNativePdb, Prefix = tmpNativePdb.split('.pdb')[0], hasPseudoGLY = hasPseudoGLY)
    mapNCOS.ReverseMap(CGPdb = tmpPdb, Prefix = tmpPdb.split('.pdb')[0], hasPseudoGLY = hasPseudoGLY)
    # fill dictionary for renderer script 
    d = {'TMPNATIVEPDB': tmpNativePdb, 'TMPPDB': tmpPdb, }
    # pymol
    if Renderer == 'pymol':    
        d['FILENAME'] = OutPrefix + '.png' if not SinglePlot else OutPrefix + '_tmp0.png'
        tmpPml = OutPrefix + '.pml'
        file(tmpPml, 'w').write(s_pymol % d)
        cmdstr1 = '%s -Qc %s' % (PYMOLEXEC, tmpPml)
        os.system(cmdstr1)
        for x in [tmpNativePdb, tmpPdb, tmpPml]: os.remove(x)
    # vmd
    elif Renderer == 'vmd':
        d['FILEPREFIX'] = OutPrefix if not SinglePlot else OutPrefix + '_tmp0'
        d['TACHYONEXEC'] = TACHYONEXEC     
        tmpTcl = OutPrefix + '.tcl'    
        file(tmpTcl, 'w').write(s_vmd % d)   
        cmdstr1 = '%s -dispdev text -eofexit -e %s > /dev/null 2>&1' % (VMDEXEC, tmpTcl)
        os.system(cmdstr1)
        for x in [tmpNativePdb, tmpPdb, tmpTcl, d['FILEPREFIX']]: os.remove(x)    
    else:
        print 'ERROR: Renderer not found'
        exit()
    # if single plot with supplied labels is requested (mostly when used from the command line)
    if SinglePlot:
        if not Label:
            rmsd = sim.geom.RMSD(pNative.Pos[BBInds], p.Pos[BBInds]) ; print rmsd
            Label = r'$RMSD = %2.2f \AA$' % rmsd
        if Renderer == 'pymol': pic0 = OutPrefix + '_tmp0.png'
        else: pic0 = OutPrefix + '_tmp0.tga'
        pic1 = OutPrefix + '_tmp1.png'
        cmdstr2 = '%s %s -trim -bordercolor white -background white -border 50x50 -quality 100 %s' % (IMAGEMAGICKEXEC, pic0, pic1)
        os.system(cmdstr2)
        pic = mpimg.imread(pic1)
        fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
        ax = fig.add_subplot(1,1,1)
        ax.imshow(pic, aspect = 'auto')
        #ax.set_title(Label, fontsize = 8)
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
        OutPrefix0 = os.path.join(OutDir, './', 'this0')
        OutPrefix1 = os.path.join(OutDir, './', 'this1')
        Overlay(NativePdb, Pdb, OutPrefix =  OutPrefix0, hasPseudoGLY = hasPseudoGLY)
        # crop pic and make it nice using imagemagick
        if Renderer == 'pymol': pic0 = OutPrefix0 + '.png'
        else: pic0 = OutPrefix0 + '.tga'
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
