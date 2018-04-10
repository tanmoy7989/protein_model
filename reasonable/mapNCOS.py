'''
/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <tanmoy.7989@gmail.com> wrote this file.  As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return. Tanmoy Sanyal
 * ----------------------------------------------------------------------------
 */
'''

#!/usr/bin/env python
import os, sys, copy, numpy as np
import protein, pdbtools
import topo

genBonds = 2 # 0 for no bonds, 1 for only BB bonds, 2 for all

NCaps = ['ACE', 'NH2']
CCaps = ['NME', 'NHE', 'NH2']
PDBFMT = "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f%6.2f"        # (atomind + 1, atomname, resname, reschainind, resnum+1, x, y, z)
BONDFMT = "%6s%5d%5d\n"                                                  # ('CONECT', a + r.StartAtom + 1, b + r.StartAtom + 1)


def findProchiralH_GLY(PosN, PosCA, PosC, PosH):
    ''' finds the prochiral H in a GLY, which can
    be used as a pseudo side chain'''
    PosH1 = PosH[0]
    PosH2 = PosH[1]
    # vectors in the N, CA, C plane using the Cahn-Ingold-Prelog rules
    vx = PosN - PosCA
    vy = PosC - PosCA
    # outward normal
    n = np.cross(vx, vy)
    # project CA-H vectors along the outward normal
    proj0 = np.dot(n, PosH[0] - PosCA)
    proj1 = np.dot(n, PosH[1] - PosCA)
    # return which hydrogen is prochiral-Si
    if proj0 > 0: return 0
    if proj1 > 0: return 1

def AddH_GLY(p):
    ''' hydrogenates the GLY residues iff they don't already have two hydrogens'''
    Filter = lambda ResName, AtomName: (ResName == 'GLY' and AtomName.__contains__('HA'))
    for i, r in enumerate(p.Seq):
        if not r == 'GLY': continue
        # does this GLY have hydrogens
        HInds = p.AtomInd(ResNum = i, UserFunc = Filter)
        if not list(HInds):
            p.TemplateAtoms(ResInd = [i], Elements = ['H'])
            p.TemplateBonds()
    return

def Map(InPdb, CGPrefix, AATraj = None, PrmTop = None, AmberEne = None, LastNFrames = 0, hasPseudoGLY = True):
    ''' Maps an all-atom pdb to a CG N-C-O-S version. If pseudo GLY side chains are requested,
    automatically hydrogenates GLY residues that don't have alpha hydrogens'''
    if hasPseudoGLY:
        print 'Using pseudo Glycines'
    # read in protein structure
    p = protein.ProteinClass(Pdb = InPdb)
    p.Update()
    # hydrogenate glycines that don't have hydrogens
    if hasPseudoGLY: AddH_GLY(p)
    # extract the all-atom pdb string
    PdbString = p.GetPdb()
    Pos = p.Pos
    # is this structure capped (# adapdted from /share/apps/pdbtools.py)
    NCap = [s for s in PdbString.split("\n") if s[0:4]=="ATOM" and s[17:20] in NCaps]
    CCap = [s for s in PdbString.split("\n") if s[0:4]=="ATOM" and s[17:20] in CCaps]
    hasNCap = len(NCap) > 0
    hasCCap = len(CCap) > 0
    # Perform some book-keeping on the sequence
    Seq = p.Seq
    if hasNCap: Seq = Seq[1:]
    if hasCCap: Seq = Seq[:-1]
    ResCount = dict( (x,Seq.count(x)) for x in set(Seq) )
    if 'GLY' in ResCount.keys() and not hasPseudoGLY:
        NCGAtoms = 3 * ResCount['GLY'] + 4 * ( sum(ResCount.values()) - ResCount['GLY'] )
    else:
        NCGAtoms = 4 * len(Seq)
    # Masks
    DecapFilter = lambda ResName, AtomName : not ( (NCaps + CCaps).__contains__(ResName) )
    SFilter_Other = lambda ResName, AtomName : not (AtomName == 'N' or AtomName == 'CA' or AtomName == 'C' or AtomName == 'O' or AtomName.__contains__('H'))
    SFilter_GLY = lambda ResName, AtomName: (ResName == 'GLY' and AtomName.__contains__('HA'))
    # Parse atom indices based on masks
    NInds = p.AtomInd(AtomName = 'N', UserFunc = DecapFilter)
    CAInds = p.AtomInd(AtomName = 'CA', UserFunc = DecapFilter)
    CInds = p.AtomInd(AtomName = 'C', UserFunc = DecapFilter)
    OInds = p.AtomInd(AtomName = 'O', UserFunc = DecapFilter)
    SInds = dict( (i, []) for i in range(len(Seq)) )
    for i, r in enumerate(Seq):
        startres = 1 if hasNCap else 0
        # sidechain hydrogens for GLY
        if r == 'GLY' and hasPseudoGLY:
            this_SInds = p.AtomInd(ResNum = startres + i, UserFunc = DecapFilter and SFilter_GLY)
            # find out which hydrogen is prochiral-Si
            ind_H = findProchiralH_GLY(PosN = Pos[NInds[i]], PosCA = Pos[CAInds[i]], PosC = Pos[CInds[i]], 
                                       PosH = ( Pos[this_SInds[0]], Pos[this_SInds[1]] ) )
            SInds[i] = [this_SInds[ind_H]]
        else:
            SInds[i] = p.AtomInd(ResNum = startres + i, UserFunc = DecapFilter and SFilter_Other)
    # Write Masks to a Map and final CG atoms to CG PDB
    s = ''
    s_bond = ''
    MapDict = dict( (i, []) for i in range(NCGAtoms) )
    n = 0
    for i, r in enumerate(Seq):
        # N lines
        if not i == 0:
            if Seq[i-1] == 'GLY' and not hasPseudoGLY: s_bond += BONDFMT % ('CONECT', n, n+1)
            else: s_bond += BONDFMT % ('CONECT', n-1, n+1)
        N_CGInd = n ; n+= 1
        MapDict[N_CGInd].append(NInds[i])
        s += PDBFMT % (N_CGInd+1, 'N  ', r, ' ', i+1, Pos[NInds[i], 0], Pos[NInds[i], 1], Pos[NInds[i], 2], 1.0, 0.0)
        s += "\n"
        # C lines
        s_bond += BONDFMT % ('CONECT', n, n+1)
        C_CGInd = n ; n+= 1
        MapDict[C_CGInd].append(CAInds[i])
        s += PDBFMT % (C_CGInd+1, 'C  ', r, ' ', i+1, Pos[CAInds[i], 0], Pos[CAInds[i] ,1], Pos[CAInds[i] ,2], 1.0, 0.0)
        s += "\n"
        # O lines
        s_bond += BONDFMT % ('CONECT', n, n+1)
        O_CGInd = n; n+=1
        MapDict[O_CGInd].extend([ CInds[i], OInds[i] ])
        COMPos = (Pos[CInds[i]] + Pos[OInds[i]]) / 2.
        s += PDBFMT % (O_CGInd+1, 'O  ', r, ' ', i+1, COMPos[0], COMPos[1], COMPos[2], 1.0, 0.0)
        s += "\n"
        # S lines
        if not r == 'GLY' or hasPseudoGLY:
            if genBonds == 2:
                s_bond += BONDFMT % ('CONECT', n-1, n+1)
            S_CGInd = n; n+= 1
            MapDict[S_CGInd].extend(SInds[i])
            COMPos = np.mean(Pos[SInds[i]], axis = 0)
            s += PDBFMT % (S_CGInd+1, 'S  ', r, ' ', i+1, COMPos[0], COMPos[1], COMPos[2], 1.0, 0.0)
            s += "\n"
    # Write CG PDB
    s += "TER\n"
    if genBonds: s += s_bond
    OutPdb = CGPrefix + '.pdb'
    with open(OutPdb, 'w') as of: of.write(s)

    # Map Traj
    if not AATraj is None:
        import sim
        simMap = sim.atommap.PosMap()
        for i in range(NCGAtoms):
            simMap += [sim.atommap.AtomMap(Atoms1 = MapDict[i], Atom2 = i)]
        AtomNames = []
        for r in Seq:
            if r == 'GLY' and not hasPseudoGLY: AtomNames.extend(['N', 'C', 'O'])
            else: AtomNames.extend(['N', 'C', 'O', 'S'])
        print 'Reading from AA Traj...'
        if PrmTop is None: Trj = sim.traj.lammps.Lammps(AATraj) # LammpsTraj
        else: Trj = sim.traj.amber.Amber(AATraj, PrmTop) # ZamTraj
        tmpinit = Trj[0]
        if Trj.FrameData.has_key('BoxL'): BoxL = Trj.FrameData['BoxL']
        else: BoxL = [0., 0., 0.]
        print 'Using Box: ', BoxL
        print 'Writing to CG Lammps Traj...'
        if LastNFrames: print 'Read %d frames, writing last %d frames' % (len(Trj), LastNFrames)
        # Note: the entire traj must be mapped to avoid file write errors
        CGTraj = CGPrefix + '.lammpstrj.gz'
        MappedTrj = sim.traj.mapped.Mapped(Trj, simMap, AtomNames = AtomNames, BoxL = BoxL)
        # now parse out only the necessary portion of the mapped traj
        if LastNFrames: MappedTrj = MappedTrj[-LastNFrames:]
        # now convert to Lammps
        sim.traj.base.Convert(MappedTrj, sim.traj.LammpsWrite, CGTraj, Verbose = True)
    
    # Write Ene File
    if not AmberEne is None:
        print 'Converting Amber Ene File...'
        CGEne = CGPrefix + '.ene.dat.gz'
        of = sim.traj.base.FileOpen(AmberEne, "rb")
        lines = of.readlines()
        start = 10 ; enefield_loc = (6,2)
        Ene = []
        for line in lines[start:]:
            l = line.split()
            if l[0] == 'L%d' % enefield_loc[0]:
                this_ene = float(l[enefield_loc[1]])
                Ene.append(this_ene)
        # parse necessary portion of Ene
        if LastNFrames: Ene = Ene[-LastNFrames:]
        np.savetxt(CGEne, Ene)
        of.close()
    return


def Map2Polymer(Pdb, PolyName, AAPdb, MappedPrefix = None, hasPseudoGLY = True, DelTmpPdb = True):
    ''' maps the given (CG) Pdb to a polymer of equivalent length 
    and returns a CG Pdb that is mapped to a poly-peptide (PolyName) of equivalent length
    Energy minimization not yet implemented'''
    # read in unmapped AA Pdb
    p_AA = protein.ProteinClass(AAPdb)
    p_AA = p_AA.Decap()
    # map this to a polymer of equivalent length
    print 'Mapping structure to a %s-%d sequence' % (PolyName, len(p_AA.Seq))
    PolySeq = [PolyName] * len(p_AA.Seq)
    p_AA = p_AA.MutateSeq(PolySeq)
    PolyAAPdb = 'polyAA.pdb'
    p_AA.WritePdb(PolyAAPdb)
    # coarse grain this pdb
    PolyCGPdb = 'polyCG.pdb'
    Map(InPdb = PolyAAPdb, CGPrefix = 'polyCG', hasPseudoGLY = hasPseudoGLY)
    # write coarse grained co-ordinates using given Pdb seq
    p_CG = protein.ProteinClass(Pdb)
    p_PolyCG = protein.ProteinClass(PolyCGPdb)
    p_CG.Pos = p_PolyCG.Pos
    if not MappedPrefix is None: MappedPdb = MappedPrefix + '.pdb'
    else:
        PdbName = Pdb.split('/')[-1].split('.pdb')[0]
        MappedPdb = os.path.join(os.getcwd(), PdbName + '_map2%s.pdb' % PolyName.lower())
    p_CG.WritePdb(MappedPdb)
    # copy over CONECT records if present
    s = ''
    with open(Pdb, 'r') as of: lines = of.readlines()
    try:
        start = [lines.index(line) for line in lines if line.startswith('CONECT')][0]
        stop = len(lines)
        s = ''.join(lines[start:stop])
    except ValueError:
        s = ''
    s0 = file(MappedPdb, 'r').read()
    if s:
        s0 += '\n'
        s0 += s
    file(MappedPdb, 'w').write(s0)
    # del temp files
    if DelTmpPdb:
        for i in [PolyAAPdb, PolyCGPdb]: os.remove(i)
    return MappedPdb


#### COMMAND LINE USAGE ####
if __name__ == '__main__':
    HelpStr = '''
python ~/protein_model/reasonable/mapNCOS.py map InPdb CGPrefix [hasPseudoGLY] [AATraj] [PrmTop] [AmberEne] [LastNFrames]
OR
python ~/protein_model/reasonable/mapNCOS.py map2poly InCGPdb Prefix AAPdb PolyName [hasPseudoGLY]
'''
    if len(sys.argv) == 1:
        print HelpStr
        exit()
    if sys.argv[1] == 'map':
        InPdb = os.path.abspath(sys.argv[2])
        CGPrefix = os.path.abspath(sys.argv[3])
        hasPseudoGLY = int(sys.argv[4]) if len(sys.argv) > 4 else 0
        AATraj = os.path.abspath(sys.argv[5]) if len(sys.argv) > 5 else None
        PrmTop = os.path.abspath(sys.argv[6]) if len(sys.argv) > 6 and sys.argv[6] else None
        AmberEne = os.path.abspath(sys.argv[7]) if len(sys.argv) > 7 and sys.argv[7] else None
        LastNFrames = int(sys.argv[8]) if len(sys.argv) > 8 else 0
        Map(InPdb = InPdb, CGPrefix = CGPrefix, AATraj = AATraj, PrmTop = PrmTop, AmberEne = AmberEne, LastNFrames = LastNFrames, hasPseudoGLY = hasPseudoGLY)

    if sys.argv[1] == 'map2poly':
        InCGPdb = os.path.abspath(sys.argv[2])
        MappedPrefix = os.path.abspath(sys.argv[3])
        AAPdb = os.path.abspath(sys.argv[4])
        PolyName = sys.argv[5]
        hasPseudoGLY = int(sys.argv[6]) if len(sys.argv) > 6 else 0
        Map2Polymer(Pdb = InCGPdb, PolyName = PolyName.upper(), AAPdb = AAPdb, MappedPrefix = MappedPrefix, hasPseudoGLY = hasPseudoGLY)

