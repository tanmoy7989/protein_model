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
import os, sys, copy, numpy as np, string, argparse
import sim, protein

genBonds = 2 # 0 for no bonds, 1 for only BB bonds, 2 for all

NCaps = ['ACE']
CCaps = ['NME', 'NHE', 'NH2']
PDBFMT = "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f%6.2f"        # (atomind + 1, atomname, resname, reschainind, resnum+1, x, y, z)
BONDFMT = "%6s%5d%5d\n"                                                  # ('CONECT', a + r.StartAtom + 1, b + r.StartAtom + 1)

# Peptide Bond Geometry 
# source: http://www.cryst.bbk.ac.uk/PPS95/course/3_geometry/peptide2.html
BondLenCN = 1.32 # A
BondLenCAC = 1.53 # A
BondLenCO = 1.24 # A
AngleCACN = 114 # degrees


def Project(PosCA, PosCGO, PosNextN):
    '''reverse maps the carbonyl C and O from the CG O site'''
    # reference vector from CA
    Vec_CAN = PosNextN - PosCA
    VecLen_CAN = np.linalg.norm(Vec_CAN)
    # characterize the peptide plane with its unit normal
    Vecx = PosCA - PosCGO
    Vecy = PosNextN - PosCGO
    VecNormal = np.cross(Vecx, Vecy)
    UnitNormal = VecNormal / np.linalg.norm(VecNormal)
    # angle of in-plane rotation about reference vector from CA
    theta = np.arcsin( (BondLenCN/VecLen_CAN) * np.sin(AngleCACN * np.pi/180.) )
    # 3D rotation matrix
    RotMat3D = sim.geom.RotMat(Vec = UnitNormal, Ang = theta)
    # projection direction from reference vector
    NewVec = np.matmul(RotMat3D, Vec_CAN)
    NewUnitVec = NewVec / np.linalg.norm(NewVec)
    # project carbonyl C
    PosC = PosCA + BondLenCAC * NewUnitVec
    # project carbonyl O from COM condition
    PosO = 2*PosCGO - PosC
    return PosC, PosO

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
    Filter = lambda ResName, AtomName: (ResName == 'GLY' and 'HA' in AtomName)
    for i, r in enumerate(p.Seq):
        if not r == 'GLY': continue
        # does this GLY have hydrogens
        HInds = p.AtomInd(ResNum = i, UserFunc = Filter)
        if not list(HInds):
            p.TemplateAtoms(ResInd = [i], Elements = ['H'])
            p.TemplateBonds()
    return


def Map(InPdb, CGPrefix, Model = None, AATraj = None, PrmTop = None, AmberEne = None, LastNFrames = 0, hasPseudoGLY = True):
    ''' Maps an all-atom pdb to a CG N-C-O-S version. If pseudo GLY side chains are requested,
    automatically hydrogenates GLY residues that don't have alpha hydrogens'''
    if hasPseudoGLY:
        print 'Using pseudo Glycines'
    # read in protein structure
    p = protein.ProteinClass(Pdb = InPdb, Model = Model)
    p.Update()
    # hydrogenate glycines that don't have hydrogens
    if hasPseudoGLY: AddH_GLY(p)
    # extract the all-atom pdb string
    PdbString = p.GetPdb()
    Pos = p.Pos
    Seq = p.Seq
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
    DecapFilter = lambda ResName, AtomName : not ( ResName in (NCaps + CCaps) )    
    SFilter_Other = lambda ResName, AtomName : not (AtomName == 'N' or AtomName == 'CA' or AtomName == 'C' or AtomName == 'O' or 'H' in AtomName)
    SFilter_GLY = lambda ResName, AtomName: (ResName == 'GLY' and 'HA' in AtomName)
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
    CurrentChain = -1    
    s = ''
    s_bond = ''
    MapDict = dict( (i, []) for i in range(NCGAtoms) )
    n = 0
    for i, r in enumerate(Seq):
        # determine if a new chain starts here
        thisChain = p.ResChain(i)
        if not thisChain == CurrentChain:
            if not i == 0: s +=  "TER\n" # inter-chain TER records
        # N lines (do not bond when starting a new chain)
        if thisChain == CurrentChain:
            if Seq[i-1] == 'GLY' and not hasPseudoGLY: s_bond += BONDFMT % ('CONECT', n, n+1)
            else: s_bond += BONDFMT % ('CONECT', n-1, n+1)
        N_CGInd = n ; n+= 1
        MapDict[N_CGInd].append(NInds[i])
        s += PDBFMT % (N_CGInd+1, 'N  ', r, string.ascii_uppercase[thisChain], i+1, Pos[NInds[i], 0], Pos[NInds[i], 1], Pos[NInds[i], 2], 1.0, 0.0)
        s += "\n"
        # now update chains
        CurrentChain = thisChain
        # C lines
        s_bond += BONDFMT % ('CONECT', n, n+1)
        C_CGInd = n ; n+= 1
        MapDict[C_CGInd].append(CAInds[i])
        s += PDBFMT % (C_CGInd+1, 'C  ', r, string.ascii_uppercase[thisChain], i+1, Pos[CAInds[i], 0], Pos[CAInds[i] ,1], Pos[CAInds[i] ,2], 1.0, 0.0)
        s += "\n"
        # O lines
        s_bond += BONDFMT % ('CONECT', n, n+1)
        O_CGInd = n; n+=1
        MapDict[O_CGInd].extend([ CInds[i], OInds[i] ])
        COMPos = (Pos[CInds[i]] + Pos[OInds[i]]) / 2.
        s += PDBFMT % (O_CGInd+1, 'O  ', r, string.ascii_uppercase[thisChain], i+1, COMPos[0], COMPos[1], COMPos[2], 1.0, 0.0)
        s += "\n"
        # S lines
        if not r == 'GLY' or hasPseudoGLY:
            if genBonds == 2:
                s_bond += BONDFMT % ('CONECT', n-1, n+1)
            S_CGInd = n; n+= 1
            MapDict[S_CGInd].extend(SInds[i])
            COMPos = np.mean(Pos[SInds[i]], axis = 0)
            s += PDBFMT % (S_CGInd+1, 'S  ', r, string.ascii_uppercase[thisChain], i+1, COMPos[0], COMPos[1], COMPos[2], 1.0, 0.0)
            s += "\n"
    # Write CG PDB
    s += "TER\n" # last TER record
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


def Map2Polymer(Pdb, AAPdb, PolyName, Model = None, MappedPrefix = None, hasPseudoGLY = True, DelTmpPdb = True):
    ''' maps the given (CG) Pdb to a polymer of equivalent length 
    and returns a CG Pdb that is mapped to a poly-peptide (PolyName) of equivalent length
    Energy minimization not yet implemented'''
    # read in unmapped AA Pdb
    p_AA = protein.ProteinClass(AAPdb, Model = Model)
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


def ReverseMap(CGPdb, Prefix, Model = None, hasPseudoGLY = False):
    # parse coarse grained pdb
    p = protein.ProteinClass(CGPdb, Model = Model)
    Seq = p.Seq
    Pos = p.Pos
    PdbString = p.GetPdb()
    # parse backbone and sidechain atoms
    NInds = p.AtomInd(AtomName = 'N')
    CInds = p.AtomInd(AtomName = 'C')
    OInds = p.AtomInd(AtomName = 'O')
    SInds = p.AtomInd(AtomName = 'S')
    # generate approx carbonyl groups for all but last residue
    s = ''
    n = 0
    CurrentChain = -1
    i_gly = 0
    for i, r in enumerate(Seq):
        # determine if a new chain starts here 
        thisChain = p.ResChain(i)
        if not thisChain == CurrentChain:
            CurrentChain = thisChain # update Chain number
            if not i == 0: s +=  "TER\n" # inter-chain TER records
        # Amide nitrogen
        PosN = Pos[NInds[i]]
        s += PDBFMT % (n+1, 'N ', r, string.ascii_uppercase[thisChain], i+1, PosN[0], PosN[1], PosN[2], 1.0, 0.0)   
        s += '\n'
        n += 1
        # Alpha carbon
        PosCA = Pos[CInds[i]]
        s += PDBFMT % (n+1, 'CA ', r, string.ascii_uppercase[thisChain], i+1, PosCA[0], PosCA[1], PosCA[2], 1.0, 0.0)
        s += '\n'
        n += 1
        # CG oxygen site
        PosCGO = Pos[OInds[i]]
        # find effective carbonyl
        if i < len(Seq) - 1:
            PosNextN = Pos[NInds[i+1]]
            PosC, PosO = Project(PosCA, PosCGO, PosNextN)
            # Carbonyl C
            s += PDBFMT % (n+1, 'C ', r, string.ascii_uppercase[thisChain], i+1, PosC[0], PosC[1], PosC[2], 1.0, 0.0)
            s += '\n'
            n += 1
            # Carbonyl O 
            s += PDBFMT % (n+1, 'O ', r, string.ascii_uppercase[thisChain], i+1, PosO[0], PosO[1], PosO[2], 1.0, 0.0)
            s += '\n'
            n += 1
        else:
            # retain CG O site for last residue
            s += PDBFMT % (n+1, 'O ', r, string.ascii_uppercase[thisChain], i+1, PosCGO[0], PosCGO[1], PosCGO[2], 1.0, 0.0)
            s += '\n'
            n += 1
        # Sidechain
        if not r == 'GLY' or hasPseudoGLY:
            PosS = Pos[SInds[i_gly]]
            s += PDBFMT % (n+1, 'S ', r, string.ascii_uppercase[thisChain], i+1, PosS[0], PosS[1], PosS[2], 1.0, 0.0)
            s += '\n'
            n += 1
            i_gly += 1
    # write AA pdb
    s += 'TER\n' # terminal record
    OutPdb = Prefix + '.pdb'
    with open(OutPdb, 'w') as of: of.write(s)



#### COMMAND LINE USAGE ####
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'AA-->CG mapper for Pdb files and LAMMPS / AMBER trajectories')
    # common required args
    parser.add_argument('InPdb', help = 'AA or CG pdb')
    parser.add_argument('OutPrefix', help = 'prefix of outputs: for trajectories, saves in current dir')
    # common optional args
    parser.add_argument('-op', '--operation', help = 'map: AA-CG map of Pdb file, backmap: CG-AA map of Pdb file, map2poly: AA-CG map of Pdb to polymer backbone')
    parser.add_argument('-gly', '--hasgly', action = 'store_true', help = 'has Pseudo GLY?')
    parser.add_argument('-m', '--model', type = int, help = 'Model number, starts from 1')
    # optional args for trajectories
    parser.add_argument('--aatraj', '--aatraj', help = 'AA trajectory, can be Lammps or Amber')
    parser.add_argument('-prmtop', '--prmtop', help = 'PrmTop file to parse Amber trajectory')
    parser.add_argument('-ene', '--ene', help = 'Amber Ene file')
    parser.add_argument('-n', '--nframes', type = int, default = 0, help = 'Work with last n frames')
    # optional args for mapping to a polymer
    parser.add_argument('-poly', '--poly', help = 'Amino acid of polymer backbone to map to')
    parser.add_argument('-aapdb', '--aapdb', help = 'AA Pdb for the CG pdb that is to be mapped to a polymer')
    
    args = parser.parse_args()
    InPdb = os.path.abspath(args.InPdb)
    OutPrefix = os.path.abspath(args.OutPrefix)
    Op = args.operation
    hasPseudoGLY = args.hasgly
    Model = args.model
    AATraj = args.aatraj
    PrmTop = args.prmtop
    AmberEne = args.ene
    LastNFrames = args.nframes
    PolyName = args.poly
    AAPdb = args.aapdb
    
    # default
    if Op is None: exit()
    
    # map
    if Op == 'map':
        Map(InPdb = InPdb, CGPrefix = OutPrefix, Model = Model, 
            AATraj = AATraj, PrmTop = PrmTop, AmberEne = AmberEne, 
            LastNFrames = LastNFrames, hasPseudoGLY = hasPseudoGLY)
    
    if Op == 'map2poly':
        Map2Polymer(Pdb = InPdb, AAPdb = AAPdb, PolyName = PolyName.upper(), 
                    Model = Model, MappedPrefix = OutPrefix, 
                    hasPseudoGLY = hasPseudoGLY) 
    
    if Op == 'backmap':
        ReverseMap(CGPdb = InPdb, Prefix = OutPrefix, Model = Model, 
                    hasPseudoGLY = hasPseudoGLY)
