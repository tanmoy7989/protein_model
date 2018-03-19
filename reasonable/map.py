'''/*
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

#### USAGE ####
HelpStr = '''python map.py InPdb [CGPrefix] [AATraj] [PrmTop] [AmberEne] [LastNFrames]'''
if len(sys.argv) == 1:
    print HelpStr
    exit()

genBonds = 2 # 0 for no bonds, 1 for only BB bonds, 2 for all

NCaps = ['ACE']
CCaps = ['NME', 'NHE']
PDBFMT = "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f%6.2f"        # (atomind + 1, atomname, resname, reschainind, resnum+1, x, y, z)
BONDFMT = "%6s%5d%5d\n"                                                  # ('CONECT', a + r.StartAtom + 1, b + r.StartAtom + 1)

# Inputs
InPdb = os.path.abspath(sys.argv[1])
CGPrefix = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else 'testcg'
AATraj = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else None
PrmTop = os.path.abspath(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] else None
AmberEne = os.path.abspath(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] else None
LastNFrames = int(sys.argv[6]) if len(sys.argv) > 6 else 0 

# read in protein structure 
p = protein.ProteinClass(Pdb = InPdb)
p.Update()
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
if not 'GLY' in ResCount.keys():
    NCGAtoms = 4 * len(Seq)
else:
    NCGAtoms = 3 * ResCount['GLY'] + 4 * ( sum(ResCount.values()) - ResCount['GLY'] )

# Masks
DecapFilter = lambda ResName, AtomName : not ( (NCaps + CCaps).__contains__(ResName) )
SFilter = lambda ResName, AtomName : not (AtomName == 'N' or AtomName == 'CA' or AtomName == 'C' or AtomName == 'O' or AtomName.__contains__('H'))
NInds = p.AtomInd(AtomName = 'N', UserFunc = DecapFilter)
CAInds = p.AtomInd(AtomName = 'CA', UserFunc = DecapFilter)
CInds = p.AtomInd(AtomName = 'C', UserFunc = DecapFilter)
OInds = p.AtomInd(AtomName = 'O', UserFunc = DecapFilter)
SInds = dict( (i, []) for i in range(len(Seq)) )
for i, r in enumerate(Seq):
	startres = 1 if hasNCap else 0
	SInds[i] = p.AtomInd(ResNum = startres + i, UserFunc = DecapFilter and SFilter)

# Write Masks to a Map and final CG atoms to CG PDB
s = ''
s_bond = ''
Map = dict( (i, []) for i in range(NCGAtoms) )
n = 0
for i, r in enumerate(Seq):
    # N lines
    if not i == 0: 
        if Seq[i-1] == 'GLY': s_bond += BONDFMT % ('CONECT', n, n+1)
        else: s_bond += BONDFMT % ('CONECT', n-1, n+1)
    N_CGInd = n ; n+= 1
    Map[N_CGInd].append(NInds[i])
    s += PDBFMT % (N_CGInd+1, 'N  ', r, ' ', i+1, Pos[NInds[i], 0], Pos[NInds[i], 1], Pos[NInds[i], 2], 1.0, 0.0)
    s += "\n"
    
    # C lines
    s_bond += BONDFMT % ('CONECT', n, n+1)
    C_CGInd = n ; n+= 1
    Map[C_CGInd].append(CAInds[i])
    s += PDBFMT % (C_CGInd+1, 'C  ', r, ' ', i+1, Pos[CAInds[i], 0], Pos[CAInds[i] ,1], Pos[CAInds[i] ,2], 1.0, 0.0)
    s += "\n"
    
    # O lines
    s_bond += BONDFMT % ('CONECT', n, n+1)
    O_CGInd = n; n+=1
    Map[O_CGInd].extend([ CInds[i], OInds[i] ])
    COMPos = (Pos[CInds[i]] + Pos[OInds[i]]) / 2.
    s += PDBFMT % (O_CGInd+1, 'O  ', r, ' ', i+1, COMPos[0], COMPos[1], COMPos[2], 1.0, 0.0)
    s += "\n"
    
    # S lines
    if not r == 'GLY':
        if genBonds == 2:
            s_bond += BONDFMT % ('CONECT', n-1, n+1)
        S_CGInd = n; n+= 1
        Map[S_CGInd].extend(SInds[i])
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
	simMap += [sim.atommap.AtomMap(Atoms1 = Map[i], Atom2 = i)]
    AtomNames = []
    for i, r in enumerate(Seq):
        if r == 'GLY': AtomNames.extend(['N', 'C', 'O'])
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
