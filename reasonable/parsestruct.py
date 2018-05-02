#!/usr/bin/env python

''' provides routines for parsing native and non-native contact information
    from supplied Pdb structure
'''

import numpy as np, copy
import protein
import reasonable.topo as topo
from const import *

Verbose = True

def ParsePdb(p, ResContactList = None):
    ''' must be supplied a cg protein object that is linked to a Pdb
    assumes that sidechains are referenced by residue number
    Can supply a ResContactList to force using those native contacts'''
    if Verbose: print 'Parsing native structure...'
    if not ResContactList is None:
        if Verbose: print 'Not calculating ResContactList. Using supplied.'
    ResPos = p.GetResPos()
    Pos = p.Pos
    if ResContactList is None: ResContactList = p.GetResContactList()
    SInds = p.GetSInds()
    c_native = []
    d_native = []
    d_ss_native = []
    c_nonnative = []
    d_nonnative = []
    d_ss_nonnative = []
    fn = open(p.Prefix + '_nativecontact.txt', 'w')
    fnn = open(p.Prefix+'_nonnativecontact.txt', 'w')
    fn.write('#Res_i Res_j d_com d_ss\n')
    fnn.write('#Res_i Res_j d_com d_ss\n')
    s = 'Used Alpha-Carbon as side chain for residues: '
    s0 = copy.copy(s)
    for i in range(len(p.Seq)-1):
        for j in range(i+1, len(p.Seq)):
            # ignore adjacent residues
            if abs(i-j) < MinCO: continue
            # check if alpha-Carbons need to be included for residues
            # with missing sidechains (usually glycine)
            res_i = p.Seq[i]
            res_j = p.Seq[j]
            m = SInds[i]
            n = SInds[j]
            if p.AtomSbyNum[i] is None:
                this_s = '%3d (%3s)' % (i, res_i)
                if not this_s in s: s += this_s + '  '
            if p.AtomSbyNum[j] is None:
                this_s = '%3d (%3s)' % (j, res_j)
                if not this_s in s: s += this_s + '  '
            # get com distance between residues
            d_com_ij = ResPos[j] - ResPos[i]
            d_com = np.sqrt(np.sum(d_com_ij * d_com_ij))
            # get sidechain distance between residues
            d_ss_mn = Pos[n] - Pos[m]
            d_ss = np.sqrt(np.sum(d_ss_mn * d_ss_mn))
            # native contacts
            if (i,j) in ResContactList:
                c_native.append((i, j))
                d_native.append(d_com)
                d_ss_native.append(d_ss)
                fn.write('%3d %3d %3.2f %3.2f\n' % (i, j, d_com, d_ss))
            # non-native contacts
            else:
                c_nonnative.append((i, j))
                d_nonnative.append(d_com)
                d_ss_nonnative.append(d_ss)
                fnn.write('%3d %3d %3.2f %3.2f\n' % (i, j, d_com, d_ss))
    fn.close()
    fnn.close()
    if s == s0: s += 'None'
    if Verbose: print s
    # convert to arrays
    d_native = np.array(d_native)
    d_ss_native = np.array(d_ss_native)
    d_nonnative = np.array(d_nonnative)
    d_ss_nonnative = np.array(d_ss_nonnative)
    # output data struct
    ContactDict = {'ResContactList' : ResContactList,
                   'c_native'       : c_native,
                   'c_nonnative'    : c_nonnative,
                   'd_native'       : d_native,
                   'd_ss_native'    : d_ss_native,
                   'd_nonnative'    : d_nonnative,
                   'd_ss_nonnative' : d_ss_nonnative
                   }
    return ContactDict

def makeMatrix(p, ContactDict, Sys):
    ''' makes a sim style matrix for native and non-native contacts
        must be supplied a cg protein object
        assumes sidechains are referenced by residue number'''
    NAID = Sys.World.NAID
    NativePairs = np.zeros([NAID, NAID], int)
    NonNativePairs = np.zeros([NAID, NAID], int)
    Topo2AID_Native = {}
    Topo2AID_NonNative = {}
    SInds = p.GetSInds()
    # native contacts
    if Verbose: print 'Generating native filters...'
    for k, (i,j) in enumerate(ContactDict['c_native']):
        # check if a residue has a missing sidechain (usually glycines)
        # such residues are ignored in making the contact matrix
        # note that they are still present in the rescontactlist
        ignoreThisPair = False
        res_i = p.Seq[i]
        res_j = p.Seq[j]
        s = ' Native contact pair: (%3d, %3d) (%3s, %3s)' % (i, j, res_i, res_j)
        if p.AtomSbyNum[i] is None or p.AtomSbyNum[j] is None:
            ignoreThisPair = True
            s+= '  Ignored'
        if Verbose: print s
        if ignoreThisPair: continue
        # populate matrix
        m = Sys.Mol[0][SInds[i]].AID
        n = Sys.Mol[0][SInds[j]].AID
        NativePairs[m, n] = 1
        NativePairs[n, m] = 1
        # store the res --> AID reference
        Topo2AID_Native[ (i,j) ] = (m, n)
    
    # non native contacts
    if Verbose: print '\nGenerating non-native filters...'
    for k, (i,j) in enumerate(ContactDict['c_nonnative']):
        # check if a residue has a missing sidechain (usually glycines)
        # such residues are ignored in making the contact matrix
        # note that they are still present in the rescontactlist
        ignoreThisPair = False
        res_i = p.Seq[i]
        res_j = p.Seq[j]
        s = ' Non-Native contact pair: (%3d, %3d) (%3s, %3s)' % (i, j, res_i, res_j)
        if p.AtomSbyNum[i] is None or p.AtomSbyNum[j] is None:
            ignoreThisPair = True
            s+= '  Ignored'
        if Verbose: print s
        if ignoreThisPair: continue
        # populate matrix
        m = Sys.Mol[0][SInds[i]].AID
        n = Sys.Mol[0][SInds[j]].AID
        NonNativePairs[m, n] = 1
        NonNativePairs[n, m] = 1
        # store the res --> AID reference
        Topo2AID_NonNative[ (i,j) ] = (m, n)
    return NativePairs, NonNativePairs, Topo2AID_Native, Topo2AID_NonNative

            

    




