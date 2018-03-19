#!/usr/bin/env python

''' provides routines for parsing native and non-native contact information
    from supplied Pdb structure
'''

import numpy as np
import protein
from const import *

Verbose = True

def ParsePdb(p):
    ''' must be supplied a cg protein object that is linked to a Pdb'''
    if Verbose: print 'Parsing native structure...'
    ResPos = p.p0.ResPos()
    Pos = p.p0.Pos
    ResContactList = p.p0.ResContactList()
    x_native = []
    d_native = []
    d_ss_native = []
    x_nonnative = []
    d_nonnative = []
    d_ss_nonnative = []
    fn = open('nativecontact.txt', 'w')
    fnn = open('nonnativecontact.txt', 'w')
    fn.write('#Res_i Res_j d_com d_ss\n')
    fnn.write('#Res_i Res_j d_com d_ss\n')
    for i in range(len(p.Seq)-1):
        for j in range(i+1, len(p.Seq)):
            # ignore adjacent residues
            if abs(i-j) < MinCO: continue
            # use alpha-carbon if glycine
            idx_i = 1 if p.Seq[i] == 'GLY' else 3
            idx_j = 1 if p.Seq[j] == 'GLY' else 3
            m = p.StartAtomInds[i] + idx_i
            n = p.StartAtomInds[j] + idx_j
            # get com distance between residues
            d_com_ij = ResPos[j] - ResPos[i]
            d_com = np.sqrt(np.sum(d_com_ij * d_com_ij))
            # get sidechain distance between residues
            d_ss_mn = Pos[n] - Pos[m]
            d_ss = np.sqrt(np.sum(d_ss_mn * d_ss_mn))
            # native contacts
            if ResContactList.__contains__( (i,j) ):
                x_native.append((i, j))
                d_native.append(d_com)
                d_ss_native.append(d_ss)
                fn.write('%3d %3d %3.2f %3.2f\n' % (i, j, d_com, d_ss))
            # non-native contacts
            else:
                x_nonnative.append((i, j))
                d_nonnative.append(d_com)
                d_ss_nonnative.append(d_ss)
                fn.write('%3d %3d %3.2f %3.2f\n' % (i, j, d_com, d_ss))
    fn.close()
    fnn.close()
    # convert to arrays
    d_native = np.array(d_native)
    d_ss_native = np.array(d_ss_native)
    d_nonnative = np.array(d_nonnative)
    d_ss_nonnative = np.array(d_ss_nonnative)
    # output data struct
    ContactDict = {'ResContactList' : ResContactList,
                   'c_native'       : x_native,
                   'c_nonnative'    : x_nonnative,
                   'd_native'       : d_native,
                   'd_ss_native'    : d_ss_native,
                   'd_nonnative'    : d_nonnative,
                   'd_ss_nonnative' : d_ss_nonnative
                   }
    return ContactDict

def makeMatrix(p, ResContactList, Sys, includeGLY = False):
    ''' makes a sim style matrix for native and non-native contacts
        must be supplied a cg protein object'''
    if Verbose: print 'Generating native and non-native filters...'
    NAID = Sys.World.NAID
    NativePairs = np.zeros([NAID, NAID], int)
    NonNativePairs = np.zeros([NAID, NAID], int)
    Topo2AID = {}
    for i in range(p.NRes-1):
        for j in range(i+1, p.NRes):
            # glycines
            if p.Seq[i] == 'GLY' or p.Seq[j] == 'GLY':
                # include glycines
                if includeGLY:
                    if Verbose: print 'Native Contact Pair: (%d, %d) Retaining C-alpha for Glycines' % (i,j)
                    idx_i = 1 if p.Seq[i] == 'GLY' else 3
                    idx_j = 1 if p.Seq[j] == 'GLY' else 3
                # don't include glycines
                else:
                    if Verbose: print 'Native Contact Pair: (%d, %d) Ignoring pairs with Glycines' % (i,j)
                    continue
            # other residues
            else:
                if Verbose: print 'Native Contact Pair: (%d, %d)' % (i,j)
                idx_i = 3
                idx_j = 3
            # populate matrix
            m = Sys.Mol[0][p.StartAtomInds[i] + idx_i].AID
            n = Sys.Mol[0][p.StartAtomInds[i] + idx_j].AID
            # native contacts
            if ResContactList.__contains__((i,j)):
                NativePairs[m, n] = 1
                NativePairs[n, m] = 1
            # non-native contacts
            else:
                NonNativePairs[m, n] = 1
                NonNativePairs[n, m] = 1
            # map AIDs
            Topo2AID[ (i,j) ] = (m, n)
    return NativePairs, NonNativePairs, Topo2AID

def ErodeNativeContact(ContactDict, Frac = 1.0):
    pass




