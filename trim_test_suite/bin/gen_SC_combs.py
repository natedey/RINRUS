#!/usr/bin/env python3

import os, sys, re
import itertools
import argparse

### Side chain subsets ###
res_atoms_sc = {
'Ala': ['CB', 'HB1', 'HB2', 'HB3'],
'Arg': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22'],
'Asn': ['CB', 'HB2', 'HB3', 'CG', 'OD1', 'ND2', 'HD21', 'HD22'],
'Asp': ['CB', 'HB2', 'HB3', 'CG', 'OD1', 'OD2'],
'Cys': ['CB', 'HB2', 'HB3', 'SG', 'HG'],
'Gln': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'NE2', 'HE21', 'HE22'],
'Glu': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'OE2'],
'Gly': [],
'His': ['CB', 'HB2', 'HB3', 'CG', 'ND1', 'CD2', 'HD2', 'CE1', 'HE1', 'NE2', 'HE2'],
'Ile': ['CB', 'HB', 'CG1', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'CD1','HD11', 'HD12', 'HD13'],
'Leu': ['CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23'],
'Lys': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ1', 'HZ2', 'HZ3'],
'Met': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'SD', 'CE', 'HE1', 'HE2', 'HE3'],
'Phe': ['CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ', 'HZ'],
'Pro': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3'],
'Ser': ['CB', 'HB2', 'HB3', 'OG', 'HG'],
'Thr': ['CB', 'HB', 'OG1', 'HG1', 'CG2', 'HG21', 'HG22', 'HG23'],
'Trp': ['CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CD2', 'NE1', 'HE1', 'CE2', 'CE3', 'HE3', 'CZ2', 'HZ2', 'HZ3', 'CZ3', 'CH2', 'HH2'],
'Tyr': ['CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ', 'OH', 'HH'],
'Val': ['CB', 'HB', 'CG1', 'HG11', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23']
}

### Function for getting all combinations ###
def SCcombs(residue):
    s = list(res_atoms_sc[residue])
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(1,len(s)+1))


if __name__ == '__main__':
    """ Usage: generate combinations file for given residue """
    parser = argparse.ArgumentParser(description='generate atom combinations')
    parser.add_argument('-res', dest='residue',default='None',help='residue')
    parser.add_argument('-file', dest='path',default='SC.dat',help='file name')

    args = parser.parse_args()
    res = args.residue
    outpath = args.path

    combs = list(SCcombs(res))

    fl = open('%s'%outpath,'w')
    for i in combs:
        fl.write(' '.join(map(str,i))+'\n')
    fl.close()

