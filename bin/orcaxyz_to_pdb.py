#!/usr/bin/env python3

import numpy as np
import sys, re, os
from read_write_pdb import *
import argparse

if __name__ == '__main__':
    """ Usage: orcaxyz_to_pdb.py -xyz orca.xyz -pdb template.pdb """
    parser = argparse.ArgumentParser(description='generate pdbfiles from orca xyz file')
    parser.add_argument('-xyz', dest='output',default='orca.xyz',help='output xyz file')
    parser.add_argument('-pdb', dest='pdbf',default='template.pdb',help='template pdb file')

    args = parser.parse_args()
    pdbf = args.pdbf
    output = args.output

    pdb, res_info, tot_charge = read_pdb(pdbf)

    lines = np.loadtxt(output, skiprows=2, usecols=range(1,4), dtype=float)
    
    xyz = []
    for atom in lines:
        xyz.append([float("%.3f"%atom[0]),float("%.3f"%atom[1]),float("%.3f"%atom[2])])

    if "/" in output:
        output = output.split('/')[-1]
    name = output.replace('.xyz','.pdb')

    sel_atom = update_xyz(pdb,xyz)
    write_pdb(name,sel_atom)

    
