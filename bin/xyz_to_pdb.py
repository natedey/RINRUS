#!/usr/bin/env python3
"""
This is a program created by Dr. Dominique Wappett and the DeYonker research group
at The University of Memphis.
Date created 10.04.2024
"""


import numpy as np
import sys, re, os
from read_write_pdb import *
import argparse

if __name__ == '__main__':
    """ Usage: xyz_to_pdb.py -xyz orca.xyz -pdb template.pdb """
    parser = argparse.ArgumentParser(description='generate pdbfiles from orca xyz file')
    parser.add_argument('-xyz', dest='xyz',default='orca.xyz',help='xyz structure file')
    parser.add_argument('-pdb', dest='pdbf',default='template.pdb',help='template pdb file')
    parser.add_argument('-frame', dest='frame',default=-1,help='select frame: -1 (final=default) or integer (count starts at 0) or \"all\"')
    parser.add_argument('-name', dest='name',default=None,help='filename for output pdb')
    # NOTE: FRAME NUMBERING STARTS AT 0 TO MATCH PYTHON AND ORCA NUMBERING #

    args = parser.parse_args()
    output = args.xyz

    ### if no name specified, use output file name
    if args.name:
        fname = args.name.replace('.pdb','')
    else:
        fname = output.split('/')[-1].replace('.xyz','')

    ### read xyz and pdb files
    pdb, res_info, tot_charge = read_pdb(args.pdbf)
    xyzf = open(output,'r').readlines()

    ### check xyz and pdb length match
    numline = xyzf[0]
    numcheck = int(numline.replace('\n',''))
    if numcheck != len(pdb):
        print('pdb and xyz lengths do not match, please check if specified pdb is correct.')
        sys.exit()

    ### extract set(s) of coords
    strucs = []
    for line in xyzf:
        # when each header line reached, save previous struc to set and then reset list to collect the next one
        if line == numline:
            try:
                strucs.append(newstruc)
                newstruc = []
            except:
                newstruc = []
            continue
        # only collect lines that contain coords, skip everything else
        try:
            line = line.replace('\n','').split()
            crd = [float("%.3f"%float(line[1])),float("%.3f"%float(line[2])),float("%.3f"%float(line[3]))]
            newstruc.append(crd)
        except:
            continue
    # add final struc to set
    strucs.append(newstruc)

    ### get chosen frame
    if args.frame == 'all':
        for i in range(len(strucs)):
            iname = f'{fname}.f{i:03d}.pdb'
            sel_atom = update_xyz(pdb,strucs[i])
            write_pdb(iname,sel_atom)
    elif int(args.frame) == -1:
        iname = fname+'.pdb'
        sel_atom = update_xyz(pdb,strucs[-1])
        write_pdb(iname,sel_atom)
    elif int(args.frame) >= 0 and int(args.frame) <= len(strucs):
        iname = f'{fname}.f{args.frame:03d}.pdb'
        sel_atom = update_xyz(pdb,strucs[int(args.frame)])
        write_pdb(iname,sel_atom)
    else:
        print("Chosen frame doesn't exist!")

