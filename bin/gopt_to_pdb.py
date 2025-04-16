#!/usr/bin/env python3
"""
This is a program written by Qianyi Cheng
at University of Memphis.
"""

from rms import *
from numpy import *
import sys, re, os
from read_write_pdb import *
import argparse
from read_gout import *

if __name__ == '__main__':
    """ Usage: gopt_to_pdb.py -gout 1.out -pdb ../template.pdb """
    parser = argparse.ArgumentParser(description='generate pdbfiles from 1.out')
    parser.add_argument('-gout', dest='output',default='1.out',help='gaussian output file')
    parser.add_argument('-pdb', dest='pdbf',default='../template.pdb',help='template pdb file')
    parser.add_argument('-frame', dest='frame',default=-1,help='select frame: -1 (final=default) or integer or \"all\"')
    parser.add_argument('-name', dest='name',default=None,help='filename for output pdb')

    args = parser.parse_args()
    pdbf = args.pdbf
    output = args.output
    if args.name != None:
        fname = args.name.replace('.pdb','')
    else:
        fname = args.name

    pdb, res_info, tot_charge = read_pdb(pdbf)
#    map, xyz_i = get_ca(pdb)
    map, xyz_i = get_fatom(pdb)
    
    natoms = len(pdb)
    
    with open(output) as f:
        lines = f.readlines()
    rot_opt = gaussian_opt_xyz(lines,natoms)
    if args.frame == 'all':
        for key in range(len(rot_opt)):
            xyz_c = array(rot_opt[key])
            (c_trans,U,ref_trans) = rms_fit(xyz_i,xyz_c[map])
            xyz_n = dot( xyz_c-c_trans, U ) + ref_trans
            sel_atom = update_xyz(pdb,xyz_n)
            if fname != None:
                name = f'{fname}.f{key:03d}.pdb'
            else:
                name = f'gopt.f{key:03d}.pdb'
            write_pdb(name,sel_atom)
    elif int(args.frame) >= -1 and int(args.frame) <= len(rot_opt):
        key = int(args.frame)
        xyz_c = array(rot_opt[key])
        (c_trans,U,ref_trans) = rms_fit(xyz_i,xyz_c[map])
        xyz_n = dot( xyz_c-c_trans, U ) + ref_trans
        sel_atom = update_xyz(pdb,xyz_n)
        if key == -1:
            if fname != None:
                name = fname+'.pdb'
            else:
                name = 'gopt.pdb'
            write_pdb(name,sel_atom)
        else:
            if fname != None:
                name = f'{fname}.f{args.frame:03d}.pdb'
            else:
                name = f'gopt.f{args.frame:03d}.pdb'
            write_pdb(name,sel_atom)
    else:
        print("The frame is not clear!")
        sys.exit()
