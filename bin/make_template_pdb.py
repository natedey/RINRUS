#!/usr/bin/env python3
"""
This is a program created by Dr. Dominique Wappett and the DeYonker research group
at The University of Memphis.
Date created 01.15.2025
"""

import os, sys, re
from numpy import *
import argparse
from read_write_pdb import *


def pdb_after_addh(tmppdb,newpdb):
    tmp_pdb, res_info, tot_charge_t = read_pdb(tmppdb)
    tmp_xyz = []
    for i in tmp_pdb:
        tmp_xyz.append([i[8],i[9],i[10]])
    new_pdb, binfo, tot_charge = read_pdb(newpdb)     #can be just xyz files from cerius or pymol
    pic_atom = []
    for line in new_pdb:
        if [line[8],line[9],line[10]] not in tmp_xyz:
            line[15] = '0 '
            pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],' H',line[15],line[16]])
        else:
            if '+' in line[15] or '-' in line[15]:
                charge = line[15]
            else:
                charge = '0 '
            idx = tmp_xyz.index([line[8],line[9],line[10]])
            line = tmp_pdb[idx]
            line[15] = charge
            if [line[2],line[5],line[6]] in res_info:
                pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'])
            else:
                pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16]])
    return pic_atom, tot_charge #, xyz, atom, hold

if __name__ == '__main__':
    ### create template pdb files for model
    ### two ways of specifying input
    ### 1) make_template_model.py -model N                            (assumes that files have default res_N.pdb and res_N_h.pdb naming format)
    ### 2) make_template_model.py -noh res_N.pdb -addh res_N_h.pdb    (specify each separately if format is different)
    parser = argparse.ArgumentParser(description='Prepare template PDB files',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-model', dest='modname', default=None, help='model number')
    parser.add_argument('-noh', dest='no_h_pdb', default=None, help='trimmed_pdb_file')
    parser.add_argument('-addh', dest='h_add_pdb', default=None, help='hadded_pdb_file')

    args = parser.parse_args()

    if args.no_h_pdb == None and args.h_add_pdb == None and args.modname != None:
        nohpdb = 'res_'+args.modname+'.pdb'
        adhpdb = 'res_'+args.modname+'_h.pdb'
    elif args.no_h_pdb != None and args.h_add_pdb != None and args.modname == None:
        nohpdb = args.no_h_pdb
        adhpdb = args.h_add_pdb
    elif args.no_h_pdb == None and args.h_add_pdb == None and args.modname == None:
        print("No model specified!")
    else:
        print("Mixed input syntax! Use either -model [num] OR -noh [pdb] -addh [pdb]")

    pdbnum = nohpdb.replace('res_','').replace('.pdb','')
    tmp_pdb = 'model_'+pdbnum+'_template.pdb'

    pic_atom, tot_charge = pdb_after_addh(nohpdb,adhpdb)
    write_pdb(tmp_pdb,pic_atom)
