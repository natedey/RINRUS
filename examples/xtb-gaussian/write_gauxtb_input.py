#!/usr/bin/env python3

### This is a script for writting xtb-gaussian input file and 1.fix file written by Qianyi Cheng ###

import os, sys, re, filecmp
from numpy import *
import argparse
from read_write_pdb import *
from glob import glob


def system_run(cmd):
    print(cmd)
    exit = os.system(cmd)
    if ( exit != 0 ):
        print('failed to run:')
        print(cmd)
        sys.exit()
 
### copy h-added pdb xyz and other information into tmppdb ###
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


def write_input(inp_name,inp_temp,charge,multiplicity,pic_atom,tot_charge,res_count):
    ### inp_name default can be 1.inp, but first is name.input such as 9.input
    ### inp_type list which includes [small/large, level_theory, basis, opt, freq,  
    ### input_template line0: size
    ### input_template line1: level
    ### input_template line2: opt + calcfc/readfc or opt + modred + info
    ### input_template line3: freq
    ### input_template line4: guess=read
    ### input_template line5: geom=checkpoint
    ### input_template line6: iop
    ### input_template line7: scf_info
    ### input_template line8: scrf_info

    print("charge= %d, tot_charge_without'-c'= %d"%(charge,tot_charge))
    with open(inp_temp) as f:
        lines = f.readlines()

    inp = open('%s'%inp_name,'w')
    for i in range(3):
        inp.write('%s'%lines[i])

    inp.write("#P ")      

    f_atom = []
    if lines[3][0] != '#' and 'modred' not in lines[3]:
        inp.write('%s '%lines[3][:-1])
    ### If modred ###
    elif lines[3][0] != '#' and 'modred' in lines[3]:
        optl = lines[3].split()
        inp.write("%s "%optl[0])
        ### get frozen atom info ###
        modred_info, modred_code = optl[1:3]
        pairs = modred_info.split(';')
        for lp in range(len(pairs)):
            lmod = []
            pair = pairs[lp]
            ar = pair.split(',')
            ### Atom name and Atom id in pair ###
            for i in range(0,len(ar),2):
                for atom in pic_atom:
                    if ar[i] == atom[2].strip() and int(ar[i+1]) == atom[6]:
                        lmod.append(pic_atom.index(atom)+1)
            f_atom.append(lmod)
    for l in range(4,8):
        if lines[l][0] != '#':
            inp.write("%s "%lines[l].strip())

    inp.write("\n\n")
    ### info_line ###
    inp.write("%s"%lines[8])
    inp.write("\n")
    inp.write("%d %d\n"%(charge+tot_charge,multiplicity))

    for atom in pic_atom:
        inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),atom[16],atom[8],atom[9],atom[10])) 
    inp.write("\n")

    if len(f_atom) != 0:
        for l in range(len(f_atom)):
            for ids in f_atom[l]:
                inp.write("%s "%ids)
            inp.write("%s\n"%modred_code)
        inp.write("\n")

    inp.close()

    ### write 1.fix file ###
    f = open("1.fix",'w')
    atom_idx = []
    for i in range(len(pic_atom)):
        if pic_atom[i][-1] == '-1':
            atom_idx.append(i+1)
    
    f.write('$hess\n')    
    for ids in atom_idx:
        f.write('    scale mass: %d, 90000000\n'%ids)
    f.close()
    ### else can be constraint ###



if __name__ == '__main__':

    #########################################################################################################################################################
    ### In working directory such as model9-ts-001 
    ### Run: write_input.py -noh nohpdb -adh haddpdb -intmp relaxh_temp -c 2 (final_pdb is saved as "template.pdb", input is saved as "1.inp", default -m 1) 
    #########################################################################################################################################################

    parser = argparse.ArgumentParser(description='Prepare template PDB files, write input files, save output PDB files in working directory',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-step', dest='step', default=0, type=int, 
            help='step0: read noh, addh pdbs, write_final_pdb and read input_template write_first_inp') 
    parser.add_argument('-wdir', dest='output_dir', default=os.path.abspath('./'), help='working dir')
    parser.add_argument('-tmp', dest='tmp_pdb', default=None, help='template_pdb_file')
    parser.add_argument('-noh', dest='no_h_pdb', default=None, help='trimmed_pdb_file')
    parser.add_argument('-adh', dest='h_add_pdb', default=None, help='hadded_pdb_file')
    parser.add_argument('-intmp', dest='input_tmp', default=None, help='template_for_write_input')
    parser.add_argument('-m', dest='multiplicity', default=1, type=int, help='multiplicity')
    parser.add_argument('-c', dest='ligand_charge', default=0, type=int, help='charge_of_ligand')
    parser.add_argument('-inpn', dest='inp_name', default='1.inp', help='input_name')
#    parser.print_help()

    args = parser.parse_args()

    step = args.step
    wdir = args.output_dir
    if args.tmp_pdb is None:
        tmp_pdb = '%s/template.pdb'%wdir
    else:
        tmp_pdb   = args.tmp_pdb

    nohpdb   = args.no_h_pdb
    adhpdb   = args.h_add_pdb
    int_tmp  = args.input_tmp
    multi    = args.multiplicity
    charge   = args.ligand_charge
    inp_name = args.inp_name
    wdir     = args.output_dir

    if step == 0:
        pic_atom, tot_charge = pdb_after_addh(nohpdb,adhpdb)
        res_count = adhpdb.split('_')[1]
        write_pdb('%s'%tmp_pdb,pic_atom,res_count)

    write_input('%s/%s'%(wdir,inp_name),int_tmp,charge,multi,pic_atom,tot_charge,res_count)
