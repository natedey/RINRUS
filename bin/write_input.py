#!/usr/bin/env python3
"""
This is a program written by Qianyi Cheng 
at university of memphis.
Date 8.10.2022
"""
import os, sys, re, filecmp
from numpy import *
import argparse
from read_write_pdb import *
from glob import glob
from input_suite import *
from pathlib import Path


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

########################################################################
### replace certain part of the pdb with provided pdb/xyz ##############
### tmppdb is the minimized structure, most of these xyz will be kept
### newpdb is the one contain transition structure/fragment
### parts is the residue name and atom name for the transition part
########################################################################
def pdb_replace(tmppdb,newpdb,parts):
    tmp_pdb, res_info, tot_charge_t = read_pdb(tmppdb)
    tmp_xyz = []
    new_xyz = []
    for i in tmp_pdb:
        tmp_xyz.append([i[6],i[2].strip()])
    new_pdb, binfo, tot_charge = read_pdb(newpdb)     #can be just xyz files from cerius or pymol
    for i in new_pdb:
        new_xyz.append([i[6],i[2].strip()])

    if parts == None:   #newpdb has the entire thing to replace the tmppdb
        for line in new_pdb:
            resatom = [line[6],line[2].strip()]
            idx = tmp_xyz.index(resatom)
            tmp_pdb[idx] = line
    else:
        for i in range(len(parts)):
            idx1 = new_xyz.index([int(parts[i][0]),parts[i][1]])
            idx2 = tmp_xyz.index([int(parts[i][0]),parts[i][1]])
            tmp_pdb[idx2] = new_pdb[idx1]

    return tmp_pdb, tot_charge_t #, xyz, atom, hold


def gen_pdbfiles(wdir,step,tmppdb):
    new_dir = '%s/step%spdbs'%(wdir,step)
    if os.path.isdir(new_dir):
        system_run('rm -r %s'%new_dir)
        system_run('mkdir %s'%new_dir)
    else:
        system_run('mkdir %s'%new_dir)
    os.chdir(new_dir)
    system_run('python3 $HOME/git/RINRUS/bin/gopt_to_pdb.py -p %s -o %s/step-%s-out -f -1'%(tmppdb,wdir,step))
    os.chdir('%s'%wdir)
    pdb_name = []
    for pdbf in glob('%s/*.pdb'%new_dir):
        m = re.search(r'(\d+).pdb',pdbf)
        pdb_name.append( int(m.group(1)) )
    return max(pdb_name), new_dir


if __name__ == '__main__':

    #########################################################################################################################################################
    ### In working directory such as model9-ts-001 
    ### Run: write_input.py -noh nohpdb -adh haddpdb -intmp relaxh_temp -c 2 (final_pdb is saved as "template.pdb", input is saved as "1.inp", default -m 1) 
    ### Run: write_input.py -step 1 -intmp modred_temp -m 1 -c 2 (will take the last pdb from 1.out and write input as "1.inp")                      
    ### Run: write_input.py -step 2 -intmp modred_temp -m 1 -c 2 -new step1pdbs/33.pdb (will take the selected pdb and write "1.inp")
    #########################################################################################################################################################

    parser = argparse.ArgumentParser(description='Prepare template PDB files, write input files, save output PDB files in working directory',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-type', dest='step', default='pdb', 
            help='hopt: read noh, addh pdbs, write_final_pdb and read input_template write_first_inp, \n' + 
            'gauout: read outputwrite_modred_inp, input_template, write_second_inp, \n' +
            'pdb: read new pdb file, input_template, write_new_inp \n' +
            'replacecoords: read pdb1 and pdb2 for replacing the fragment in pdb1 with pdb2 coordinates and write new input\n' +
            'fsapt: F-SAPT0 calculation (psi4 only, ignores format specification)')
    # general options
    parser.add_argument('-m', dest='multiplicity', default=1, type=int, help='multiplicity')
    parser.add_argument('-c', dest='ligand_charge', default=0, type=int, help='charge_of_ligand')
    parser.add_argument('-format',dest='fmat',default=None,help="input_file_format eg.'gaussian','qchem','gau-xtb'")
    parser.add_argument('-intmp', dest='input_tmp', default=None, help='template_for_write_input')
    parser.add_argument('-inpn', dest='inp_name', default='1.inp', help='input_name')
    parser.add_argument('-basisinfo',dest='basisinfo',default=None,help="'intmp' if in template file, use dictionary otherwise")
    parser.add_argument('-wdir', dest='output_dir', default=os.path.abspath('./'), help='working dir')
    # type = pdb
    parser.add_argument('-pdb', dest='new_pdb', default=None, help='new_pdb_file')
    # type = hopt
    parser.add_argument('-noh', dest='no_h_pdb', default=None, help='trimmed_pdb_file')
    parser.add_argument('-adh', dest='h_add_pdb', default=None, help='hadded_pdb_file')
    parser.add_argument('-tmp', dest='tmp_pdb', default=None, help='template_pdb_file')
    # type = gauout
    parser.add_argument('-outf', dest='gau_out', default='1.out', help='output_name')
    parser.add_argument('-ckp', dest='check_point', default=None, help='check_file_frame')
    # type = replacecoords
    parser.add_argument('-pdb1', dest='pdb1', default=None, help='minima_pdb_file')
    parser.add_argument('-pdb2', dest='pdb2', default=None, help='ts_pdb_file')
    parser.add_argument('-parts', dest='parts', default=None, help='ts_frag_indo')
    # type = fsapt
    parser.add_argument('-seed', dest='seed', default=None, help='seed')
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
    gauout   = args.gau_out
    inp_name = args.inp_name
    multi    = args.multiplicity
    charge   = args.ligand_charge
    wdir     = args.output_dir
    ifmat    = args.fmat
    basisinfo = args.basisinfo

    hopt = 0

    if step == 'hopt':
        #pic_atom, tot_charge = pdb_after_addh(nohpdb,adhpdb)
        #res_count = adhpdb.split('_')[1]
        #write_pdb('%s'%tmp_pdb,pic_atom,res_count)
        newpdb = args.new_pdb
        pic_atom, res_info, tot_charge = read_pdb(newpdb)
        res_count = args.new_pdb
        hopt = 1
    elif step == 'gauout':
        i_name = []
        for gau_input in glob('%s/*inp'%wdir):
            m = re.search(r'-(\d+)-inp', gau_input)
            if m:
                i_name.append( int(m.group(1)) )
        if len(i_name) == 0:
            i_step = 1
            if os.path.isfile('%s/1.out'%wdir):
                system_run( 'cp 1.inp step-%s-inp'%(i_step) )
                system_run( 'cp 1.out step-%s-out'%(i_step) )
                system_run( 'cp 1.chk step-%s-chk'%(i_step) )
        else:
            i_step = max(i_name)
            if filecmp.cmp('1.out','%s/step-%s-out'%(wdir,i_step)) is False:
                i_step += 1
                system_run( 'cp 1.inp step-%s-inp'%(i_step) )
                system_run( 'cp 1.out step-%s-out'%(i_step) )
                system_run( 'cp 1.chk step-%s-chk'%(i_step) )
            elif filecmp.cmp('1.inp','%s/step-%s-inp'%(wdir,i_step)) is True and filecmp.cmp('1.out','%s/step-%s-out'%(wdir,i_step)) is True:
                i_step = i_step
                print("1.inp is the same as step-%s-inp, and 1.out is the same as step-%s-out, please check!")
#                sys.exit()
            else:
                print("check if the files are propagated correctly!")
#                sys.exit()
        res_count = 'step-%s'%i_step
        if step == 'gauout':
            pdb_file, new_dir = gen_pdbfiles(wdir,i_step,tmp_pdb)
            pic_atom, res_info, tot_charge = read_pdb('%s/%s.pdb'%(new_dir,pdb_file))
    elif step == 'pdb':
        newpdb = args.new_pdb
        pic_atom, res_info, tot_charge = read_pdb(newpdb)
        res_count = args.new_pdb
    elif step == 'replacecoords':
        pdb1 = args.pdb1
        pdb2 = args.pdb2
        parts = args.parts
        if parts != None:
           with open(parts) as f:
               plines = f.readlines()
           parts = []
           for line in plines:
               parts.append(line.split())
        pic_atom, tot_charge = pdb_replace(pdb1,pdb2,parts)    
        res_count = args.pdb1
        # write the modified pdb so it can be checked
        pdb1name = pdb1.split('/')[-1]
        pdb1name = pdb1name.split('.')[0]
        write_pdb('%s/input_parts_replaced.pdb'%wdir,pic_atom)

    if int_tmp == None or int_tmp == '':
        tmpltdir = Path.home() / 'git' / 'RINRUS' / 'template_files'
        int_tmp = tmpltdir / f'{ifmat}_input_template.txt'
        print("Using default opt+freq input template: "+str(int_tmp))

    if ifmat == "gaussian":
        write_gau_input('%s/%s'%(wdir,inp_name),int_tmp,charge,multi,pic_atom,tot_charge,res_count,basisinfo,hopt)
    elif ifmat == "qchem":
        write_qchem_input('%s/%s'%(wdir,inp_name),int_tmp,charge,multi,pic_atom,tot_charge,res_count)
    elif ifmat == "gau-xtb":
        write_xtb_input('%s/%s'%(wdir,inp_name),int_tmp,charge,multi,pic_atom,tot_charge,res_count)
    elif ifmat == "orca":
        write_orca_input('%s/%s'%(wdir,inp_name),int_tmp,charge,multi,pic_atom,tot_charge,res_count,hopt)
    elif ifmat == "psi4-fsapt":
        if inp_name == "1.inp":
            inp_name = "input.dat"
        seed = args.seed
        write_psi4_fsapt_input('%s/%s'%(wdir,inp_name),int_tmp,charge,multi,pic_atom,tot_charge,res_count,seed)
    else:
        print("ERROR: 'ifmat' not set. Please provide an input format for your calculations!")

