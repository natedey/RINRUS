#!/usr/bin/env python3
"""
This is a program written by Qianyi Cheng 
at university of memphis.
"""
from numpy import *
import os, sys, re, argparse
from res_atoms import *

def get_fg_sapt(fgfile):
    with open(fgfile,'r') as f:
        lines = f.readlines()

    res_list = {}
    for i in range(1,len(lines)):
        v = lines[i].split()
        key = v[0]
        res_list[key] = []
        for value in v[1:]:
            res_list[key].append(float(value))
    return res_list

#res_list = get_fg_sapt('/home/qcheng1/projects/gnmt/sapt-rank/fsapt/gly/FG-SAPT.dat')

def gen_res_atom(res_list,old_res_atom):
    new_res_list = sorted(res_list.items(),key=lambda x:abs(x[1][-1]),reverse=True)
    with open(old_res_atom,'r') as f:
        lines = f.readlines()
    old_atoms = {}
    for i in range(len(lines)):
        v = lines[i].split()
        key = (v[0],v[1])
        atoms = []
        for j in range(3,len(v)):
            atoms.append(v[j])
        old_atoms[key] = atoms

    res_part = []
    res_value = []
    for i in new_res_list:
        res_part.append(i[0])
        res_value.append(i[-1])

    key_list = []
    res_atom = []
    res_score = []
    for lx in range(len(res_part)):
        key = res_part[lx]
        chain, res_id, part = key.split(':')
        nkey = (chain,res_id)
        if part == 'SC':
            if nkey not in key_list and nkey in old_atoms.keys():
                atom_list = []
                for atom in old_atoms[nkey]:
                    if atom not in ['CA','N','C','O','HA','H','HA1','HA2','HA3']:
                        atom_list.append(atom)
                if len(atom_list) > 0:
                    res_atom.append(atom_list)
                    key_list.append(nkey)
                    res_score.append(res_value[lx])
            elif nkey in key_list and nkey in old_atoms.keys():
                atom_list = res_atom[key_list.index(nkey)]
                cb_idx = []
                for atom in old_atoms[nkey]:
                    if atom not in ['CA','N','C','O','HA','H','HA1','HA2','HA3']:
                        cb_idx.append(old_atoms[nkey].index(atom))
                if len(cb_idx) > 0:
                    res_atom.append(old_atoms[nkey])
                    key_list.append(nkey)
                    res_score.append(res_value[lx])
        if part == 'MC':
            if nkey not in key_list and nkey in old_atoms.keys():
                atom_list = []
                for atom in old_atoms[nkey]:
                    if atom in ['CA','N','C','O','HA','H','HA1','HA2','HA3']:
                        atom_list.append(atom)
                if len(atom_list) > 0:
                    key_list.append(nkey)
                    res_atom.append(atom_list)
                    res_score.append(res_value[lx])
            elif nkey in key_list and nkey in old_atoms.keys():
                atom_list = res_atom[key_list.index(nkey)]
                cb_idx = []
                for atom in old_atoms[nkey]:
                    if atom in ['CA','N','C','O','HA','H','HA1','HA2','HA3']:
                        cb_idx.append(old_atoms[nkey].index(atom))
                if len(cb_idx) > 0:
                    res_atom.append(old_atoms[nkey])
                    key_list.append(nkey)
                    res_score.append(res_value[lx])

    return res_part, key_list, res_atom, res_score

#res_part, key_list, res_atom = gen_res_atom(res_list,'/home/qcheng1/projects/gnmt/arpeggio/sam-gly/contact_counts.dat')

def write_res_atom(res_list,res_part,key_list,res_atom,seed,res_score):
    f = open('res_atoms_fsapt.dat','w')
    seeds_in = seed.split(',')
    seeds = []
    idx_list = []
    for sel in seeds_in:
        v = sel.split(':')
        seeds.append((v[0],v[1]))
    for sd in seeds:
        if sd not in key_list:
            f.write('%8s%8s\n'%(sd[0],sd[1]))
        else:
            idx = key_list.index(sd)
            idx_list.append(idx)
            f.write('%8s%8s'%(key_list[idx][0],key_list[idx][1]))
            f.write('%12.3f'%res_score[idx][-1])
            for atom in res_atom[idx]:
                f.write('%6s'%atom)
            f.write('\n')
    for i in range(len(key_list)):
        if i not in idx_list:
            f.write('%8s%8s'%(key_list[i][0],key_list[i][1]))
            f.write('%12.3f'%res_score[i][-1])
            for atom in res_atom[i]:
                f.write('%6s'%atom)
            f.write('\n')
    f.close()

if __name__ == '__main__':
    ################################################################################
    ### python3 sapt2rin.py -f FG-sapt.dat -c contact_counts.dat
    ################################################################################

    parser = argparse.ArgumentParser(description='Get contact information from SAPT energy')
    parser.add_argument('-p', dest='fg',help='FG-sapt energy dat')
    parser.add_argument('-c', dest='cf',help='contact data from other programs')
    parser.add_argument('-s', dest='seed',help='seed residues, in the format "A:293,A;294"')

    args = parser.parse_args()
    fg_file = args.fg
    ct_file = args.cf
    seed = args.seed
    
    res_list = get_fg_sapt(fg_file)
    res_part, key_list, res_atom, res_score = gen_res_atom(res_list,ct_file)


    write_res_atom(res_list,res_part,key_list,res_atom,seed,res_score)
