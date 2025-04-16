#!/usr/bin/env python
"""
This is a program created by Dr. Dominique Wappett and the DeYonker research group
at The University of Memphis.
Date created 04.09.2025
"""

import sys, os, argparse
import pickle
from read_write_pdb import *

def addtodict(dct,key,num):
    try:
        dct[key].append(num)
    except:
        dct[key] = [num]

elmts = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 
        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 
        'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
        'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create fA.dat and fB.dat files for FSAPT analysis of all FGs at once')
    parser.add_argument('-pdb', dest='pdbf', default=None, help='model pdb')
    parser.add_argument('-s','-seed', dest='seed', default=None, help='seed, examples: A:300,A:301,A:302')
    parser.add_argument('-inp', dest='inpf', default='../input.dat', help='input file for fsapt calculation')
    
    args = parser.parse_args()
    pdbf = args.pdbf
    seed = args.seed.split(',')

    ### read input, get atom numbers and fragment A/B/C from coord lines
    inp = open(args.inpf, 'r').read().split('--')
    num = 1
    atoms = {}
    frags = ['A','B','C']
    fnum = 0
    for block in inp:
        hadats = 0
        for line in block.split('\n'):
            line = line.split()
            if line and line[0] in elmts:
                key = (float(line[1]),float(line[2]),float(line[3]))
                atoms[key] = {'num': num, 'fragment': frags[fnum]}
                num +=1
                hadats = 1
        if hadats == 1:
            fnum += 1
    
    ### get info for each atom from pdb by matching coords
    pdb, res_info, tot_charge = read_pdb(pdbf)
    for line in pdb:
        key = (line[8],line[9],line[10])
        info = {'chain': line[5], 'id': line[6], 'resn': line[4], 'name': line[2].strip(), 'pdbnum': line[1]}
        atoms[key].update(info)
        if f'{line[5]}:{line[6]}' in seed:
            atoms[key]['seed'] = True
        else:
            atoms[key]['seed'] = False

    # atom dict contents for each atom: 'fragment'/'num'/'chain'/'id'/'resn'/'name'/'pdbnum'/'seed'

    ### get fA/fB number lists and seed/FG atom lists for res_atoms
    fAnum = {}
    fBnum = {}
    fBatom = {}
    seedinfo = {}
    fAinfo = {}
    FGinfo = {}
    for a in atoms.keys():
        resid = f"{atoms[a]['chain']}:{atoms[a]['id']}"
        if atoms[a]['fragment']=='A':
            addtodict(fAnum,resid,atoms[a]['num'])
            if atoms[a]["seed"]:
                addtodict(seedinfo,resid,atoms[a]['name'])
            else:
                addtodict(fAinfo,resid,atoms[a]['name'])
        elif atoms[a]['fragment']=='B':
            # get FG name
            if atoms[a]['resn'] in ['WAT','HOH']:
                key = f'{atoms[a]["chain"]}:{atoms[a]["id"]}:{atoms[a]["resn"]}'
            elif atoms[a]["name"] in ['C','O']:
                key = f'{atoms[a]["chain"]}:{atoms[a]["id"]}:MC'
            elif atoms[a]["name"] in ['N','H']:
                key = f'{atoms[a]["chain"]}:{atoms[a]["id"]-1}:MC'
            else:
                key = f'{atoms[a]["chain"]}:{atoms[a]["id"]}:{atoms[a]["resn"]}'
            # add num to fBnum and atom name to FGinfo
            addtodict(fBnum,key,atoms[a]['num'])
            addtodict(fBatom,key,atoms[a]["name"])
            # add atom name to FG info and also seed info if 
            if atoms[a]["seed"]:
                addtodict(seedinfo,resid,atoms[a]['name'])
            else:
                addtodict(FGinfo,key,atoms[a]["name"])
        elif atoms[a]['fragment']=='C':
            if atoms[a]["seed"]:
                addtodict(seedinfo,resid,atoms[a]['name'])

    ### write fA.dat
    fa1 = open('fA.dat','w')
    for key in fAnum.keys():
        fa1.write(f'{key} {" ".join([str(i) for i in fAnum[key]])}\n')
    fa1.close()

    ### write fB.dat
    fb1 = open('fB.dat','w')
    # filter out CA/HA atoms only present for capping/connecting from "proper" FGs
    caps = []
    for key in fBnum.keys():
        if key.split(':')[-1] in ['WAT','HOH','MC','GLY']:
            fb1.write(f'{key} {" ".join([str(i) for i in fBnum[key]])}\n')
        elif not set(fBatom[key]).issubset({'CA','HA','HA2','HA3','H01','H02','H03'}):
            fb1.write(f'{key} {" ".join([str(i) for i in fBnum[key]])}\n')
        else:
            caps = caps + fBnum[key]
            del FGinfo[key]
    if caps:
        fb1.write(f'cap_ats {" ".join([str(i) for i in caps])}\n')
    fb1.close()

    ### save atom lists for res_atoms later
    atom_info = {'seed': seedinfo, 'other_fA': fAinfo, 'FG': FGinfo}
    with open('fdict.pkl', 'wb') as fp:
        pickle.dump(atom_info, fp)

    ### human readable file for checking atom numbers
    with open('frag_atom_check.dat','w') as f:
        f.write('frag Ninp ch   ID     res  atom Npdb\n')
        lastfrag = ''
        for a in atoms.keys():
            if lastfrag != atoms[a]["fragment"]:
                f.write('------------------------------------\n')
            f.write(f'{atoms[a]["fragment"]:<4} {atoms[a]["num"]:<4} {atoms[a]["chain"]:<4} {atoms[a]["id"]:<6} {atoms[a]["resn"]:<4} {atoms[a]["name"]:<4} {atoms[a]["pdbnum"]:<4}\n')
            lastfrag = atoms[a]["fragment"]
