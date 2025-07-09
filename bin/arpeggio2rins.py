#!/usr/bin/env python3
"""
This is a program created by Dr. Dominique Wappett 
(based on the old program by Dr. Qianyi Cheng)
and the DeYonker research group at The University of Memphis.
Date created 2025-06-18
"""


import os, sys, re
from copy import *
import argparse
import numpy as np
import pandas as pd
from res_atoms import *
from read_write_pdb import *


def contact_dict_add(atomdict,contatom,conts):
### add identified contacts to dictionary ###
    key = (contatom[0],int(contatom[1]))
    if key not in atomdict.keys():
        atomdict[key] = {}
    if contatom[2] not in atomdict[key].keys():
        atomdict[key][contatom[2]] = conts
    else:
        #atomdict[key][contatom[2]] = atomdict[key][contatom[2]] + conts
        atomdict[key][contatom[2]] = [sum(x) for x in zip(atomdict[key][contatom[2]], conts)]

def arpeggio_atom(arpeggiofile,seedlist,noprox):
### get contact counts between any atom and seed ###
    atomdict = {}
    seeddict = {s: {'atoms': [], 'counts': [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]} for s in seedlist}
    lines = open(arpeggiofile,'r').readlines()
    for line in lines:
        line = line.split()
        at1 = line[0].split('/')
        at2 = line[1].split('/')
        # only look for protein-seed interactions, ignore protein-protein or seed-seed
        if at1[0]+':'+at1[1] in seedlist and at2[0]+':'+at2[1] not in seedlist:
            seedatom = at1
            contatom = at2
        elif at1[0]+':'+at1[1] not in seedlist and at2[0]+':'+at2[1] in seedlist:
            seedatom = at2
            contatom = at1
        else: 
            continue
        conts = [int(x) for x in line[2:-1]]
        # avoid adding atoms with only proximal interactions to dictionary if not wanted
        if noprox:
            conts[4] = int(0)
        if sum(conts) > 0:
            contact_dict_add(atomdict,contatom,conts)
            seedkey = str(seedatom[0])+':'+str(seedatom[1])
            seeddict[seedkey]['counts'] = [sum(x) for x in zip(seeddict[seedkey]['counts'], conts)]
            if seedatom[2] not in seeddict[seedkey]['atoms']:
                seeddict[seedkey]['atoms'].append(seedatom[2])
    return atomdict, seeddict

def arpeggio_fg(atomdict,pdb_res_name,seeddict):
### collect atoms and contact counts by functional groups ###
    FGcontacts = {}
    FGatomlists = {}
    for key in atomdict.keys():
        # if known amino acid res type, split into side and main chains
        if pdb_res_name[key] in res_atoms_all.keys():
            for atom in atomdict[key].keys():
                # assign FG
                if atom in ['C','O']:
                    fg = key[0]+':'+str(key[1])+':MC'
                elif atom in ['N','H']:
                    fg = key[0]+':'+str(key[1]-1)+':MC'
                else:
                    fg = key[0]+':'+str(key[1])+':'+pdb_res_name[key]+'_SC'
                # make sure FG is in dict
                if fg not in FGcontacts.keys():
                    FGcontacts[fg] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                # add atom counts to overall FG counts
                FGcontacts[fg] = [sum(x) for x in zip(FGcontacts[fg], atomdict[key][atom])]
                # add atom name to atom list
                if fg not in FGatomlists.keys():
                    FGatomlists[fg] = []
                FGatomlists[fg].append(atom)
       # otherwise whole res is FG (ligands, cofactors, water, etc)
        else:
            for atom in atomdict[key].keys():
                fg = key[0]+':'+str(key[1])+':'+pdb_res_name[key]
                # make sure FG is in dict
                if fg not in FGcontacts.keys():
                    FGcontacts[fg] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                # add atom counts to overall FG counts
                FGcontacts[fg] = [sum(x) for x in zip(FGcontacts[fg], atomdict[key][atom])]
                # add atom name to atom list
                if fg not in FGatomlists.keys():
                    FGatomlists[fg] = []
                FGatomlists[fg].append(atom)
    allcont = FGcontacts
    allats = FGatomlists
    for s in seeddict.keys():
        allcont[s] = seeddict[s]['counts']
        allats[s] = seeddict[s]['atoms']

    return allcont, allats


def rinrus_arpeggio_outputs(allcont,allats,seeddict,noprox,arpfile,sel_res):
### make table of contacts, sort by total and types, save to file, write res_atoms files ###
    contypes = ['clash','covalent','vdw_clash','vdw','proximal','hbond','weak_hbond','halogen_bond','ionic','metal_complex','aromatic','hydrophobic','carbonyl','polar','weak_polar']
    #shortnames = ['clash','cov','vdw_cl','vdw','prox','hbond','wk_hb','halog','ionic','metal','arom','hydrphb','crbnyl','polar','wk_pol']
    df = pd.DataFrame.from_dict(allcont,orient='index',columns=contypes)
    if noprox:
        df = df.drop(columns='proximal')
    tc = df.apply(lambda row: len([c for c in row if c >0]), axis='columns')
    df.insert(loc=0,column='total',value=df[list(df.columns)].sum(axis=1))
    df.insert(loc=1,column='types',value=tc.values)
    df = df.query('total > 0')
    df = df.sort_values(by=['total','types'], key=abs, ascending=False)
    # add seed label to seed frags just for printed contact table
    df.index = [i+' (seed)' if i in seeddict.keys() else i for i in df.index] 
    df.to_string(buf=f'FG_arpeggio_counts.dat')
    df.index = [i.replace(' (seed)','') for i in df.index]

    ### make counts res atoms file ###
    f = open('res_atoms.dat','w')
    if noprox:
        proxinf = 'proximal contacts ignored'
    else:
        proxinf = 'proximal contacts included'
    f.write(f'# selected/ranked functional groups by arpeggio contact counts: arpeggio file {arpfile}; seed {sel_res}; {proxinf}\n')
    for g in df.index.values:
        # if MC, make sure atoms correspond to correct res id
        if g.split(':')[-1] == 'MC': 
            idx = int(g.split(":")[1])
            if set(allats[g]).issubset({'C','O'}):
                f.write(f'{g.split(":")[0]:<4} {idx:<8} {df.loc[g,"total"]:<8}')
                for a in allats[g]:
                    f.write(f' {a:<4}')
                f.write('\n')
            elif set(allats[g]).issubset({'N','H'}):
                f.write(f'{g.split(":")[0]:<4} {idx+1:<8} {df.loc[g,"total"]:<8}')
                for a in allats[g]:
                    f.write(f' {a:<4}')
                f.write('\n')
            else:
                f.write(f'{g.split(":")[0]:<4} {idx:<8} {df.loc[g,"total"]:<8}')
                for a in [a for a in allats[g] if a in ['C','O']]:
                    f.write(f' {a:<4}')
                f.write('\n')
                f.write(f'#{g.split(":")[0]:<3} {idx+1:<8} {df.loc[g,"total"]:<8}')
                for a in [a for a in allats[g] if a in ['N','H']]:
                    f.write(f' {a:<4}')
                f.write('\n')
        else:
            f.write(f'{g.split(":")[0]:<4} {g.split(":")[1]:<8} {df.loc[g,"total"]:<8}')
            for a in allats[g]:
                f.write(f' {a:<4}')
            f.write('\n')
    f.close()

    ### make types res atoms file ###
    df = df.sort_values(by=['types','total'], key=abs, ascending=False)
    f = open('res_atoms_types.dat','w')
    f.write(f'# selected/ranked functional groups by arpeggio contact types: arpeggio file {arpfile}; seed {sel_res}; {proxinf}\n')
    for g in df.index.values:
        # if MC, make sure atoms correspond to correct res id
        if g.split(':')[-1] == 'MC':
            idx = int(g.split(":")[1])
            if set(allats[g]).issubset({'C','O'}):
                f.write(f'{g.split(":")[0]:<4} {idx:<8} {df.loc[g,"types"]:<8}')
                for a in allats[g]:
                    f.write(f' {a:<4}')
                f.write('\n')
            elif set(allats[g]).issubset({'N','H'}):
                f.write(f'{g.split(":")[0]:<4} {idx+1:<8} {df.loc[g,"types"]:<8}')
                for a in allats[g]:
                    f.write(f' {a:<4}')
                f.write('\n')
            else:
                f.write(f'{g.split(":")[0]:<4} {idx:<8} {df.loc[g,"types"]:<8}')
                for a in [a for a in allats[g] if a in ['C','O']]:
                    f.write(f' {a:<4}')
                f.write('\n')
                f.write(f'#{g.split(":")[0]:<3} {idx+1:<8} {df.loc[g,"types"]:<8}')
                for a in [a for a in allats[g] if a in ['N','H']]:
                    f.write(f' {a:<4}')
                f.write('\n')
        else:
            f.write(f'{g.split(":")[0]:<4} {g.split(":")[1]:<8} {df.loc[g,"types"]:<8}')
            for a in allats[g]:
                f.write(f' {a:<4}')
            f.write('\n')
    f.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate interaction information from arpeggio file.\
            Usage: arpeggio2rins.py -pdb pdb -f probe_file -s seed [-prox]')
    parser.add_argument('-pdb',dest='pdb',default=None,help='pdb file')
    parser.add_argument('-f',dest='arpfile',default=None,help='arpeggio contacts file')
    parser.add_argument('-s','-seed',dest='seed',default=None,help='seed for selecting RIN, in the format of "A:300,A:301,A:302"')
    parser.add_argument('-prox','-proximal',dest='prox',action='store_true',help='include proximal interactions')

    ### parse arguments ###
    args = parser.parse_args()
    arpfile = args.arpfile
    pdbf = args.pdb
    sel_res = args.seed
    seedlist = sel_res.split(',')
    noprox = not args.prox

    ### process pdb to get res names ###
    pdb, res_info, tot_charge = read_pdb(pdbf)
    pdb_res_name = {}
    for line in pdb:
        key = (line[5], int(line[6]))
        pdb_res_name[key] = line[4].strip()

    ### get and sort contacts ###
    atomdict, seeddict = arpeggio_atom(arpfile,seedlist,noprox)
    allcont, allats = arpeggio_fg(atomdict,pdb_res_name,seeddict)

    ### make contact count table and res_atoms files ###
    rinrus_arpeggio_outputs(allcont,allats,seeddict,noprox,arpfile,sel_res)
