#!/usr/bin/env python
"""
This is a program created by Dr. Dominique Wappett
(based on the old program by Dr. Qianyi Cheng)
and the DeYonker research group at The University of Memphis.
Date created 2025-06-18
"""


import os, sys, re
from copy import *
from numpy import *
import argparse
import pandas as pd
from res_atoms import *


def atom_split(c):
### separate out atom labels from probe file ###
    cha    = c[:2].strip()          # chain
    res_id = c[2:6].strip()         # res_id
    res_nm = c[6:10].strip()        # res_name
    name   = c[10:-1].strip()       # atom_name
    atom = [cha, res_id, res_nm, name]
    return atom

def contact_dict_add(atomdict,contatom,ctype):
### organize and count identified contacts ###
    key = (contatom[0],int(contatom[1]))
    if key not in atomdict.keys():
        atomdict[key] = {'res': contatom[2]}
    if contatom[3] not in atomdict[key].keys():
        atomdict[key][contatom[3]]={'wc': 0, 'cc': 0, 'so': 0, 'bo': 0, 'hb': 0}
    atomdict[key][contatom[3]][ctype] += 1

def probe_atom(probefile,seedlist):
### get contact counts between any atom and seed ###
    atomdict = {}
    seedcontact = {s: {'wc': 0, 'cc': 0, 'so': 0, 'bo': 0, 'hb': 0} for s in seedlist}
    seedatoms = {s: [] for s in seedlist}
    lines = open(probefile,'r').readlines()
    for line in lines:
        line = line.split(':')
        ctype = line[2]
        at1 = atom_split(line[3])
        at2 = atom_split(line[4])
        # ignore if neither atom in line is in seed or both in seed
        if at1[0]+':'+at1[1] not in seedlist and at2[0]+':'+at2[1] not in seedlist:
            continue
        if at1[0]+':'+at1[1] in seedlist and at2[0]+':'+at2[1] in seedlist:
            continue
        # pick which atom to add to dict
        if at1[0]+':'+at1[1] in seedlist:
            contatom = at2
            seedatom = at1
        else:
            contatom = at1
            seedatom = at2
        contact_dict_add(atomdict,contatom,ctype)
        seedkey = seedatom[0]+':'+seedatom[1]
        seedcontact[seedkey][ctype] += 1
        seedatoms[seedkey].append(seedatom[3])
    return atomdict, seedcontact, seedatoms

def probe_fg(atomdict,seedcontact,seedatoms):
### collect atoms and contact counts by functional groups ###
    FGcontacts = {}
    FGatomlists = {}
    for key in atomdict.keys():
        # if known amino acid res type, split into side and main chains
        if atomdict[key]['res'] in res_atoms_all.keys():
            for atom in [a for a in atomdict[key].keys() if a != 'res']:
                # assign FG
                if atom in ['C','O']:
                    fg = key[0]+':'+str(key[1])+':MC'
                elif atom in ['N','H']:
                    fg = key[0]+':'+str(key[1]-1)+':MC'
                else:
                    fg = key[0]+':'+str(key[1])+':'+atomdict[key]['res']+'_SC'
                # make sure FG is in dict
                if fg not in FGcontacts.keys():
                    FGcontacts[fg] = {'wc': 0, 'cc': 0, 'so': 0, 'bo': 0, 'hb': 0}
                # add atom counts to overall FG counts
                for ct in ['wc','cc','so','bo','hb']:
                    FGcontacts[fg][ct] += atomdict[key][atom][ct]
                # add atom name to atom list
                if fg not in FGatomlists.keys():
                    FGatomlists[fg] = []
                FGatomlists[fg].append(atom)
        # otherwise whole res is FG (ligands, cofactors, water, etc)
        else:
            for atom in [a for a in atomdict[key].keys() if a != 'res']:
                fg = key[0]+':'+str(key[1])+':'+atomdict[key]['res']
                # make sure FG is in dict
                if fg not in FGcontacts.keys():
                    FGcontacts[fg] = {'wc': 0, 'cc': 0, 'so': 0, 'bo': 0, 'hb': 0}
                # add atom counts to overall FG counts
                for ct in ['wc','cc','so','bo','hb']:
                    FGcontacts[fg][ct] += atomdict[key][atom][ct]
                # add atom name to atom list
                if fg not in FGatomlists.keys():
                    FGatomlists[fg] = []
                FGatomlists[fg].append(atom)
    allcont = FGcontacts
    allats = FGatomlists
    for s in seedcontact.keys():
        allcont[s] = seedcontact[s]
        allats[s] = unique(seedatoms[s]).tolist()

    #return FGcontacts, FGatomlists
    return allcont, allats

def rinrus_probe_outputs(allcont,allats,seedcontact,probefile,sel_res):
    ### make contact table and save ###
    df = pd.DataFrame.from_dict(allcont,orient='index')
    df.insert(loc=0,column='total',value=df[list(df.columns)].sum(axis=1))
    df = df.sort_values(by='total', key=abs, ascending=False)
    df.index = [i+' (seed)' if i in seedcontact.keys() else i for i in df.index]
    df.to_string(buf=f'FG_probe_counts.dat')
    df.index = [i.replace(' (seed)','') for i in df.index]

    ### make res atoms file ###
    f = open('res_atoms.dat','w')
    f.write(f'# selected/ranked functional groups by probe: probe file {probefile}; seed {sel_res}\n')
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
                f.write(f'{g.split(":")[0]:<4} {idx+1:<8} {df.loc[g,"total"]:<8}')
                for a in [a for a in allats[g] if a in ['N','H']]:
                    f.write(f' {a:<4}')
                f.write('\n')
        else:
            f.write(f'{g.split(":")[0]:<4} {g.split(":")[1]:<8} {df.loc[g,"total"]:<8}')
            for a in allats[g]:
                f.write(f' {a:<4}')
            f.write('\n')
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate interaction information from probe file.\
            Usage: probe2rins.py -f probe_file -s seed')
    parser.add_argument('-f',dest='probefile',default=None,help='probe_file')
    parser.add_argument('-s','-seed',dest='seed',default=None,help='seed for select RIN, in the format of "A:300,A:301,A:302"')

    ### parse arguments ###
    args = parser.parse_args()
    probefile = args.probefile
    sel_res = args.seed
    seedlist = sel_res.split(',')

    ### get and sort contacts ###
    atomdict, seedcontact, seedatoms = probe_atom(probefile,seedlist)
    allcont, allats = probe_fg(atomdict, seedcontact, seedatoms)

    ### make rinrus outputs ###
    rinrus_probe_outputs(allcont,allats,seedcontact,probefile,sel_res)
