#!/usr/bin/env python
"""
This is a program created by Dr. Qianyi Cheng and the DeYonker research group
at The University of Memphis.
Date created 08.12.2022
"""

import sys, os, argparse

parser = argparse.ArgumentParser(description='Generates file documenting interactions of functional groups of arpeggio contacts file')
parser.add_argument('-c', help='name of arpeggio contacts file to analyze')
parser.add_argument('-p', help='name of PDB file to analyze')
parser.add_argument('-s', nargs="+", default=['A/128'], help='residue(s) of seed in format chain/residue e.g. A/128')
parser.add_argument('-save', default='FG-arpeggio.dat', help='name of savefile')
args = parser.parse_args()

#ID waters from PDB file
waters = set()
with open(args.p, 'r') as pdbfile:
    for line in pdbfile:
        if line[0:3] == "END": continue
        if line[17:20] == "WAT" or line[17:20] == "HOH":
            waters.add(line[21]+"/"+line[22:26].strip())

#Grab interaction info
interactions = {}

def nameFG(atom, waterset):
    chain = atom.split("/")[0]
    res = int(atom.split("/")[1])
    a = atom.split("/")[2]

    if atom.rsplit("/",1)[0] in waterset:
        return chain+":"+str(res)+":WAT"
    elif a in ['C', 'O']:
        return chain+":"+str(res)+":MC"
    elif a in ['N', 'H']:
        return chain+":"+str(res-1)+":MC"
    else:
        return chain+":"+str(res)+":SC"
        

with open(args.c, 'r') as readfile:
    for line in readfile:
        res1 = line.split()[0].rsplit("/",1)[0]
        res2 = line.split()[1].rsplit("/",1)[0]

        if res1 in args.s and res2 in args.s: continue
        if res1 not in args.s and res2 not in args.s: continue

        if res1 in args.s:
            FG = nameFG(line.split()[1], waters)
        else:
            FG = nameFG(line.split()[0], waters)

        values = [int(x) for x in line.split()[2:17]]

        if FG not in interactions.keys():
            interactions[FG] = values
        else:
            interactions[FG] = [x+y for x,y in zip(values, interactions[FG])]

names = list(interactions.keys())
names.sort()
with open(args.save, 'w') as savefile:
    for item in names:
        savefile.write(item+" "+" ".join(map(str,interactions[item]))+"\n")

        
