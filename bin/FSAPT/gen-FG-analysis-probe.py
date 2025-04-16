#!/usr/bin/env python
"""
This is a program created by Dr. Qianyi Cheng and the DeYonker research group
at The University of Memphis.
Date created 08.12.2022
"""

import sys, os, argparse
import fnmatch
from collections import defaultdict

parser = argparse.ArgumentParser(description='Generates file documenting contacts of functional groups of PDB file')
parser.add_argument('-p', help='name of probe file to analyze')
parser.add_argument('-s', nargs="+", default=['A:128'], help='residue(s) of seed in format chain/residue e.g. A:128')
parser.add_argument('-save', default='FG-probe.dat', help='name of savefile')
args = parser.parse_args()

contacts = defaultdict(int)
MCatoms = ['C', 'N', 'O', 'H']

#Grab seed info
with open(args.p, 'r') as readfile:
    for line in readfile:
        if line[10] + ":" + line[11:15].strip() in args.s and line[27] + ":" + line[28:32].strip() in args.s: continue

        elif line[10] + ":" + line[11:15].strip() in args.s:
            atom = line[37:41].strip()
            atomid = int(line[28:32])
            atomname = line[33:36]
            if atomname == "HOH" or atomname == "WAT":
                FG = line[27]+":"+str(atomid)+":WAT" 
            elif atom in MCatoms:
                if atom == "C" or atom == "O":
                    FG = line[27]+":"+str(atomid)+":MC"
                else:
                    FG = line[27]+":"+str(atomid-1)+":MC"
            else:
                FG = line[27]+":"+str(atomid)+":SC"
            contacts[FG] +=1

        elif line[27] + ":" + line[28:32].strip() in args.s:
            atom = line[20:24].strip()
            atomid = int(line[11:15])
            atomname = line[16:19]
            if atomname == "HOH" or atomname == "WAT":
                FG = line[10]+":"+str(atomid)+":WAT"
            elif atom in MCatoms:
                if atom == "C" or atom == "O":
                    FG = line[10]+":"+str(atomid)+":MC"
                else:
                    FG = line[10]+":"+str(atomid-1)+":MC"
            else:
                FG = line[10]+":"+str(atomid)+":SC"
            contacts[FG] +=1

names = list(contacts.keys())
names.sort()

with open(args.save, 'w') as savefile:
    for item in names:
        savefile.write(item+" "+str(contacts[item])+"\n")


