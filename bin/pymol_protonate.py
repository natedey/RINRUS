#!/usr/bin/env python3
"""
This is a program created by Dr. Qianyi Cheng and the DeYonker research group
at The University of Memphis.
Date created 06.12.2019
"""
import os, sys
import subprocess


def system_run(cmd):
    print(cmd)
    exit = os.system(cmd)
    if exit != 0:
        print("failed to run:")
        print("pymol may be set as an alias in your shell. Please run 'pymol -qc log.pml'")
        print(cmd)
        sys.exit()


import argparse

##########
# SYNTAX examples:
# -pdb res_12.pdb
# -ignore_ids A:25,A:26
# -ignore_atoms A:25:C+O,A:26:N
# -ignore_atnames NH1,NH2
##########

parser = argparse.ArgumentParser()
parser.add_argument("-pdb", nargs="+")
parser.add_argument("-ignore_ids")
parser.add_argument("-ignore_atoms")
parser.add_argument("-ignore_atnames")
args = parser.parse_args()
ignoreid = args.ignore_ids
ignoreatom = args.ignore_atoms
ignoreatname = args.ignore_atnames

# Process ignore ids and atoms
notres = ""
if ignoreid != None and ignoreid != '' and ignoreid != ' ':
    ignoreid = ignoreid.split(',')
    for i in ignoreid:
        res_id = i.split(':')
        res = f" and not (chain {res_id[0]} and resi {int(res_id[1])})"
        notres += res

notats = ""
if ignoreatom != None and ignoreatom != '' and ignoreatom != ' ':
    notats = ""
    ignoreatom = ignoreatom.split(',')
    for i in ignoreatom:
        res_id = i.split(':')
        atlist = res_id[2].split('+')
        ats = "name " + str(atlist[0])
        for a in atlist[1:]:
            ats += f" or name {a}"
        res = f" and not (chain {res_id[0]} and resi {int(res_id[1])} and ({ats}))"
        notats += res
        
notatn = "name NH1 or name NH2"
if ignoreatname != None and ignoreatname != '' and ignoreatname != ' ':
    ignoreatname = ignoreatname.split(',')
    for i in ignoreatname:
        atom = f" or name {i}"
        notatn += atom

with open("log.pml", "w") as logf:
    for pdbfilename in args.pdb:
        name = os.path.splitext(pdbfilename)[0]
        outputfilename = f"{name}_h.pdb"
        logf.write(f"load {pdbfilename}\n")
        logf.write(f'cmd.select("sel","{name}{notres}{notats} and not ({notatn})")\n')
        logf.write('cmd.h_add("sel")\n')
        logf.write(f'cmd.save("./{outputfilename}")\n')

cmd = "pymol -qc log.pml"
system_run(cmd)
