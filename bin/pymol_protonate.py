#!/usr/bin/env python3
"""
This is a program written by qianyi cheng in deyonker research group
at university of memphis.
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

parser = argparse.ArgumentParser()
parser.add_argument("-pdb", nargs="+")
parser.add_argument("-ignore_ids")
parser.add_argument("-ignore_ats")
args = parser.parse_args()
ignoreid = args.ignore_ids
ignoreat = args.ignore_ats

# Process ignore ids and atoms
ignoreid = ignoreid.split(',')
notres = ""
for i in ignoreid:
    res_id = i.split(':')
    res = f" and not (chain {res_id[0]} and resi {int(res_id[1])})"
    notres += res
notat = "name NH1 or name NH2"
if ignoreat is not None:
    ignoreat = ignoreat.split(',')
    for i in ignoreat:
        atom = f" or name {i}"
        notat += atom

with open("log.pml", "w") as logf:
    for pdbfilename in args.pdb:
        name = os.path.splitext(pdbfilename)[0]
        outputfilename = f"{name}_h.pdb"
        logf.write(f"load {pdbfilename}\n")
        if args.ignore_ids is not None:
            logf.write(
                f'cmd.select("sel","{name}{notres} and not ({notat})")\n'
            )
            logf.write('cmd.h_add("sel")\n')
        else:
            logf.write(f'cmd.h_add("{name}")\n')
        logf.write(f'cmd.save("./{outputfilename}")\n')

cmd = "pymol -qc log.pml"
system_run(cmd)
