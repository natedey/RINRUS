#!/usr/bin/env python
"""
This is a program created by Dr. Dominique Wappett and the DeYonker research group
at The University of Memphis.
Date created 04.09.2025
"""

import sys, os, argparse
import subprocess
from pathlib import Path
from read_write_pdb import *
import pandas as pd
import pickle

def addtodict(dct,key,val):
    try:
        dct[key].append(val)
    except:
        dct[key] = [val]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='FSAPT analysis, summary and res_atoms creation')
    parser.add_argument('-path', dest='path', default=None, help='path to PSI4 fsapt.py script')
    parser.add_argument('-save', dest='savepath', default='../', help='location for saving sorted summary file and res_atoms.dat')
    parser.add_argument('-rank', dest='rank', default='whole', help='seed part to use for Eint ranking ("whole" or one of the labels used in fA.dat)')

    args = parser.parse_args()
    savepath = args.savepath
    rank = args.rank
    if rank == 'whole':
        rank = 'Whole fA'

    ### make sure psi4 path formatted correctly for script to work, if not specified then define from psi4 location
    if args.path:
        psipath = os.path.expanduser(args.path)
        if psipath[-1] == '/':
            psipath = psipath[:-1]
    else:
        psipath = subprocess.run(['which', 'psi4'], capture_output=True, text=True).stdout.replace('\n','')
        if psipath:
            psipath = Path(psipath).resolve().parents[1]
            psipath = psipath / "share" / "psi4" / "fsapt"
            psipath = str(psipath)
        else:
            print("Can't find psi4 to get fsapt.py analysis script! Check conda env or specify path with -path flag")
            sys.exit()

    ### run fsapt analysis
    subprocess.run(["python",psipath+"/fsapt.py"])

    ### collect total values and sort
    output = open('fsapt.dat', 'r').read()
    output = [sect for sect in output.split('\n\n') if sect != '']
    output = output[-1].split('\n')[1:]
    
    vals = {}
    fgid = []
    for line in output:
        l = line.split()
        if l[0] not in vals.keys():
            vals[l[0]]=[]
        if l[1].lower() != 'all' and l[1].lower() != 'cap_ats' :
            vals[l[0]].append(float(l[-1]))
            if l[1] not in fgid:
                fgid.append(l[1])
    df = pd.DataFrame.from_dict(vals)
    df.index = fgid
    df = df.rename(columns={'All': 'Whole fA'})
    df = df.sort_values(by=rank, key=abs, ascending=False)
    fg_sorted = df.index.values

    ### load atom dictionary
    with open('fdict.pkl', 'rb') as fp:
        atoms = pickle.load(fp)

    ### write res atoms
    with open(f'{savepath}res_atoms_fsapt.dat','w') as f:
        f.write(f'# ranked functional groups by F-SAPT |Eint| for seed part: {rank}\n')
        for group in atoms['seed'].keys():
            gsp = group.split(':')
            ats = [f'{i:<4}' for i in atoms['seed'][group]]
            f.write(f"{gsp[0]:<4} {gsp[1]:<8}     seed  {' '.join(ats)}\n")
        for group in atoms['other_fA'].keys():
            gsp = group.split(':')
            ats = [f'{i:<4}' for i in atoms['other_fA'][group]]
            f.write(f"{gsp[0]:<4} {gsp[1]:<8} other fA  {' '.join(ats)}\n")
        for group in fg_sorted:
            if group in atoms['FG'].keys():
                gsp = group.split(':')
                if gsp[-1] == 'MC':
                    ats = ['C   ','O   ']
                else:
                    ats = [f'{i:<4}' for i in atoms['FG'][group]]
                f.write(f"{gsp[0]:<4} {gsp[1]:<8} {df.loc[group,rank]:>8}  {' '.join(ats)}\n")

    ### save sorted summary file
    df.to_string(buf=f'{savepath}FG-SAPT-ranked.dat')
