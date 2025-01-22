#!/usr/bin/env python3
import sys, os, re
from read_write_pdb import *
import pandas as pd

def get_model_FGs(f,seedlist):
    pdb, res_info, tot_charge_t = read_pdb(f)
    pdbres={}
    pdbats={}
    FGs = []
    for line in pdb:
        reskey=(line[5],int(line[6]))
        if reskey not in pdbres.keys():
            pdbres[reskey]=line[4].strip()
        atom=line[2].strip()
        if reskey not in pdbats.keys():
            pdbats[reskey]=[atom]
        else:
            pdbats[reskey].append(atom)
    for res in sorted(pdbats.keys()):
        resn=pdbres[res]
        if res in seedlist or resn == 'HOH' or resn == 'WAT':
            FGs.append((f'{res[0]}:{res[1]}',resn))
        elif resn == 'GLY':
            FGs.append((f'{res[0]}:{res[1]}',f'SC_{resn}'))
            if set(['C','O']).issubset(pdbats[res]):
                FGs.append((f'{res[0]}:{res[1]}','MC'))
        else:
            if 'CB' in pdbats[res]:
                FGs.append((f'{res[0]}:{res[1]}',f'SC_{resn}'))
            if set(['C','O']).issubset(pdbats[res]):
                FGs.append((f'{res[0]}:{res[1]}','MC'))
    return FGs

def write_model_building(seedlist,mlist,seednamed):
    modelfgs=[]
    allfgs=[]
   
    for i in mlist:
        model = f"res_{i}.pdb"
        fgs=get_model_FGs(model,seedlist)
        modelfgs.append(fgs)
        for group in fgs:
            if group not in allfgs and group not in seednamed:
                allfgs.append(group)

    ### uncomment line below to have table ordered by chain sequence rather than when added to model
    #allfgs=sorted(sorted(allfgs, key = lambda x : x[1], reverse = True), key = lambda x : x[0])
    allfgs = seednamed + allfgs

    tab=[]
    for k in allfgs:
        fgline=[]
        for m in modelfgs:
            if k in m:
                fgline.append('x ')
            else:
                fgline.append(' ')
        tab.append(fgline)

    tab=pd.DataFrame(tab)
    tab.index=[f'{fg[0]}:{fg[1]}' for fg in allfgs]
    #tab.columns=[f'mod{m}' for m in mlist]
    tab.columns=mlist

    with open('seq_model_contents.dat', 'w') as f:
        f.write(tab.to_string())
        f.write('\n\n')

        for i in range(len(tab.columns)):
            for j in range(i+1,len(tab.columns)):
                if tab.iloc[:,i].equals(tab.iloc[:,j]):
                    f.write(f'Models {mlist[i]} and {mlist[j]} are the same!\n')

    return

