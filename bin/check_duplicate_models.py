#!/usr/bin/env python3
"""
This is a program created by Dr. Dominique Wappett and the DeYonker research group
at The University of Memphis.
Date created 01.15.2025
"""
import sys, os, re
from read_write_pdb import *
from numpy import *
from itertools import *

if __name__ == '__main__':
    ### Compare pdbs to each other
    ### usage:  check_duplicate_models.py model1.pdb model2.pdb ...
    
    pdbfiles=sys.argv[1:]
    pdblengths=[]
    pdbdictsall=[]
    pdbdictsheavy=[]

    for f in pdbfiles:
        pdb, res_info, tot_charge_t = read_pdb(f)
        pdblengths.append(len(pdb))
        pdbats={}
        pdbheavy={}
        for i in range(len(pdb)):
            reskey=(pdb[i][5],int(pdb[i][6]))
            atom=pdb[i][2].strip()
            if reskey not in pdbats.keys():
                pdbats[reskey]=[atom]
            else:
                pdbats[reskey].append(atom)
            if atom[0] != 'H':
                if reskey not in pdbheavy.keys():
                    pdbheavy[reskey]=[atom]
                else:
                    pdbheavy[reskey].append(atom)
        pdbdictsall.append(pdbats)
        pdbdictsheavy.append(pdbheavy)

    matches=0
    maybes=0
    for i in range(len(pdbfiles)):
        pdb1=pdbfiles[i]
        for j in range(i+1,len(pdbfiles)):
            pdb2=pdbfiles[j]
            checklen=0
            checkres=0
            checkheavy=0
            checkall=0
            #compare lengths
            if pdblengths[i] != pdblengths[j]:
                continue
            #compare res ids
            if sorted(pdbdictsall[i].keys()) != sorted(pdbdictsall[j].keys()) or sorted(pdbdictsheavy[i].keys()) != sorted(pdbdictsheavy[j].keys()):
                continue
            #compare heavy atoms for each res id
            mismatch=0
            for key in pdbdictsheavy[i].keys():
                if sorted(pdbdictsheavy[i][key]) != sorted(pdbdictsheavy[j][key]):
                    mismatch+1
            if mismatch > 0:
               continue 
            #compare all atoms for each res id
            mismatch=0
            for key in pdbdictsall[i].keys():
                if sorted(pdbdictsall[i][key]) != sorted(pdbdictsall[j][key]):
                    mismatch+=1
            if mismatch == 0:
                print(f"{pdb1} and {pdb2} are the same model")
                matches+=1
            else:
                print(f"{pdb1} and {pdb2} only differ in the hydrogens")
                maybes+=1

    if matches == 0 and maybes == 0:
        print("all these models are distinct")

