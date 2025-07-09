"""
This is a program written by Qianyi Cheng in DeYonker Research Group
at University of Memphis.
Version 1.1
Date 6.18.2020
"""

import sys, os, re
from numpy import *
from res_atoms import *


def check_mc(res,value):
    case1 = ['N','H']
    case2 = ['C','O']
    case3 = ['N','CA','C','O','H','HA','HA2','HA3']
    case4 = ['CA','HA','HA2','HA3']
    ### if res name is PRO, then include entire residue ###
    ### DAW 2025-06-16: unless only C-terminus, that can be included on its own ###
    if res == 'PRO':
        pro_sc_nt = set(res_atoms_sc['PRO']+['N'])
        if not bool(set(value)&pro_sc_nt) and ('C' in value or 'O' in value):
            for i in ['CA','C','O','HA','HA2','HA3']:
                if i not in value:
                    value.append(i)
        else:
            value = res_atoms_all[res]

    if bool(set(value)&set(case1)) and 'C' not in value and 'O' not in value:
        for i in ['N','CA','H','HA','HA2','HA3']:
            if i not in value:
                value.append(i)
    elif bool(set(value)&set(case2)) and 'N' not in value and 'H' not in value:
        for i in ['CA','C','O','HA','HA2','HA3']:
            if i not in value:
                value.append(i)
    elif 'N' not in value and 'H' not in value and 'C' not in value and 'O' not in value:
        for i in ['CA','HA','HA2','HA3']:
            if i not in value:
                value.append(i)
    elif bool(set(value)&set(case1)) and bool(set(value)&set(case2)):
        for i in case3:
            if i not in value:
                value.append(i)
    else:
        for i in case3:
            value.append(i)
    return value

def check_sc(res,value,ncres_atoms_sc):
    if res == 'PRO':
        if bool(set(value)&set(res_atoms_sc[res])):
            value = res_atoms_all[res]
    elif res in res_atoms_sc.keys():
        if bool(set(value)&set(res_atoms_sc[res])):
            for i in res_atoms_sc[res]:
                if i not in value:
                    value.append(i)
    elif res in ncres_atoms_sc.keys():
        if bool(set(value)&set(ncres_atoms_sc[res])):
            for i in ncres_atoms_sc[res]:
                if i not in value:
                    value.append(i)
    else:   
        print("Residue %s not recognized so it might be trimmed incorrectly, please check!"%res)
    return value

def get_noncanonical_resinfo(ncres):
    with open(ncres) as f:
        lines = f.readlines()
    ncres_atoms_all = {}
    ncres_atoms_sc = {}
    for line in lines:
        c = line.split()
        if c[0] not in ncres_atoms_all.keys():
            ncres_atoms_all[c[0]] = []
            ncres_atoms_sc[c[0]] = []
        for at in c[1:]:
            ncres_atoms_all[c[0]].append(at)
            if at not in ['C','O','N','H','CA','HA','HA2','HA3']:
                ncres_atoms_sc[c[0]].append(at)
    return ncres_atoms_all, ncres_atoms_sc

def final_pick2(pdb,res_atom,res_info,sel_key):
    res_pick = []
    for line in pdb:
        key = (line[5],line[6])
        if key in res_atom.keys() and line[2].strip() in res_atom[key]:
            if line[2].strip() in res_info[key]:
                res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'] )
            else:
                res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],' 0'] )
    return res_pick, res_info

def get_sel_keys(seed_list):  
    seeds = seed_list.split(',')
    sel_keys = []
    for seed in seeds:
        res_id = seed.split(':')
        sel_keys.append((res_id[0],int(res_id[1])))
    sel_keys.sort(key = lambda x: x[1]) #by res id
    sel_keys.sort(key = lambda x: x[0]) #by chain id
    return sel_keys

