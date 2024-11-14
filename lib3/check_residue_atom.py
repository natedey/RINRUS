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
    if res == 'PRO':
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

def check_sc(res,value,cres_atoms_sc):
    if res == 'PRO':
        value = res_atoms_all[res]
    elif res != 'PRO' and res in res_atoms_sc.keys():
        if bool(set(value)&set(res_atoms_sc[res])):
            for i in res_atoms_sc[res]:
                if i not in value:
                    value.append(i)
    elif res != 'PRO' and res in cres_atoms_sc.keys():
        if bool(set(value)&set(cres_atoms_sc[res])):
            for i in cres_atoms_sc[res]:
                if i not in value:
                    value.append(i)
    else:   
        print("Residue %s is either canonical or nocanonical, please check!"%res)
    return value

def get_noncanonical_resinfo(cres):
    with open(cres) as f:
        lines = f.readlines()
    cres_atoms_all = {}
    cres_atoms_sc = {}
    for i in range(len(lines)):
        c = lines[i].split()
        if c[0] not in cres_atoms_all.keys():
            cres_atoms_all[c[0]] = []
            for j in range(1,len(c)):
                cres_atoms_all[c[0]].append(c[j])
        else:
            if c[0] not in cres_atoms_sc.keys():
                cres_atoms_sc[c[0]] = []
                for j in range(1,len(c)):
                    cres_atoms_sc[c[0]].append(c[j])
    return cres_atoms_all, cres_atoms_sc

def final_pick(pdb,res_atom,res_info,sel_key):
    list_cb = ['ARG','LYS','GLU','GLN','MET','TRP','TYR','PHE']
    res_pick = []
    for line in pdb:
        if (line[5],line[6]) in res_atom.keys() and line[2].strip() in res_atom[(line[5],line[6])]:
            if line[2].strip() == 'CA':
                res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'] )
                if (line[5],line[6]) not in res_info.keys():
                    res_info[(line[5],line[6])] = ['CA']
                else:
                    if 'CA' not in res_info[(line[5],line[6])]:
                        res_info[(line[5],line[6])].append('CA')
            elif line[2].strip() == 'CB':
                if line[4].strip() in list_cb and 'CG' in res_atom[(line[5],line[6])]:
                    res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'] )
                    if 'CB' not in res_info[(line[5],line[6])]:
                        res_info[(line[5],line[6])].append('CB')
                else:
                    res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],' 0'] )
            else:
                res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],' 0'] )
    return res_pick,res_info

def final_pick2(pdb,res_atom,res_info,sel_key):
    list_cb = ['ARG','LYS','GLU','GLN','MET','TRP','TYR','PHE']
    res_pick = []
    for line in pdb:
#        if (line[5],line[6]) in res_atom.keys() and line[4].strip() in list_cb:
#            if line[2].strip() in res_atom[(line[5],line[6])]:
#                if line[2].strip() == 'CB':
#                    res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'] )
#                    if 'CB' not in res_info[(line[5],line[6])]:
#                        res_info[(line[5],line[6])].append('CB')
#                else:
#                    res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],' 0'] )
        if (line[5],line[6]) in res_atom.keys() and line[2].strip() in res_atom[(line[5],line[6])]:
            if line[2].strip() == 'CB':
                if line[4].strip() in list_cb and (line[5],line[6]) not in sel_key:
                    res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'] )
                    if 'CB' not in res_info[(line[5],line[6])]:
                        res_info[(line[5],line[6])].append('CB')
                else:
                    if 'CB' in res_info[(line[5],line[6])]:
                        res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'] )
                    else:
                        res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],' 0'] )

            else:
                if line[2].strip() in res_info[(line[5],line[6])]:
                    res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'] )
#            elif line[2].strip() not in res_info[(line[5],line[6])]:
#                if line[2].strip() == 'CA':
#                    res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'] )
#                    res_info[(line[5],line[6])] = ['CA']
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

