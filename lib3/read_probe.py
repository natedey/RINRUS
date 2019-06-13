"""
This is a program written by Qianyi Cheng in DeYonker Research Group
at University of Memphis.
Version 1.0
"""

import sys, os
from numpy import *
from res_atoms import *

def get_sel_atoms(pdb,res_list,res_atom,res_name,res_info,pdb_res_name):   # res_list in format ('A',400,'A',62)
    #pdb [ record, serial, atomname, altloc, resname, chain, resnum, achar, x, y, z, occ, tfactor, segid, elsymbol, charge.strip(), fix ] 
    # 1      0       1       2       3       4           5   6       7      8  9  10 11   12         13      14      15              16
    for i in range(0,len(res_list),2):
        chain = res_list[i]
        resnum = res_list[i+1]
        key = (chain,resnum)
        res_atom[key] = []
        for j in pdb:
            pdb_res_name[(j[5],j[6])] = j[4].strip()
            if j[5] == chain and j[6] == resnum:
                res_atom[key].append(j[2].strip())
                resname = j[4].strip()
                res_name[key] = resname
                if resname in res_atoms_all.keys():
                    res_info[key] = ['CA']
                else:
                    res_info[key] = []

    return res_atom, res_name, res_info, pdb_res_name


def get_side(column):
    ### case1 : A1138 HOH  O   : ###
    ### case2 : A 137 TYR  CE2B: ###
    ###        0123456789012345
    side = []
    side.append(column[1])                  # chain
    side.append(int(column[2:6]))           # res_id
    side.append(column[6:10].strip())       # res_name
    side.append(column[10:-1].strip())      # atom_name
    if column[-1] != ' ':                   # uncertainty
        side.append(column[-1])
    else:
        side.append('A')
    return side


def analyse_side(side,res_name,res_count,res_atom):
    key = (side[0],side[1])
    if key not in res_name.keys():
        res_name[key] = side[2]
        res_count[key] = 1          # count based on chain+res_id
        res_atom[key] = [side[3]]
    else:
        res_count[key] += 1
        if side[3] not in res_atom[key]:
            res_atom[key].append(side[3])
    return res_name, res_count, res_atom


def get_probe_atoms(probefile,res_name,res_atom,res_cout):     # # res_list in format ('A',400,'A',62)
    sel_keys = res_name.keys()
    with open(probefile) as f:
        lines = f.readlines()
    for i in range(len(lines)):
        c = lines[i].split(':')
        #c[1] : '1->1'
        #c[2] : 'wc':wide contact,'cc': close contact,'so':small overlap,'bo':big overlap,'hb':hydrogen bond
        #c[3] : ' A  38 TRP  C   ' or 'A1138 HOH  O   '
        #c[4] : ' A  38 TRP  CD2 '
        #c[5],c[6]: '0.025', '0.335'
        #c[7]-c[9]: '-10.370', '-17.062', '-19.661'
        #c[10]: '0.000'   score
        #c[11]: '0.0104'  raw_score in rinalyzer
        #c[12]: 'C'
        #c[13]: 'C'
        #c[14]-c[16]:'-10.370', '-17.062', '-19.661'
        #c[17],c[18]: '30.32', '26.49
        
        l_side = get_side(c[3])
        r_side = get_side(c[4])
        if (l_side[0],l_side[1]) in sel_keys and (r_side[0],r_side[1]) not in sel_keys:
            res_name, res_cout, res_atom = analyse_side(r_side,res_name,res_cout,res_atom)
        if (l_side[0],l_side[1]) not in sel_keys and (r_side[0],r_side[1]) in sel_keys:
            res_name, res_cout, res_atom = analyse_side(l_side,res_name,res_cout,res_atom)
    return res_atom, res_name, res_cout


def get_res_list(res_atom):
    res_list = {}
    for key in sorted(res_atom.keys()):
        if key[0] not in res_list.keys():
            res_list[key[0]] = [key[1]]
        else:
            res_list[key[0]].append(key[1])
    return res_list


def check_mc(res,value,atoms):
    case1 = ['N','H']
    case2 = ['C','O']
    case3 = ['N','CA','C','O','H','HA','HA2','HA3']
    case4 = ['CA','HA','HA2','HA3']
    if bool(set(value)&set(case1)) and 'C' not in value and 'O' not in value:
        for i in ['N','CA','H','HA','HA2','HA3']:
            atoms.append(i)
    elif bool(set(value)&set(case2)) and 'N' not in value and 'H' not in value:
        for i in ['CA','C','O','HA','HA2','HA3']:
            atoms.append(i)
    elif 'N' not in value and 'H' not in value and 'C' not in value and 'O' not in value:
        for i in ['CA','HA','HA2','HA3']:
            atoms.append(i)
    elif bool(set(value)&set(case1)) and bool(set(value)&set(case2)):
        for i in case3:
            atoms.append(i)
    else:
        for i in case3:
            atoms.append(i)
    return atoms

def check_sc(res,value,atoms):
    for i in res_atoms_sc[res]:
        atoms.append(i)
    return atoms

def get_res_parts(res,value):   # res is res_name[key], value is res_atom[key] check mainchain sidechain
    case1 = ['N','H']
    case2 = ['C','O']
    case3 = ['N','CA','C','O','H','HA','HA2','HA3']

    #if res in ['PRO','GLY']:
    #    print res
    #    atoms = res_atoms_all[res]
    if res == 'PRO':
        atoms = res_atoms_all[res]
    else:
        atoms = []
        for i in value:
            atoms.append(i)
#        print res, bool(set(value)&set(case3)), bool(set(value)&set(res_atoms_sc[res]))
        if bool(set(value)&set(case3)) and bool(set(value)&set(res_atoms_sc[res])):
            atoms = check_mc(res,value,atoms)
            atoms = check_sc(res,value,atoms)
        elif bool(set(value)&set(case3)) and not bool(set(value)&set(res_atoms_sc[res])):
            atoms = check_mc(res,value,atoms)
        elif bool(set(value)&set(res_atoms_sc[res])) and not bool(set(value)&set(case3)):
            atoms = check_sc(res,value,atoms)
        else:
            print("Something is wrong!")

    return atoms


def check_s(chain,id,atom,res_info,nres_atom):
    if (chain,id) not in nres_atom.keys():
        nres_atom[(chain,id)] = []
        for i in atom:
            nres_atom[(chain,id)].append(i)
    else:
        for i in atom:
            nres_atom[(chain,id)].append(i)
            
    if 'CA' in atom:
        if (chain,id) in res_info.keys():
            res_info[(chain,id)].append('CA')
        else:
            res_info[(chain,id)] = ['CA']
    else:
        if 'CB' in atom:
            nres_atom[(chain,id)].append('CA')
            if (chain,id) in res_info.keys():
                res_info[(chain,id)].append('CA')
                res_info[(chain,id)].append('CB')
            else:
                res_info[(chain,id)] = ['CA','CB']
    return nres_atom, res_info

def check_a(chain,id,atom,nres_info,nres_atom,res_name,pdb_res_name):
    if 'C' in atom and (chain,id+1) in pdb_res_name.keys():
        if (chain,id+1) not in nres_atom.keys(): 
            nres_atom[(chain,id+1)] = ['N','CA']
        else:
            nres_atom[(chain,id+1)].append('N')
            nres_atom[(chain,id+1)].append('CA')
        if (chain,id+1) not in nres_info.keys():
            nres_info[(chain,id+1)] = ['CA']
        else:
            nres_info[(chain,id+1)].append('CA')
        if (chain,id+1) not in res_name.keys():
            res_name[(chain,id+1)] = pdb_res_name[(chain,id+1)]
    return nres_atom, nres_info, res_name

def check_b(chain,id,atom,nres_info,nres_atom,res_name,pdb_res_name):
    if 'N' in atom and (chain,id-1) in pdb_res_name.keys():
        if (chain,id-1) not in nres_atom.keys():
            nres_atom[(chain,id-1)] = ['CA','C','O']
        else:
            nres_atom[(chain,id-1)].append('CA')
            nres_atom[(chain,id-1)].append('C')
            nres_atom[(chain,id-1)].append('O')
        if (chain,id-1) not in nres_info.keys():
            nres_info[(chain,id-1)] = ['CA']
        else:
            nres_info[(chain,id-1)].append('CA')
        if (chain,id-1) not in res_name.keys():
            res_name[(chain,id-1)] = pdb_res_name[(chain,id-1)]
    return nres_atom, nres_info, res_name


def check_bb(chain,id,nres_atom):
    if 'CA' in nres_atom[(chain,id)]:
        if (chain,id-1) in nres_atom.keys() and 'CA' in nres_atom[(chain,id-1)]:
            nres_atom[(chain,id)].append('N') 
            nres_atom[(chain,id)].append('H') 
        if (chain,id+1) in nres_atom.keys() and 'CA' in nres_atom[(chain,id+1)]:
            nres_atom[(chain,id)].append('C')
            nres_atom[(chain,id)].append('O')
    return nres_atom

def final_pick(pdb,nres_atom,nres_info):
    res_pick = []
    for line in pdb:
        if (line[5],line[6]) in nres_atom.keys() and line[2].strip() in nres_atom[(line[5],line[6])]:
            if (line[5],line[6]) in nres_info.keys() and line[2].strip() in nres_info[(line[5],line[6])]:
                res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'] )
            else:
                res_pick.append( [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],' 0'] )
    return res_pick

