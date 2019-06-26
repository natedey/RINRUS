#!/usr/bin/env python

#This code is done by tsmmers1 and is a modification of qcheng1 code.
#Its function is to extract the PDB coordinates for the combinatoric cluster models
#Usage: ./scriptname.py 

import os, sys, re
from PDB import read_pdb
from read_probe import *
from copy import *
from string import ascii_lowercase

pdb = read_pdb('/home/tsmmers1/chem/comt/3bwm_h_mg_ts_noGlu199.ent')
idx_list = ['A',300,'A',301,'A',302,'A',411]
probe = '/home/tsmmers1/chem/comt/3bwm_h.probe'
freqf = '/home/tsmmers1/chem/comt/pdbs/3bwm_wPropCoord_noGlu199.txt'


res_name = {}
res_atom = {}
res_cout = {}
res_info = {}
pdb_res_name = {}
cres_atom = {}

### get key residues ###
res_atom, res_name, res_info, pdb_res_name = get_sel_atoms(pdb,idx_list,res_atom,res_name,res_info,pdb_res_name)
sel_key = res_name.keys()   ### sel_key = [('A',300),('A',301),('A',302)]

### Sort residue by frequency ###
qf = {}
qf[len(sel_key)] = {}
for i in sel_key:
    try:
        qf[len(sel_key)][i[0]].append(i[1])
    except:
        qf[len(sel_key)][i[0]] = [i[1]]

with open(freqf) as f:
    lines = f.readlines()
sm = len(lines)
j = len(sel_key)
for i in range(sm):
    c = lines[i].split(':')
    c = re.split(r'[(,)]', c[0])
    c = list(filter(None, c))
    c = [int(k) for k in c]
    cha = 'A'
    j += 1
    qf[j] = deepcopy(qf[j-1])
    qf[j][cha] = c

### read in probe file ###
res_atom, res_name, res_cout = get_probe_atoms(probe,res_name,res_atom,res_cout)

#res_list = get_res_list(res_atom)
#nres_list = get_res_list(res_atom)
from os.path import isfile, join
mypath = os.getcwd()
names_list = [f for f in os.listdir(mypath) if isfile(join(mypath, f))]

keyfile = open("fileskey.txt", "a")

for nm_res in sorted(qf.keys()):
    res_list = qf[nm_res]
    print nm_res, qf[nm_res]
    for key in res_list.keys():
        for res in sorted(res_list[key]):
            if (key,res) in sel_key or res_name[(key,res)] == 'HOH':
                cres_atom[(key,res)] = res_atom[(key,res)]
            else:
                cres_atom[(key,res)] = get_res_parts(res_name[(key,res)],res_atom[(key,res)])


    nres_atom = {}
    for key in res_list.keys():
        for res in sorted(res_list[key]):
            key1 = (key,res)
            if res_name[key1] == 'HOH': 
                res_info[key1] = []
                nres_atom[key1] = cres_atom[key1]
            else:
                nres_atom, res_info, res_name = check_b(key,res,cres_atom[key1],res_info,nres_atom,res_name,pdb_res_name)
                nres_atom, res_info = check_s(key,res,cres_atom[key1],res_info,nres_atom)
                nres_atom, res_info, res_name = check_a(key,res,cres_atom[key1],res_info,nres_atom,res_name,pdb_res_name)
    
    
    for key in sorted(nres_atom.keys()):
        nres_atom = check_bb(key[0],key[1],nres_atom)
    
    res_pick = final_pick(pdb,nres_atom,res_info)
    setsize = str(res_list).split('[')
    setsize = setsize[1].split(']')
    setsize = str(len(setsize[0].split(',')))
    setname = ""
    for let1 in ascii_lowercase:
        mybool = False
        for let2 in ascii_lowercase:
            tryname = 'res_' + setsize + '-' + let1 + let2 + '.pdb'
            if tryname not in names_list:
                 setname = tryname
                 names_list.append(tryname)
                 mybool = True
                 break
        if mybool == True:
            break
    outf = setname
    print outf
    write_pdb(outf,res_pick)
    keyfile.write(outf + " : " + str(qf[nm_res]) + "\n")

keyfile.close()
