"""
This is a program written by qianyi cheng in deyonker research group
at university of memphis.
Version 1.0
"""
import os, sys, re
from read_write_pdb import *
from read_probe import *
from copy import *

##############   Example   ########################################
#pdb = read_pdb('../comt2/3bwm_h_mg.ent')
#probe = '../comt2/3bwm_h.probe'
#freqf = 'freq_per_res.dat'
#idx_list = ['A',300,'A',301,'A',302]
###################################################################


pdb, res_info, tot_charge = read_pdb(sys.argv[1])
probe = sys.argv[2]
freqf = sys.argv[3]
list4 = sys.argv[4]

idx_list = []
c = list4.split(',')
for i in range(0,len(c),2):
    idx_list.append(c[i])
    idx_list.append(int(c[i+1]))
print(idx_list)

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
    print(lines[i].split())
    c = lines[i].split()
    cha = c[0]
    res = int(c[1])
    freq = int(c[2])
    if (cha,res) in sel_key: continue
    j += 1
    qf[j] = deepcopy(qf[j-1])
    qf[j][cha].append(res)

### read in probe file ###
res_atom, res_name, res_cout = get_probe_atoms(probe,res_name,res_atom,res_cout)

#res_list = get_res_list(res_atom)
#nres_list = get_res_list(res_atom)
for nm_res in sorted(qf.keys()):
    res_list = qf[nm_res]
    print(nm_res, qf[nm_res])
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
    outf = 'res_%s.pdb'%str(nm_res)
    write_pdb(outf,res_pick)
