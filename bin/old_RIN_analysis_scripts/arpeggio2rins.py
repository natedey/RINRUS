#!/usr/bin/env python3
"""
This is a program created by Dr. Qianyi Cheng and the DeYonker research group
at The University of Memphis.
Date created 08.09.2022
"""
from numpy import *
import os, sys, re, argparse
from res_atoms import *

def check_atom(key,atom,dic):
    if key not in dic.keys():
        dic[key] = [atom]
    else:
        if atom not in dic[key]:
            dic[key].append(atom)
    return dic

def de_contacts(contactsf,seed,opt):
    ### get seed information ###
    seeds_in = seed.split(',')
    seeds = []
    for sel in seeds_in:
        selstr = sel.replace(':','/')
        selstr = selstr + '/'
        seeds.append(selstr)
    ### Arpeggio interaction types 16 total ###
    rnames = ['clash','covalent','vdw_clash','vdw','proximal','hbond','weak_hbond','halogen_bond','ionic','metal_complex','aromatic','hydrophobic','carbonyl','polar','weak_polar']

    rinrus = {}  ### format dict[(chain,resid)] = [all atoms interact with seed]
    resinfo = {} ### format dict[(chain,resid)] = [each interaction' counts]
    inter_detail = {} ### format dict[(chain/resid/atom, chain/resid/atom)] = [each interaction' counts]
    #Interacting entities    string from (INTER, INTRA_NON_SELECTION, INTRA_SELECTION, SELECTION_WATER, NON_SELECTION_WATER, WATER_WATER)    Distinguishes how this atom pair relates to the selected atoms: see below

    dat = genfromtxt(contactsf,dtype=str)
    dat_gp1 = dat[:,0]
    dat_gp2 = dat[:,1]
    dat_int = dat[:,2:-1].astype(int)
    dat_end = dat[:,-1]
    for i in range(len(dat_gp1)):
        chain1,resid1,aname1 = dat_gp1[i].split('/')
        chain2,resid2,aname2 = dat_gp2[i].split('/')
        st1 = chain1 + '/' + resid1 + '/'
        st2 = chain2 + '/' + resid2 + '/'
        key1 = (chain1,int(resid1))
        key2 = (chain2,int(resid2))
        ### get residue/atom interactions ###
        if st1 in seeds and st2 not in seeds:
            key_i = (dat_gp1[i],dat_gp2[i]) 
            if key_i not in inter_detail.keys():
                inter_detail[key_i] = dat_int[i]
            else:
                inter_detail[key_i] += dat_int[i]
        #elif st2 in seeds and st1 not in seeds:
        elif st1 not in seeds and st2 in seeds:
            key_i = (dat_gp2[i],dat_gp1[i]) 
            if key_i not in inter_detail.keys():
                inter_detail[key_i] = dat_int[i]
            else:
                inter_detail[key_i] += dat_int[i]
        ### get residue interactions ###
        if st1 in seeds and st2 not in seeds:
            rinrus = check_atom(key1,aname1,rinrus)
            key = key1 
            if key not in resinfo.keys():
                resinfo[key] = dat_int[i]
            else:
                resinfo[key] = resinfo[key] + dat_int[i]
            rinrus = check_atom(key2,aname2,rinrus)
            key = key2 
            if key not in resinfo.keys():
                resinfo[key] = dat_int[i]
            else:
                resinfo[key] = resinfo[key] + dat_int[i]
        elif st1 not in seeds and st2 in seeds:
            rinrus = check_atom(key2,aname2,rinrus)
            key = key2 
            if key not in resinfo.keys():
                resinfo[key] = dat_int[i]
            else:
                resinfo[key] = resinfo[key] + dat_int[i]
            rinrus = check_atom(key1,aname1,rinrus)
            key = key1 
            if key not in resinfo.keys():
                resinfo[key] = dat_int[i]
            else:
                resinfo[key] = resinfo[key] + dat_int[i]
        
    f1 = open('node_info.dat','w')
    for key in sorted(inter_detail.keys()):
        f1.write('%s %s '%(key[0],key[1]))
        if opt == 0:
            for l in range(len(inter_detail[key])):
                if inter_detail[key][l] != 0:
                    f1.write('%s %d '%(rnames[l],inter_detail[key][l]))
        elif opt == 1:
            for l in range(len(inter_detail[key])):
                if l != 4 and inter_detail[key][l] != 0:
                    f1.write('%s %d '%(rnames[l],inter_detail[key][l]))
        f1.write('\n')

    f1.close()

    f2 = open('contype_counts.dat','w')
    f2.write(f'# selected/ranked residues by arpeggio: contact file {cf}; seed {seeds_in}; contact types\n')
    f3 = open('contact_counts.dat','w')
    f3.write(f'# selected/ranked residues by arpeggio: contact file {cf}; seed {seeds_in}; contact counts\n')

    if opt == 0:
        newlist = sorted(resinfo.items(),key=lambda x:sum(x[1]),reverse=True)
        newlist_t = sorted(resinfo.items(),key=lambda x:count_nonzero(x[1]),reverse=True)
    elif opt == 1:
        for key in resinfo.keys():
            resinfo[key][4] = 0
        newlist = sorted(resinfo.items(),key=lambda x:sum(x[1]),reverse=True)
        newlist_t = sorted(resinfo.items(),key=lambda x:count_nonzero(x[1]),reverse=True)

    sort_r = dict(newlist)
    sort_t = dict(newlist_t)
    for key in sort_t.keys():
        for it in key:
            f2.write('%8s'%it)
        f2.write('%8d'%count_nonzero(resinfo[key]))
        for it in rinrus[key]:
            f2.write('%6s'%it)
        f2.write('\n')    
    for key in sort_r.keys():
        for it in key:
            f3.write('%8s'%it)
        f3.write('%8d'%sum(resinfo[key]))
        for it in rinrus[key]:
            f3.write('%6s'%it)
        f3.write('\n')    
    f2.close()
    f3.close()

if __name__ == '__main__':

    ################################################################################
    ### python3 /home/qcheng1/projects/p450-md/farpeggio.py -f /home/qcheng1/projects/p450-md/syringol/arpeggio-run/md3/md3_55_align25ncb.contacts -s (:404,):403 -p 1 ###
    ################################################################################

    parser = argparse.ArgumentParser(description='Get contact information from Arpeggio output')
    parser.add_argument('-p', dest='opt', default=1, type=int, help='ignor proximal or not')
    parser.add_argument('-f', dest='cf', help='Arpeggio contact output file')
    parser.add_argument('-s', dest='seed', help='Seed residues, in the format "A:601,A:602"')

    args = parser.parse_args()
    opt = args.opt
    cf = args.cf
    seed = args.seed

    if opt == 0:
        de_contacts(cf,seed,0)        
    elif opt == 1:
        de_contacts(cf,seed,1)
    else:
        print("If you want to not ignor proximal please add '-p 0' when running the script!")

