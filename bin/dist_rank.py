#!/usr/bin/env python3
"""
This is a program created by Dr. Dominique Wappett and the DeYonker research group
at The University of Memphis.
Date created 03.14.2025
"""
import os, sys, re
from numpy import *
from read_write_pdb import *
from dist_tools import *
from copy import *
import argparse
import pandas as pd
import numpy as np
import pickle
from res_atoms import *


def getseeddist(cres,pdb,noH):
    # get list of seed atoms and their coords
    center_res = cres.split(',')
    if len(center_res[0].split(':')) == 2:
        idx_list = []
        for i in range(len(center_res)):
            a,b = center_res[i].split(':')
            idx_list.append((a,int(b)))
        sel_atoms = []
        sel_coord = []
        seeddict = {}
        for atom in pdb:
            if atom[14].strip() == 'H' and noH:
                continue
            elif (atom[5], atom[6]) in idx_list:
                sel_atoms.append(Atoms[atom[14].strip()][1])
                sel_coord.append([atom[8],atom[9],atom[10]])
                seeddict[(atom[5], atom[6], atom[2].strip())]=[atom[8], atom[9], atom[10]]
    elif len(center_res[0].split(':')) == 3:
        idx_list = {}
        for i in range(len(center_res)):
            a,b,c = center_res[i].split(':')
            if (a,int(b)) not in idx_list.keys():
                idx_list[(a,int(b))] = [c]
            else:
                idx_list[(a,int(b))].append(c)
        sel_atoms = []
        sel_coord = []
        seeddict = {}
        for atom in pdb:
            if atom[14].strip() == 'H' and noH:
                continue
            elif (atom[5], atom[6]) in idx_list.keys() and atom[2].strip() in idx_list[(atom[5], atom[6])]:
                sel_atoms.append(Atoms[atom[14].strip()][1])
                sel_coord.append([atom[8],atom[9],atom[10]])
                seeddict[(atom[5], atom[6], atom[2].strip())]=[atom[8], atom[9], atom[10]]
    sel_coord = array(sel_coord)
    sel_atoms = array(sel_atoms)

    # get COM and avg xyz centres as well
    com = center_of_mass(sel_atoms,sel_coord)
    avg = avg_coord(sel_coord)
    return seeddict, com, avg


def getatomdists(pdb,com,avg,seeddict,noH):
    allinfo = []
    allats = {}
    alldclos = {}
    for atom in pdb:
        # get all distance values
        if atom[14].strip() == 'H' and noH == True:
            continue
        coords = array([atom[8], atom[9], atom[10]])
        distcom = calc_dist(coords, com)
        distavg = calc_dist(coords, avg)
        distatom = {}
        for seedatom in seeddict.keys():
            distatom[seedatom] = calc_dist(coords,array(seeddict[seedatom]))
        distclosest = min(list(distatom.values()))
        
        # categorise atoms
        # atom[5] is chain id, atom[6] is resid, atom[4] is res name, atom[2] is atom name
        key = atom[5], int(atom[6])
        atname = f'{atom[2].strip():<4}'
        if atom[4].strip()=="HOH" or atom[4].strip()=="WAT":
            FG = "solv"
            dkey = f'{atom[5]}:{atom[6]}:{atom[4].strip()}'
        elif (atom[5],atom[6],atom[2].strip()) in seeddict.keys():
            FG = "seed"
            dkey = f'{atom[5]}:{atom[6]}:{atom[4].strip()}'
        elif atom[4].strip() not in res_atoms_all.keys():
            FG = "unrc"
            dkey = f'{atom[5]}:{atom[6]}:{atom[4].strip()}'
        else:
            if atom[2].strip()=="N" or atom[2].strip()=="H":
                FG = "Ntrm"
                dkey = f'{atom[5]}:{int(atom[6])-1}:MC'
            elif atom[2].strip()=="C" or atom[2].strip()=="O":
                FG = "Ctrm"
                dkey = f'{atom[5]}:{int(atom[6])}:MC'
            elif atom[2].strip()=="CA" or atom[2].strip()=="HA" or atom[2].strip()=="HA2" or atom[2].strip()=="HA3":
                FG = "alph"
                dkey = f'{atom[5]}:{int(atom[6])}:{atom[4].strip()}_SC'
            else:
                FG = "side"
                dkey = f'{atom[5]}:{int(atom[6])}:{atom[4].strip()}_SC'

        allinfo.append([atom[5], int(atom[6]), atom[4], FG, atom[2].strip(), round(distcom,3), round(distavg,3), round(distclosest,3)])
        if key not in allats.keys():
            allats[key] = {'name': atom[4], 'all': []}
        allats[key]['all'].append(atname)
        if FG not in allats[key].keys():
            allats[key][FG] = [atname]
        else:
            allats[key][FG].append(atname)
        
        if dkey not in alldclos.keys():
            alldclos[dkey] = {atom[2].strip(): distclosest}
        else:
            alldclos[dkey][atom[2].strip()] = distclosest
    
    return allinfo,allats,alldclos
    

def apply_cutoff(allinfo,dtype,cut):
    inradius = []
    selats = {}
    for line in allinfo:
        #line = ['chID','resI','resN','FGrp','atom','dCOM','dAvg','dClo']
        if (dtype=="mass" and line[5] <= cut) or (dtype=="avg" and line[6] <= cut) or (dtype=="closest" and line[7] <= cut):
            key = (line[0],line[1])
            atname = f'{line[4]:<4}'
            FG = line[3]
            inradius.append(line)
            if key not in selats.keys():
                selats[key] = {'name': line[2], 'all': []}
            selats[key]['all'].append(atname)
            if FG not in selats[key].keys():
                selats[key][FG] = [atname]
            else:
                selats[key][FG].append(atname)
    atomdf = pd.DataFrame(inradius, columns=['chID','resI','resN','FGrp','atom','dCOM','dAvg','dClo'])
    return atomdf, selats


def res_grouping(atomdf,seldist,selats):
    resdf = atomdf.sort_values(seldist)
    resdf = resdf.drop_duplicates(subset=['chID','resI','resN'], keep='first')
    removedist = [d for d in ['dCOM','dAvg','dClo'] if d != seldist]
    resdf = resdf.drop(columns=['FGrp',removedist[0],removedist[1]])
    #resdf = resdf.assign(atomlist = [' '.join(selats[(x[0],x[1])]['all']) for x in np.array(resdf)])
    return resdf


def fg_grouping(atomdf,seldist,selats):
    FGdf = atomdf.sort_values(seldist)
    FGdf = FGdf.drop_duplicates(subset=['chID','resI','resN','FGrp'], keep='first')
    removedist = [d for d in ['dCOM','dAvg','dClo'] if d != seldist]
    FGdf = FGdf.drop(columns=removedist)
    peps = []
    for i in range(FGdf.shape[0]):
        # only keep one half of each peptide bond
        if FGdf.iloc[i,3]=='Ctrm':
            if [FGdf.iloc[i,0],FGdf.iloc[i,1]] in peps:
                FGdf.iloc[i,3] = 'DUPL'
            else:
                FGdf.iloc[i,2] = 'PEP'
                peps.append([FGdf.iloc[i,0],FGdf.iloc[i,1]])
        elif FGdf.iloc[i,3]=='Ntrm':
            if [FGdf.iloc[i,0],FGdf.iloc[i,1]-1] in peps:
                FGdf.iloc[i,3] = 'DUPL'
            else:
                FGdf.iloc[i,2] = 'PEP'
                peps.append([FGdf.iloc[i,0],FGdf.iloc[i,1]-1])
        elif FGdf.iloc[i,3]=='alph':
            # CA/HA should be considered part of side chain if it's going to be added or res is glycine
            if 'side' in selats[(FGdf.iloc[i,0],FGdf.iloc[i,1])].keys() or selats[(FGdf.iloc[i,0],FGdf.iloc[i,1])]['name'] == 'GLY':
                FGdf.iloc[i,3] = 'side'
                try:
                    for atom in selats[(FGdf.iloc[i,0],FGdf.iloc[i,1])]['alph']:
                        selats[(FGdf.iloc[i,0],FGdf.iloc[i,1])]['side'].append(atom)
                except:
                    selats[(FGdf.iloc[i,0],FGdf.iloc[i,1])]['side']=selats[(FGdf.iloc[i,0],FGdf.iloc[i,1])]['alph']
            # otherwise ignore alphas - will be added if MCs are added and if there's nothing else in radius prob too far out to be important?
            else:
                FGdf.iloc[i,3] = 'DUPL'
    #remove any remaining duplicates
    FGdf = FGdf[FGdf.FGrp != 'DUPL']
    FGdf = FGdf.drop_duplicates(subset=['chID','resI','resN','FGrp'], keep='first')
    return FGdf


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Distance based selection scheme with both residue and functional group partitioning')
    parser.add_argument('-pdb', dest='pdbf', default=None, help='pdb file to use')
    parser.add_argument('-s','-seed', dest='seed', default=None, help='seed, examples: A:300,A:301,A:302')
    parser.add_argument('-satom', dest='seedatom', default=None, help='atoms to use as center, examples: A:300:CA')
    parser.add_argument('-type', dest='dtypef', default=None, help='"closest" or "mass" or "avg"')
    parser.add_argument('-max', dest='coff', default=5, help='cutoff/max distance from ligand, default 5A')
    parser.add_argument('-noH','-noh', dest='noH', action='store_true', help='ignore hydrogens')
    parser.add_argument('-store', dest='storeinfo', action='store_true', help='store distances and pdb/seed/noH info in all_dists.pkl')
    parser.add_argument('-recut', dest='recut', action='store_true', help='apply new cutoff to atom distances read from all_dists.pkl')    
    parser.add_argument('-byres', dest='byres', action='store_true', help='use residue grouping instead of FG grouping')

    args = parser.parse_args()
    if args.seed is not None and args.seedatom is not None:
        print("Please only use one of seed and satom to specify distance center")
        sys.exit(1)
    elif args.seed is not None:
        cres = args.seed
    elif args.seedatom is not None:
        cres = args.seedatom

    pdbf = args.pdbf
    pdb, res_info, tot_charge = read_pdb(pdbf)
    cut = float(args.coff)
    noH = args.noH
 
    dtype = args.dtypef.lower()
    if dtype == "mass":
        seldist = 'dCOM'
    elif dtype == "avg":
        seldist = 'dAvg'
    elif dtype == "closest":
        seldist = 'dClo'
 
    if args.recut:
        with open('all_dists.pkl','rb') as fp:
            loaded = pickle.load(fp)
        pdbf = loaded['pdbf']
        cres = loaded['seed']
        noH = loaded['noH']
        allinfo = loaded['allinfo']
        print(f'Distance calculation settings used in alldists.pkl: pdb = {pdbf}, seed = {cres}, noH = {noH}')
    else:
        seeddict, com, avg = getseeddist(cres,pdb,noH)    
        allinfo,allats,alldclos = getatomdists(pdb,com,avg,seeddict,noH)
        allatomdf = pd.DataFrame(allinfo, columns=['chID','resI','resN','FGrp','atom','dCOM','dAvg','dClo'])
        with open('all_dists_table.dat','w') as f:
            f.write(f'# Distances from seed {cres} in pdb {pdbf}, noH = {noH}\n')
        allatomdf.to_csv('all_dists_table.dat', sep='\t', index=False, mode='a')
        if args.storeinfo:
            savedists = {'pdbf': pdbf, 'seed': cres, 'noH': noH, 'allinfo': allinfo}
            with open('all_dists.pkl', 'wb') as fp:
                pickle.dump(savedists, fp)
        
    atomdf,selats = apply_cutoff(allinfo,dtype,cut)

    if args.byres:
        resdf = res_grouping(atomdf,seldist,selats)
        resdf.to_csv('sorted_residues_%s_%.1f_A.dat'%(seldist,cut), sep='\t', index=False)   
        f1 = open('res_atoms.dat','w')
        f1.write(f'# selected/ranked residues by distance: pdb {pdbf}; seed {cres}; disttype {dtype}; cutoff {cut}; ignore H {str(noH)}\n')
        distcol = resdf.columns.get_loc(seldist)
        for i in range(resdf.shape[0]):
            all_ats = ' '.join(selats[(FGdf.iloc[i,0],FGdf.iloc[i,1])]['all'])
            f1.write(f'{resdf.iloc[i,0]:<4} {resdf.iloc[i,1]:<8} {resdf.iloc[i,distcol]:<8} {all_ats}\n')
        f1.close()        
    else:
        FGdf = fg_grouping(atomdf,seldist,selats)
        FGdf.to_csv('sorted_FGs_%s_%.1f_A.dat'%(seldist,cut), sep='\t', index=False) 
        f2 = open('res_atoms.dat','w')
        f2.write(f'# selected/ranked functional groups by distance: pdb {pdbf}; seed {cres}; disttype {dtype}; cutoff {cut}; ignore H {str(noH)}\n')
        distcol = FGdf.columns.get_loc(seldist)
        for i in range(FGdf.shape[0]):
            fg_ats = ' '.join(selats[(FGdf.iloc[i,0],FGdf.iloc[i,1])][FGdf.iloc[i,3]])
            f2.write(f'{FGdf.iloc[i,0]:<4} {FGdf.iloc[i,1]:<8} {FGdf.iloc[i,distcol]:<8} {fg_ats}\n')
        f2.close()

