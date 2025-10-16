#!/usr/bin/env python3
"""
This is a program created by Dr. Dominique Wappett and the DeYonker research group
at The University of Memphis.
Prototype created 2025-07-07
"""

import os, sys, re
import argparse
import subprocess
from subprocess import Popen, PIPE,STDOUT
from pathlib import Path
import pandas as pd
import pickle
import glob
#biopython funcs for getting nearest waters
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB import NeighborSearch
#RINRUS functions
from probe2rins import *
from arpeggio2rins import *
from dist_rank import *

pd.options.mode.copy_on_write = True

### contact types from probe and arpeggio ###
probe_types = ['wc','cc','so','bo','hb']
arp_types = ['a_cl','a_cov','a_vdwcl','a_vdw','a_prx','a_hb','a_whb','a_hal','a_ion','a_met','a_arom','a_hp','a_CO','a_polar','a_wpol'] 

### probe analysis. runs probe first if no [name].probe file found ###
def run_and_process_probe(path_to_scripts,pdb,fname,seedlist,batch_rinfo,atom_master,MD_wats):
    f = pdb.replace('.pdb','.probe')
    # if no probe file, run probe. otherwise use existing
    if not os.path.isfile(f):
        print(f'Running probe on {pdb}')
        args = [path_to_scripts+'probe -unformated -MC -self "all" -Quiet '+ pdb +' > '+ f]
        out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    # do probe analysis stuff
    atomdict, seedcontact, seedatoms = probe_atom(f,seedlist)
    allcont, allats = probe_fg(atomdict, seedcontact, seedatoms)
    
    # frame analysis/res_atoms
    f = f.rsplit('/',1)
    rinrus_probe_outputs(allcont,allats,seedcontact,f[-1],seedlist)
    f[-1] = f[-1].replace('.probe','')
    if len(f) > 1:
        os.rename('res_atoms.dat',f'{f[0]}/{fname}.res_atoms..dat')
        os.rename('FG_probe_counts.dat',f'{f[0]}/{fname}.FG_probe_counts.dat')
    else:
        os.rename('res_atoms.dat',f'{fname}.res_atoms.dat')
        os.rename('FG_probe_counts.dat',f'{fname}.FG_probe_counts.dat')

    # add to master dict by fragment/filename
    for key in allcont.keys():
        if key in seedlist:
            continue
        renamed_probecont = {f'p_{c}': allcont[key][c] for c in probe_types}
        if key not in batch_rinfo.keys():
            batch_rinfo[key] = {fname: renamed_probecont}
        elif fname not in batch_rinfo[key].keys():
            batch_rinfo[key][fname] = renamed_probecont
        else:
            batch_rinfo[key][fname].update(renamed_probecont)
        batch_rinfo[key][fname]['p_tot'] = sum([allcont[key][c] for c in probe_types])
        if key.endswith('WAT'):
            if fname not in MD_wats.keys():
                MD_wats[fname] = {key: {'probe': sum([allcont[key][c] for c in probe_types])}}
            elif key not in MD_wats[fname].keys():
                MD_wats[fname][key] = {'probe': sum([allcont[key][c] for c in probe_types])}
            else:
                MD_wats[fname][key]['probe'] = sum([allcont[key][c] for c in probe_types])
    for key in allats.keys():
        if key not in atom_master.keys():
            atom_master[key] = allats[key]
        else:
            for atom in allats[key]:
                if atom not in atom_master[key]:
                    atom_master[key].append(atom)
    return batch_rinfo, atom_master, MD_wats

### arpeggio analysis. runs arpeggio first if no [name].contacts file found ###
def run_and_process_arpeggio(path_to_scripts,pdb,fname,pdb_res_name,seedlist,batch_rinfo,atom_master,MD_wats):
    f = str(pdb).replace('pdb','contacts')
    if not os.path.isfile(f):
        print(f'Running arpeggio on {pdb}')
        path = os.path.expanduser(path_to_scripts+'arpeggio/arpeggio.py')
        arg = [sys.executable,path,str(pdb)]
        result = subprocess.run(arg,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    if not os.path.isfile(f):
        print('Error running arpeggio!')
        if 'WARNING: Problems reading a PDB file' in stdout:
            print(f'Open Babel could not read PDB {pdb}. Check formatting and try again')
        sys.exit()
    # do arpeggio stuff
    # arpeggio_atom third arg is noprox. setting to true and removing proximal
    atomdict, seeddict = arpeggio_atom(f,seedlist,True)
    allcont, allats = arpeggio_fg(atomdict,pdb_res_name,seeddict)

    # frame analysis/res_atoms
    f = f.rsplit('/',1)
    rinrus_arpeggio_outputs(allcont,allats,seeddict,True,f[-1],seedlist)
    if len(f) > 1:
        os.rename('res_atoms.dat',f'{f[0]}/{fname}.res_atoms.dat')
        os.rename('res_atoms_types.dat',f'{f[0]}/{fname}.res_atoms_types..dat')
        os.rename('FG_arpeggio_counts.dat',f'{f[0]}/{fname}.FG_arpeggio_counts.dat')
    else:
        os.rename('res_atoms.dat',f'{fname}.res_atoms.dat')
        os.rename('res_atoms_types.dat',f'{fname}.res_atoms_types.dat')
        os.rename('FG_arpeggio_counts.dat',f'{fname}.FG_arpeggio_counts.dat')

    # add to master dict by fragment/filename
    for key in allcont.keys():
        if key in seedlist:
            continue
        contdict = {arp_types[i]: allcont[key][i] for i in range(len(arp_types))}
        if key not in batch_rinfo.keys():
            batch_rinfo[key] = {fname: contdict}
        elif fname not in batch_rinfo[key].keys():
            batch_rinfo[key][fname] = contdict
        else:
            batch_rinfo[key][fname].update(contdict)
        batch_rinfo[key][fname]['a_tot'] = sum(allcont[key])
        if key.endswith('WAT'):
            if fname not in MD_wats.keys():
                MD_wats[fname] = {key: {'arpeggio': sum(allcont[key])}}
            elif key not in MD_wats[fname].keys():
                MD_wats[fname][key] = {'arpeggio': sum(allcont[key])}
            else:
                MD_wats[fname][key]['arpeggio'] = sum(allcont[key])
    for key in allats.keys():
        if key not in atom_master.keys():
            atom_master[key] = allats[key]
        else:
            for atom in allats[key]:
                if atom not in atom_master[key]:
                    atom_master[key].append(atom)
    return batch_rinfo, atom_master, MD_wats

### get seed-water distances for adding waters to batch RIN ###
def nearestwaters(seedlist,pdbf,fname,MD_wats):
    pdb_parser = PDBParser(QUIET=True)
    s = pdb_parser.get_structure('structure', pdbf)
    s_atoms = list(s.get_atoms())    
    # extract seed and water atoms only
    seedatoms = []
    for a in s_atoms:
        chain = a.get_parent().get_parent().get_id()
        resid = [i for i in a.get_parent().get_id() if i and not str(i).isspace()][0]
        if f'{chain}:{resid}' in seedlist:
            seedatoms.append(a)
    watatoms = [a for a in s_atoms if a.get_parent().get_resname() == 'WAT']
    # do neighbour search on this subset and extract seed-water pairs
    entity = list(seedatoms + watatoms)
    ns = NeighborSearch(entity)
    sele = ns.search_all(10)
    sele = [pair for pair in sele if len(set(pair).intersection(seedatoms)) == 1]
    if fname not in MD_wats.keys():
        MD_wats[fname] = {}
    for pair in sele:
        wat = [i for i in pair if i.get_parent().get_resname() == 'WAT'][0]
        chain = wat.get_parent().get_parent().get_id()
        resid = [i for i in wat.get_parent().get_id() if i and not str(i).isspace()][0]
        key = f'{chain}:{resid}:WAT'
        if key not in MD_wats[fname].keys():
            MD_wats[fname][key] = {'d_clo': pair[0] - pair[1]}
        elif 'd_clo' not in MD_wats[fname][key].keys():
            MD_wats[fname][key]['d_clo'] = pair[0] - pair[1]
        elif pair[0] - pair[1] < MD_wats[fname][key]['d_clo']:
            MD_wats[fname][key]['d_clo'] = pair[0] - pair[1]
    return MD_wats 

### run all processing and generate dataframe of all info ###
def batch_processing(path_to_scripts,pdbfiles,seedlist,seedinp,whichtests,mdwat):
    batch_rinfo = {}
    atom_master = {}
    MD_wats = {}

    for pdbf in pdbfiles:
        fname = pdbf.split('/')[-1].replace('.pdb','')
        pdb, res_info, tot_charge = read_pdb(pdbf)
        pdb_res_name = {}
        for line in pdb:
            key = (line[5].strip(), int(line[6]))
            pdb_res_name[key] = line[4].strip()

        ### do contact analysis ###
        if whichtests['probe']:
            batch_rinfo, atom_master, MD_wats = run_and_process_probe(path_to_scripts,pdbf,fname,seedlist,batch_rinfo,atom_master,MD_wats)
        if whichtests['arpeggio']:
            batch_rinfo, atom_master, MD_wats = run_and_process_arpeggio(path_to_scripts,pdbf,fname,pdb_res_name,seedlist,batch_rinfo,atom_master,MD_wats)
        if mdwat:
            MD_wats = nearestwaters(seedlist,pdbf,fname,MD_wats)
                

    ### turn collected results into dataframe ###
    reorg = {(key,frame): batch_rinfo[key][frame] for key in batch_rinfo.keys() for frame in batch_rinfo[key].keys()}
    rawdf = pd.DataFrame.from_dict(reorg,orient='index')
    if 'a_prx' in rawdf.columns:
        rawdf = rawdf.drop(columns='a_prx')

    ### MD water analysis ###
    watcount = {}
    for frame in MD_wats.keys():
        MD_wats[frame] = pd.DataFrame.from_dict(MD_wats[frame],orient='index')
        MD_wats[frame].columns = ['contacts' if c in ['probe','arpeggio'] else c for c in MD_wats[frame].columns]
        watcount[frame] = {'contacts': MD_wats[frame]['contacts'].count()}
    watcount = pd.DataFrame.from_dict(watcount,orient='index')
    MD_wats['counts'] = pd.concat([watcount, watcount.apply(['mean','median','max','min']).round(2)])

    return rawdf, atom_master, MD_wats


### all data from every frame collected in rawdf: now condense, sort and rank here as required ###
def process_raw_df(rawdf):
    ### remove waters
    wats = [res for res in rawdf.index.levels[0].values if res.endswith('WAT')]
    rawdf = rawdf.drop(index=wats,level=0)
    
    ### do summary statistics across set of frames, extract specific info, rank based on info
    aggdict = {col: ['count','sum','median','mean','std','var','max','min',np.ptp] for col in rawdf.columns.values if col.startswith('p_') or col.startswith('a_')} #np.ptp for range
    resdf = rawdf.groupby(level=0).agg(aggdict)

    cols = {'p_tot': 'probe', 'a_tot': 'arpeggio'}
    selcols = [c for c in cols.keys() if c in resdf.columns]
    summdf = resdf[selcols].round(2)
    summdf.columns = summdf.columns.remove_unused_levels()
    selcols = [cols[c] for c in selcols]
    summdf.columns = summdf.columns.set_levels([selcols, ['nf' if c == 'count' else 'range' if c == 'ptp' else c for c in summdf.columns.levels[1]]])
    summdf.insert(1,(selcols[0],'prop_f'),round(summdf[(selcols[0],'nf')]/len(pdbfiles),2))
    summdf.insert(summdf.columns.get_loc((selcols[0],'mean'))+1,(selcols[0],'wt_mean'),round(summdf[(selcols[0],'prop_f')]*summdf[(selcols[0],'mean')],2))
    ### sort
    summdf = summdf.sort_values(by=[(selcols[0],'nf'),(selcols[0],'mean')],ascending=False)
    rankcol = summdf[(selcols[0],'prop_f')].copy()

    return rawdf, resdf, summdf, rankcol


### make res atoms file ###
def write_batch_res_atoms(seedinp,seedlist,atom_master,df,resatsource,batchfile):
    f = open(batchfile,'w')
    f.write(f'# master res_atoms file from batch processing info. seed: {seedinp}\n# processed directory: {resatsource}\n')
    for s in seedlist:
        f.write(f'{s.split(":")[0]:<4} {s.split(":")[1]:<8} seed      ')
        for a in atom_master[s]:
            f.write(f' {a:<4}')
        f.write('\n')
    for g in df.index.values:
        # if MC, make sure atoms correspond to correct res id
        if g.split(':')[-1] == 'MC':
            idx = int(g.split(":")[1])
            if set(atom_master[g]).issubset({'C','O'}):
                f.write(f'{g.split(":")[0]:<4} {idx:<8} pf={df.loc[g]:<6} ')
                for a in atom_master[g]:
                    f.write(f' {a:<4}')
                f.write('\n')
            elif set(atom_master[g]).issubset({'N','H'}):
                f.write(f'{g.split(":")[0]:<4} {idx+1:<8} pf={df.loc[g]:<6} ')
                for a in atom_master[g]:
                    f.write(f' {a:<4}')
                f.write('\n')
            else:
                f.write(f'{g.split(":")[0]:<4} {idx:<8} pf={df.loc[g]:<6} ')
                for a in [a for a in atom_master[g] if a in ['C','O']]:
                    f.write(f' {a:<4}')
                f.write('\n')
                f.write(f'{g.split(":")[0]:<4} {idx+1:<8} pf={df.loc[g]:<6} ')
                for a in [a for a in atom_master[g] if a in ['N','H']]:
                    f.write(f' {a:<4}')
                f.write('\n')
        else:
            f.write(f'{g.split(":")[0]:<4} {g.split(":")[1]:<8} pf={df.loc[g]:<6} ')
            for a in atom_master[g]:
                f.write(f' {a:<4}')
            f.write('\n')
    f.close()

def res_atoms_add_waters(seedlist,MD_wats,batchfile):
    with open(batchfile) as bf:
        batchres = bf.readlines()[2:]
    nwat = int(MD_wats['counts'].loc['median','contacts'])
    for frame in [key for key in MD_wats.keys() if key != 'counts']:
        MD_wats[frame] = MD_wats[frame].sort_values(by=['contacts','d_clo'],ascending=[False,True])
        sel_wat = list(MD_wats[frame].index.values[0:nwat])
        #filename = f'res_atoms_batch.{frame}.dat'
        filename = batchfile.replace('.dat',f'.{frame}.dat')
        f = open(filename,'w')
        f.write(f'# batch res_atoms file plus first {nwat} waters (by contacts then closest atom-atom distance) for frame {frame}\n')
        f.write(''.join(batchres))
        for wat in sel_wat:
            if np.isnan(MD_wats[frame].loc[wat,'contacts']):
                val = np.format_float_positional(MD_wats[frame].loc[wat,'d_clo'], precision=2)
                f.write(f'{wat.split(":")[0]:<4} {wat.split(":")[1]:<8}  d={val:<6}  O    H1   H2   \n')
            else:
                val = MD_wats[frame].loc[wat,'contacts']
                f.write(f'{wat.split(":")[0]:<4} {wat.split(":")[1]:<8}  c={val:<6}  O    H1   H2   \n')
        f.close()



if __name__ == '__main__':
    """ 
    Usage: batch_rin.py -s [seed] -dir [dir] -type [probe/arpeggio]
    Decisions hard-coded:
    - arpeggio proximal ignored
    - res_atoms_batch.dat has no waters since they differ between frames (median no. waters can be added with -addwat flag)
    - sorting of summary stats/res_atoms_batch by number of frames then mean contact count
    """
    parser = argparse.ArgumentParser(description='batch contact analysis of MD frames')
    parser.add_argument('-s','-seed',dest='seed',help='seed ch:ID[,ch:ID,ch:ID]')
    parser.add_argument('-dir',dest='dir',default='./',help='directory to process')
    parser.add_argument('-type',dest='rintype',default='probe',help='type of processing: probe/arpeggio')
    parser.add_argument('-addwat',dest='mdwat',action='store_true',help='add waters to batch rin for each frame')

    args = parser.parse_args()
    seedinp = args.seed
    seedlist = seedinp.split(',')
    startdir = os.getcwd()
    print(f'Processing all PDB files in directory: {args.dir}')
    if args.dir.endswith('/'):
        args.dir = args.dir[0:-1]
    os.chdir(f'{args.dir}')
    pdbfiles = glob.glob(f'*.pdb')
    pdbfiles.sort()
    whichtests = {'probe': False, 'arpeggio': False}
    whichtests[args.rintype] = True

    # location of this file is used to find probe/arpeggio binaries if they need to be run
    path_to_scripts = str(Path(__file__).parent)+'/'

    ### run all processing ###
    rawdf, atom_master, MD_wats = batch_processing(path_to_scripts,pdbfiles,seedlist,seedinp,whichtests,args.mdwat)
    rawdf, resdf, summdf, rankcol = process_raw_df(rawdf)
    
    ### print the summary dataframe to the console and save as text file
    print(summdf)
    summdf.to_string(f'batch_stats.dat')
    ### also print numbers of waters
    MD_wats['counts'].to_string(buf='waters_per_frame.dat',header=['N_waters_with_contacts'])
    
    ### save the raw dataframe to a pkl file for merging directory batches
    with open('batch_rawdf.pkl', 'wb') as fp:
        pickle.dump(rawdf, fp)
    with open(f'batch_atominfo.pkl', 'wb') as fp:
        pickle.dump(atom_master, fp)

    ### write res_atoms based on the rankings of summdf
    write_batch_res_atoms(seedinp,seedlist,atom_master,rankcol,args.dir,'batch_res_atoms.dat')
    if args.mdwat:
        res_atoms_add_waters(seedlist,MD_wats,'batch_res_atoms.dat')
        del MD_wats['counts']
        with open(f'batch_waterinfo.pkl', 'wb') as fp:
            pickle.dump(MD_wats, fp)
    
    os.chdir(startdir)
