#!/usr/bin/env python3
"""
This is a program created by Dr. Dominique Wappett and the DeYonker research group
at The University of Memphis.
Prototype created 2025-07-08
"""

import os, sys, re
import argparse
import pandas as pd
import pickle
import glob
#RINRUS functions
from MD_batch_rin import *

pd.options.mode.copy_on_write = True
pd.set_option('display.max_rows', 500)

if __name__ == '__main__':
    """ 
    Usage: merge_batches.py -dirs [dirs]
    """
    parser = argparse.ArgumentParser(description='merge data from batch processed dirs')
    parser.add_argument('-s','-seed',dest='seed',help='seed ch:ID[,ch:ID,ch:ID]')
    parser.add_argument('-dirs',dest='dirs',nargs='*',help='directories to merge data from')
    parser.add_argument('-addwat',dest='mdwat',action='store_true',help='add waters to batch rin for each frame')

    args = parser.parse_args()
    seedinp = args.seed
    seedlist = seedinp.split(',')
    dirlist = [d[0:-1] if d.endswith('/') else d for d in args.dirs]

    atom_master = {}
    MD_wats = {}
    Ntot = 0
    watperframe = {}
    for d in dirlist:
        # count no. frames in dir
        Ntot += len(glob.glob(f'{d}/*.pdb'))
        # load info from dir
        df = pd.read_pickle(f'{d}/batch_rawdf.pkl')
        with open(f'{d}/batch_atominfo.pkl','rb') as fp:
            batch_atoms = pickle.load(fp)
        # accumulate
        if dirlist.index(d) == 0:
            rawdf = df
        else:
            rawdf = pd.concat([rawdf,df])
        for key in batch_atoms.keys():
            if key not in atom_master.keys():
                atom_master[key] = batch_atoms[key]
            else:
                for atom in batch_atoms[key]:
                    if atom not in atom_master[key]:
                        atom_master[key].append(atom)
        if args.mdwat:
            with open(f'{d}/batch_waterinfo.pkl','rb') as fp:
                batch_wat = pickle.load(fp)
            MD_wats.update(batch_wat)
        
        ### get n waters per frame info even if waters not included
        wats = open(f'{d}/waters_per_frame.dat','r').readlines()[1:-4]
        for w in wats:
            w = w.strip().split()
            watperframe[w[0]] = float(w[1])


    ### all data from probe/arpeggio/distance collected in rawdf ###
    ### condense, sort and rank below this point as required     ###
    
    ### do summary statistics across set of frames, extract specific info, rank based on info
    aggdict = {col: ['count','sum','median','mean','std','var','max','min',np.ptp] for col in rawdf.columns.values if col.startswith('p_') or col.startswith('a_')} #np.ptp for range
    resdf = rawdf.groupby(level=0).agg(aggdict)

    roundto = min([4,len(str(Ntot))])
    cols = {'p_tot': 'probe', 'a_tot': 'arpeggio'}
    selcols = [c for c in cols.keys() if c in resdf.columns]
    summdf = resdf[selcols].round(roundto)
    summdf.columns = summdf.columns.remove_unused_levels()
    selcols = [cols[c] for c in selcols]
    summdf.columns = summdf.columns.set_levels([selcols, ['nf' if c == 'count' else 'range' if c == 'ptp' else c for c in summdf.columns.levels[1]]])
    summdf.insert(1,(selcols[0],'prop_f'),round(summdf[(selcols[0],'nf')]/Ntot,roundto))
    summdf.insert(summdf.columns.get_loc((selcols[0],'mean'))+1,(selcols[0],'wt_mean'),round(summdf[(selcols[0],'prop_f')]*summdf[(selcols[0],'mean')],roundto))
    # sort
    summdf = summdf.sort_values(by=[(selcols[0],'nf'),(selcols[0],'mean')],ascending=False)
    rankcol = summdf[(selcols[0],'prop_f')].copy()


    ### print the summary dataframe to the console and save as text file
    print(summdf)
    summdf.to_string(f'batch_stats.dat')

    ### save the raw dataframe to a pkl file to be reopened and analysed if needed
    with open('batch_rawdf.pkl', 'wb') as fp:
        pickle.dump(rawdf, fp)

    ### write a res_atoms file based on the rankings of summdf
    batchfile = 'batch_res_atoms.dat'
    write_batch_res_atoms(seedinp,seedlist,atom_master,rankcol,f'directories: {dirlist}','batch_res_atoms.dat')
    if args.mdwat:
        watcount = {}
        for frame in MD_wats.keys():
            MD_wats[frame].columns = ['contacts' if c in ['probe','arpeggio'] else c for c in MD_wats[frame].columns]
            watcount[frame] = {'contacts': MD_wats[frame]['contacts'].count()}
        watcount = pd.DataFrame.from_dict(watcount,orient='index')
        MD_wats['counts'] = pd.concat([watcount, watcount.apply(['mean','median','max','min']).round(2)])
        res_atoms_add_waters(seedlist,MD_wats,batchfile)
        MD_wats['counts'].to_string(buf='waters_per_frame.dat',header=['N_waters_with_contacts'])
    else:
        #still write waters per frame even without full water analysis
        watperframe = pd.DataFrame.from_dict(watperframe,orient='index',columns=['N_waters_with_contacts'])
        watperframe = pd.concat([watperframe, watperframe.apply(['mean','median','max','min']).round(2)])
        watperframe.to_string(buf='waters_per_frame.dat',header=['N_waters_with_contacts'])

