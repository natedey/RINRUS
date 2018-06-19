#!/usr/bin/env python

import sys, os, re
from PDB import *
from numpy import *
from glob import glob

def compare_2(idx1,idx2,pdb1,pdb2):
    if len(pdb1) == len(pdb2):
        heavy = 0
        count = 0
        for i in range(len(pdb1)):
            if 'H' not in pdb1[i][2]:
                heavy += 1
                if pdb1[i][2] == pdb2[i][2]:
                    count += 1
        if heavy == count:
            print "final_%d.pdb and final_%d.pdb probably are the same, please double check!"%(idx1,idx2)


if __name__ == '__main__':
    files = []
    pdb_dic = {}
    for file in glob('%s/final_*.pdb'%sys.argv[1]):
        count = file.split('/')[-1][:-4].split('_')[-1]
        files.append(int(count))
        pdb_dic[int(count)] = read_pdb(file)
    
    for i in sorted(files):
        print i
#        if i-1 in files:
#            compare_2(i-1,i,pdb_dic[i-1],pdb_dic[i])
        if i+1 in files:
            compare_2(i,i+1,pdb_dic[i],pdb_dic[i+1])
