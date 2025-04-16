#! /usr/bin/env python
"""
This is a program created by Dr. Tejaskumar Suhagia and the DeYonker research group
at The University of Memphis.
Date created 05.08.2024
"""

import argparse
import sys

residuedict = {}
lastreesidue = []
used_chains = set()
used_newids = set()
used_chainID= set()
lastOldNewchainID = {}

if __name__ == '__main__':
    """ Usage: duplicate_atomfinder.py -p template.pdb """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p', dest='pdbf', default='template.pdb', help='pdb file')
    args = parser.parse_args()
    pdbf = args.pdbf
    output_file_path = 'ChainCorrected_'+pdbf  
    def reskey1(line1):
        res2, chain1, resid1 = line1[17:20].strip(), line1[20:22].strip(), line1[22:26].strip()
        if chain1.strip() == "":
            chain1 = "A"
        return res2 + resid1       
    def reschainkey1(line1):
        res2, chain1, resid1 = line1[17:20].strip(), line1[20:22].strip(), line1[22:26].strip()
        if chain1.strip() == "":
            chain1 = "A"
        return res2 + chain1 + resid1      
    def resid1(line1):
        res2, chain1, resid1 = line1[17:20].strip(), line1[20:22].strip(), line1[22:26].strip()
        if chain1.strip() == "":
            chain1 = "A"
        return resid1      
    def reschain1(line1):
        res2, chain1, resid1 = line1[17:20].strip(), line1[20:22].strip(), line1[22:26].strip()
        if chain1.strip() == "":
            chain1 = "A"
        return chain1
    def reschainid1(line1):
        res2, chain1, resid1 = line1[17:20].strip(), line1[20:22].strip(), line1[22:26].strip()
        if chain1.strip() == "":
            chain1 = "A"
        return chain1 + resid1 
        
    with open(f"{pdbf}", "r") as readstart:
        lines = readstart.readlines()
        
        for line in lines[:]:
            if "ATOM" in line or "HETATM" in line:
                if len(lastreesidue) == 0:
                    lastreesidue.append(reschainid1(line))
                    used_chains.add(reschain1(line))
                    residuedict[reschainid1(line)] = +1
                elif len(lastreesidue) == 1:    
                    reskey = reskey1(line)
                    reschain = reschain1(line)
                    resid=resid1(line)
                    res2, chain2, resid2 = line[17:20].strip(), line[20:22].strip(), line[22:26].strip()
                    used_chains.add(reschain)
                    if chain2.strip() == "":
                        chain2 = "A"
                    reschainid = chain2+resid2
                    if reschainid in lastreesidue:
                        pass
                    elif reschainid is not lastreesidue:
                        if reschainid in residuedict:
                            residuedict[reschainid1(line)] = residuedict[reschainid1(line)]+1
                            lastreesidue.pop()
                            lastreesidue.append(reschainid)
                        elif reschainid not in residuedict:
                            residuedict[reschainid1(line)] = 1
                            lastreesidue.pop()
                            lastreesidue.append(reschainid)
                                                                    
duplicate_chainid=0
duplicates=0
for key, value in residuedict.items():
    if value > 1:
        print(key, value, "duplicate found")
        duplicates=duplicates+value
        duplicate_chainid=duplicate_chainid+1
print("total duplicate chainid found =", duplicate_chainid)
print("total duplicates found =", duplicates)
