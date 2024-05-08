#! /usr/bin/env python
import argparse
import sys
import os

"""
This is a program written by Tejaskumar Suhgaia in Deyonker research group
at university of memphis.
"""
residuedict = {}
lastreesidue = []
used_chains = set()
used_newids = set()
used_chainID= set()
lastOldNewchainID = {}

if __name__ == '__main__':
    """ Usage: res_extract.py -p template.pdb """
    parser = argparse.ArgumentParser(description='renumber residue_id for hydrogen water from 0 to sorresponding oxygen from pdb file')
    parser.add_argument('-p', dest='pdbf', default='83REAC_QChem1.pdb', help='template pdb file')
    args = parser.parse_args()
    pdbf = args.pdbf
    output_file_path = 'waterid_ocrrected_'+pdbf  
    
    def resid1(line1):
        res2, chain1, resid1 = line1[17:20].strip(), line1[20:22].strip(), line1[22:26].strip()
        return resid1

    def resname1(line1):
        res2, chain1, resid1 = line1[17:20].strip(), line1[20:22].strip(), line1[22:26].strip()
        return res2            
        
    with open(f"{pdbf}", "r") as readstart:
        lines = readstart.readlines()
        
        for i, line in enumerate(lines):
            if "ATOM" in line or "HETATM" in line:
                if resname1(line) == "HOH": 
                    if resid1(line) == "0":
                        if i + 1 < len(lines):
                            next_line = lines[i + 1]
                            previous_line = lines[i - 1]
                            if resid1(next_line) != "0" and resname1(next_line) == "HOH":                          
                                print(f"{line[:22]}{resid1(next_line):>4}{line[26:]}", end="")
                            elif resid1(next_line) == "0" and resname1(next_line) == "HOH":
                                print(f"{line[:22]}{resid1(previous_line):>4}{line[26:]}", end="")
                            elif resid1(previous_line) != "0" and resname1(previous_line) == "HOH":
                                print(f"{line[:22]}{resid1(previous_line):>4}{line[26:]}", end="")
                    elif resid1(line) != "0":
                        print(f"{line[:22]}{resid1(line):>4}{line[26:]}", end="")
                elif resname1(line) != "HOH":
                    print(line)
            else:
                print(line,end="")