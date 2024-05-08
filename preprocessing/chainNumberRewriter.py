#! /usr/bin/env python
import argparse
import sys

"""
This is a program written by Tejaskumar Suhgaia in Deyonker research group
at university of memphis.
"""

################################################################################################
## This script will check the chain ID of the each residue in the pdb file
## 1. if it is blank, it will assign it to "A" chain.
## 2. if it is not blank, it will assign it to the same chain.
## 3. if it is not blank, but it is not the same as the previous residue,
##    it will assign it to the next chain ID alphabet.
## 4. AMBER(with no chain number) and RCSB protein databank pdb's with no chain ID assigned or 
##    duplicate chains for same residue will be corrected (eg. 1uhe.pdb with multiple models 
##    containing same resID and chainID).
## 5. it will create a new pdb file named "ChainCorrected_'givenPDBNAME'".pdb
## 6. residueID with hex number(vmd) will be assigned to new chainID with new 4 digit number.
###############################################################################################


def isNumber(s):    
    for i in range(len(s)):
        if s[i].isdigit() != True:
            return False
    return True
def is_hex(s):
    try:
        int(s, 16)
        return True
    except ValueError:
        return False
def is_hex_or_int(value):
    try:
        if isNumber(value) == True:
            return "int"
        elif is_hex(value) == True:
            return "hex"        
    except ValueError:
        return "neither"
def convert_hex_to_four_digit(hex_value):
    return int(str(int(hex_value, 16)).zfill(4))
def assign_new_chain(used_chains1,used_chainID1):
    new_chain = chr(ord("A") + 5)
    while new_chain+"1" in used_chainID1:
        new_chain = chr(ord(new_chain) + 1)
    new_id=1
    while new_chain+str(new_id) in used_chainID1:
        new_id+=1
    used_chainID1.add(new_chain+str(new_id))
    used_chains1.add(new_chain)
    return new_chain,new_id
def generate_chain_id(uniquiechains1,used_chains1, used_chainID1):
    alphabet = 'FGHIJKLMNOPQRSTUVWXYZ'
    
    for letter in alphabet:
        for i in range(1, 10000):
            chain_id = f'{letter}{str(i).rjust(4)}'
            
            if chain_id not in used_chainID1 and letter not in uniquiechains1:
                used_chainID.add(chain_id)
                used_chains.add(letter)
                return letter,i,chain_id
def find_unique_chains(pdbf):
    unique_chains = set()

    with open(pdbf, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21]
                if chain_id.isalpha():
                    unique_chains.add(chain_id)

    return list(unique_chains)

residuedict = {}
lastreesidue = []
used_chains = set()
used_newids = set()
used_chainID= set()
lastOldNewchainID = {}

if __name__ == '__main__':
    """ Usage: chainNumberRewriter.py -p template.pdb """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p', dest='pdbf', default='template.pdb', help='pdb file')
    args = parser.parse_args()
    pdbf = args.pdbf
    uniquiechains = find_unique_chains(pdbf)
    output_file_path = 'ChainCorrected_'+pdbf  
    with open(output_file_path, 'w') as output_file:
        sys.stdout = output_file
        
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
        
        with open(f"{pdbf}", "r") as readstart:
            lines = readstart.readlines()
            
            for line in lines[:]:
                if "ATOM" in line or "HETATM" in line:
                    if len(lastreesidue) == 0:
                        lastreesidue.append(reschainkey1(line))
                        used_chains.add(reschain1(line))
                        residuedict[reskey1(line)] = reschain1(line) 
                        if is_hex_or_int(resid1(line))== "int":
                            print(line[:20] + " " + (residuedict[reskey1(line)]) + line[22:], end="")
                        elif is_hex_or_int(resid1(line))== "hex":
                            resid=convert_hex_to_four_digit(resid1(line))
                            if resid>9999:
                                reschain,resid,reschainid=generate_chain_id(uniquiechains,used_chains,used_chainID)
                                lastOldNewchainID[reschainkey1(line)]=reschain.rjust(2)+str(resid).rjust(4)
                                print(line[:20] +reschain.rjust(2) +str(resid).rjust(4)+  line[26:], end="")

                    elif len(lastreesidue) == 1:    
                        reskey = reskey1(line)
                        reschain = reschain1(line)
                        resid=resid1(line)
                        res2, chain2, resid2 = line[17:20].strip(), line[20:22].strip(), line[22:26].strip()
                        used_chains.add(reschain)
                        if chain2.strip() == "":
                            chain2 = "A"
                        reschainkey = res2+chain2+resid2
                        if reschainkey in lastreesidue:
                            if is_hex_or_int(resid1(line))== "int":
                                print(line[:20] + " " + (residuedict[reskey1(line)]) + line[22:], end="")
                            elif is_hex_or_int(resid1(line))== "hex":
                                resid=convert_hex_to_four_digit(resid1(line))
                                if resid>9999:
                                    reschainidkey=lastOldNewchainID.get(reschainkey)
                                    print(line[:20] +reschainidkey+  line[26:], end="")
                        elif reschainkey not in lastreesidue:
                            if is_hex_or_int(resid1(line))== "int":
                                if reskey not in residuedict:
                                    lastreesidue.pop()
                                    lastreesidue.append(reschainkey)
                                    residuedict[reskey] = reschain
                                    print(line[:20] + " " + (residuedict[reskey]) + line[22:], end="")
                                elif reskey in residuedict:
                                    res2, chain2, resid2 = line[17:20].strip(), line[20:22].strip(), line[21:26].strip()
                                    if chain2.strip() == "":
                                        chain2 = "A"
                                    newchain = (chr(ord(residuedict[reskey]) + 1))
                                    reschainkeynew = res2 + newchain + resid2
                                    lastreesidue.pop()
                                    lastreesidue.append(reschainkey)
                                    residuedict[reskey] = newchain
                                    print(line[:20] + " " + (residuedict[reskey]) + line[22:], end="")
                            elif is_hex_or_int(resid1(line))== "hex":
                                reskey = reskey1(line)
                                reschain = reschain1(line)
                                resid=resid1(line)
                                res2, chain2, resid2 = line[17:20].strip(), line[20:22].strip(), line[22:26].strip()
                                used_chains.add(reschain)
                                if chain2.strip() == "":
                                    chain2 = "A"
                                reschainkey = res2+chain2+resid2
                                resid=convert_hex_to_four_digit(resid1(line))
                                if resid>9999:
                                    if reschainkey not in lastOldNewchainID:
                                        reschain,resid,reschainid=generate_chain_id(uniquiechains,used_chains,used_chainID)
                                        lastOldNewchainID[reschainkey1(line)]=reschain.rjust(2)+str(resid).rjust(4)
                                        print(line[:20] +reschain.rjust(2) +str(resid).rjust(4)+  line[26:], end="")

                                    elif reschainkey in lastOldNewchainID:
                                        reschainkey = res2+chain2+resid2
                                        reschainidkey=lastOldNewchainID.get(reschainkey)
                                        print(line[:20] +reschainidkey+  line[26:], end="")
                             
    sys.stdout = sys.__stdout__
