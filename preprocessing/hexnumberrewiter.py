#! /usr/bin/env python
import argparse
import sys

"""
This is a program written by Tejaskumar Suhgaia in Deyonker research group
at university of memphis.
"""

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
            return "integer"
        elif is_hex(value) == True:
            return "hexadecimal"        
    except ValueError:
        return "neither"
def convert_hex_to_four_digit(hex_value):
    return str(int(hex_value, 16)).zfill(4)

def assign_new_chain(old_chain, used_chains):
    new_chain = chr(ord(old_chain) + 1)
    while new_chain in used_chains:
        new_chain = chr(ord(new_chain) + 1)
    used_chains.append(new_chain)
    return new_chain
def process_pdb_file(file_path):
    new_pdb_lines = []
    used_chains = set()


    # print(convert_hex_to_four_digit("2e5f"))
    # print(convert_hex_to_four_digit("2e5e"))
    # print(convert_hex_to_four_digit("2e60"))
    # integ1="42"
    # print(is_hex("2e60"))
    # print(isNumber("2160"))
    
if __name__ == '__main__':
    """ Usage: hexnumberrewiter.py -p template.pdb """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p', dest='pdbf', default='template.pdb', help='pdb file')
    args = parser.parse_args()
    pdbf = args.pdbf
    output_file_path = 'HexCorrected_'+pdbf  
    new_pdb_lines=[]
    with open(pdbf, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract relevant information
                old_chain = line[21]
                hex_residue_number = line[22:26]
                used_chains=[]

                # Convert hexadecimal residue number to four-digit number
                new_residue_number = convert_hex_to_four_digit(hex_residue_number)

                # Assign a new chain alphabet
                new_chain = assign_new_chain(old_chain, used_chains)

                # Replace old chain and residue number with new values
                new_line = line[:21] + new_chain + line[22:26] + new_residue_number + line[26:]
                new_pdb_lines.append(new_line)
            else:
                new_pdb_lines.append(line)

    # Write the modified PDB file
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(new_pdb_lines)

if __name__ == "__main__":
    pdb_file_path = "your_pdb_file.pdb"
    process_pdb_file(pdb_file_path)
