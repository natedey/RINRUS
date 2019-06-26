#!/usr/bin/env python

###################################################
### Usage: gen_methyl-shift_input.py argv[1] 
### argv[1] pdbfile of structure to shift methyl on
###################################################

#Script written by tsmmers1
#Function: generate a new pdb file with the COMT methyl in a TS location

import sys, re, os

pdbfilename = sys.argv[1]
savefilename = str(sys.argv[1])[:-8] + "-mod.pdb"
pdbfile = open(pdbfilename, "r")
savefile = open(savefilename, "w")

for line in pdbfile:
    res = line[17:20]
    atom = line[12:16].strip()
    if res == "SAM" and atom == "CE":
        savefile.write("HETATM  218  CE  SAM A 301      -6.462  -9.262 -16.490  1.00 28.78           C0       0\n")
    elif res == "SAM" and atom == "HE1":
        savefile.write("HETATM  245  HE1 SAM A 301      -6.316  -9.343 -17.555  1.00 28.78           H0       0\n")
    elif res == "SAM" and atom == "HE2":
        savefile.write("HETATM  246  HE2 SAM A 301      -6.175  -8.359 -15.968  1.00 28.78           H0       0\n")
    elif res == "SAM" and atom == "HE3":
        savefile.write("HETATM  247  HE3 SAM A 301      -7.118  -9.953 -15.979  1.00 28.78           H0       0\n")
    else:
        savefile.write(line)

pdbfile.close()
savefile.close()
print "Generated: " + savefilename
