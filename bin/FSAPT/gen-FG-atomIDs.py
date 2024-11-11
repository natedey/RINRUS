#!/usr/bin/env python

import sys, os, argparse
import fnmatch
from collections import defaultdict

parser = argparse.ArgumentParser(description='Generates file documenting atom IDs of functional groups of PDB file')
parser.add_argument('-p', help='name of PDB file to analyze')
parser.add_argument('-s', help='residue(s) of seed in format chain:resID,chain:resID,...')
parser.add_argument('-save', default='pdbFG.dat', help='name of savefile')
args = parser.parse_args()

modelsize = 0
seedatoms = []
watatomnames = defaultdict(list)
resatomnames = defaultdict(list)
GLYlist = []
FGlist = []

#Grab seed info
with open(args.p, 'r') as readfile:
    for line in readfile:
        if line == "END\n" or line[0:4] == "TER " or line[0:6] == "CONECT": continue
        if line[21]+":"+line[22:26].strip() in args.s.split(","):
            seedatoms.append(line[21]+"/"+line[22:26].strip()+"/"+line[12:16].strip())
        if line[0:4] == "ATOM" or line[0:6] == "HETATM": modelsize +=1

#Grab residue and water info
with open(args.p, 'r') as readfile:
    atomID = len(seedatoms)+1
    for line in readfile:
        if line == "END\n": continue
        if line[21]+":"+line[22:26].strip() in args.s.split(","): continue
        if line[17:20] == "WAT" or line[17:20] == "HOH":
            watatomnames[line[21]+":"+line[22:26].strip()].append(line[12:16].strip()+":"+str(atomID))
            atomID +=1
        else:
            resatomnames[line[21]+":"+line[22:26].strip()].append(line[12:16].strip()+":"+str(atomID))
            atomID +=1
            if line[17:20] == "GLY":
                GLYlist.append(line[21]+":"+line[22:26].strip())

for res in resatomnames.keys():
    #account for residue SC not bound to adjacent residues
    if res.split(":")[0]+":"+str(int(res.split(":")[1])-1) not in resatomnames.keys() and res.split(":")[0]+":"+str(int(res.split(":")[1])+1) not in resatomnames.keys():
        FGlist.append([res+":SC"]+[x.split(":")[1] for x in resatomnames[res]])
        continue
    #account for SC functional groups and GLY SCs
    if fnmatch.filter(resatomnames[res], "CB:*") or res in GLYlist:
        SCnums = []
        for atom in resatomnames[res]:
            if atom.split(":")[0] not in ["C", "O", "N", "H"]:
                SCnums.append(atom.split(":")[1])
        FGlist.append([res+":SC"]+SCnums)
    #account for MC functional groups
    if res.split(":")[0]+":"+str(int(res.split(":")[1])+1) in resatomnames.keys():
        if fnmatch.filter(resatomnames[res], "C:*") and fnmatch.filter(resatomnames[res], "O:*") and fnmatch.filter(resatomnames[res.split(":")[0]+":"+str(int(res.split(":")[1])+1)], "N:*"):
            if fnmatch.filter(resatomnames[res.split(":")[0]+":"+str(int(res.split(":")[1])+1)], "H:*"):
                FGlist.append([res+":MC", fnmatch.filter(resatomnames[res], "C:*")[0].split(":")[1], fnmatch.filter(resatomnames[res], "O:*")[0].split(":")[1], 
                    fnmatch.filter(resatomnames[res.split(":")[0]+":"+str(int(res.split(":")[1])+1)], "N:*")[0].split(":")[1], fnmatch.filter(resatomnames[res.split(":")[0]+":"+str(int(res.split(":")[1])+1)], "H:*")[0].split(":")[1]])
            else:
                FGlist.append([res+":MC", fnmatch.filter(resatomnames[res], "C:*")[0].split(":")[1], fnmatch.filter(resatomnames[res], "O:*")[0].split(":")[1], 
                    fnmatch.filter(resatomnames[res.split(":")[0]+":"+str(int(res.split(":")[1])+1)], "N:*")[0].split(":")[1]])

#Save information to file
with open(args.save, 'w') as savefile:
    savefile.write(str(modelsize)+"\n")
    savefile.write(" ".join(seedatoms)+"\n")
    for fgs in FGlist:
        savefile.write(" ".join(fgs)+"\n")
    for wat in watatomnames.keys():
        savefile.write(" ".join([wat+":WAT"]+[x.split(":")[1] for x in watatomnames[wat]])+"\n")




