#!/usr/bin/env python

import sys, os, argparse
import subprocess

parser = argparse.ArgumentParser(description='RUN COMMAND WITHIN PSI4 FSAPT DIRECTORY. Takes a PDB FG analysis file, generates fA/fB.dat files analyzes FSAPT results, saves compilation')
parser.add_argument('-p', default='../pdbFG.dat', help='name/location of PDB FG datafile to analyze, default = ../pdbFG.dat')
parser.add_argument('-path', help='path to PSI4 fsapt.py script; example: /home/ndyonker/git/psi4/objdir/stage/share/psi4/fsapt/')
parser.add_argument('-a', nargs = '+', help='optional: specify specific FG atoms of seed for analysis in format inside pdbFG.dat file (e.g. A/1/CA A/1/CB ...)')
parser.add_argument('-save', default='../FG-SAPT.dat', help='name/location of savefile, default = ../FG-SAPT.dat')
args = parser.parse_args()

def genFA(argA, seedA):
    if argA:
        partA = [seedA.index(x)+1 for x in argA]
        partB = list(set(range(1, len(seedA)+1))-set(partA))
    else:
        partA = list(range(1, len(seedA)+1))
        partB = []

    with open('fA.dat', 'w') as savefile:
        savefile.write("seedA "+" ".join(map(str, partA))+"\n")
        savefile.write("seedB "+" ".join(map(str, partB))+"\n")
    return

def genFB(seedatoms, totatoms, fatoms):
    partA = [int(x) for x in fatoms]
    partB = list(set(range(seedatoms+1, totatoms+1))-set(partA))

    with open('fB.dat', 'w') as savefile:
        savefile.write("enzyA "+" ".join(map(str, partA))+"\n")
        savefile.write("enzyB "+" ".join(map(str, partB))+"\n")
    return

def grabFSAPTenergy():
    datafile = open('fsapt.dat', 'r').readlines()
    for i in range(len(datafile)):
        line = datafile[i]
        if "F-ISAPT: Links 50-50" in line and "No exact dispersion present" in datafile[i+2]:
            index = datafile.index(line)
            data = datafile[index+7].split()[2:]
        elif "F-ISAPT: Links 50-50" in line and "Full Analysis" in datafile[i+2]:
            index = datafile.index(line)
            data = datafile[index+5].split()[2:]
    return data


if __name__ == '__main__':
    
    if not os.path.isfile(args.p):
        sys.exit("Error: unable to find -p (default ../pdbFG.dat) file")

    fginfo = open(args.p, 'r').readlines()
    natoms = int(fginfo[0])
    
    #Generate fA.dat file
    genFA(args.a, fginfo[1].split())

    #For each functional group
    fsaptdata = dict()

    for fgline in fginfo[2:]:
        fg = fgline.split()[0]
        fgatoms = fgline.split()[1:]

        #Generate fB.dat file
        genFB(len(fginfo[1].split()), natoms, fgatoms)

        #Run PSI4 fsapt.py script
        subprocess.call(["python",args.path+"fsapt.py"])
        
        #Store relevant information from generated fsapt.dat file
        fsaptdata[fg] = grabFSAPTenergy()

    #Save compiled data
    with open(args.save, 'w') as savefile:
        savefile.write("FG Elst Exch IndAB IndBA Disp EDisp Total\n")
        for item in fsaptdata.keys():
            savefile.write(" ".join([item]+fsaptdata[item])+"\n")
    print("Completion of FSAPT FG data compilation: "+os.getcwd())
