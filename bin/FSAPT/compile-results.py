#!/usr/bin/env python
"""
This is a program created by Dr. Qianyi Cheng and the DeYonker research group
at The University of Memphis.
Date created 08.12.2022
"""

import sys, os, argparse

parser = argparse.ArgumentParser(description='Compiles FG-SAPT.dat data from multiple directories into one .csv file')
parser.add_argument('-l', help='name of listfile of PDBs containing FG-SAPT.dat, FG-probe.dat and FG-arpeggio.dat files to analyze')
parser.add_argument('-save', default='compiled-FG-SAPT.csv', help='name of .csv savefile, default = compiled-FG-SAPT.csv')
args = parser.parse_args()

if __name__ == '__main__':
    
    dirlist = open(args.l, 'r').readlines()
    for directory in dirlist:
        if not os.path.isfile(directory.strip()+"/FG-SAPT.dat"):
            sys.exit("Error: unable to find FG-SAPT.dat file in "+directory.strip())
        if not os.path.isfile(directory.strip()+"/FG-probe.dat"):
            sys.exit("Error: unable to find FG-probe.dat file in "+directory.strip())
        if not os.path.isfile(directory.strip()+"/FG-arpeggio.dat"):
            sys.exit("Error: unable to find FG-arpeggio.dat file in "+directory.strip())

    savefile = open(args.save, 'w')
    savefile.write("Model,FG,Elst,Exch,IndAB,IndBA,Disp,Total,Contacts,Clash,Covalent,VdWClash,VdW,Proximal,HBond,weakHBond,Halogen,Ionic,Metal,Aromatic,Hydrophobic,Carbonyl,Polar,weakPolar\n")

    for directory in dirlist:

        #Grab Probe info
        contacts = {}
        with open(directory.strip()+"/FG-probe.dat", 'r') as readfile:
            for line in readfile:
                contacts[line.split()[0]] = line.split()[1]

        #Grab Arpeggio
        arpeggio = {}
        with open(directory.strip()+"/FG-arpeggio.dat", 'r') as readfile:
            for line in readfile:
                arpeggio[line.split()[0]] = line.split()[1:]

        #Grab SAPT and save to compiled file
        with open(directory.strip()+"/FG-SAPT.dat", 'r') as readfile:
            for line in readfile:
                if line.split()[0] == "FG": continue
                else:
                    FG = line.split()[0]

                    model = directory.strip()
                    saptstr = ",".join(line.split())

                    if FG in contacts.keys():
                        probestr = contacts[FG]
                    else:
                        probestr = "0"

                    if FG in arpeggio.keys():
                        arpstr = ",".join(arpeggio[FG])
                    else:
                        arpstr = ",".join([0]*15)

                    savefile.write(",".join([model,saptstr,probestr,arpstr])+"\n")

    savefile.close()



