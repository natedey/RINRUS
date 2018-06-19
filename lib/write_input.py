import os, sys, re
from numpy import *


def write_gassian_input(head,foot,pic_atom,count):
    with open(head) as f:
        hlines = f.readlines()
    with open(foot) as f:
        flines = f.readlines()

    input = open('%d.input'%count,'w')
    for i in range(len(hlines)):
        input.write(hlines[i])
    
    for i in range(len(pic_atom)):
        input.write("%4s %6s %8.3f %8.3f %8.3f\n"%(atom[i],hold[i],xyz[i][0],xyz[i][1],xyz[i][2])) 
    
    input.write("\n")

    for i in range(len(flines)):
        input.write(flines[i])

    input.write("\n")
    input.close()
    
#        if "%" in hlines[i] or "#" in hlines[i]:
#            input.write(hlines[i])
#
#    input = open('%d.input'%count,'w')
#    input.write("%chk=1.chk\n")     #write check file into 1.chk
#    input.write("%nproc=24\n")      #use 24 processors
#    input.write("%mem=14GB\n")      #use memory
#    
#    input.write("#P b3lyp/gen opt freq scf=(xqc,maxcon=128,maxcyc=128)\n")      #many things can be changed to user input
#    input.write("\n")
#    title = raw_input("Input the descriptive name of this calculation: ")
#    input.write("%s\n"%title)
#    input.write("\n")
#    input.write("1 1\n")        #charge and multiplicity, also need to be changed case by case
#    
#    
#    for i in range(len(pic_atom)):
#        input.write("%4s %6s %8.3f %8.3f %8.3f\n"%(atom[i],hold[i],xyz[i][0],xyz[i][1],xyz[i][2])) 
#    
#    input.write("\n")
#    input.write("O N S\n")
#    input.write("6-31G(d')\n")
#    input.write("****\n")
#    input.write("C H\n")
#    input.write("6-31G\n")
#    input.write("****\n")
#    input.write("\n")
#    input.close()
