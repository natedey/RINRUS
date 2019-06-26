#!/usr/bin/env python

#This script was written by tsmmers1
#Warning: this script is not designed for application beyond the initial analytic function
#Usage: ./scriptname.py
#Function: This script examines the lists of residue sets among four different model construction methods
#          and determines overlap among the different set compositions

import os
import sys
import re

file1 = open("distance-based","r")
file1lines = file1.readlines()
file1.close()

file2 = open("trang-based","r")
file2lines = file2.readlines()
file2.close()

file3 = open("freq-based","r")
file3lines = file3.readlines(0)
file3.close()

file4 = open("random-based","r")
file4lines = file4.readlines(0)
file4.close()

distset = []
trangset = []
freqset = []
randset = []
kulikset = [[300,301,302],[141,169,170,300,301,302,411],[42,66,90,141,169,170,199,300,301,302,411],[40,41,42,64,67,68,71,72,90,91,119,141,142,169,170,199,300,301,302,411],[38,40,41,42,64,66,67,68,71,72,89,90,91,92,95,118,119,120,139,141,142,143,144,146,169,170,198,199,300,301,302,411]]
martinset = [[300,301,302],[42,144,199,300,301,302,411],[42,90,119,141,144,170,199,300,301,302,411],[42,66,68,71,72,90,91,119,141,142,144,170,199,300,301,302,411],[40,42,66,68,71,72,89,90,91,119,120,139,141,142,144,146,170,198,199,300,301,302,411],[38,40,42,66,67,68,71,72,89,90,91,119,120,139,141,142,143,144,146,170,198,199,300,301,302,411],[38,40,41,42,66,67,68,71,72,89,90,91,92,95,118,119,120,139,141,142,143,144,146,170,198,199,300,301,302,411],[6,38,40,41,42,46,64,66,67,68,70,71,72,89,90,91,92,95,117,118,119,120,139,141,142,143,144,146,169,170,173,174,198,199,300,301,302,411],[5,6,7,10,38,40,41,42,43,46,64,65,66,67,68,69,70,71,72,73,189,90,91,92,95,117,118,119,120,139,141,142,143,144,145,146,147,150,169,170,173,174,198,199,205,300,301,302,411],[4,5,6,7,8,9,10,11,38,39,40,41,42,43,46,64,65,66,67,68,69,70,71,72,73,76,88,89,90,91,92,95,96,98,99,115,117,118,119,120,121,139,140,141,142,143,144,145,146,147,150,169,170,173,174,198,199,205,300,301,302,411],[5,6,7,10,38,40,41,42,43,46,64,65,66,67,68,69,70,71,72,73,89,90,91,92,95,117,118,119,120,139,141,142,143,144,145,146,147,169,170,173,174,198,199,300,301,302,411],[4,5,6,7,8,9,10,11,38,39,40,41,42,43,46,64,65,66,67,68,69,70,71,72,73,76,88,89,90,91,92,95,96,98,99,115,117,118,119,120,121,139,140,141,142,143,144,145,146,147,169,170,173,174,198,199,300,301,302,411]]

for line in file1lines:
    resset = line.split('[')[1].split(']')[0] 
    resset = [int(x) for x in resset.split(',')]
    resset.sort()
    distset.append(resset)

for line in file2lines:
    resset = line.split('[')[1].split(']')[0]
    resset = [int(x) for x in resset.split(',')]
    resset.sort()
    trangset.append(resset)

for line in file3lines:
    resset = line.split('[')[1].split(']')[0]
    resset = [int(x) for x in resset.split(',')]
    resset.sort()
    freqset.append(resset)

for line in file4lines:
    resset = line.split('(')[1].split(')')[0]
    resset = [int(x) for x in resset.split(',')]
    resset.sort()
    randset.append(resset)

counter = 0
counter2 = 0
for res in trangset:
    counter += 1
    if res not in randset:
        print res
        print " " + str(len(res))
        counter2 += 1
    else:
        print('\t',res)
print(str(counter) + " " + str(counter2) + " " + str(counter2/counter*100))



