#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    print "\nUsage: ./sort-residues.py arg1    \n\t\t arg1 = filename containing list of residues to be sorted\n"
    sys.exit()

filename = sys.argv[1]
corrected = []

with open(filename) as fp:
    for line in fp:
        header = line.split("[")[0]
        reslist = line.split("[")[1].split("]")[0]
        footer = line.split("]")[1]
        reslist = [int(x) for x in reslist.split(",")]
        reslist.sort()
        newstring = header + str(reslist) + footer
        corrected.append(newstring)

savefile = open("sorted-" + filename, "w")
for item in corrected:
    savefile.write(item)
savefile.close()

