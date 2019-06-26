#!/usr/bin/env python

#Usage: ./scriptname.py
#This script was written by tsmmers1
#Warning: this script is not designed to function beyond the limited scope of this work
#Function: Analyze the presence/frequency of specific residues in the combinatoric set of COMT models.

import os
import sys
import re

filename = "unique_sets_freq.txt"

myfile = open(filename)
mylines = myfile.readlines()
myfile.close()

#Determine frequency of residues missing
list141 = []
list169 = []
list170 = []
list300 = []
list301 = []
list302 = []
list411 = []
list199 = []

total141 = 0
total169 = 0
total170 = 0
total300 = 0
total301 = 0
total302 = 0
total411 = 0
total199 = 0

for i in range(len(mylines)):
    linesplit = mylines[i].split(":")
    if "141" not in linesplit[0]:
        list141.append(i)
        total141 += int(linesplit[1])
    if "169" not in linesplit[0]:
        list169.append(i)
        total169 += int(linesplit[1])
    if "170" not in linesplit[0]:
        list170.append(i)
        total170 += int(linesplit[1])
    if "300" not in linesplit[0]:
        list300.append(i)
        total300 += int(linesplit[1])
    if "301" not in linesplit[0]:
        list301.append(i)
        total301 += int(linesplit[1])
        print(linesplit)
    if "302" not in linesplit[0]:
        list302.append(i)
        total302 += int(linesplit[1])
    if "411" not in linesplit[0]:
        list411.append(i)
        total411 += int(linesplit[1])
    if "199" not in linesplit[0]:
        list199.append(i)
        total199 += int(linesplit[1])

print("\nInfo on models missing coordinating species:\n")
print("\t" + str(len(list141)) + " (" + str(len(list141)*100/2325) + "%) unique sets do not have ASP141. These account for " + str(total141) + " (" + str(total141*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list169)) + " (" + str(len(list169)*100/2325) + "%) unique sets do not have ASP169. These account for " + str(total169) + " (" + str(total169*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list170)) + " (" + str(len(list170)*100/2325) + "%) unique sets do not have ASN170. These account for " + str(total170) + " (" + str(total170*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list300)) + " (" + str(len(list300)*100/2325) + "%) unique sets do not have Mg300. These account for " + str(total300) + " (" + str(total300*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list301)) + " (" + str(len(list301)*100/2325) + "%) unique sets do not have SAM301. These account for " + str(total301) + " (" + str(total301*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list302)) + " (" + str(len(list302)*100/2325) + "%) unique sets do not have DNC302. These account for " + str(total302) + " (" + str(total302*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list411)) + " (" + str(len(list411)*100/2325) + "%) unique sets do not have HOH411. These account for " + str(total411) + " (" + str(total411*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list199)) + " (" + str(len(list199)*100/2325) + "%) unique sets do not have GLU199. These account for " + str(total199) + " (" + str(total199*100/524287) +  "%) of the total sets.") 

compiled = set(list300) | set(list301) | set(list302)
compiledtot = 0
for i in compiled:
    compiledtot += int(mylines[i].split(":")[1])
print("\n\t" + str(len(compiled)) + " (" + str(len(compiled)*100/2325) + "%) unique sets are missing one or more primary ligands (Res 300,301,302). These account for " + str(compiledtot) + " (" + str(compiledtot*100/524287) +  "%) of the total sets.")

compiled = set(list300) | set(list301) | set(list302) | set(list141) | set(list169) | set(list170) | set(list411)
compiledtot = 0
for i in compiled:
    compiledtot += int(mylines[i].split(":")[1])
print("\n\t" + str(len(compiled)) + " (" + str(len(compiled)*100/2325) + "%) unique sets are missing one or more Mg-associated residues. These account for " + str(compiledtot) + " (" + str(compiledtot*100/524287) +  "%) of the total sets.")

compiled = set(list300) | set(list301) | set(list302) | set(list141) | set(list169) | set(list170) | set(list411) | set(list199)
compiledtot = 0
for i in compiled:
    compiledtot += int(mylines[i].split(":")[1])
print("\n\t" + str(len(compiled)) + " (" + str(len(compiled)*100/2325) + "%) unique sets are missing one or more Mg-associated residues and Glu199. These account for " + str(compiledtot) + " (" + str(compiledtot*100/524287) +  "%) of the total sets.")


correctedsets = []
correctedfreq = []

for i in range(len(mylines)):
    linesplit = mylines[i].split(":")
    linesplit2 = linesplit[0].split("(")
    linesplit2 = linesplit2[1].split(")")
    linesplit2 = set(map(int, linesplit2[0].split(",")))
    if linesplit2 in correctedsets:
        index = correctedsets.index(linesplit2)
        correctedfreq[index] += int(linesplit[1])
    elif linesplit2 not in correctedsets:
        correctedsets.append(linesplit2)
        correctedfreq.append(int(linesplit[1]))

myzip = zip(correctedsets, correctedfreq)
myzip.sort(key = lambda t: (len(t[0]),t[1]), reverse = True)

file1 = open("OrigResList.txt", "w")
for i in range(len(myzip)):
   file1.write(str(myzip[i]) + "\n")
file1.close()

file1 = open("OrigResList.txt", "r")
file1lines = file1.readlines()
file1.close()

file1 = open("OrigResList.txt", "w")
linelength = 0
for line in file1lines:
    linesplit = re.split('\[([^[\]]*)\]',line)
    freq = linesplit[2].replace(")","")
    freq = int(freq.replace(",",""))
    splitline = linesplit[1].split(",")
    charge = 0
    splitline = list(map(int, splitline))
    splitline.sort()
    for i in range(len(splitline)):
        if int(splitline[i])==199:
            charge += -1
        if int(splitline[i])==141:
            charge += -1
        if int(splitline[i])==90:
            charge += -1
        if int(splitline[i])==169:
            charge += -1
        if int(splitline[i])==146:
            charge += 1
        if int(splitline[i])==144:
            charge += 1
        if int(splitline[i])==300:
            charge += 2
        if int(splitline[i])==301:
            charge += -1	     
    if len(splitline) == linelength:
        file1.write("(" + ",".join(map(str,splitline)) + "):" + str(freq) + " (" + str(charge) + ")" + "\n")
    else:
#        file1.write("\n#" + str(len(splitline)) + " Residues\n")
        file1.write("(" + ",".join(map(str,splitline)) + "):" + str(freq) + " (" + str(charge) + ")" +  "\n")
        linelength = len(splitline)
file1.close()    


##############

correctedsets = []
correctedfreq = []

for i in range(len(mylines)):
    linesplit = mylines[i].split(":")
    linesplit2 = linesplit[0].split("(")
    linesplit2 = linesplit2[1].split(")")
    linesplit2 = set(map(int, linesplit2[0].split(",")))
    linesplit2.add(300)
    linesplit2.add(301)
    linesplit2.add(302)
    if linesplit2 in correctedsets:
        index = correctedsets.index(linesplit2)
        correctedfreq[index] += int(linesplit[1])
    elif linesplit2 not in correctedsets:
        correctedsets.append(linesplit2)
        correctedfreq.append(int(linesplit[1]))
print("\nInfo on Corrected Models:\n")
print("\tFinal number of unique corrected residue sets when (Res 300,301,302) added: " + str(len(correctedsets)) + " (" + str(len(correctedsets)*100/2325) + "%)")

list141 = []
list169 = []
list170 = []
list411 = []
list199 = []

total141 = 0
total169 = 0
total170 = 0
total411 = 0
total199 = 0

for i in range(len(correctedsets)):
    if 141 not in correctedsets[i]:
        list141.append(i)
        total141 += int(correctedfreq[i])
    if 169 not in correctedsets[i]:
        list169.append(i)
        total169 += int(correctedfreq[i])
    if 170 not in correctedsets[i]:
        list170.append(i)
        total170 += int(correctedfreq[i])
    if 411 not in correctedsets[i]:
        list411.append(i)
        total411 += int(correctedfreq[i])
    if 199 not in correctedsets[i]:
        list199.append(i)
        total199 += int(correctedfreq[i])
print("\t" + str(len(list141)) + " (" + str(len(list141)*100/len(correctedsets)) + "%) unique sets do not have ASP141. These account for " + str(total141) + " (" + str(total141*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list169)) + " (" + str(len(list169)*100/len(correctedsets)) + "%) unique sets do not have ASP169. These account for " + str(total169) + " (" + str(total169*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list170)) + " (" + str(len(list170)*100/len(correctedsets)) + "%) unique sets do not have ASN170. These account for " + str(total170) + " (" + str(total170*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list411)) + " (" + str(len(list411)*100/len(correctedsets)) + "%) unique sets do not have HOH411. These account for " + str(total411) + " (" + str(total411*100/524287) +  "%) of the total sets.")
print("\t" + str(len(list199)) + " (" + str(len(list199)*100/len(correctedsets)) + "%) unique sets do not have GLU199. These account for " + str(total199) + " (" + str(total199*100/524287) +  "%) of the total sets.")

compiled = set(list141) | set(list169) | set(list170) | set(list411)
compiledtot = 0
for i in compiled:
    compiledtot += int(mylines[i].split(":")[1])
print("\n\t" + str(len(compiled)) + " (" + str(len(compiled)*100/2325) + "%) unique sets are missing one or more Mg-associated residues. These account for " + str(compiledtot) + " (" + str(compiledtot*100/524287) +  "%) of the total sets.")

compiled = set(list141) | set(list169) | set(list170) | set(list411) | set(list199)
compiledtot = 0
for i in compiled:
    compiledtot += int(mylines[i].split(":")[1])
print("\n\t" + str(len(compiled)) + " (" + str(len(compiled)*100/2325) + "%) unique sets are missing one or more Mg-associated residues and Glu199. These account for " + str(compiledtot) + " (" + str(compiledtot*100/524287) +  "%) of the total sets.")


myzip = zip(correctedsets, correctedfreq)
myzip.sort(key = lambda t: (len(t[0]),t[1]), reverse = True)

file1 = open("3bwm_w300-301-302.txt", "w")
for i in range(len(myzip)):
   file1.write(str(myzip[i]) + "\n")
file1.close()

file1 = open("3bwm_w300-301-302.txt", "r")
file1lines = file1.readlines()
file1.close()

file1 = open("3bwm_w300-301-302.txt", "w")
linelength = 0
for line in file1lines:
    linesplit = re.split('\[([^[\]]*)\]',line)
    freq = linesplit[2].replace(")","")
    freq = int(freq.replace(",",""))
    splitline = linesplit[1].split(",")
    charge = 0
    splitline = list(map(int, splitline))
    splitline.sort()
    for i in range(len(splitline)):
        if int(splitline[i])==199:
            charge += -1
        if int(splitline[i])==141:
            charge += -1
        if int(splitline[i])==90:
            charge += -1
        if int(splitline[i])==169:
            charge += -1
        if int(splitline[i])==146:
            charge += 1
        if int(splitline[i])==144:
            charge += 1
        if int(splitline[i])==300:
            charge += 2
        if int(splitline[i])==301:
            charge += -1	    
    if len(splitline) == linelength:
        file1.write("(" + ",".join(map(str,splitline)) + "):" + str(freq) + " (" + str(charge) + ")" + "\n")
    else:
#        file1.write("\n#" + str(len(splitline)) + " Residues\n")
        file1.write("(" + ",".join(map(str,splitline)) + "):" + str(freq) + " (" + str(charge) + ")" +  "\n")
        linelength = len(splitline)
file1.close()    

#########

correctedsets = []
correctedfreq = []

for i in range(len(mylines)):
    linesplit = mylines[i].split(":")
    linesplit2 = linesplit[0].split("(")
    linesplit2 = linesplit2[1].split(")")
    linesplit2 = set(map(int, linesplit2[0].split(",")))
    linesplit2.add(300)
    linesplit2.add(301)
    linesplit2.add(302)
    linesplit2.add(141)
    linesplit2.add(169)
    linesplit2.add(170)
    linesplit2.add(411)
    if linesplit2 in correctedsets:
        index = correctedsets.index(linesplit2)
        correctedfreq[index] += int(linesplit[1])
    elif linesplit2 not in correctedsets:
        correctedsets.append(linesplit2)
        correctedfreq.append(int(linesplit[1]))
print("\nInfo on Corrected Models:\n")
print("\tFinal number of unique corrected residue sets when Mg-coordinated species added: " + str(len(correctedsets)) + " (" + str(len(correctedsets)*100/2325) + "%)")

list199 = []
total199 = 0
for i in range(len(correctedsets)):
    if 199 not in correctedsets[i]:
        list199.append(i)
        total199 += int(correctedfreq[i])
print("\t" + str(len(list199)) + " (" + str(len(list199)*100/len(correctedsets)) + "%) unique sets do not have GLU199. These account for " + str(total199) + " (" + str(total199*100/524287) +  "%) of the total sets.")

myzip = zip(correctedsets, correctedfreq)
myzip.sort(key = lambda t: (len(t[0]),t[1]), reverse = True)

file1 = open("3bwm_wPropCoord.txt", "w")
for i in range(len(myzip)):
   file1.write(str(myzip[i]) + "\n")
file1.close()

file1 = open("3bwm_wPropCoord.txt", "r")
file1lines = file1.readlines()
file1.close()

file1 = open("3bwm_wPropCoord.txt", "w")
linelength = 0
for line in file1lines:
    linesplit = re.split('\[([^[\]]*)\]',line)
    freq = linesplit[2].replace(")","")
    freq = int(freq.replace(",",""))
    splitline = linesplit[1].split(",")
    charge = 0
    splitline = list(map(int, splitline))
    splitline.sort()
    for i in range(len(splitline)):
        if int(splitline[i])==199:
            charge += -1
        if int(splitline[i])==141:
            charge += -1
        if int(splitline[i])==90:
            charge += -1
        if int(splitline[i])==169:
            charge += -1
        if int(splitline[i])==146:
            charge += 1
        if int(splitline[i])==144:
            charge += 1
        if int(splitline[i])==300:
            charge += 2
        if int(splitline[i])==301:
            charge += -1	     
    if len(splitline) == linelength:
        file1.write("(" + ",".join(map(str,splitline)) + "):" + str(freq) + " (" + str(charge) + ")" + "\n")
    else:
#        file1.write("\n#" + str(len(splitline)) + " Residues\n")
        file1.write("(" + ",".join(map(str,splitline)) + "):" + str(freq) + " (" + str(charge) + ")" +  "\n")
        linelength = len(splitline)
file1.close()    


##################
correctedsets = []
correctedfreq = []

for i in range(len(mylines)):
    linesplit = mylines[i].split(":")
    linesplit2 = linesplit[0].split("(")
    linesplit2 = linesplit2[1].split(")")
    linesplit2 = set(map(int, linesplit2[0].split(",")))
    linesplit2.add(300)
    linesplit2.add(301)
    linesplit2.add(302)
    linesplit2.add(141)
    linesplit2.add(169)
    linesplit2.add(170)
    linesplit2.add(411)
    linesplit2.add(199)
    if linesplit2 in correctedsets:
        index = correctedsets.index(linesplit2)
        correctedfreq[index] += int(linesplit[1])
    elif linesplit2 not in correctedsets:
        correctedsets.append(linesplit2)
        correctedfreq.append(int(linesplit[1]))
print("\nInfo on Corrected Models:\n")
print("\tFinal number of unique corrected residue sets when Mg-coordinated species and Glu199 added: " + str(len(correctedsets)) + " (" + str(len(correctedsets)*100/2325) + "%)")

myzip = zip(correctedsets, correctedfreq)
myzip.sort(key = lambda t: (len(t[0]),t[1]), reverse = True)

file1 = open("3bwm_wPropCoord_w199.txt", "w")
for i in range(len(myzip)):
   file1.write(str(myzip[i]) + "\n")
file1.close()

file1 = open("3bwm_wPropCoord_w199.txt", "r")
file1lines = file1.readlines()
file1.close()

file1 = open("3bwm_wPropCoord_w199.txt", "w")
linelength = 0
for line in file1lines:
    linesplit = re.split('\[([^[\]]*)\]',line)
    freq = linesplit[2].replace(")","")
    freq = int(freq.replace(",",""))
    splitline = linesplit[1].split(",")
    charge = 0
    splitline = list(map(int, splitline))
    splitline.sort()
    for i in range(len(splitline)):
        if int(splitline[i])==199:
            charge += -1
        if int(splitline[i])==141:
            charge += -1
        if int(splitline[i])==90:
            charge += -1
        if int(splitline[i])==169:
            charge += -1
        if int(splitline[i])==146:
            charge += 1
        if int(splitline[i])==144:
            charge += 1
        if int(splitline[i])==300:
            charge += 2
        if int(splitline[i])==301:
            charge += -1	    	    
    if len(splitline) == linelength:
        file1.write("(" + ",".join(map(str,splitline)) + "):" + str(freq) + " (" + str(charge) + ")" + "\n")
    else:
#        file1.write("\n#" + str(len(splitline)) + " Residues\n")
        file1.write("(" + ",".join(map(str,splitline)) + "):" + str(freq) + " (" + str(charge) + ")" +  "\n")
        linelength = len(splitline)
file1.close()    

#################











