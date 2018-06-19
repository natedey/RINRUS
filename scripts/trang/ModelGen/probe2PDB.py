import sys
import os

if (len(sys.argv) < 4):
	print ("Usage: probe_file pdb_file residue_ID_list")
	print ("Example: 1RFF_h.probe 1RFF.pdb 699")
	print ("Example: 1RFF_h.probe 1RFF.pdb 699,263,493")
	exit()

inputfile1 = sys.argv[1] # probe file
inputfile2 = sys.argv[2] # pdb file
inputRes = sys.argv[3]

inresId = inputRes.split(",")

resId = []
resName = []
if (len(inresId) == 0):
	if (":" in inputRes):
		parts = inputRes.split(":")
		resId.append(parts[1])
		resName.append(parts[3])
	else:
		resId.append(inputRes.strip())
		resName.append("")
else:
	for res in inresId:
		if (":" in res):
			parts = res.split(":")
			resId.append(parts[1])
			resName.append(parts[3])
		else:
			resId.append(res.strip())
			resName.append("")


logfile = open("logfile.txt", "w")
resList = []
for i in range(len(resId)):
	with open(inputfile1, "r") as probefile:
		query = resId[i]+" "+resName[i]
		# print ("query ", query)
		for line in probefile:
			# print (line[12:20], line[30:37])
			if (line[10] == 'A') and (line[27] == 'A'):
				if (int(resId[i]) == int(line[11:15])) or (int(resId[i]) == int(line[28:32])):
					logfile.write(line)

					if (line[11:19] not in resList):
						resList.append(line[11:19])
					if (line[28:36] not in resList):
						resList.append(line[28:36])

logfile.close()

# print (sorted(resList))

fname = os.path.basename(inputfile2)
fname = fname.replace(".pdb", "-") + "-".join(resId) + ".pdb"
# print (fname)
outfile = open(fname, "w")
for k in sorted(resList):
	key = k.strip().split(" ")
	# print (key[0])
	with open(inputfile2, "r") as pdbfile:
		for line in pdbfile:
			if (line[21] == 'A'):
				if (line[0:4] == "ATOM") or (line[0:6] == "HETATM"):
					if (int(key[0]) == int(line[23:26])):
						outfile.write(line)

outfile.close()

print ("Result is saved in ", os.path.join(os.getcwd(),fname))



