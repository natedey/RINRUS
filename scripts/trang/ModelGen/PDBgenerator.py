import sys
import os

def makePDB(inputfile1, inputfile2):
	# read frequency txt file
	key_list = []
	key_freq = {}
	with open(inputfile1, "r") as txtfile:
		dir_name = os.path.dirname(txtfile.name)
		for line in txtfile:
			if (">" not in line):
				parts = line.strip().split(":")
				parts[0] = parts[0].replace("(","").replace(")","").replace(" ","")
				keys = parts[0].split(",")
				key_list.append(tuple(keys))
				key_freq[tuple(keys)] = parts[1]


	# create result directory
	result_dir = dir_name+"-"+os.path.basename(inputfile2).replace(".pdb","")
	path = os.path.join(os.getcwd(),result_dir)
	if (os.path.exists(path) == False):
		os.mkdir(path)

	print ("Results are saved in ", result_dir)

	idx = 0
	idx_file = open(os.path.join(result_dir,"idx_file.csv"), "w")
	idx_file.write("Filename,Residue list,Frequency\n")
	for keys in key_list:
		idx = idx + 1

		fname = os.path.basename(inputfile2).replace(".pdb","")
		fname = fname+"-"+str(idx)+".pdb"
		idx_file.write(fname+","+"-".join(keys)+","+key_freq[keys]+"\n")
		# open output pdb file
		outfile = open(os.path.join(result_dir,fname), "w")
		for k in keys:
			# open pdb file
			with open(inputfile2, "r") as pdbfile:
				for line in pdbfile:
					if (line[0:4] == "ATOM") or (line[0:6] == "HETATM"):
						if (int(k) == int(line[23:26])):
							outfile.write(line)
		
		outfile.close()

	idx_file.close()

if __name__ == "__main__":
	if (len(sys.argv) < 3):
		print ("Usage: file_contains_unique_models pdb_to_trim")
		exit()

	inputfile1 = sys.argv[1] # set frequency txt file
	inputfile2 = sys.argv[2] # pdb file

	makePDB(inputfile1, inputfile2)