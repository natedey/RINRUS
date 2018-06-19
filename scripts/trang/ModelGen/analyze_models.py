import sys, os
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict

def Strip_analyze(inputfile):

	# read input file
	# input file has format:
	# edge 1, edge 2,....edge n | res 1, res 2, ..., res n
	print ("Reading ", inputfile)
	combo = {}
	with open(inputfile,"r") as txtfile:
		for line in txtfile:
			line = line.strip()
			if ("|" in line):
				parts = line.split("|")
				key = parts[0].split(",")
				residue_list = parts[1].split(",")
				combo[tuple(key)] = residue_list

	# print (combo)

	# create result directory
	result_dir = inputfile.replace(".txt","")
	path = os.path.join(os.getcwd(),result_dir)
	if (os.path.exists(path) == False):
		os.mkdir(path)

	# creat ID list
	print("Creating ID list\n")
	id_list = {}
	combo_id = {}
	for key in combo:
		combo_id[key] = []
		for i in combo[key]:
			parts = i.split(":")
			combo_id[key].append(int(parts[1]))
			if parts[1] not in id_list:
				id_list[parts[1]] = ""
				id_list[parts[1]] = parts[3]

	csv_id = open(os.path.join(path, "id_list.csv"), "w")
	for key in id_list:
		csv_id.write(key+","+id_list[key]+"\n")

	csv_id.close()
	
	###################################
	# histogram of all sets
	###################################
	val = []
	max_key = 0
	for key in combo:
		val.append(len(key))
		if (max_key < len(key)):
			max_key = len(key)

	plt.clf()
	plt.hist(val,bins=max_key,rwidth=0.9)
	plt.xlabel("Number of edges")
	plt.ylabel("Number of sets")
	plt.title("Total # of sets : "+str(len(combo)))
	plt.plot()
	plt.savefig(os.path.join(path,'all_sets.png'))

	# get neighbors by ID
	neighbors = {}
	for key in combo_id:
		for i in combo_id[key]:
			if i not in neighbors:
				neighbors[i] = 1
			else:
				neighbors[i] = neighbors[i] + 1

	# print (neighbors)
	###################################
	# bar chart of neighbors by ID
	###################################
	keys = []
	values = []
	for key in neighbors:
		keys.append(int(key))
		values.append(neighbors[key])

	y_pos = np.arange(len(keys))
	 
	plt.clf()
	plt.bar(y_pos, values, align='center', alpha=0.5)
	plt.xticks(y_pos, keys, rotation=90)
	plt.ylabel('Count')
	plt.xlabel('Residue ID')
	plt.title("Residue frequency in all sets\nTotal # of sets : "+str(len(combo)))
	plt.tight_layout() # make rooms for the x-axis labels
	plt.plot()
	plt.savefig(os.path.join(path,'neighbor.png'))

	###################################
	# count unique sets
	# edges by unique set
	###################################
	unique_sets = {}
	edges_by_unique_sets = {}
	for key in combo_id:
		l = tuple(sorted(combo_id[key]))
		if (l not in unique_sets):
			unique_sets[l] = 1
			
		else:
			unique_sets[l] = unique_sets[l] + 1

		if (l not in edges_by_unique_sets):
			edges_by_unique_sets[l] = []

		edges_by_unique_sets[l].append(key)

	# histogram of unique residue groups

	keys = []
	values = []
	short_keys = []
	for key in sorted(unique_sets, key=unique_sets.get, reverse=True):
		keys.append(key)
		if (len(key)>3):
			short_keys.append("")
		values.append(unique_sets[key])

	y_pos = np.arange(len(keys))

	plt.clf() 
	plt.bar(y_pos, values, align='center')
	plt.xticks([], [], rotation=90)
	plt.ylabel('Count')
	plt.xlabel('Unique set')
	plt.title("Frequency of unique sets\nTotal # of unique sets : "+str(len(unique_sets)))
	plt.tight_layout() # make rooms for the x-axis labels
	plt.plot()
	plt.savefig(os.path.join(path,'unique_sets_freq.png'), dpi=200)

	# print unique_sets_freq to file
	txt_unique_sets = open(os.path.join(path,"unique_sets_freq.txt"), "w")
	# txt_unique_sets.write("Total "+str(len(unique_sets))+"\n")
	for key in sorted(unique_sets, key=unique_sets.get, reverse=True):
		txt_unique_sets.write(str(key)+":"+str(unique_sets[key])+"\n")

	txt_unique_sets.close()

	# print edges by unique sets
	txt_edges_by_unique_sets = open(os.path.join(path, "edges_by_unique_sets.txt"),"w")
	for key in edges_by_unique_sets:
		txt_edges_by_unique_sets.write(str(key)+"|"+str(edges_by_unique_sets[key])+"\n")

	txt_edges_by_unique_sets.close()

	###################################
	# frequency of residue in each unique sets
	###################################
	residue_in_unique_sets = {}
	for key in unique_sets:
		for i in key:
			if (i not in residue_in_unique_sets):
				residue_in_unique_sets[i] = 1
			else:
				residue_in_unique_sets[i] = residue_in_unique_sets[i] + 1

	# print (residue_in_unique_sets)

	# sorted by residue #
	keys = []
	values = []
	for key in residue_in_unique_sets:
		keys.append(int(key))
		values.append(residue_in_unique_sets[key])

	y_pos = np.arange(len(keys))

	plt.clf() 
	plt.bar(y_pos, values, align='center', alpha=0.5)
	plt.xticks(y_pos, keys, rotation=90)
	plt.ylabel('Count')
	plt.xlabel('Residue ID')
	plt.title("Residue frequency in unique sets\nTotal # of unique sets : "+str(len(unique_sets)))
	plt.tight_layout() # make rooms for the x-axis labels
	plt.plot()
	plt.savefig(os.path.join(path,'residue_in_unique_sets.png'))

	# sorted by frequency
	keys = []
	values = []
	for key in sorted(residue_in_unique_sets, key=residue_in_unique_sets.get, reverse=True):
		keys.append(int(key))
		values.append(residue_in_unique_sets[key])

	y_pos = np.arange(len(keys))

	plt.clf() 
	plt.bar(y_pos, values, align='center', alpha=0.5)
	plt.xticks(y_pos, keys, rotation=90)
	plt.ylabel('Count')
	plt.xlabel('Residue ID')
	plt.title("Residue frequency in unique sets\nTotal # of unique sets : "+str(len(unique_sets)))
	plt.tight_layout() # make rooms for the x-axis labels
	plt.plot()
	plt.savefig(os.path.join(path,'residue_in_unique_sets_sorted.png'))

	print ("Results are saved in "+path)

if __name__ == "__main__":
	if (len(sys.argv) < 2):
		print ("Usage: file_contains_models")
		exit()

	inputfile = sys.argv[1]

	Strip_analyze(inputfile)

