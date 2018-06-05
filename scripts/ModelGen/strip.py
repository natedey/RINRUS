import os
import sys
import itertools
import collections

# combination of all elements in an input list
# remove empty set
def combine(inlist):
	combo = {}

	for L in range(0, len(inlist)+1):
		for subset in itertools.combinations(inlist, L):
			# print (subset)
			if (len(subset) > 0): # remove empty set
				s = tuple(sorted(subset))
				combo[s] = []

	return combo

# save output to csv file
def saveToCSV(combo, output, inputfile, selected_nodes):
	# padding combo
	padded_combo = {}
	padded_output = {}
	padding = ('',)
	for i in output:
		padded_output[i] = []
		max_length = 0
		for j in range(len(output[i])):
			len_element = len(combo[output[i][j]])
			if (max_length < len_element):
				max_length = len_element

		for j in output[i]:
			s = j
			if (len(j) < max_length):
				s = s + padding*(max_length-len(j))
			
			padded_output[i].append(s)

			padded_combo[s] = combo[j]
			if (len(combo[j]) < max_length):
				padded_combo[s] = padded_combo[s] + list(padding)*(max_length-len(combo[j]))

	# name = "-".join(selected_nodes)
	filename = inputfile.replace(".sif", "-")+"-".join(selected_nodes)+"_py.csv"
	outfile = os.path.join(os.getcwd(), filename)
	csv_file = open(outfile, "w")

	for i in padded_output:
		for it2 in range(len(padded_output[i][0])):
			for it in range(len(padded_output[i])):
				# print (padded_output[i][it][it2], ",",padded_combo[padded_output[i][it]][it2], ",,", end="")
				csv_file.write(padded_output[i][it][it2]+","+padded_combo[padded_output[i][it]][it2]+",,")
			csv_file.write("\n")
		csv_file.write("\n")
	csv_file.close()

# save output to txt file
def saveToTXT(combo, output, inputfile, selected_nodes):
	filename = inputfile.replace(".sif", "-")+"-".join(selected_nodes)+"_py.txt"
	outfile = os.path.join(os.getcwd(), filename)
	txt_file = open(outfile, "w")

	for key in output:
		txt_file.write(str(key)+"\n")
		for k in output[key]:
			for i in range(len(k)):
				txt_file.write(k[i])
				if (i < (len(k)-1)):
					txt_file.write(",")
			txt_file.write("|")

			for i in range(len(combo[k])):
				txt_file.write(combo[k][i])
				if (i < (len(combo[k])-1)):
					txt_file.write(",")
			txt_file.write("\n")
			
	txt_file.close()

def Strip(inputfile, input_nodes):
	input_nodes = input_nodes.strip()
	select_nodes = input_nodes.split(",")

	print ("\n--------------------------------------------------------------------------------")
	print ("Processing file "+inputfile+" for nodes ("+input_nodes+") :")
	print ("--------------------------------------------------------------------------------")

	nodes_set = []
	open(inputfile, "r")
	with open(inputfile, "r") as sif_file:
		for line in sif_file:
			for i in range(len(select_nodes)):
				# if the line read from file does not contain
				# combi:all_all proceed
				if ("combi:all_all" not in line):
					# if the line does contain the node from selected nodes
					if (select_nodes[i] in line):
						line = line.strip()
						nodes_set.append(line)

	selected_nodes = []
	for i in range(len(select_nodes)):
		nodes = select_nodes[i].split(":")
		selected_nodes.append(nodes[1])

	print ("\n--------------------------------------------------------------------------------")
	print ("List of first neighbors for nodes ("+input_nodes+"):")
	print ("--------------------------------------------------------------------------------")

	#--------------------------------------------
	# Get the different types of edges and loads
	# into an aray

	# Ignores combi:all_all
	#--------------------------------------------

	edges_list = []
	counter = 0
	for i in range(len(nodes_set)):
		print (nodes_set[i])
		components = nodes_set[i].split(" ")
		# components[0] is source
		# components[1] is edge
		# components[2] is target
		if (components[1] not in "combi:all_all"):
			if (components[1] not in edges_list):
				edges_list.append(components[1])
				counter = counter + 1

	print ("\n--------------------------------------------------------------------------------")
	print ("List of distinct edge types for nodes ("+input_nodes+"):")
	print ("--------------------------------------------------------------------------------")

	for i in range(len(edges_list)):
		print (edges_list[i])

	print ("\n--------------------------------------------------------------------------------")
	print ("Generating the output file for selected nodes ("+input_nodes+"):")
	print ("--------------------------------------------------------------------------------")

	combo = combine(edges_list)

	for node in nodes_set:
		components = node.split(" ")
		source = components[0]
		edge = components[1]
		target = components[2]

		for c in combo:
			if (edge in c):
				# check if the source node is in our selection list
				if (source not in combo[c]):
					combo[c].append(source)

				# check if the target node is in our selection list
				if (target not in combo[c]):
					combo[c].append(target)

	# sort 
	for c in combo:
		combo[c].sort()

	# # save to file
	output = {}
	for key in combo:
		l = len(key)
		if (l not in output):
			output[l] = []
		
		output[l].append(key)

	# save to txt file
	saveToTXT(combo, output, inputfile, selected_nodes)

	# save to csv file
	# saveToCSV(combo, output, inputfile, selected_nodes)

	print ("\nSuccess!")

if __name__ == "__main__":
	if (len(sys.argv) < 3):
		print ("Usage: filename_to_process node1,node1,...,node_n")
		print ("Example: pdb1v0y_h.sif A:340:_:SER")
		print ("Example: pdb1v0y_h.sif A:340:_:SER,A:334:_:LYS")
		exit()


	# get input parameter from command line
	inputfile = sys.argv[1]
	input_nodes = sys.argv[2]

	Strip(inputfile, input_nodes)
