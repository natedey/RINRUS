import os
import sys
from strip import *
from analyze_models import *
from makePDB import *

if (len(sys.argv) < 4):
	print ("Usage: filename_to_process node1,node1,...,node_n pdb_to_trim")
	print ("Example: pdb1v0y_h.sif A:340:_:SER 1v0y_h.pdb")
	print ("Example: pdb1v0y_h.sif A:340:_:SER,A:334:_:LYS 1v0y_h.pdb")
	exit()


# get input parameter from command line
inputfile = sys.argv[1]
input_nodes = sys.argv[2]
inputfile2 = sys.argv[3]

# strip
Strip(inputfile, input_nodes)

select_nodes = input_nodes.strip().split(",")
selected_nodes = []
for i in range(len(select_nodes)):
	nodes = select_nodes[i].split(":")
	selected_nodes.append(nodes[1])

filename = inputfile.replace(".sif", "-")+"-".join(selected_nodes)+"_py.txt"

# strip analyze
Strip_analyze(filename)

dir_unique_models = filename.replace(".txt", "")
unique_models = os.path.join(dir_unique_models,"unique_sets_freq.txt")

makePDB(unique_models, inputfile2)
