import sys
from trimProbe import *
from probe2Sif import *

if __name__ == "__main__":
	sel_id = "3bwm_lig"
	pdb_file_name = '3bwm_prot.pdb'
	pdb_path = '3bwm'
	probe_file_name = "3bwm/3bwm_prot_h.probe"
	output_probe = "new_3bwm_probes-prot.probe"

	component1 = ['3BWM_h', 'protein', [[0,'A',[' ','_',' '],[' ','_',' ']]]]
	component2 = ['HOH', 'water', [[0,'A',['W',401,' '],['W',510,' ']]]]
	component3 = ['SAM', 'ligand', [[0,'A',['H_SAM',301,' '],['H_DNC',302,' ']]]]

	selection_lst = [component1, component2, component3]

	stsel_obj = set_str_sel(sel_id,pdb_file_name,pdb_path,selection_lst)

	# read probe file
	with open(probe_file_name, "r") as f:
		lines = f.readlines()

	probes = trim_Probe(lines,stsel_obj)

	# write trimmed probe
	with open(output_probe, "w") as pf:
		for i in probes:
			pf.write(i)

	print("Trimmed probe: ", output_probe)

	print("Generating SIF file...\n")
	probe2Sif(output_probe)

	print("SIF file: ", output_probe.replace(".probe", ".sif"))