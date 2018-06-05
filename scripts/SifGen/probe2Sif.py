import os, sys, re
from res_mcatoms_dic import *

def get_side(resID, atom):
	if resID in aa_trans_dic.keys():
		if atom in mc_atoms_dic.keys():
			side = 'mc'
		else:
			side = 'sc'
	else:
		if resID == 'HOH':
			side = 'solvent'
		else:
			side = 'ligand'

	return side

def get_inttype(c):
	# 'wc':wide contact,'cc': close contact,'so':small overlap
	# 'bo':big overlap
	# 'hb':hydrogen bond
	if c in ['wc','cc']:
		action = 'cnt'
	elif c in ['bo','so']:
		action = 'ovl'
	elif c == 'hb':
		action = 'hbond'
	else:
		print ('Cannot find interaction type!')
		return "None"

	return action

def set_mc_sc_ligand(side1, side2):
	# order: mc sc ligand solvent
	if side1==None or side2==None:
		return None, None
	elif side1 == 'mc':
		return side1, side2
	elif side1 == 'sc':
		if side2 == 'mc':
			return side2, side1
		else:
			return side1, side2
	elif side1 == 'ligand':
		if (side2 in ['mc', 'sc']):
			return side2, side1
		else:
			return side1, side2
	elif side1 == 'solvent':
		return side2, side1

def probe2Sif(probefile):

	with open(probefile,'r') as f:
		lines = f.readlines()

	res_parts = {}
	res_parts1 = {}

	for i in range(len(lines)):
	    c = lines[i].split(':')
	    #c[1] : '1->1'
	    #c[2] : 'wc':wide contact,'cc': close contact,'so':small overlap,'bo':big overlap,'hb':hydrogen bond
	    #c[3] : ' A  38 TRP  C   '
	    #c[4] : ' A  38 TRP  CD2 '
	    #c[5],c[6]: '0.025', '0.335'
	    #c[7]-c[9]: '-10.370', '-17.062', '-19.661'
	    #c[10]: '0.000'   score
	    #c[11]: '0.0104'  raw_score in rinalyzer
	    #c[12]: 'C'
	    #c[13]: 'C'
	    #c[14]-c[16]:'-10.370', '-17.062', '-19.661'
	    #c[17],c[18]: '30.32', '26.49

	    if c[3][1:2] != 'A' or c[4][1:2] != 'A': continue
	    if int(c[3][2:6]) == int(c[4][2:6]): continue

	    res1  = c[3][7:10]
	    atom1 = c[3][11:15].strip()
	    side1 = get_side(res1, atom1)
	    # altloc1 = c[3][15:16]
	    part1 = c[3][1:2]+':'+c[3][2:6].strip()+':'+'_'+':'+c[3][7:10] 

	    res2 = c[4][7:10]
	    atom2 = c[4][11:15].strip()
	    side2 = get_side(res2, atom2)
	    altloc2 = c[4][15:16]
	    part2 = c[4][1:2]+':'+c[4][2:6].strip()+':'+'_'+':'+c[4][7:10] 

	    action = get_inttype(c[2])

	    sideL, sideR = set_mc_sc_ligand(side1, side2)
	    act = action+':'+sideL+'_'+sideR
	    
	    if (part1,act,part2) not in res_parts.keys() and (part2,act,part1) in res_parts.keys():
	    	key = (part2,act,part1)

	    	if side1 not in res_parts[key][1]:
	    		res_parts[key][1].append(side1)
	    	if side2 not in res_parts[key][0]:
	    		res_parts[key][0].append(side2)
	    	key1 = (part2,part1)

	    elif (part1,act,part2) not in res_parts.keys() and (part2,act,part1) not in res_parts.keys():
	    	key = (part1,act,part2)
	    	res_parts[key] = [[side1],[side2]]
	    	key1 = (part1,part2)

	    else:
	    	if side1 not in res_parts[key][0]:
	    		res_parts[key][0].append(side1)
	    	if side2 not in res_parts[key][1]:
	    		res_parts[key][1].append(side2)
	    	key1 = (part1,part2)

	
	interaction = {}
	combi_interaction = []

	for key in res_parts.keys():
		# if 'solvent' in key[1]: continue
		# if "hb" in key[1]: continue

		if key not in interaction:
			interaction[key] = res_parts[key]

			
		k1 = [key[0], key[2]]
		k2 = [key[2], key[0]]
		if (k1 not in combi_interaction) and (k2 not in combi_interaction):
			key1 = (key[0],'combi:all_all',key[2])
			interaction[key1] = []
			combi_interaction.append(k1)

	
	f = open(probefile.replace(".probe", ".sif"), "w")
	for key in interaction.keys():
		if ('combi' not in key[1]):
			act = key[0]+' '+key[1]+' '+key[2]
			f.write(act+"\n")
			# print (act)

	for key in interaction.keys():
		if ('combi' in key[1]):
			act = key[0]+' '+key[1]+' '+key[2]
			f.write(act+"\n")
			# print (act)

	f.close()

if __name__ == "__main__":
	if (len(sys.argv) < 2):
		print ("Usage: probe_file")
		exit()

	probe2Sif(sys.argv[1])

	

	
