"""
This is a program written by qianyi cheng in deyonker research group
at university of memphis.
"""
import os, sys, re
from copy import *
from numpy import *
import argparse


aa_trans_dic = {
'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M', 'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H', 'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G',
'2AS': 'D', '3AH': 'H', '5HP': 'E', 'ACL': 'R', 'AIB': 'A', 'ALM': 'A', 'ALO': 'T', 'ALY': 'K', 'ARM': 'R', 'ASA': 'D',
'ASB': 'D', 'ASK': 'D', 'ASL': 'D', 'ASQ': 'D', 'AYA': 'A', 'BCS': 'C', 'BHD': 'D', 'BMT': 'T', 'BNN': 'A', 'BUC': 'C',
'BUG': 'L', 'C5C': 'C', 'C6C': 'C', 'CCS': 'C', 'CEA': 'C', 'CHG': 'A', 'CLE': 'L', 'CME': 'C', 'CSD': 'A', 'CSO': 'C',
'CSP': 'C', 'CSS': 'C', 'CSW': 'C', 'CXM': 'M', 'CY1': 'C', 'CY3': 'C', 'CYG': 'C', 'CYM': 'C', 'CYQ': 'C', 'DAH': 'F',
'DAL': 'A', 'DAR': 'R', 'DAS': 'D', 'DCY': 'C', 'DGL': 'E', 'DGN': 'Q', 'DHA': 'A', 'DHI': 'H', 'DIL': 'I', 'DIV': 'V',
'DLE': 'L', 'DLY': 'K', 'DNP': 'A', 'DPN': 'F', 'DPR': 'P', 'DSN': 'S', 'DSP': 'D', 'DTH': 'T', 'DTR': 'W', 'DTY': 'Y',
'DVA': 'V', 'EFC': 'C', 'FLA': 'A', 'FME': 'M', 'GGL': 'E', 'GLZ': 'G', 'GMA': 'E', 'GSC': 'G', 'HAC': 'A', 'HAR': 'R',
'HIC': 'H', 'HIP': 'H', 'HMR': 'R', 'HPQ': 'F', 'HTR': 'W', 'HYP': 'P', 'IIL': 'I', 'IYR': 'Y', 'KCX': 'K', 'LLP': 'K',
'LLY': 'K', 'LTR': 'W', 'LYM': 'K', 'LYZ': 'K', 'MAA': 'A', 'MEN': 'N', 'MHS': 'H', 'MIS': 'S', 'MLE': 'L', 'MPQ': 'G',
'MSA': 'G', 'MSE': 'M', 'MVA': 'V', 'NEM': 'H', 'NEP': 'H', 'NLE': 'L', 'NLN': 'L', 'NLP': 'L', 'NMC': 'G', 'OAS': 'S',
'OCS': 'C', 'OMT': 'M', 'PAQ': 'Y', 'PCA': 'E', 'PEC': 'C', 'PHI': 'F', 'PHL': 'F', 'PR3': 'C', 'PRR': 'A', 'PTR': 'Y',
'SAC': 'S', 'SAR': 'G', 'SCH': 'C', 'SCS': 'C', 'SCY': 'C', 'SEL': 'S', 'SEP': 'S', 'SET': 'S', 'SHC': 'C', 'SHR': 'K',
'SOC': 'C', 'STY': 'Y', 'SVA': 'S', 'TIH': 'A', 'TPL': 'W', 'TPO': 'T', 'TPQ': 'A', 'TRG': 'K', 'TRO': 'W', 'TYB': 'Y',
'TYQ': 'Y', 'TYS': 'Y', 'TYY': 'Y', 'AGM': 'R', 'GL3': 'G', 'SMC': 'C', 'ASX': 'B', 'CGU': 'E', 'CSX': 'C', 'GLX': 'Z'
}

mc_atoms_dic = {'N': '', 'CA': '', 'C': '', 'O': '', 'H': '', 'HA': '', 'OXT': '', 'HA2': '', 'HA3': '', 'H?': '', 'W': ''}

def get_res_type(resID, atom):
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
        print('Cannot find interaction type!')
        return "None"
    return action

def set_mc_sc_ligand(side1, side2):
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

def get_side(c):
    cha    = c[:2].strip()          # chain 
    res_id = c[2:6].strip()         # res_id 
    res_nm = c[6:10].strip()        # res_name
    atom   = c[10:-1].strip()       # atom_name
    if c[-1] == ' ':
        spe_c = 'A'
    else:
        spe_c = c[-1]
    return cha, res_id, res_nm, atom, spe_c

def check_repeat(a, a_list):
    if a not in a_list:
        a_list.append(a)
    return a_list

def check_dict_repeat(key,a,dict):
    try:
        dict[key].append(a)
    except:
        dict[key] = [a]
    return dict

def probe_analysis(probefile,sel_res):
    ### sel_res format A:300,A:301,A:302 ###
    sel_res_list = sel_res.split(',')
    res_list = deepcopy(sel_res_list)

    res_acts = {}
    actions  = []
    siflines = {}
    res_atoms = {}

    with open(probefile,'r') as f:
        lines = f.readlines()

    for line in lines:
        c = line.split(':')
        obj1 = c[1]
        acts = c[2]

        cha1, res_id1, res_nm1, atom1, spe_c1 = get_side(c[3])           
        cha2, res_id2, res_nm2, atom2, spe_c2 = get_side(c[4])           

        if cha1 == cha2 and res_id1 == res_id2: continue
        key1 = cha1+':'+res_id1
        key2 = cha2+':'+res_id2
        if key1 in sel_res_list or key2 in sel_res_list:
            res_list = check_repeat(key1,res_list)
            res_list = check_repeat(key2,res_list)
            side1 = get_res_type(res_nm1,atom1)
            side2 = get_res_type(res_nm2,atom2)
            acts1 = acts+':'+side1+'_'+side2
            acts2 = acts+':'+side2+'_'+side1
            actions = check_repeat(acts1,actions)
            actions = check_repeat(acts2,actions)
            if key1 not in res_acts.keys():
                res_acts[key1] = [acts1]
            else:
                res_acts[key1] = check_repeat(acts1,res_acts[key1])
            if key2 not in res_acts.keys():
                res_acts[key2] = [acts2]
            else:
                res_acts[key2] = check_repeat(acts2,res_acts[key2])
            if cha1 < cha2 or ( cha1 == cha2 and int(res_id1) < int(res_id2) ):
                order = key1+' '+acts1+' '+acts2+' '+key2
            else:
                order = key2+' '+acts2+' '+acts1+' '+key1
            if order not in siflines.keys():
                siflines[order] = 1
            else:
                siflines[order] += 1
            res_atoms = check_dict_repeat(key1,atom1,res_atoms)
            res_atoms = check_dict_repeat(key2,atom2,res_atoms)
    return sorted(res_list), res_acts, actions, siflines, res_atoms            


def order_reslist(res_list):
    res_dic = {}
    for i in res_list:
        c = i.split(':')
        try:
            res_dic[c[0]].append(int(c[1]))
        except:
            res_dic[c[0]] = [int(c[1])]
    return res_dic

def print_list(res_list):
    for i in res_list:
        print(i)


def print_dict(dict):
    for i in sorted(dict.keys()):
        print(i, dict[i])


def write_res_freq(res_list, res_acts, res_atoms):
    f_res = open('freq_per_res.dat','w')
    for k in sorted(res_atoms, key=lambda k:len(res_atoms[k]),reverse=True):
        cha,res = k.split(':')
        f_res.write('%-4s %-8s %-8d'%(cha,res,len(res_atoms[k])))
        for act in res_acts[k]:
            f_res.write(' %-20s'%act)
        f_res.write('\n')
    f_res.close()


def write_res_atom(res_dict,res_atoms):
    f_res = open('res_atoms.dat','w')
    for key in sorted(res_dict.keys()):
        for i in sorted(res_dict[key]):
            f_res.write('%-4s %-8d'%(key,i))
            akey = key+':'+str(i)
            for act in unique(res_atoms[akey]):
                f_res.write(' %-4s'%act)
            f_res.write('\n')
    f_res.close()

def write_rin(actions):
    f_act = open('rin_list.dat','w')
    for i in actions:
        f_act.write('%-30s\n'%i)
    f_act.close()

def write_sif(siflines):
    probef_n = probefile.split('/')[-1] 
    f_sif = open(probef_n.replace(".probe",".sif"),'w')
    for key in sorted(siflines.keys()):
        f_sif.write('%-60s %10d\n'%(key,siflines[key]))
    f_sif.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate interaction information from probe file.\
            Usage: probe2rins.py -p probe_file -s seed')
    parser.add_argument('-f',dest='probefile',default=None,help='probe_file')
    parser.add_argument('-s',dest='seed',default=None,help='seed for select RIN, in the format of "A:300,A:301,A:302"')

    args = parser.parse_args()
    probefile = args.probefile
    sel_res = args.seed

    res_list, res_acts, actions, siflines, res_atoms = probe_analysis(probefile,sel_res)
    res_dict = order_reslist(res_list)

    write_sif(siflines)
    write_rin(actions)
    write_res_atom(res_dict,res_atoms)
    write_res_freq(res_list,res_acts,res_atoms)
