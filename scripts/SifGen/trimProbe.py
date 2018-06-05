'''
Most parts of this script are reused from RINerator program
http://rinalyzer.de/rinerator.php
'''

import Bio.PDB
import os.path
import sys 

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
'TYQ': 'Y', 'TYS': 'Y', 'TYY': 'Y', 'AGM': 'R', 'GL3': 'G', 'SMC': 'C', 'ASX': 'B', 'CGU': 'E', 'CSX': 'C', 'GLX': 'Z',
'O': 'W'}

mc_atoms_dic = {'N': '', 'CA': '', 'C': '', 'O': '', 'H': '', 'HA': '', 'OXT': '', 'HA2': '', 'HA3': '', 'H?': '', 'W': ''}

def set_str_sel(sel_id,pdb_file_bname, pdb_path, selection_lst):
    """
    Set the structure selection
    selection_lst = [[comp_name, comp_type, res_ranges],...]
    """
    pdb_file = os.path.join(pdb_path,pdb_file_bname)
    stsel_obj = StSel(sel_id,pdb_file)
    stsel_obj.get_pdb()
    for component_desc in selection_lst:
        [comp_name, comp_type, res_ranges] = component_desc
        stsel_obj.set_component(comp_name, comp_type, res_ranges)
    stsel_obj.set_res_id_lst()
    
    
    stsel_obj.set_subid_dic()
    print ("INFO: Read PDB %s for selection %s" % (pdb_file_bname, sel_id))
    return stsel_obj


class StSel:
    """
    Structure Selection
    """
    
    
    def __init__(self, sel_id,pdb_file):
        self.sel_id = sel_id
        self.pdb_file = pdb_file
        self.st_obj = None   # Bio.PDB structure object
        self.seq_obj = None  # Sequence object
        self.protein_names_lst = [] # list of names of selected proteins
        self.protein_sel_dic = {}   # dict of selected proteins: dic[comp_name] = mcomp_dic
        self.ligand_names_lst = [] # list of names of selected ligands
        self.ligand_sel_dic = {}  # dict of selected ligands: dic[comp_name] = mcomp_dic
        # mcomp_dic['comp_name'] = comp_name ,label
        # mcomp_dic['comp_type'] = comp_type , protein, ligand, metal...
        # mcomp_dic['res_ranges'] = res_ranges , res_ranges = [range1,range2,...], range = [model, chain, [st_res_id_st, end_res_id]]
        # mcomp_dic['comp_res_id_lst'] = [res_id, res_id,...]
        self.res_id_dic = {} # dict of res_id dic[res_id] = [comp_name,comp_type]
        self.res_id_lst = [] # sorted lsit of res_id, according to order of components in protein_names_lst, and then in ligand_names_lst
        self.res_subid_dic = {} # dict of subset of residue identifiers, to interpret whatif data, ignores model and hetero flag
        #                      dic[(chain_id,res_seq,icode,resname)] = res_id
        self.atom_label_dic = {}  # dict of atom_label dic[atom_label] = [comp_name, comp_type] atom_label = (atom_id,atom_altloc)
        self.atom_subid_dic = {} #  dict of subset of atom identifiers, to interpret whatif data, ignores model hetero flag alternative_loc
        #                      dic[(chain_id,res_seq,icode,resname,atom_name)] = atom_label
    
    
    def get_pdb(self):
        parser_obj = Bio.PDB.PDBParser()
        self.st_obj = parser_obj.get_structure(self.sel_id, self.pdb_file)
    
    
    def set_component(self,comp_name, comp_type,res_ranges):
        mcomp_dic = self.set_mcomp(comp_name,comp_type)
        mcomp_dic = self.set_res(mcomp_dic,res_ranges)
        if comp_type == 'protein':
            self.protein_names_lst.append(comp_name)
            self.protein_sel_dic[comp_name] = mcomp_dic
        else:
            self.ligand_names_lst.append(comp_name)
            self.ligand_sel_dic[comp_name] = mcomp_dic
    
    
    def set_res_id_lst(self):
        """
        set a ordered list of all residues in selection
        first the protein components ordered in protein_names_lst
        then the ligand components ordered in ligand_names_lst
        """
        self.res_id_lst = []
        comp_names_lst = []
        for protein_name in self.protein_names_lst:
            mcomp_dic = self.protein_sel_dic[protein_name]
            comp_res_id_lst = self.get_comp_res_id_lst(mcomp_dic)
            for res_id in comp_res_id_lst:
                self.res_id_lst.append(res_id)
        for ligand_name in self.ligand_names_lst:
            mcomp_dic = self.ligand_sel_dic[ligand_name]
            comp_res_id_lst = self.get_comp_res_id_lst(mcomp_dic)
            for res_id in comp_res_id_lst:
                self.res_id_lst.append(res_id)
    
    
    def set_subid_dic(self ):
        """
        set a dict of subset of residue identifiers and a dict of subset of atom identifiers
        only needed to interpret data from whatif
        """
        self.res_subid_dic={}
        self.atom_subid_dic = {}

        # print "res_is_lst\n"
        # for i in self.res_id_lst:
        #     print i

        for res_id in self.res_id_lst:
            chain_id = self.get_chain_from_res_id(res_id)
            res_seq = self.get_res_seq_from_res_id(res_id)
            icode = self.get_icode_from_res_id(res_id)
            resname = self.get_resname_from_res_id(res_id)
            dic_key = self.get_res_id_key(chain_id,res_seq,icode,resname)
            # res_subid_dic
            if dic_key in self.res_subid_dic:
                print ("WARNING: Cannot identify residues. Repeated residue identifier: chain %s res_seq %s icode %s resname %s" % (chain_id,res_seq,icode,resname) )
            else:
                self.res_subid_dic[dic_key] = res_id
            # atom_subid_dic
            res_obj = self.res_obj_from_res_id(res_id)
            atom_lst = res_obj.get_list()
            # print ("atom list ", res_id, atom_lst)
            for atom_obj in atom_lst:
                atom_name = atom_obj.get_name()
                altloc=atom_obj.get_altloc()
                # print ("atom_name, altloc ",atom_name, altloc)
                dic_key = self.get_atom_subid_key(chain_id,res_seq,icode,resname,atom_name,altloc)
                if dic_key in self.atom_subid_dic:
                    print ("WARNING: Cannot identify atoms. Repeated atom identifier: chain %s res_seq %s icode %s resname %s atom_name %s altloc" % (chain_id,res_seq,icode,resname,atom_name,altloc))
                else:
                    atom_label=self.get_atom_label(atom_obj)
                    self.atom_subid_dic[dic_key] =atom_label

        # print "atom_subid_dic\n"
        # for i in self.atom_subid_dic:
        #     print (i, self.atom_subid_dic[i])
    
    
    def set_mcomp(self,comp_name,comp_type):
        """
        Initialise molecular component dictionary
        Can be a protein or ligand...
        """
        mcomp_dic = {}
        mcomp_dic['comp_name'] = comp_name # label
        mcomp_dic['comp_type'] = comp_type # protein, ligand, metal...
        mcomp_dic['res_ranges'] = {}  # res_ranges , res_ranges = [range1,range2,...], range = [model, chain, [st_res_id_st, end_res_id]]
        mcomp_dic['comp_res_id_lst'] = {} # [res_id, res_id,...]
        return mcomp_dic
    
    def set_mcomp_res_ranges(self,mcomp_dic,res_ranges):
        mcomp_dic['res_ranges'] = res_ranges
        return mcomp_dic
    
    def set_mcomp_comp_res_id_lst(self,mcomp_dic,comp_res_id_lst):
        mcomp_dic['comp_res_id_lst'] = comp_res_id_lst
        return mcomp_dic
    
    def get_comp_name(self,mcomp_dic):
        return mcomp_dic['comp_name']
    
    def get_comp_type(self,mcomp_dic):
        return mcomp_dic['comp_type']
    
    def get_comp_res_id_lst(self,mcomp_dic):
        return mcomp_dic['comp_res_id_lst']
    
    def set_res(self,mcomp_dic,res_ranges):
        comp_res_id_lst = []
        # use the model of first segment, models of other segments ignored
        model_id = self.get_model_id_from_res_ranges(res_ranges)
        for seg_id in res_ranges:
            chain_id = self.get_chain_id_from_seg_id(seg_id)
            chain_obj = self.st_obj[model_id][chain_id]
            st_flag = False
            res_lst = chain_obj.get_list()
            for res_obj in res_lst:
                if self.check_aa(res_obj) == False and self.get_comp_type(mcomp_dic) == 'protein':
                    continue
                res_id = res_obj.get_full_id()
                if st_flag == False:
                    if self.check_st_res(res_id, seg_id) == True:
                        comp_res_id_lst.append(res_id)
                        self.store_id_dic(res_id,mcomp_dic)
                        if self.get_end_res_seq_from_seg_id(seg_id) == None:
                            break
                        st_flag = True
                else:
                    comp_res_id_lst.append(res_id)
                    self.store_id_dic(res_id,mcomp_dic)
                    if self.check_end_res(res_id, seg_id) == True:
                        break
        mcomp_dic = self.set_mcomp_res_ranges(mcomp_dic,res_ranges)
        mcomp_dic = self.set_mcomp_comp_res_id_lst(mcomp_dic,comp_res_id_lst)
        return mcomp_dic
    
    
    def check_aa(self, res_obj):
        hetf = res_obj.get_id()[0]
        if hetf == ' ' or hetf[:2] == 'H_':
            res_name = res_obj.get_resname()
            if res_name in aa_trans_dic:
                return True
        print ("INFO: Non amino acid: %s" % str(res_obj.get_full_id()))
        return False
    
    
    def check_st_res(self,res_id, seg_id):
        hetf = self.get_hetf_from_res_id(res_id)
        res_seq = self.get_res_seq_from_res_id(res_id)
        icode = self.get_icode_from_res_id(res_id)
        st_hetf = self.get_st_hetf_from_seg_id(seg_id)
        st_res_seq = self.get_st_res_seq_from_seg_id(seg_id)
        st_icode = self.get_st_icode_from_seg_id(seg_id)
        if st_res_seq == '_':
            return True
        if res_seq == st_res_seq:
            if icode == st_icode:
                if hetf == st_hetf:
                    return True
        return False
    
    
    def check_end_res(self,res_id, seg_id):
        hetf = self.get_hetf_from_res_id(res_id)
        res_seq = self.get_res_seq_from_res_id(res_id)
        icode = self.get_icode_from_res_id(res_id)
        end_res_seq = self.get_end_res_seq_from_seg_id(seg_id)
        end_icode = self.get_end_icode_from_seg_id(seg_id)
        end_hetf = self.get_end_hetf_from_seg_id(seg_id)
        if end_res_seq == '_':
            return False
        if res_seq == end_res_seq:
            if icode == end_icode:
                if hetf == end_hetf:
                    return True
        return False
    
    
    def store_id_dic(self,res_id, mcomp_dic):
        comp_name=self.get_comp_name(mcomp_dic)
        comp_type=self.get_comp_type(mcomp_dic)
        self.store_res_id_dic(res_id,comp_name,comp_type)
        self.store_atom_label_dic(res_id,comp_name,comp_type)
    
    
    def store_res_id_dic(self,res_id,comp_name,comp_type):
        self.res_id_dic[res_id]=[comp_name,comp_type]
    
    
    def store_atom_label_dic(self,res_id,comp_name,comp_type):
        res_obj = self.res_obj_from_res_id(res_id)
        atom_lst = res_obj.get_list()
        for atom_obj in atom_lst:
            atom_label = self.get_atom_label(atom_obj)
            self.atom_label_dic[atom_label]=[comp_name,comp_type]
    
    
    def res_obj_from_res_id(self,res_id):
        res_obj = self.st_obj[res_id[1]][res_id[2]][res_id[3]]
        return res_obj
    
    
    def get_model_id_from_res_ranges(self,res_ranges):
        """
        get model id from ranges
        use the model of first segment, models of other segments ignored
        """
        model_id = res_ranges[0][0]
        return model_id
    
    
    def get_res_id_key(self,chain_id,res_seq,icode,resname):
        dic_key = (chain_id,res_seq,icode,resname)
        return dic_key
    
    
    def get_atom_subid_key(self,chain_id,res_seq,icode,resname,atom_name,altloc):
        subid= (chain_id,res_seq,icode,resname,atom_name,altloc)
        # print ("atom name", atom_name, "subid", subid)
        return subid
    
    
    def get_atom_label_from_subid(self,chain_id,res_seq,icode,resname,atom_name,altloc):
        atom_key = self.get_atom_subid_key(chain_id,res_seq,icode,resname,atom_name,altloc)
        if atom_key in self.atom_subid_dic:
            atom_label = self.atom_subid_dic[atom_key]
            # print (atom_label)
        else:
            # print (atom_key)
            # print ("ERROR: No atom with chain_id %s res_seq %s icode %s resname %s atom_name %s" % (chain_id,str(res_seq),icode,resname,atom_name))
            # raw_input()
            atom_label = None
        return atom_label
    
    
    def get_chain_id_from_seg_id(self,seg_id):
        chain_id = seg_id[1]
        return chain_id
    
    
    def get_st_hetf_from_seg_id(self,seg_id):
        return seg_id[2][0]
    
    def get_st_res_seq_from_seg_id(self,seg_id):
        return seg_id[2][1]
    
    def get_st_icode_from_seg_id(self,seg_id):
        return seg_id[2][2]
    
    
    def get_end_hetf_from_seg_id(self,seg_id):
        return seg_id[3][0]
    
    def get_end_res_seq_from_seg_id(self,seg_id):
        return seg_id[3][1]
    
    def get_end_icode_from_seg_id(self,seg_id):
        return seg_id[3][2]
    
    
    def get_chain_from_res_id(self,res_id):
        return res_id[2][0]
    
    def get_hetf_from_res_id(self,res_id):
        return res_id[3][0]

    def get_res_seq_from_res_id(self,res_id):
        return res_id[3][1]
    
    def get_icode_from_res_id(self,res_id):
        return res_id[3][2]
    
    def get_resname_from_res_id(self,res_id):
        res_obj = self.st_obj[res_id[1]][res_id[2]][res_id[3]]
        resname = res_obj.get_resname()
        return resname
    
    
    def get_atom_label(self,atom_obj):
        atom_id = atom_obj.get_full_id()
        atom_altloc = atom_obj.get_altloc()
        return (atom_id, atom_altloc)
    
    
    def get_atom_id_from_label(self,atom_label):
        return atom_label[0]
    
    
    def get_atom_name_from_label(self,atom_label):
        return atom_label[0][4][0]
    
    
    def check_mc_sc_ligand_from_label(self,atom_label):
        """
        return:
        mc
        sc
        ligand
        """
        comp_type = self.get_comp_type_from_label(atom_label)
        if comp_type == None:
            return None
        #njd
        if comp_type == 'water':
            return 'solvent'
        elif comp_type != 'protein':
            return 'ligand'
        atom_name = self.get_atom_name_from_label(atom_label)
        if atom_name in mc_atoms_dic:
            return 'mc'
        return 'sc'
    
    
    def get_comp_type_from_label(self,atom_label):
        if atom_label in self.atom_label_dic:
            return self.atom_label_dic[atom_label][1]
        else:
            return None
    
    
    def atom_label_to_res_id(self,atom_label):
        atom_id = self.get_atom_id_from_label(atom_label)
        res_id = self.get_res_id_from_atom_id(atom_id)
        return res_id
    
    
    def get_res_id_from_atom_id(self,atom_id):
        return atom_id[:-1]
    
    
    # def get_seq(self):
    #     char_lst = []
    #     for res_id in self.res_id_lst:
    #         res_obj = self.res_obj_from_res_id(res_id)
    #         resname = res_obj.get_resname()
    #         res_code = aa_trans_dic[resname]
    #         char_lst.append(res_code)
    #     sep = ''
    #     res_seq_str = sep.join(char_lst)
    #     seq_id = self.sel_id
    #     seq_obj = sequence.Sequence()
    #     seq_obj.set_seq(seq_id,res_seq_str)
    #     self.seq_obj = seq_obj
    
    
    def get_res_label_from_seq_idx(self,seq_idx):
        res_id = self.res_id_lst[seq_idx]
        chain_id=self.get_chain_from_res_id(res_id)
        res_seq=self.get_res_seq_from_res_id(res_id)
        icode = self.get_icode_from_res_id(res_id)
        if chain_id == ' ':
            chain_id = '_'
        if icode == ' ':
            icode = '_'
        res_label = "%s:%s:%s" % (chain_id,str(res_seq),icode)
        return res_label

def atom_desc_to_atom_label(atom_desc,stsel_obj):
    chain_id = atom_desc[1:2]
    res_seq = int(atom_desc[2:6])
    icode = atom_desc[6:7]
    resname = atom_desc[7:10]
    atom_name = atom_desc[11:15].strip()
    altloc = atom_desc[15:16]
    atom_label=stsel_obj.get_atom_label_from_subid(chain_id,res_seq,icode,resname,atom_name,altloc)
    return atom_label

def trim_Probe(lines, stsel_obj):
    probes = []
    for line in lines:
        line_lst = line.split(':')
        atom_desc1 = line_lst[3]
        atom_desc2 = line_lst[4]
        atom_label1 = atom_desc_to_atom_label(atom_desc1,stsel_obj)
        atom_label2 = atom_desc_to_atom_label(atom_desc2,stsel_obj)
        if (atom_label1 == None) or (atom_label2 == None):
            continue
        else:
            probes.append(line)

    return probes

# if __name__ == "__main__":
# 	sel_id = "3bwm_lig"
# 	pdb_file_name = '3bwm_prot.pdb'
# 	pdb_path = 'new'

# 	component1 = ['3BWM_h', 'protein', [[0,'A',[' ','_',' '],[' ','_',' ']]]]
# 	component2 = ['HOH', 'water', [[0,'A',['W',401,' '],['W',510,' ']]]]
# 	component3 = ['SAM', 'ligand', [[0,'A',['H_SAM',301,' '],['H_DNC',302,' ']]]]

# 	selection_lst = [component1, component2, component3]

# 	stsel_obj = set_str_sel(sel_id,pdb_file_name,pdb_path,selection_lst)

# 	probe_file_name = "new/3bwm_prot_h.probe"
# 	with open(probe_file_name, "r") as f:
# 		lines = f.readlines()

# 	probe = trim_Probe(lines)

# 	with open("new_3bwm_probes-prot.probe", "w") as pf:
# 		for i in probes:
# 			pf.write(i)

