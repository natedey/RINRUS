from numpy import *

Res31 ={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
        'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
        'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
        'TRP':'W','TYR':'Y','ASX':'N','GLX':'Q','UNK':'X','INI':'K',
        'AAR':'R','ACE':'X','ACY':'G','AEI':'T','AGM':'R','ASQ':'D',
        'AYA':'A','BHD':'D','CAS':'C','CAY':'C','CEA':'C','CGU':'E',
        'CME':'C','CMT':'C','CSB':'C','CSD':'C','CSE':'C','CSO':'C',
        'CSP':'C','CSS':'C','CSW':'C','CSX':'C','CXM':'M','CYG':'C',
        'CYM':'C','DOH':'D','EHP':'F','FME':'M','FTR':'W','GL3':'G',
        'H2P':'H','HIC':'H','HIP':'H','HTR':'W','HYP':'P','KCX':'K',
        'LLP':'K','LLY':'K','LYZ':'K','M3L':'K','MEN':'N','MGN':'Q',
        'MHO':'M','MHS':'H','MIS':'S','MLY':'K','MLZ':'K','MSE':'M',
        'NEP':'H','NPH':'C','OCS':'C','OCY':'C','OMT':'M','OPR':'R',
        'PAQ':'Y','PCA':'Q','PHD':'D','PRS':'P','PTH':'Y','PYX':'C',
        'SEP':'S','SMC':'C','SME':'M','SNC':'C','SNN':'D','SVA':'S',
        'TPO':'T','TPQ':'Y','TRF':'W','TRN':'W','TRO':'W','TYI':'Y',
        'TYN':'Y','TYQ':'Y','TYS':'Y','TYY':'Y','YOF':'Y','FOR':'X',
        '---':'-','PTR':'Y','LCX':'K','SEC':'D','MCL':'K','LDH':'K',
        'HIE':'H',}

Res13 = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS',
         'Q':'GLN', 'E':'GLU', 'G':'GLY', 'H':'HIS', 'I':'ILE',
         'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO',
         'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL','X':'UNK'}


def read_pdb(pdbfile,TER=False):
    f = open(pdbfile,'r')
    pdb = []
    res_info = []
    for line in f:
        record = line[:6]
        if ( record != 'ATOM  ' and record != 'HETATM' ): continue
        serial = int( line[6:11] )
        atomname = line[12:16]
        altloc = line[16]
        resname = line[17:20]
        chain = line[21]
        resnum = int( line[22:26] )
        achar = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        try:
            occ = float( line[54:60] )
        except:
            occ = 1.0
        try:
            tfactor = float( line[60:66] )
        except:
            tfactor = 1.0
        try:
            segid = line[72:76]
        except:
            segid = ''
        try:
            elsymbol = line[76:78]
        except:
            elsymbol = ''
        try:
            charge = line[78:80]
        except:
            charge = '0 '
        # 0:  record
        # 1:  serial
        # 2:  atomname
        # 3:  altloc
        # 4:  resname
        # 5:  chain
        # 6:  resnum
        # 7:  achar
        # 8:  x
        # 9:  y
        # 10: z
        # 11: occ
        # 12: tfactor
        # 13: segid
        # 14: elsymbol
        # 15: charge
        try:
            fix = line[85:87]
            if int(fix) == -1:
                res_info.append([atomname, chain, resnum])
        except:
            fix = " 0"
            
        pdb.append( [ record, serial, atomname, altloc, resname, chain, resnum, achar, x, y, z, occ, tfactor, segid, elsymbol, charge, fix ] )
        #     1          0       1       2       3       4           5   6       7      8  9  10 11   12         13      14      15      16
    f.close()
    return pdb, res_info



def write_pdb(filename,pdb,renum_atom=True,hydrogen=True,renum_res=False):
    if ( isinstance(filename,file) ):
        f = filename
        file_opened = False
    else:
        f = open(filename,'w')
        file_opened = True
    serial = 1
    serial_res = 1
    prev_res = pdb[0][6]
    for p in pdb:
        t = p[:] # copy
        if p[0] == 'TER':
            f.write('TER\n')
            continue
        if ( hydrogen == False and p[2].strip()[0] == 'H' ): continue
        if ( renum_atom ):
            t = [t[0],serial] + t[2:17]
            serial += 1
        else:
            t = t[:17]
        if ( renum_res ):
            if ( prev_res != p[6] ):
                serial_res += 1
                prev_res = p[6]
            t = t[0:6]+[serial_res]+t[7:17]
        f.write('%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%3s%6s\n'%tuple(t))
    if ( file_opened ):
        f.close()


def get_ca(pdb,unique_residue=True):
    pdb_ca = []
    seen_residue = set()
    for p in pdb:
        if len(p) < 2:
            continue
        if ( p[2].strip() == 'CA' ):
            if ( unique_residue ):
                reskey = (p[6],p[7])
                if ( reskey in seen_residue ): continue
                seen_residue.add(reskey) 
            pdb_ca.append(p)
    return pdb_ca

def get_coord(pdb):
    coord = []
    for p in pdb:
        coord.append((p[8],p[9],p[10]))
    return coord

