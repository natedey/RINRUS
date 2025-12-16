"""
This is a program written by Qianyi Cheng 
at university of memphis.
Date 8.10.2022
"""
import os, sys, re, filecmp
from numpy import *
import argparse
from read_write_pdb import *
from glob import glob
from gaussian_basis_dict import *

"""
write gaussian input files using gaussian template file
"""
def write_gau_input(inp_name,inp_temp,charge,multiplicity,pic_atom,tot_charge,res_count,basisinfo,hopt):
    ### inp_name default is 1.inp, but first is name.input such as 9.input
    ### inp_type list which includes [small/large, level_theory, basis, opt, freq,  
    ### input_template line0: nprocshared and mem 
    ### input_template line1: level
    ### input_template line2: opt + calcfc/readfc or opt + modred + info
    ### input_template line3: freq
    ### input_template line4: guess=read
    ### input_template line5: geom=checkpoint
    ### input_template line6: iop
    ### input_template line7: scf_info
    ### input_template line8: scrf_info
    ### after the basis sets there needs to be an empty line in the end    

    print("charge= %d, tot_charge_without'-c'= %d"%(charge,tot_charge))
    with open(inp_temp) as f:
        lines = f.readlines()

    inp = open('%s'%inp_name,'w')
    ### Start from checkpoint file on line 1 ###
    inp.write("%chk=1.chk\n")           #write check file into 1.chk
    v = lines[0].split()                #In line[0] of input_template file are number for nprocshared and mem
    pro_line = '%nprocshared='
    inp.write("%s%d\n"%(pro_line,int(v[0])))
    mem_line = '%mem='
    inp.write("%s%dGB\n"%(mem_line,int(v[1])))
    inp.write("#P ")      

    ### get atom list and work out if ECPs needed
    atomlist = []
    for atom in pic_atom:
        if atom[14].strip() not in atomlist:
            atomlist.append(atom[14].strip())
    ecplist = []
    for i in atomlist:
        if i not in ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar']:
            ecplist.append(i)

    if lines[1][0] != '#':
        method = lines[1].strip().split('/')
        if len(ecplist) > 0:
            inp.write("%s/genecp "%method[0])
        elif len(method) == 1 or method[1] == 'gen':
            inp.write("%s/gen "%method[0])
        else:
            inp.write("%s "%lines[1].strip())
        
    if lines[2][0] != '#':
        optl = lines[2].split()
        inp.write("%s "%optl[0])
        if 'modred' in lines[2]:
            f_atom = []
            modred_info, modred_code = optl[1:3]
            pairs = modred_info.split(';')
            lmod = []
            for pair in pairs:
                ar = pair.split(',')
                lmod.append(int(len(ar)/2))
                for i in range(0,len(ar),2):
                    for atom in pic_atom:
                        if ar[i] == atom[2].strip() and int(ar[i+1]) == atom[6]:
                            f_atom.append(pic_atom.index(atom)+1)
    for l in range(3,9):
        if lines[l][0] != '#':
            inp.write("%s "%lines[l].strip())

    inp.write("\n\n")
    inp.write("%s\n"%res_count)
    inp.write("\n")
    inp.write("%d %d\n"%(charge+tot_charge,multiplicity))

    if lines[5][0] == '#':
        ### For pre-optimization where all heavy atoms kept frozen with "-1" except H atoms
        if hopt == 1:
            for atom in pic_atom:
                if atom[14].strip() == 'H':
                    inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),'0',atom[8],atom[9],atom[10])) 
                else:
                    inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),'-1',atom[8],atom[9],atom[10])) 
        ### Normal optimization
        else:
            for atom in pic_atom:
                inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),atom[16],atom[8],atom[9],atom[10])) 
    inp.write("\n")

    count = 0
    if lines[2][0] != '#' and 'modred' in lines[2]:
        for l in lmod:
            for i in range(l):
                inp.write("%d "%f_atom[count+i])
            inp.write("%s\n"%modred_code)
            count += l
        inp.write("\n")

    ### use basis info from template if specified otherwise use basis dictionary
    if basisinfo == 'intmp':
        if lines[10][0] != '#' and 'basis' in lines[10]:
            for l in lines[11:]:
                inp.write("%s"%l)
            if len(lines) <= 10:
                inp.write('\n')
    else:
        for atom in atomlist:
            basis = basis_dict[atom]
            inp.write(basis+'\n')
        if len(ecplist) > 0:
            inp.write('\n')
            for atom in ecplist:
                ecp = ecp_dict[atom]
                inp.write(ecp+'\n')
        inp.write('\n')

    if lines[8][0] != '#' and 'scrf' in lines[8]:
        inp.write('radii=uff\nalpha=1.2\neps=4.0\n\n')
    
    inp.close()


"""
write qchem input files using qchem template file
"""
def write_qchem_input(inp_name,inp_temp,charge,multiplicity,pic_atom,tot_charge,res_count):
    ### inp_name default can be 1.inp, but first is name.input such as 9.input
    ### From input_temple file ###
    ### jobtype opt/freq ###
    ### method ###
    ### basis    ###
    ### ecp      ###
    ### solvent_method ###

    print("charge= %d, tot_charge_without'-c'= %d"%(charge,tot_charge))
    with open(inp_temp) as f:
        lines = f.readlines()
#    if line.startswith('#'): continue
    for line in lines:
        if line.startswith('#'): continue
        if line.startswith('jobtype'):
            v = line.split(':')
            if "/" in v[1]:
                jobtype = [x.strip() for x in v[1].split('/')]
            else:
                jobtype = v[1].strip()
        if line.startswith('method'):
            method = line.split(':')[1].strip()
        if line.startswith('basis'):
            basis = line.split(':')[1].strip()
        if line.startswith('ecp'):
            ecp = line.split(':')[1].strip()
        if line.startswith('solvent_method'):
            sol = line.split(':')[1].strip()
        if line.startswith('def_basis'):
            start = lines.index(line)
            def_basis = []
            for l in range(start+1,len(lines)):
                def_basis.append(lines[l])
         
    frozen = []
    inp = open('%s'%inp_name,'w')
    inp.write('$molecule\n')
    inp.write("%d %d\n"%(charge+tot_charge,multiplicity))
    for i in range(len(pic_atom)):
        atom = pic_atom[i]
        inp.write("%6s %8.3f %8.3f %8.3f\n"%(atom[14].strip(),atom[8],atom[9],atom[10])) 
        if atom[16] == '-1':
            frozen.append(i+1)
    inp.write("$end\n\n")

    inp.write("$rem\n")
    inp.write("jobtype %s\n"%jobtype[0])
    inp.write("method %s\n" % method)
    inp.write("basis %s\n"%basis)
    inp.write("ecp %s\n"%ecp)
    inp.write("solvent_method pcm\n")
    inp.write("symmetry false\n")
    inp.write("sym_ignore true\n")
    inp.write("no_reorient true\n")
    inp.write("geom_opt_max_cycles 1500\n")
    inp.write("mem_total 120000\n")
    inp.write("DFT_D D3_BJ\n")
    inp.write("$end\n\n")

    inp.write("$opt\nFIXED\n")
    for i in frozen:
        inp.write("%s XYZ\n"%i)
    inp.write("ENDFIXED\n")
    inp.write("$end\n\n")
    
    if len(def_basis) > 0:
        for l in def_basis:
            inp.write("%s"%l)
        inp.write('\n')

    if len(jobtype) > 1:
        inp.write("@@@@\n\n")
        for j in range(1,len(jobtype)):
            inp.write("$molecule\n")
            inp.write("read\n")
            inp.write("$end\n\n")
            inp.write("$rem\n")
            inp.write("jobtype %s\n"%jobtype[j])
            inp.write("method %s\n" % method)
            inp.write("basis %s\n"%basis)
            inp.write("ecp %s\n"%ecp)
            inp.write("solvent_method pcm\n")
            inp.write("symmetry false\n")
            inp.write("no_reorient true\n")
            inp.write("mem_total 120000\n")
            inp.write("DFT_D D3_BJ\n")
            inp.write("SCF_GUESS read\n")
            inp.write("$end\n\n")
        if len(def_basis) > 0:
            for l in def_basis:
                inp.write("%s"%l)
            inp.write('\n')


    inp.close()


"""
write gau-xtb input use gau-xtb template file
"""
def write_xtb_input(inp_name,inp_temp,charge,multiplicity,pic_atom,tot_charge,res_count):
    ### inp_name default can be 1.inp, but first is name.input such as 9.input
    ### inp_type list which includes [small/large, level_theory, basis, opt, freq,  
    ### input_template line0: size
    ### input_template line1: level
    ### input_template line2: opt + calcfc/readfc or opt + modred + info
    ### input_template line3: freq
    ### input_template line4: guess=read
    ### input_template line5: geom=checkpoint
    ### input_template line6: iop
    ### input_template line7: scf_info
    ### input_template line8: scrf_info

    print("charge= %d, tot_charge_without'-c'= %d"%(charge,tot_charge))
    with open(inp_temp) as f:
        lines = f.readlines()

    inp = open('%s'%inp_name,'w')
    for i in range(3):
        inp.write('%s'%lines[i])

    inp.write("#P ")      

    f_atom = []
    if lines[3][0] != '#' and 'modred' not in lines[3]:
        inp.write('%s '%lines[3][:-1])
    ### If modred ###
    elif lines[3][0] != '#' and 'modred' in lines[3]:
        optl = lines[3].split()
        inp.write("%s "%optl[0])
        ### get frozen atom info ###
        modred_info, modred_code = optl[1:3]
        pairs = modred_info.split(';')
        for lp in range(len(pairs)):
            lmod = []
            pair = pairs[lp]
            ar = pair.split(',')
            ### Atom name and Atom id in pair ###
            for i in range(0,len(ar),2):
                for atom in pic_atom:
                    if ar[i] == atom[2].strip() and int(ar[i+1]) == atom[6]:
                        lmod.append(pic_atom.index(atom)+1)
            f_atom.append(lmod)
    for l in range(4,8):
        if lines[l][0] != '#':
            inp.write("%s "%lines[l].strip())

    inp.write("\n\n")
    ### info_line ###
    inp.write("%s"%res_count)
    inp.write("\n\n")
    inp.write("%d %d\n"%(charge+tot_charge,multiplicity))

    for atom in pic_atom:
        inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),atom[16],atom[8],atom[9],atom[10])) 
    inp.write("\n")

    if len(f_atom) != 0:
        for l in range(len(f_atom)):
            for ids in f_atom[l]:
                inp.write("%s "%ids)
            inp.write("%s\n"%modred_code)
        inp.write("\n")

    inp.close()

    ### write 1.fix file ###
    f = open("1.fix",'w')
    atom_idx = []
    for i in range(len(pic_atom)):
        if pic_atom[i][-1] == '-1':
            atom_idx.append(i+1)
    
    f.write('$hess\n')    
    for ids in atom_idx:
        f.write('    scale mass: %d, 90000000\n'%ids)
    f.close()
    ### else can be constraint ###



"""
write orca input files using orca template file
"""
def write_orca_input(inp_name,inp_temp,charge,multiplicity,pic_atom,tot_charge,res_count,hopt):
    ### From input_template file ###
    ### nprocs
    ### maxcore
    ### jobtype opt/freq
    ### method
    ### disp
    ### basis
    ### auxbasis
    ### solvent_method

    with open(inp_temp) as f:
        lines = f.readlines()
    hline = ["!"]
    for line in lines:
        if line.startswith('#'): continue
        if line.startswith('jobtype'):
            v = line.split(':')
            if "/" in v[1]:
                jobtype = [x.strip() for x in v[1].split('/')]
            else:
                jobtype = v[1].strip()
            hline.append(jobtype)
        if line.startswith('nprocs'):
            nprocs = line.split(':')[1].strip()
        if line.startswith('maxcore'):
            maxcore = line.split(':')[1].strip()
        if line.startswith('method'):
            method = line.split(':')[1].strip()
            hline.append(method)
        if line.startswith('disp'):
            disp = line.split(':')[1].strip()
            hline.append(disp)
        if line.startswith('basis'):
            basis = line.split(':')[1].strip()
            hline.append(basis)
        if line.startswith('auxbasis'):
            auxbasis = line.split(':')[1].strip()
            hline.append(auxbasis)
        if line.startswith('ecp'):
            ecp = line.split(':')[1].strip()
            hline.append(ecp)
        if line.startswith('solvent'):
            solv = line.split(':')[1].strip()
            hline.append(solv)
        if line.startswith('cpcm_epsilon'):
            soleps = line.split(':')[1].strip()
        # DAW note: idea was that custom basis could be put at end of input template so read from keyword to end. 
        # but I haven't yet got around to adding code to format and add the custom basis into the input file...
        #if line.startswith('def_basis'):
        #    start = lines.index(line)
        #    def_basis = []
        #    for l in range(start+1,len(lines)):
        #        def_basis.append(lines[l])

    ## write input
    inp = open('%s'%inp_name,'w')
    #inp.write("! %s %s %s %s %s CPCM(water)\n"%(jobtype,method,disp,basis,auxbasis))
    inp.write(" ".join(hline)+"\n")
    inp.write("%%pal nprocs %s end\n"%nprocs)
    inp.write("%%maxcore %s000\n"%maxcore)
    if 'soleps' in locals() and soleps:
        inp.write("%%cpcm epsilon %s end\n"%soleps)
    inp.write("%geom\n  constraints\n")
    for i in range(len(pic_atom)):
        atom = pic_atom[i]
        if atom[16] == '-1':
            atn = int(atom[1]) - 1
            inp.write("  { C %d C }\n"%atn)
    inp.write("  end\n")
    if hopt == 1:
        inp.write("  optimizeHydrogens true\n")
    inp.write("end\n")
    inp.write("*xyz %d %d\n"%(charge+tot_charge,multiplicity))
    for i in range(len(pic_atom)):
        atom = pic_atom[i]
        if atom[16] == '-1':
            inp.write("%6s %8.3f %8.3f %8.3f  M = 999999999999.9\n"%(atom[14].strip(),atom[8],atom[9],atom[10]))
        else:
            inp.write("%6s %8.3f %8.3f %8.3f\n"%(atom[14].strip(),atom[8],atom[9],atom[10]))
    inp.write("*\n")

    inp.close()

"""
write psi4 fsapt input file
"""
def write_psi4_fsapt_input(inp_name,inp_temp,charge,multiplicity,pic_atom,tot_charge,res_count,seed):
    ### From input_template file ###
    ### memory
    ### set_num_threads
    ### molecule
    ### set_globals
    ### set_sapt
    ### energy

    with open(inp_temp) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('#'): continue
        if line.startswith('memory'):
            memory = line.split(':')[1].strip()
        if line.startswith('set_num_threads'):
            threads = line.split(':')[1].strip()
        if line.startswith('molecule'):
            mol = line.split(':')[1].split(',')
        if line.startswith('set_globals'):
            setglob = line.split(':')[1].split(',')
        if line.startswith('energy'):
            energy = line.split(':')[1].strip()

    # parse seed, separate fragments
    seedlist = []
    seed = seed.split(',')
    for i in seed:
        i = i.split(':')
        seedlist.append((i[0],int(i[1])))
    seedatoms = []
    enzatoms = []
    for i in range(len(pic_atom)):
        atom = pic_atom[i]
        if (atom[5],int(atom[6])) in seedlist:
            seedatoms.append(atom)
        else:
            enzatoms.append(atom)

    ## write input
    inp = open('%s'%inp_name,'w')
    inp.write("memory %s\n\n"%memory)
    inp.write("set_num_threads(%s)\n\n"%threads)
    inp.write("molecule frame {\n")
    inp.write("%s %s\n"%(charge,multiplicity))
    for i in range(len(seedatoms)):
        atom = seedatoms[i]
        inp.write("%1s %8.3f %8.3f %8.3f\n"%(atom[14].strip(),atom[8],atom[9],atom[10]))
    inp.write("--\n")
    inp.write("%s 1\n"%tot_charge)
    for i in range(len(enzatoms)):
        atom = enzatoms[i]
        inp.write("%1s %8.3f %8.3f %8.3f\n"%(atom[14].strip(),atom[8],atom[9],atom[10]))
    inp.write("\n")
    for m in mol:
        inp.write("\t%s\n"%m)
    inp.write("} \n \n")
    inp.write("set globals {\n")
    for g in setglob:
        inp.write("\t%s\n"%g)
    inp.write("} \n \n")
    inp.write("set sapt print 1\n\n")
    inp.write("energy('%s')\n"%energy)

    inp.close()

