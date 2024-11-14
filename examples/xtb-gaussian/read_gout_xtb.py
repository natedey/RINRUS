"""
This is a program written by Qianyi Cheng
at University of Memphis.
"""

from rms import *
from numpy import *
import sys, re, os
from read_write_pdb import *
from sys_tools import *
import argparse

def get_scf(lines):
    scfe = []
    for l in lines:
        if "SCF Done" in l:
            scfe.append(float(l.split()[4]))
#            scfe.append(lines.index(l))
    return scfe


def gaussian_num(gfile):
    with open(gfile) as f:
        fline = f.read().replace('\n ','')
    m = re.search('NImag=(\w+)',fline)
    try:
        nimag = int(m.group(0).split('=')[-1])
    except:
        nimag = None
    m = re.search('NAtoms=(\s+)(\d+)',fline)
    natoms = int(m.group(0).split('=')[-1])
#    m = re.search('NBasis=(\s+)(\d+)',fline)
#    nbasis = int(m.group(0).split('=')[-1])
    m = re.search('Charge =(\s+)[+-]?(\d+) Multiplicity =(\s+)(\d+)',fline)
    v = m.group(0).split()
    charge = int(v[2])
    multip = int(v[5])

    return nimag, natoms, charge, multip


def get_optl(lines):
    p_start = []
    for i in range(len(lines)):
        if 'Standard orientation' in lines[i]:
            p_start.append(i+5)
    return p_start


def gaussian_opt_xyz(lines,natoms):
    p_start = get_optl(lines)

    ### Read in stardard orientation geom
    opt = []
    for j in range(len(p_start)):
        xyz = []
        for i in range(p_start[j],p_start[j]+natoms):
            v = lines[i].split()
            xyz.append([float(v[3]),float(v[4]),float(v[5])])
        opt.append(xyz)
    return opt


def gaussian_energy(lines):
    for i in range(len(lines)):
        l = lines[i]
        #Zero-point correction=                           2.285998 (Hartree/Particle)
        #Thermal correction to Energy=                    2.407590
        #Thermal correction to Enthalpy=                  2.408534
        #Thermal correction to Gibbs Free Energy=         2.131712
        #Sum of electronic and zero-point Energies=          -6682.107609
        #Sum of electronic and thermal Energies=             -6681.986017
        #Sum of electronic and thermal Enthalpies=           -6681.985073
        #Sum of electronic and thermal Free Energies=        -6682.261895
        #
        #                    E (Thermal)             CV                S
        #                     KCal/Mol        Cal/Mol-Kelvin    Cal/Mol-Kelvin
        #Total                 1510.786            480.805            582.621
        if "Zero-point correction="  in l:
            zero = float(l.split()[-2])
        if "Thermal correction to Energy=" in l:
            ther = float(l.split()[-1])
        if "Thermal correction to Enthalpy=" in l:
            he = float(l.split()[-1])
        if "Thermal correction to Gibbs Free Energy=" in l:
            ge = float(l.split()[-1])
        if "E (Thermal)             CV                S" in l:
            se = float(lines[i+2].split()[-1])
        if "Sum of electronic and zero-point Energies" in l:
            szero = float(l.split()[-1])
        if "Sum of electronic and thermal Energies" in l:
            sele = float(l.split()[-1])
        if "Sum of electronic and thermal Enthalpies" in l:
            enthalpy = float(l.split()[-1])
        if "Sum of electronic and thermal Free Energies" in l:
            freeg = float(l.split()[-1])
    #return get_scf(lines)[-1], zero, ther, he, ge, se, szero, sele, enthalpy, freeg   
    return szero-zero, zero, ther, he, ge, se, szero, sele, enthalpy, freeg   


def gaussian_freq(lines):
    f_lines = []
    for i in range(len(lines)):
        if 'Frequencies' in lines[i][1:12]:
            f_lines.append(i)
    ### Read in imaginary frequency
    freq_xyz  = {}
    freq_info = {}
    keys = []
    num_atoms = f_lines[1]-f_lines[0]-7
    atom_idx = {}
    for i in range(len(f_lines)):
        start = f_lines[i] - 2
        end = f_lines[i]+4+num_atoms
        modes = lines[start].split()
        sym   = lines[start+1].split()
        mag   = lines[start+2].split()[2:]
        red   = lines[start+3].split()[3:]
        force = lines[start+4].split()[3:]
        ir_it = lines[start+5].split()[3:]
        for j in range(3):
            key = int(modes[j])
            keys.append(key)
            freq_info[key] = (sym[j],float(mag[j]),float(red[j]),float(force[j]),float(ir_it[j]))
            freq_xyz[key]  = []
        atom_idx[keys[-1]] = []
        for j in range(start+7,end+1):
            v = lines[j].split()
            atom_idx[keys[-1]].append(int(v[0])-1)
            for k in range(3):
                freq_xyz[keys[i*3+k]].append([float(v[3*k+2]),float(v[3*k+3]),float(v[3*k+4])])
                
    return atom_idx[3], freq_xyz, freq_info


def gaussian_atom_names(lines,natoms):
    atom_name = []
    for i in range(len(lines)):
        if "Redundant internal coordinates found in file." in lines[i]:
            if "Recover coordinates from OT file" in lines[i+1]:
                i += 1
            for idx in range(i+1,i+1+natoms):
                if ',' in lines[idx]:
                    atom_name.append(lines[idx][:-1].split(','))
                else:
                    va = lines[idx][:-1].split()
                    if len(va) == 4:
                        atom_name.append([va[0],'0'])
                    else:
                        atom_name.append(lines[idx][:-1].split()[:2])
    return atom_name                        

def gaussian_eigen(lines):
    f = open('etrack.dat','w')
    for l in lines:
       m = re.search('Done|YES|NO |Step number|     Eigenvalues ---   -',l)
       if m:
           print(l[:-1])
           f.write(l)
#       print( m.span())
    f.close()

#if __name__ == '__main__':
#    """ Usage: gopt_to_pdb.py -o ../1.out -s -1 -p ../template.pdb """
#    parser = argparse.ArgumentParser(description='generate pdbfiles from 1.out')
#    parser.add_argument('-f', dest='frame',default=None,help='select frame/range')
#    parser.add_argument('-o', dest='outputs',default='1.out',help='output files')
#    parser.add_argument('-s', dest='start',default=0,type=int,help='output steps')
#    parser.add_argument('-p', dest='pdbf',default='template.pdb',help='template pdb file')
#
#    args = parser.parse_args()
#    steps = []
#    files = []
#    if ',' in args.outputs:
#        files = args.outputs.split(',')
#        steps = args.start.split(',')
#    else:
#        files = [args.outputs]
#        steps = [args.start]
#    pdbf = args.pdbf
#    for gfile in files:
#        with open(gfile) as f:
#            lines = f.readlines()
#        gaussian_eigen(lines)
#        scfe, zero,ther, he, ge, se, szero, sele, enthalpy, freeg = gaussian_energy(lines)
#        print(scfe)
#        stand_opt = gaussian_opt_xyz(lines,natoms)
#        atom_idx, freq, info = gaussian_freq(lines)
