"""
This is a program written by Qianyi Cheng
at University of Memphis.
"""


from read_write_pdb import *
from rms import *
from numpy import *
import sys, re, os
import argparse
from read_gout_xtb import *

def form_irc_xyz(opt,atom_idx,xyz,scale):
    irc1 = []
    irc2 = []
    for i in range(opt.shape[0]):
        if i not in atom_idx:
            irc1.append(opt[i,:])
            irc2.append(opt[i,:])
        else:
            irc1.append(opt[i,:]+scale*array(xyz[atom_idx.index(i)]))
            irc2.append(opt[i,:]-scale*array(xyz[atom_idx.index(i)]))
    return array(irc1), array(irc2)

def write_input(inp_f,dir1,dir2,irc1,irc2,atom_name,charge,multip,scale):
    f1_list = []
    f2_list = []
    with open(inp_f) as f:
        lines = f.readlines()
    lines[3] = re.sub(r' iop\S+','',lines[3])
    #lines[3] = re.sub(r' \S+checkpoint\S+','',lines[3])
    lines[3] = re.sub(r' geom\S+','',lines[3])
    lines[3] = re.sub(r' guess\S+','',lines[3])
    lines[3] = re.sub(r'\(ts\S+','',lines[3])
    if 'freq' not in lines[3]:
        index = lines[3].find('scf=')
        lines[3] = lines[3][:index]+'freq '+lines[3][index:]
    if 'opt' not in lines[3]:
        index = lines[3].find('freq')
        lines[3] = lines[3][:index]+'opt(nomicro) '+lines[3][index:]
    elif 'opt' in lines[3]:
        index1 = lines[3].find('opt')
        index2 = lines[3].find('freq')
        lines[3] = lines[3][:index1]+'opt(nomicro) '+lines[3][index2:]
    for i in range(5):
        f1_list.append(lines[i])
        f2_list.append(lines[i])
    f1_list.append("the positive pertubation structure with scale +%.2f\n"%scale)
    f2_list.append("the negative pertubation structure with scale -%.2f\n"%scale)
    for i in range(6,8):
        f1_list.append(lines[i])
        f2_list.append(lines[i])
    for atom in range(len(irc1)):
        f1_list.append('%6s%6s%12.6f%12.6f%12.6f\n'%(atom_name[atom][0],atom_name[atom][1],irc1[atom,0],irc1[atom,1],irc1[atom,2]))
        f2_list.append('%6s%6s%12.6f%12.6f%12.6f\n'%(atom_name[atom][0],atom_name[atom][1],irc2[atom,0],irc2[atom,1],irc2[atom,2]))
    f1_list.append("\n")
    f2_list.append("\n")
    
    for i in range(9,len(lines)):
        v =lines[i].split()
        if len(v) == 5 and v[1] in ['0','-1']:
            continue
        elif len(v) > 0 and v[-1] in ['A','F']:
            continue
        else:
            if lines[i].strip() == '' and f1_list[-1] == '\n':
                continue
            else:
                f1_list.append(lines[i])
                f2_list.append(lines[i])

    f1 = open('%s/1.inp'%dir1,'w')
    f2 = open('%s/1.inp'%dir2,'w')
    for l in range(len(f1_list)):
        f1.write("%s"%f1_list[l])
        f2.write("%s"%f2_list[l])
    
#            if len(v) == 2 and v[0] == str(charge) and v[1] == str(multip):
    f1.close()
    f2.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Write irc inputs')
    parser.add_argument('-wdir', dest='output_dir', default=os.path.abspath('./'), help='working dir')
    parser.add_argument('-dir1', dest='dir1', default=os.path.abspath('./irc1'), help='irc1 dir')
    parser.add_argument('-dir2', dest='dir2', default=os.path.abspath('./irc2'), help='irc2 dir')
    parser.add_argument('-o', dest='gau_out', default=None, help='output_name')
    parser.add_argument('-i', dest='gau_inp', default=None, help='input_name')
    parser.add_argument('-s', dest='scale', type=float,default=0.1, help='scale_factor')
    parser.add_argument('-n', dest='num_freq', type=int,default=1, help='scale_factor')


    args = parser.parse_args()
    wdir = args.output_dir
    scale = args.scale
    num_freq = args.num_freq
    dir1 = args.dir1
    dir2 = args.dir2
    os.mkdir(dir1)
    os.mkdir(dir2)
    if args.gau_out is None:
        out_f = '%s/1.out'%wdir
    else:
        out_f = args.gau_out

    if args.gau_inp is None:
        inp_f = '%s/1.inp'%wdir
    else:
        inp_f = args.gau_inp

    nimag, natoms, charge, multip = gaussian_num(out_f)
    with open(out_f) as f:
        lines = f.readlines()

    opt = array(gaussian_opt_xyz(lines,natoms)[-1])
    atom_name = gaussian_atom_names(lines,natoms)
    atom_idx, freq_xyz, freq_info = gaussian_freq(lines)


    xyz = freq_xyz[num_freq]
    irc1, irc2 = form_irc_xyz(opt,atom_idx,xyz,scale)
    write_input(inp_f,dir1,dir2,irc1,irc2,atom_name,charge,multip,scale)
