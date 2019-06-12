import os, sys, re, filecmp
from numpy import *
import argparse
from read_write_pdb import *
from glob import glob


def system_run(cmd):
    print(cmd)
    exit = os.system(cmd)
    if ( exit != 0 ):
        print('failed to run:')
        print(cmd)
        sys.exit()
 
### copy h-added pdb xyz and other information into tmppdb ###
def pdb_after_addh(tmppdb,newpdb):
    tmp_pdb, res_info, tot_charge_t = read_pdb(tmppdb)
    tmp_xyz = []
    for i in tmp_pdb:
        tmp_xyz.append([i[8],i[9],i[10]])
    new_pdb, binfo, tot_charge = read_pdb(newpdb)     #can be just xyz files from cerius or pymol
    pic_atom = []
    for line in new_pdb:
        if [line[8],line[9],line[10]] not in tmp_xyz:
            line[15] = '0 '
            pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],' H',line[15],line[16]])
        else:
            if '+' in line[15] or '-' in line[15]:
                charge = line[15]
            else:
                charge = '0 '
            idx = tmp_xyz.index([line[8],line[9],line[10]])
            line = tmp_pdb[idx]
            line[15] = charge
            if [line[2],line[5],line[6]] in res_info:
                pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],'-1'])
            else:
                pic_atom.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16]])
    return pic_atom, tot_charge #, xyz, atom, hold

########################################################################
### replace certain part of the pdb with provided pdb/xyz ##############
### tmppdb is the minimized structure, most of these xyz will be kept
### newpdb is the one contain transition structure/fragment
### parts is the residue name and atom name for the transition part
########################################################################
def pdb_replace(tmppdb,newpdb,parts):
    tmp_pdb, res_info, tot_charge_t = read_pdb(tmppdb)
    tmp_xyz = []
    new_xyz = []
    for i in tmp_pdb:
        tmp_xyz.append([i[4].strip(),i[2].strip()])
    new_pdb, binfo, tot_charge = read_pdb(newpdb)     #can be just xyz files from cerius or pymol
    for i in new_pdb:
        new_xyz.append([i[4].strip(),i[2].strip()])
    print(new_xyz[0])

    if parts == None:   #newpdb has the entire thing to replace the tmppdb
        for line in new_pdb:
            resatom = [line[4].strip(),line[2].strip()]
            idx = tmp_xyz.index(resatom)
            tmp_pdb[idx] = line
    else:
        for i in range(len(parts)):
            idx1 = new_xyz.index([parts[i][0],parts[i][1]])
            idx2 = tmp_xyz.index([parts[i][0],parts[i][1]])
            tmp_pdb[idx2] = new_pdb[idx1]

    return tmp_pdb, tot_charge_t #, xyz, atom, hold

def write_input(inp_name,inp_temp,charge,multiplicity,pic_atom,tot_charge):
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

    print("charge= %d, tot_charge= %d"%(charge,tot_charge))
    with open(inp_temp) as f:
        lines = f.readlines()

    inp = open('%s'%inp_name,'w')
    inp.write("%chk=1.chk\n")           #write check file into 1.chk
    inp.write("%nprocshared=10\n")      #use 10 processors

    if 'small' in lines[0]:
        inp.write("%mem=20GB\n")        #use memory if it is large job use 80GB
    else:
        inp.write("%mem=80GB\n")
    inp.write("#P ")      

    if lines[1][0] != '#':
        inp.write("%s "%lines[1].strip())      

#    if bool(lines[2].strip()): line is empty
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
                lmod.append(len(ar)/2)
                for i in range(0,len(ar),2):
                    for atom in pic_atom:
                        if ar[i] in atom[2] and ar[i+1] in atom[4]:
                            f_atom.append(pic_atom.index(atom)+1)
    for l in range(3,9):
#        if bool(lines[l].strip()):
        if lines[l][0] != '#':
            inp.write("%s "%lines[l].strip())

    inp.write("\n\n")
    inp.write("info_line\n")
    inp.write("\n")
    inp.write("%d %d\n"%(charge+tot_charge,multiplicity))

    if lines[5][0] == '#':
        ### pm7 with opt only will be relax h step/ pm7 with opt(modred) otherwise
        if lines[1][0] != '#' and 'pm' in lines[1] and lines[2].strip() == 'opt':
            for atom in pic_atom:
                if atom[14].strip() == 'H':
                    inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),'0',atom[8],atom[9],atom[10])) 
                else:
                    inp.write("%4s %6s         %8.3f %8.3f %8.3f\n"%(atom[14].strip(),'-1',atom[8],atom[9],atom[10])) 
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

    if lines[10][0] != '#' and 'basis' in lines[10]:
        for l in lines[11:]:
            inp.write("%s"%l)
        inp.write('\n')

    if lines[8][0] != '#' and 'scrf' in lines[8]:
        inp.write('radii=uff\nalpha=1.2\neps=4.0\n\n')
    
    inp.close()


def gen_pdbfiles(wdir,step,tmppdb):
    new_dir = '%s/step%spdbs'%(wdir,step)
    if os.path.isdir(new_dir):
        system_run('rm -r %s'%new_dir)
        system_run('mkdir %s'%new_dir)
    else:
        system_run('mkdir %s'%new_dir)
    os.chdir(new_dir)
#    system_run('gopt_to_pdb.py %s %s/step-%s-out 0'%(tmppdb,wdir,step))
    system_run('gopt_to_pdb.py -p %s -o %s/step-%s-out -f -1'%(tmppdb,wdir,step))
    os.chdir('%s'%wdir)
    pdb_name = []
    for pdbf in glob('%s/*.pdb'%new_dir):
        m = re.search('(\d+).pdb',pdbf)
        pdb_name.append( int(m.group(1)) )
    return max(pdb_name), new_dir


if __name__ == '__main__':

    #########################################################################################################################################################
    ### In working directory such as model9-ts-001 
    ### Run: write_input.py -noh nohpdb -adh haddpdb -intmp relaxh_temp -c 2 (final_pdb is saved as "template.pdb", input is saved as "1.inp", default -m 1) 
    ### Run: write_input.py -step 1 -intmp modred_temp -m 1 -c 2 (will take the last pdb from 1.out and write input as "1.inp")                      
    ### Run: write_input.py -step 2 -intmp modred_temp -m 1 -c 2 -new step1pdbs/33.pdb (will take the selected pdb and write "1.inp")
    #########################################################################################################################################################

    parser = argparse.ArgumentParser(description='Prepare template PDB files, write input files, save output PDB files in working directory')
    parser.add_argument('-step', dest='step', default=0, type=int, 
            help='step0: read noh, addh pdbs, write_final_pdb and read input_template write_first_inp,\
          step1: read outputwrite_modred_inp, input_template, write_second_inp,\
          step2: read outputwrite_modred_inp, input_template, write_new_inp')
    parser.add_argument('-wdir', dest='output_dir', default=os.path.abspath('./'), help='working dir')
    parser.add_argument('-tmp', dest='tmp_pdb', default=None, help='template_pdb_file')
    parser.add_argument('-noh', dest='no_h_pdb', default=None, help='trimmed_pdb_file')
    parser.add_argument('-adh', dest='h_add_pdb', default=None, help='hadded_pdb_file')
    parser.add_argument('-new', dest='new_pdb', default=None, help='new_pdb_file')
    parser.add_argument('-intmp', dest='input_tmp', default=None, help='template_for_write_input')
    parser.add_argument('-outf', dest='gau_out', default='1.out', help='output_name')
    parser.add_argument('-inpn', dest='inp_name', default='1.inp', help='input_name')
    parser.add_argument('-m', dest='multiplicity', default=1, type=int, help='multiplicity')
    parser.add_argument('-c', dest='ligand_charge', default=0, type=int, help='charge_of_ligand')
    parser.add_argument('-pdb1', dest='pdb1', default=None, help='minima_pdb_file')
    parser.add_argument('-pdb2', dest='pdb2', default=None, help='ts_pdb_file')
    parser.add_argument('-parts', dest='parts', default=None, help='ts_frag_indo')

    args = parser.parse_args()

    step = args.step
    wdir = args.output_dir
    if args.tmp_pdb is None:
        tmp_pdb = '%s/template.pdb'%wdir
    else:
        tmp_pdb   = args.tmp_pdb

    nohpdb   = args.no_h_pdb
    adhpdb   = args.h_add_pdb
#    newpdb   = args.new_pdb
    int_tmp  = args.input_tmp
    gauout   = args.gau_out
    inp_name = args.inp_name
    multi    = args.multiplicity
    charge   = args.ligand_charge
    wdir     = args.output_dir

    if step == 0:
        pic_atom, tot_charge = pdb_after_addh(nohpdb,adhpdb)
        write_pdb('%s'%tmp_pdb,pic_atom)
    elif step == 1:
        i_name = []
        for gau_input in glob('%s/*inp'%wdir):
            m = re.search(r'-(\d+)-inp', gau_input)
            if m:
                i_name.append( int(m.group(1)) )
        if len(i_name) == 0:
            i_step = 1
            if os.path.isfile('%s/1.out'%wdir):
                system_run( 'cp 1.inp step-%s-inp'%(i_step) )
                system_run( 'cp 1.out step-%s-out'%(i_step) )
                system_run( 'cp 1.chk step-%s-chk'%(i_step) )
        else:
            i_step = max(i_name)
            if filecmp.cmp('1.inp','%s/step-%s-inp'%(wdir,i_step)) is False and filecmp.cmp('1.out','%s/step-%s-out'%(wdir,i_step)) is False:
                i_step += 1
                system_run( 'cp 1.inp step-%s-inp'%(i_step) )
                system_run( 'cp 1.out step-%s-out'%(i_step) )
                system_run( 'cp 1.chk step-%s-chk'%(i_step) )
            elif filecmp.cmp('1.inp','%s/step-%s-inp'%(wdir,i_step)) is True and filecmp.cmp('1.out','%s/step-%s-out'%(wdir,i_step)) is True:
                i_step = i_step
                print("1.inp is the same as step-%s-inp, and 1.out is the same as step-%s-out, please check!")
#                sys.exit()
            else:
                print("check if the files are propagated correctly!")
                sys.exit()
        if step == 1:
            pdb_file, new_dir = gen_pdbfiles(wdir,i_step,tmp_pdb)
            pic_atom, res_info, tot_charge = read_pdb('%s/%s.pdb'%(new_dir,pdb_file))
    elif step == 2:
        newpdb = args.new_pdb
        pic_atom, res_info, tot_charge = read_pdb(newpdb)
    elif step == 3:
        pdb1 = args.pdb1
        pdb2 = args.pdb2
        parts = args.parts
        if parts != None:
           with open(parts) as f:
               plines = f.readlines()
           parts = []
           for line in plines:
               parts.append(line.split())
        pic_atom, tot_charge = pdb_replace(pdb1,pdb2,parts)    
    write_input('%s/%s'%(wdir,inp_name),int_tmp,charge,multi,pic_atom,tot_charge)
        
