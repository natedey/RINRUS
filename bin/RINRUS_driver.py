#!/usr/bin/env python3

import os
import argparse
import shutil
import shlex,subprocess
from pathlib import Path
import sys
from subprocess import Popen, PIPE,STDOUT
import logging
import io
from datetime import datetime

#log header info, takes details of commit
def log_header():
    gitpath = str(Path(__file__).resolve().parents[1])
    pwd = os.getcwd()
    #gitver = subprocess.run(f"cd {gitpath}; git show -s --pretty='format:%h'; cd {pwd}",shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    gitver = subprocess.run(f"cd {gitpath}; git show -s --pretty='format:%h %cd' --date=format-local:'%Y-%m-%d %H:%M'; cd {pwd}",shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    gitver = gitver.stdout.split()
    #fetch = subprocess.run(f'stat -c %y {gitpath}/.git/FETCH_HEAD',shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    #fetchdate = fetch.stdout[0:10]
    #fetchtime = fetch.stdout[11:16]
    headtxt = ('--------------------------------------------------------------------------------------\n'
    '              RINRUS: The Residue Interaction Network ResidUe Selector                \n'
    '--------------------------------------------------------------------------------------\n'
    f'(C) 2018-2024. Using version {gitver[0]}, published on github {gitver[1]} at {gitver[2]}.        \n'
    'Developed in the group of Prof. Nathan DeYonker at the University of Memphis, TN USA. \n'
    'Contributors: Q. Cheng, N. DeYonker, D. Wappett, T. Summers, D. Agbaglo, T. Suhagia,  \n'
    '    T. Santaloci, J. Bachega.                                                         \n'
    'Acknowledge RINRUS by citing: github.com/natedey/RINRUS, DOI:10.1016/j.bpj.2021.07.029\n'
    '    and DOI:10.1039/D3CP06100K                                                        \n'
    '--------------------------------------------------------------------------------------' )

    clbanner = ('--------------------------------------------------------------------------------------\n'
    '           Running RINRUS: The Residue Interaction Network ResidUe Selector           \n'
    'Developed in the group of Prof. Nathan DeYonker at the University of Memphis, TN USA. \n'
    f'(C) 2018-2024. Using version {gitver[0]}, published on github {gitver[1]} at {gitver[2]}.        \n'
    '--------------------------------------------------------------------------------------\n')

    return headtxt,clbanner

def driver_file_reader(file,logger):
    path_to_RIN = ''
    red = 'True'
    pdb = ''
    Seedlist = []
    seed = ''
    must_include = ''
    RIN_program = ''
    selatom_file = ''
  #  Histidine = ''
    model_num = ''
    charge = ''
    multi = ''
  #  residues_not_protonated = ''
    Computational_program = ''
    template_path = ''
    basis_set_library = ''
    with open(file,'r') as fp:
        data = fp.readlines()

    # print input verbatim
    rawinp = ''
    for line in data:
        rawinp += '> ' + line
    logger.info('Inputs being read from: '+file+'\n'+f'----------raw contents----------\n'+rawinp+'--------------------------------')

    # process inputs
    for line in data:
        if '#' in line:
            line
        else:
            if 'Path_to_scripts:' in line:
                path_to_RIN = line.replace('Path_to_scripts:', '').replace('\n','').replace(' ','')
            if 'PDB' in line:
                pdb = line.replace('PDB:','').replace('\n','').replace(' ','')
            if 'Protonate_initial:' in line:
                red = line.replace('Protonate_initial:','').replace('\n','').replace(' ','')
            if 'Seed:' in line:
                seed = line.replace('Seed:','').replace('\n','').replace(' ','')
                if ',' in seed:
                    line = seed.split(',')
                    for i in line:
                        if i != '':
                            Seedlist.append(i)
                else:
                    Seedlist.append(seed)
            if 'Must_include' in line:
                must_include = line.replace('Must_include:','').replace('\n','').replace(' ','')
            if 'RIN_program:' in line:
                RIN_program = line.replace('RIN_program:','').replace('\n','').replace(' ','')
            if 'RIN_info_file:' in line:
                selatom_file = line.replace('RIN_info_file:','').replace('\n','').replace(' ','')
            if 'Seed_charge:' in line:
                charge = line.replace('Seed_charge:','').replace('\n','').replace(' ','')
            if 'Multiplicity:' in line:
                multi = line.replace('Multiplicity:','').replace('\n','').replace(' ','')
            if 'Computational_program:' in line:
                Computational_program = line.replace('Computational_program:','').replace('\n','').replace(' ','')
            if 'Input_template_path:' in line:
                template_path = line.replace('Input_template_path:','').replace('\n','').replace(' ','')
            if 'Gaussian_basis_intmp:' in line:
                gbtmp = line.replace('Gaussian_basis_intmp:','').replace('\n','').replace(' ','')
                if gbtmp.lower() in ['true', 't', '1', 'y']:
                    basis_set_library = 'intmp'
                else:
                    basis_set_library = 'default_dict'
            #if 'Basisset_library:' in line:
            #    basis_set_library = line.replace('basisset_library:','').replace('\n','').replace(' ','')
            if 'Model(s):' in line:
                model_num = line.replace('Model(s):','').replace('\n','').replace(' ','')
    
    
    if red.lower() in ['false', 'f', '0', 'n']:
        red = 'False'
    elif red.lower() in ['true', 't', '1', 'y']:
        red = 'True'

    logger.info('Path to the RINRUS scripts bin directory: ' + str(path_to_RIN))
    logger.info('PDB Name: ' + pdb)
    logger.info('Protonate initial PDB: ' + red)
    logger.info('RIN program: ' + RIN_program)
    if RIN_program.lower() == 'manual':
        logger.info('Manual RIN info file: ' + selatom_file)
    logger.info('Seed: ' + str(seed))
    logger.info('Seed Charge: ' + str(charge))
    if must_include == '':
        logger.info('Fragments that must be included: none')
    else:
        logger.info('Fragments that must be included: ' + str(must_include))
    logger.info('Multiplicity: ' + str(multi))
    logger.info('Model number selection: ' + str(model_num))
    logger.info('Computational Program: ' + str(Computational_program))
    if Computational_program.lower() != 'none':
        if template_path == '':
            logger.info('Path to the input template file (default): ~/git/RINRUS/template_files/'+str(Computational_program)+'_input_template.txt')
        else:
            logger.info('Path to the input template file (specified): ' + str(template_path))
        if Computational_program.lower() == 'gaussian':
            logger.info('Source of basis sets for Gaussian: ' + str(basis_set_library)) 
        #logger.info('Path to the basis set library: ' + str(basis_set_library))
    
    return pdb,red,Seedlist,must_include,RIN_program,selatom_file,charge,multi,Computational_program,template_path,basis_set_library,seed,model_num,path_to_RIN


def res_atom_count(filename,must_include,Seedlist):
    seednum = len(Seedlist)
    # residues in res_atoms
    resnum = 0
    with open(filename,'r') as fp:
        data = fp.readlines()
        for i in data:
            key = i.split()[0] + ':' + i.split()[1]
            if i != '' and not i.startswith('#'):
                resnum += 1
    # must_include fragments
    if must_include != '':
        addnum = len(must_include.split(','))
    else:
        addnum = 0
    totnum = resnum + addnum
    return seednum,resnum,addnum,totnum


def run_reduce(pdb,logger,path_to_RIN):
    """
    runs reduce with -NOFLIP flag and uses logger to log output
    """
    print('Protonating input PDB file')
    path = os.path.expanduser(path_to_RIN+'/reduce')
    pdb_2 = pdb.replace('.pdb','')
    args = path_to_RIN + '/reduce -NOFLIP -Quiet  '+ str(pdb)+ ' > '+ str(pdb_2)+'_h.pdb'
    io.StringIO(initial_value='', newline='\r')
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    #print(out.stdout)
    logger.info('Reduce run as: '+ str(out.args))
    logger.info('Return code: '+ str(out.returncode))
    logger.info('Output from Reduce: '+ out.stdout)
    #shutil.copy(str(pdb_2)+'_h.pdb',str(pdb_2)+'_h_modify.pdb')
    #mod_pdb = str(pdb_2)+'_h_modify.pdb'
    mod_pdb = str(pdb_2)+'_h.pdb'
    return mod_pdb


def select_by_probe(pdb,seed,logger,path_to_RIN):
    """_summary_
    Hardcoded Flags:
        -unformated
        -MC
        -self "all"
        -Quiet
    Args:
        pdb (_type_): _h modified pdb
        seed (_type_): seed fragments
        logger (_type_): shows the probe input, return code and output from the probe program
        path_to_RIN (_type_): path to where file is

    Returns:
        _type_: .probe file
    """
    
    print('Generating probe RIN')
    probe = pdb.replace('.pdb','')
    args = [path_to_RIN+'/probe -unformated -MC -self "all" -Quiet '+ pdb +' > '+ probe + '.probe']
    probe = probe + '.probe'
    logger.info('Probe run as: '+ str(' '.join(args)))
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Return code: '+ str(out.returncode))
    logger.info('Output from Probe: '+ out.stdout)

    print('Analyzing probe RIN')
    path = os.path.expanduser(path_to_RIN+'/probe2rins.py')
    args =  sys.executable + ' ~/git/RINRUS/bin/probe2rins.py -f '+ str(probe)+ ' -s ' + seed
    logger.info('RIN analysis run as: '+ str(args))
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Return code= '+ str(out.returncode))
    #logger.info('Output: \n'+ out.stdout)
    return probe
    

def select_by_distance(calc_type,hydro,pdb,seed,cut,logger,path_to_RIN):
    path = os.path.expanduser(path_to_RIN+'/dist_rank.py')
    arg = [sys.executable, path ,'-pdb',str(pdb),'-s',str(seed),'-max',cut,'-type',calc_type]
    if hydro.lower() == "no":
        arg.append('-noH')
    logger.info('Distance selection run as: '+ str(' '.join(arg)))
    result = subprocess.run(arg)
    logger.info('Return code: '+ str(result.returncode))
    logger.info('Output:\n'+ str(result.stdout))
    return


def select_by_arpeggio(pdb,seed,path_to_RIN,logger):
    print('Generating arpeggio RIN')
    path = os.path.expanduser(path_to_RIN+'/arpeggio/arpeggio.py')
    arg = [sys.executable,path,str(pdb)]
    logger.info('Arpeggio run as: ' + str(' '.join(arg)))
    result = subprocess.run(arg)
    #logger.info('arpeggio.py output' + result.args)
    print('Analyzing arpeggio RIN')
    path = os.path.expanduser(path_to_RIN+'/arpeggio2rins.py')
    arg = [sys.executable,path,'-f',str(pdb).replace('pdb','contacts'),'-s',seed]
    result = subprocess.run(arg)
    logger.info('RIN analysis run as: '+ str(' '.join(arg)))
    #logger.info('arpeggio2rins.py output' + result.args)
    return
    

def trim_model(seed,must_include,pdb,selfile,model_num,path_to_RIN,RIN_program,logger):
    path = os.path.expanduser(path_to_RIN+'/rinrus_trim2_pdb.py')
    args = [sys.executable,path, '-s',str(seed).replace('\n',''), '-pdb',str(pdb),'-c',str(selfile), '-model', str(model_num), '-mustadd', str(must_include)]
    result = subprocess.run(args)
    logger.info('Model trimming run as: ' + str(' '.join(args)))
    return


def protonate_model(freeze,model_num,path_to_RIN,logger):
    path = os.path.expanduser(path_to_RIN+'/pymol_protonate.py')
    name = 'res_' + str(model_num)+'.pdb'
    arg= [sys.executable,path,'-pdb', name]
    if freeze[0] != '':
        arg.append('-ignore_ids')
        arg.append(str(freeze[0]))
    if freeze[1] != '':
        arg.append('-ignore_atoms')
        arg.append(str(freeze[1]))
    if freeze[2] != '':
        arg.append('-ignore_atnames')
        arg.append(str(freeze[2]))
    out = subprocess.run(arg)
    logger.info('Model protonation run as: '+ str(' '.join(arg)))
    return

def make_temp_pdb(model_num,path_to_RIN,logger):
    path = os.path.expanduser(path_to_RIN+'/make_template_pdb.py')
    name = 'res_'+str(model_num)
    arg = [sys.executable,path,'-name',name]
    out = subprocess.run(arg)
    logger.info('Creation of template pdb run as: '+ str(' '.join(arg)))
    return

def create_input_file(template,format,basisinfo,charge,model_num,path_to_RIN,logger):
    path = os.path.expanduser(path_to_RIN+'/write_input.py')
    modpdb = f'model_{model_num}_template.pdb'
    inpn = f'{model_num}.inp'
    arg = [sys.executable,path,'-format',str(format),'-c',str(charge),'-pdb',modpdb,'-inpn',inpn]
    if template != '' and template != None:
        arg.append('-intmp')
        arg.append(str(template))
    if format.lower() == 'gaussian' and basisinfo == 'intmp':
        arg.append('-basisinfo')
        arg.append(basisinfo)
    result = subprocess.run(arg)
    logger.info('Write_input run as: '+ str(' '.join(arg)))
    return



def run_rinrus_driver(file):
    ### SET UP LOGGING
    dt = datetime.now().strftime("%Y-%m-%d")
    lf = f"rinrus_log_{dt}.out"
    logging.basicConfig(level=logging.DEBUG,
                    filename=lf,
                    format='%(message)s',
                    filemode='w')
    logger = logging.getLogger()

    #write header 
    header,clbanner = log_header()
    logger.info(header+'\n\n')
    print(clbanner)

    #change to date + message format for main part of log
    logger.handlers[0].setFormatter(logging.Formatter(fmt='%(asctime)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S'))
    logger.info('RUNNING RINRUS DRIVER:')
    logger.info('Python executable being used: ' + sys.executable + '\n\n')

    
    ### PARSE INPUTS
    logger.info('READING DRIVER INPUT FILE:')
    driver_input_file = file
    #logger.info('Inputs from: ' + file)
    pdb,red,Seedlist,must_include,RIN_program,selatom_file,charge,multi,Computational_program,template_path,basis_set_library,seed,model_num,path_to_RIN = driver_file_reader(file,logger)
    RIN_program = RIN_program.lower()
    model_num = model_num.strip()
    logger.info('section done\n\n')

    
    ### INITIAL PROTONATION
    if red == 'True':
        logger.info('PROTONATION OF INITIAL PDB:')
        mod_pdb = run_reduce(pdb,logger,path_to_RIN)
        logger.info('section done\n\n')
    elif red == 'False':
        mod_pdb = pdb


    ### SELECTING ACTIVE SITE
    logger.info('SELECTING ACTIVE SITE:')
    if RIN_program.lower() == 'probe':
        probe = select_by_probe(mod_pdb,seed,logger,path_to_RIN)
        selfile = 'res_atoms.dat'
    elif RIN_program.lower() == 'arpeggio':
        select_by_arpeggio(pdb,seed,path_to_RIN,logger)
        selfile = 'contact_counts.dat'
    elif RIN_program.lower() == 'distance':
        print('Distance based selection scheme needs some more details')
        logger.info('Distance based selection scheme selected. Getting user input.')
        #calc_type = input("Do you want distance based calc to use average Cartesian coordinates or center of mass of the seed? (avg or mass): \n")
        calc_type = input("What type of distance? closest or avg or mass: \n")
        logger.info('Selected calc_type: '+ str(calc_type))
        hydro = input("Do you want to include protons in the distance calc? yes or no \n")
        logger.info('Protons included in distance calc: '+ str(hydro))
        cut = input("What is the cutoff distance in angstroms? \n")
        logger.info('Selected cut off distance: '+ str(cut))
        select_by_distance(calc_type,hydro,pdb,seed,cut,logger,path_to_RIN)
        selfile = 'res_atoms_by_FG.dat'
    elif RIN_program.lower() == 'manual':
        if selatom_file != '':
            selfile = selatom_file
        else:
            selfile = 'res_atoms.dat'
    print('Using ' + str(selfile) + ' as active site selection for trimming procedure')
    logger.info('Using ' + str(selfile) + ' as active site selection for trimming procedure')
    logger.info('section done\n\n')
        
        
    ### TRIMMING AND CAPPING MODEL, WRITING INPUT FILE
    logger.info('CREATING MODEL AND WRITING INPUT FILE:')
    seednum,resnum,addnum,totnum = res_atom_count(selfile,must_include,Seedlist)
    if addnum == 0:
        minsize = seednum+1
    else:
        minsize = seednum+addnum

    logger.info('Size of maximal model is: ' + str(totnum))
    print('Size of maximal model is: ' + str(totnum))
    
    option = list(range(minsize,totnum+1))
    option.append('all')
    x = False
    logger.info('Valid options for building models ' + str(option))
    logger.info('The user selected the model option in driver_input: '+ str(model_num))
    print('Model selected in driver input: ' + str(model_num))
    
    if (model_num.strip().isnumeric() == True and int(model_num) in option) or model_num.lower() in ['all','max','maximal']:
        if model_num.lower() in ['max','maximal']:
            model_num = totnum
        x = True
    else:
        print('Model selected in input file is not a valid option. Valid options are:')
        print(option)
        logger.info("The user did not input a correct model number in driver_input file. Requesting input.")    
    while x != True:
        model_num = input('Which model would you like (max = ' + str(totnum) + ')?\n')
        if model_num.isnumeric() == True and int(model_num) in option:
            logger.info('The user selected the model: '+ str(model_num))
            x = True
        elif model_num.lower() == 'all':
            logger.info('The user selected all models')
            x = True
        else:
            print('Model selected is not a valid option. Please try again.')
            logger.info("The user did not input a correct model number. Trying again.")

    seed_name = ''
    for i in Seedlist:
        seed_name+=i + ','
    
    print('Should anything be avoided in the capping protonation step? Typically the seed (' + seed_name[0:-1] + ') is avoided')
    freezeinp = input('Exclusion options: nothing/seed/manual \n')
    if freezeinp.lower() == 'nothing':
        freeze = ['','','']
    elif freezeinp.lower() == 'seed':
        freeze = [seed_name[0:-1],'','']
    else:
        print('Manual input of exclusions (see documentation for pymol_protonate.py for help):')
        freezeid = input('-ignore_ids (whole residues to ignore): ')
        freezeat = input('-ignore_atoms (specific atoms to ignore): ')
        freezename = input('-ignore_atnames (atom types to ignore): ')
        freeze = [freezeid, freezeat, freezename]
    logger.info('Residues excluded from protonation: '+ str(freeze[0]))
    if freeze[1] != '':
        logger.info('Specific atoms excluded from protonation: '+ str(freeze[1]))
    if freeze[2] != '':
        logger.info('Atom types excluded from protonation: '+ str(freeze[2]))

    if model_num=='all':
        #tot = []
        for num in range(minsize,totnum+1):
            #tot.append(num)
            trim_model(seed,must_include,mod_pdb,selfile,num,path_to_RIN,RIN_program,logger)
            protonate_model(freeze,num,path_to_RIN,logger)
            make_temp_pdb(num,path_to_RIN,logger)
            if Computational_program.lower() != 'none':
                create_input_file(template_path,Computational_program,basis_set_library,charge,str(num),path_to_RIN,logger)
                #shutil.copy('1.inp',str(num)+'.inp')
                #shutil.copy('template.pdb','template_'+str(num)+'.pdb')
    else:
        trim_model(seed,must_include,mod_pdb,selfile,model_num,path_to_RIN,RIN_program,logger)
        protonate_model(freeze,model_num,path_to_RIN,logger)
        make_temp_pdb(model_num,path_to_RIN,logger)
        if Computational_program.lower() != 'none':
            create_input_file(template_path,Computational_program,basis_set_library,charge,str(model_num),path_to_RIN,logger)       
        
    logger.info('section done\n')
         
    return

if __name__ == '__main__':
    """ Usage: RINRUS_driver.py -i driver_input """
    parser = argparse.ArgumentParser(description='generate output containing script lines from driver_input')
    parser.add_argument('-i',
                        dest='driver_input',
                        default='driver_input',
                        help='This is the driver input file')

    args = parser.parse_args()
    
    run_rinrus_driver(args.driver_input)
