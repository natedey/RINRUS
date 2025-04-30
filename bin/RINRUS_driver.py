#!/usr/bin/env python3
"""
This is a program created by Dr. Taylor Santaloci, Dr. Dominique Wappett and the DeYonker research group
at The University of Memphis.
Initial version created 04.03.2023
Production version created 10.24.2024
"""

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

### LIST OF CURRENTLY RECOGNIZED OPTIONS           ###
### DICTIONARY STRUCTURE:                          ###
### option: [description for log file, value help] ###
opts = {'path_to_scripts': ['Path to the RINRUS scripts bin directory','dir path'],
        'pdb': ['Input PDB', 'filename'],
        'protonate_initial': ['Protonate input PDB','true OR false'],
        'seed': ['Seed', 'ch:ID[,ch:ID,...]'],
        'rin_program': ['RIN program','probe OR arpeggio OR distance OR manual'],
        'res_atoms_file': ['Specified res_atoms file','filename    (only used if rin_program = manual)'], 
        'arpeggio_rank': ['Arpeggio ranking','contacts OR types    (only used if rin_program = arpeggio)'], 
        'dist_type': ['Distance type', 'avg OR com OR closest    (only used if rin_program = distance)'], 
        'dist_satom': ['Distance satom', 'ch:ID:atom[,ch:ID:atom,...]    (only used if rin_program = distance)'], 
        'dist_max': ['Distance max radius', 'number    (only used if rin_program = distance)'],
        'dist_noh': ['Calculate distance excluding hydrogens','true OR false    (only used if rin_program = distance)'], 
        'must_add': ['Fragments that must be included', 'ch:ID[:S/:N/:C] etc'],
        'model': ['Selected model(s)', 'all OR maximal OR max OR number'],
        'model_prot_ignore_ids': ['Residues avoided in model protonation','ch:ID[,ch:ID,...]'],
        'model_prot_ignore_atoms': ['Specific atoms avoided in model protonation','ch:ID:atom[,ch:ID:atom,...]'],
        'model_prot_ignore_atnames': ['Atom types avoided in model protonation','atom[,atom,...]'],
        'qm_input_format': ['Format for QM input file(s)', 'gaussian OR orca OR qchem OR gau-xtb OR psi4-fsapt'],
        'qm_input_template': ['Specified QM input template','filename    (only used if qm_input_format defined)'],
	'gaussian_basis_intmp': ['Source of basis sets for Gaussian input', 'true OR false    (only used if qm_input_format = gaussian)'], 
        'qm_calc_hopt': ['Freeze all heavy atoms in QM input', 'true OR false   (only used if qm_input_format defined and not psi4-fsapt)'],
        'seed_charge': ['Seed charge','integer    (only used if qm_input_format defined)'],
        'multiplicity': ['Multiplicity','integer    (only used if qm_input_format defined)'],
        'fsapt_fa': ['F-SAPT fragment A','ch:ID[,ch:ID,...]    (only used if qm_input_format = psi4-fsapt)']}
# list of which opts have true/false values and need to be converted to booleans
tfopts = ['protonate_initial','dist_noh','gaussian_basis_intmp','qm_calc_hopt']

### define log header info, gets details of commit ###
def log_header(year):
    gitpath = str(Path(__file__).resolve().parents[1])
    pwd = os.getcwd()
    gitver = subprocess.run(f"cd {gitpath}; git show -s --pretty='format:%h %cd' --date=format-local:'%Y-%m-%d %H:%M'; cd {pwd}",shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    gitver = gitver.stdout.split()
    headtxt = ('--------------------------------------------------------------------------------------\n'
    '              RINRUS: The Residue Interaction Network ResidUe Selector                \n'
    '--------------------------------------------------------------------------------------\n'
    f'(C) 2018-{year}. Using version {gitver[0]}, published on github {gitver[1]} at {gitver[2]}.        \n'
    'Developed in the group of Prof. Nathan DeYonker at the University of Memphis, TN USA. \n'
    'Contributors: Q. Cheng, N. DeYonker, D. Wappett, T. Summers, D. Agbaglo, T. Suhagia,  \n'
    '    T. Santaloci, J. Bachega.                                                         \n'
    'Acknowledge RINRUS by citing: github.com/natedey/RINRUS, DOI:10.1016/j.bpj.2021.07.029\n'
    '    and DOI:10.1039/D3CP06100K                                                        \n'
    '--------------------------------------------------------------------------------------' )

    clbanner = ('--------------------------------------------------------------------------------------\n'
    '           Running RINRUS: The Residue Interaction Network ResidUe Selector           \n'
    'Developed in the group of Prof. Nathan DeYonker at the University of Memphis, TN USA. \n'
    f'(C) 2018-{year}. Using version {gitver[0]}, published on github {gitver[1]} at {gitver[2]}.        \n'
    '--------------------------------------------------------------------------------------\n')
    return headtxt,clbanner


### parse driver input file into dict of options ###
def driver_file_reader(inpfile,logger,scriptpath):
    data = open(inpfile,'r').readlines()
    
    # print input verbatim
    rawinp = ''
    for line in data:
        rawinp += '> ' + line
    logger.info('Inputs being read from: '+inpfile+'\n'+f'----------raw contents----------\n'+rawinp+'--------------------------------')

    # remove comments and empty lines from input and turn into dict
    data = [line for line in data if not line.startswith('#') and line.strip()]
    inpdict = {line.split(':',1)[0].lower(): line.split(':',1)[1].strip() for line in data}

    # check dict
    checked_dict = {}
    for key in inpdict.keys():
        # ignore if no value given for option
        if not inpdict[key]:
            continue
        # convert relevant options to bool values
        if key in tfopts:
            if inpdict[key].lower() in ['true', 't', '1']:
                checked_dict[key] = True
            else:
                checked_dict[key] = False
        # add other valid opts to dict
        elif key in opts.keys():
            checked_dict[key] = inpdict[key]
        # print warnings for unrecognized options
        else:
            logger.info(f'Unrecognized option {key} in driver input will be ignored.')
    # warn user and quit if any of the 4 necessary options have not been specified
    if 'pdb' not in checked_dict.keys() or 'seed' not in checked_dict.keys() or 'rin_program' not in checked_dict.keys() or 'model' not in checked_dict.keys():
        print('Input needs to contain PDB, seed, RIN_program and model options at minimum. Please check input!')
        logger.info('Input needs to contain PDB, seed, RIN_program and model options at minimum. Quitting RINRUS')
        sys.exit()
    # add split up seed list to dict
    checked_dict['seedlist'] = [s for s in checked_dict['seed'].split(',') if s]    
    # remove any suboptions that don't match parent option
    if checked_dict['rin_program'].lower() != 'arpeggio' and 'arpeggio_rank' in checked_dict.keys():
        del checked_dict['arpeggio_rank']
    if checked_dict['rin_program'].lower() != 'distance':
        for i in ['dist_type','dist_satom','dist_max','dist_noh']:
            if i in checked_dict.keys():
                del checked_dict[i]
    if checked_dict['rin_program'].lower() != 'manual' and 'res_atoms_file' in checked_dict.keys():
        del checked_dict['res_atoms_file']
    if 'qm_input_format' in checked_dict.keys():
        if checked_dict['qm_input_format'] != 'psi4-fsapt' and 'fsapt_fa' in checked_dict.keys():
            del checked_dict['fsapt_fa']
        if checked_dict['qm_input_format'].lower() != 'gaussian' and 'gaussian_basis_intmp' in checked_dict.keys():
            del checked_dict['gaussian_basis_intmp']
    else:
        for i in ['qm_input_template','gaussian_basis_intmp','qm_calc_hopt','seed_charge','multiplicity','fsapt_fA']:
            if i in checked_dict.keys():
                del checked_dict[i]
    # flag to use gaussian basis info in template file or not
    if 'qm_input_format' in checked_dict.keys() and checked_dict['qm_input_format'].lower() == 'gaussian':
        if 'gaussian_basis_intmp' in checked_dict.keys() and checked_dict['gaussian_basis_intmp']:
            checked_dict['gaussian_basis_intmp'] = 'basis in template'
        else:
            checked_dict['gaussian_basis_intmp'] = 'default basis dict'
    # clean up script path and/or define from RINRUS_driver.py file location 
    if 'path_to_scripts' in checked_dict.keys():
        if checked_dict['path_to_scripts'].startswith('~'):
            checked_dict['path_to_scripts'] = os.path.expanduser(checked_dict['path_to_scripts'])
        elif checked_dict['path_to_scripts'].startswith('$'):
            checked_dict['path_to_scripts'] = os.path.expandvars(checked_dict['path_to_scripts'])
    else:
        checked_dict['path_to_scripts'] = str(scriptpath)
    if not checked_dict['path_to_scripts'].endswith('/'):
        checked_dict['path_to_scripts'] += '/'

    # logging:
    logger.info('Inputs parsed and checked. Using options:')
    for key in checked_dict.keys():
        if key in opts.keys():
            logger.info(f'{opts[key][0]}: {checked_dict[key]}')

    return checked_dict


### INTERNAL DRIVER FUNCTIONS ###

def run_reduce(pdb,logger,path_to_scripts):
    """
    runs reduce with -NOFLIP flag and uses logger to log output
    """
    print('Protonating input PDB file')
    path = os.path.expanduser(path_to_scripts+'reduce')
    pdb_2 = pdb.replace('.pdb','')
    args = path_to_scripts + '/reduce -NOFLIP -Quiet  '+ str(pdb)+ ' > '+ str(pdb_2)+'_h.pdb'
    io.StringIO(initial_value='', newline='\r')
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Reduce run as: '+ str(out.args))
    if out.stdout:
        logger.info('Output from Reduce: \n'+ out.stdout)
    if out.returncode == 1:
        logger.info('Return code = 1')
        logger.info('Quitting RINRUS')
        print('Error in Reduce\n'+str(out.stdout))
        print('Quitting RINRUS')
        sys.exit()
    mod_pdb = str(pdb_2)+'_h.pdb'
    return mod_pdb

def select_by_probe(pdb,seed,logger,path_to_scripts):
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
        path_to_scripts (_type_): path to where file is

    Returns:
        _type_: .probe file
    """    
    print('Generating probe RIN')
    probe = pdb.replace('.pdb','.probe')
    args = [path_to_scripts+'probe -unformated -MC -self "all" -Quiet '+ pdb +' > '+ probe]
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Probe run as: '+ str(' '.join(args)))
    if out.stdout:
        logger.info('Output from probe: \n'+ str(out.stdout))
    if out.returncode == 1:
        logger.info('Return code = 1')
        logger.info('Quitting RINRUS')
        print('Error in probe\n'+str(out.stdout))
        print('Quitting RINRUS')
        sys.exit()
    print('Analyzing probe RIN')
    path = os.path.expanduser(path_to_scripts+'probe2rins.py')
    args =  [sys.executable,path,'-f',str(probe),'-s',seed]
    out = subprocess.run(args,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('RIN analysis run as: '+ str(' '.join(args[1:])))
    return probe
    
def select_by_arpeggio(pdb,seed,path_to_scripts,logger):
    print('Generating arpeggio RIN')
    path = os.path.expanduser(path_to_scripts+'arpeggio/arpeggio.py')
    arg = [sys.executable,path,str(pdb)]
    result = subprocess.run(arg,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Arpeggio run as: ' + str(' '.join(arg[1:])))
    if result.stdout:
        logger.info('Output from arpeggio: \n'+ result.stdout)
    if result.returncode == 1:
        logger.info('Return code = 1')
        logger.info('Quitting RINRUS')
        print('Error in arpeggio\n'+str(result.stdout))
        print('Quitting RINRUS')
        sys.exit()
    print('Analyzing arpeggio RIN')
    path = os.path.expanduser(path_to_scripts+'arpeggio2rins.py')
    arg = [sys.executable,path,'-f',str(pdb).replace('pdb','contacts'),'-s',seed]
    result = subprocess.run(arg,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('RIN analysis run as: '+ str(' '.join(arg[1:])))
    return

def select_by_distance(checked_dict,logger):
    path = os.path.expanduser(checked_dict['path_to_scripts']+'dist_rank.py')
    arg = [sys.executable, path ,'-pdb',checked_dict['pdb'],'-type',checked_dict['dist_type']]
    if 'dist_satom' in checked_dict.keys():
        arg.append('-satom')
        arg.append(checked_dict['dist_satom'])
    else:
        arg.append('-seed')
        arg.append(checked_dict['seed'])
    if 'dist_max' in checked_dict.keys():
        arg.append('-max')
        arg.append(checked_dict['dist_max'])
    if 'dist_noh' in checked_dict.keys():
        arg.append('-noH')
    result = subprocess.run(arg,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Distance selection run as: '+ str(' '.join(arg[1:])))
    return
    
def res_atom_count(filename,must_add,Seedlist):
    seednum = len(Seedlist)
    # residues in res_atoms
    resnum = 0
    with open(filename,'r') as fp:
        data = fp.readlines()
        for i in data:
            key = i.split()[0] + ':' + i.split()[1]
            if i != '' and not i.startswith('#') and key not in Seedlist:
                resnum += 1
    # must_add fragments
    if must_add != '':
        addnum = len(must_add.split(','))
    else:
        addnum = 0
    totnum = seednum + resnum + addnum
    return seednum,resnum,addnum,totnum

def trim_model(checked_dict,model_num,selfile,logger):
    path = os.path.expanduser(checked_dict['path_to_scripts']+'rinrus_trim2_pdb.py')
    args = [sys.executable,path, '-s',checked_dict['seed'], '-pdb',checked_dict['pdb'],'-ra',str(selfile), '-model', str(model_num)]
    if 'must_add' in checked_dict.keys():
        args.append('-mustadd')
        args.append(checked_dict['must_add'])
    result = subprocess.run(args,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Model trimming run as: ' + str(' '.join(args[1:])))
    return

def protonate_model(checked_dict,model_num,logger):
    path = os.path.expanduser(checked_dict['path_to_scripts']+'pymol_protonate.py')
    name = 'res_' + str(model_num)+'.pdb'
    arg= [sys.executable,path,'-pdb', name]
    if 'model_prot_ignore_ids' in checked_dict.keys():
        arg.append('-ignore_ids')
        arg.append(checked_dict['model_prot_ignore_ids'])
    if 'model_prot_ignore_atoms' in checked_dict.keys():
        arg.append('-ignore_atoms')
        arg.append(checked_dict['model_prot_ignore_atoms'])
    if 'model_prot_ignore_atnames' in checked_dict.keys():
        arg.append('-ignore_atnames')
        arg.append(checked_dict['model_prot_ignore_atnames'])
    out = subprocess.run(arg,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Model protonation run as: '+ str(' '.join(arg[1:])))
    return

def make_temp_pdb(model_num,path_to_scripts,logger):
    path = os.path.expanduser(path_to_scripts+'make_template_pdb.py')
    arg = [sys.executable,path,'-model',str(model_num)]
    out = subprocess.run(arg,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Creation of template pdb run as: '+ str(' '.join(arg[1:])))
    return

def create_input_file(checked_dict,model_num,logger):
    path = os.path.expanduser(checked_dict['path_to_scripts']+'write_input.py')
    modpdb = f'model_{model_num}_template.pdb'
    inpn = f'{model_num}.inp'
    arg = [sys.executable,path,'-format',checked_dict['qm_input_format'],'-c',checked_dict['seed_charge'],'-pdb',modpdb,'-inpn',inpn]
    if 'qm_input_template' in checked_dict.keys():
        arg.append('-intmp')
        arg.append(checked_dict['qm_input_template'])
    if checked_dict['qm_input_format'].lower() == 'gaussian' and checked_dict['gaussian_basis_intmp'] == 'basis in template':
        arg.append('-basisinfo')
        arg.append('intmp')
    if checked_dict['qm_input_format'].lower() == 'psi4-fsapt':
        arg.append('-fA')
        arg.append(checked_dict['fsapt_fa'])
    if 'qm_input_hopt' in checked_dict.keys() and checked_dict['qm_input_hopt']:
        arg.append('-type')
        arg.append('hopt')
    result = subprocess.run(arg,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Write_input run as: '+ str(' '.join(arg[1:])))
    return


### DRIVER ITSELF ###
def run_rinrus_driver(inpfile,scriptpath):
    ### SET UP LOGGING
    dt = datetime.now().strftime("%Y-%m-%d")
    lf = f"rinrus_log_{dt}.out"
    logging.basicConfig(level=logging.DEBUG,
                    filename=lf,
                    format='%(message)s',
                    filemode='w')
    logger = logging.getLogger()

    #write header 
    header,clbanner = log_header(dt[0:4])
    logger.info(header+'\n\n')
    print(clbanner)

    #change to date + message format for main part of log
    logger.handlers[0].setFormatter(logging.Formatter(fmt='%(asctime)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S'))
    logger.info('RUNNING RINRUS DRIVER:')
    logger.info('Python executable being used: ' + sys.executable + '\n\n')

    
    ### PARSE INPUTS
    logger.info('READING DRIVER INPUT FILE:')
    checked_dict = driver_file_reader(inpfile,logger,scriptpath)
    logger.info('section done\n\n')

    
    ### INITIAL PROTONATION
    if checked_dict['protonate_initial']:
        logger.info('PROTONATION OF INITIAL PDB:')
        mod_pdb = run_reduce(checked_dict['pdb'],logger,checked_dict['path_to_scripts'])
        checked_dict['pdb'] = modpdb
        logger.info('section done\n\n')


    ### SELECTING ACTIVE SITE
    logger.info('SELECTING ACTIVE SITE:')
    logger.info('Using: '+checked_dict['rin_program'])
    if checked_dict['rin_program'].lower() == 'probe':
        probe = select_by_probe(checked_dict['pdb'],checked_dict['seed'],logger,checked_dict['path_to_scripts'])
        selfile = 'res_atoms.dat'
    elif checked_dict['rin_program'].lower() == 'arpeggio':
        select_by_arpeggio(checked_dict['pdb'],checked_dict['seed'],checked_dict['path_to_scripts'],logger)
        if 'arpeggio_rank' in checked_dict.keys() and checked_dict['arpeggio_rank'] == 'counts':
            selfile = 'contact_counts.dat'
        elif 'arpeggio_rank' in checked_dict.keys() and checked_dict['arpeggio_rank'] == 'types':
            selfile = 'contype_counts.dat'
        else: 
            arprank = input('Build models with arpeggio contact count or type ranking? (counts/types)\n')
            if arprank == 'counts':
                selfile = 'contact_counts.dat'
            elif arprank == 'types':
                selfile = 'contype_counts.dat'
                logger.info('Command line arpeggio_rank input: ' + arprank)
    elif checked_dict['rin_program'].lower() == 'distance':
        if not checked_dict['dist_type']:
            checked_dict['dist_type'] = input("Distance type not specified in input file: select closest or avg or mass \n")
            logger.info('Command line dist_type input: '+ str(checked_dict['dist_type']))
        select_by_distance(checked_dict,logger)
        selfile = 'res_atoms_by_FG.dat'
    elif RIN_program.lower() == 'manual':
        if checked_dict['res_atoms_file']:
            selfile = checked_dict['res_atoms_file']
        else:
            selfile = 'res_atoms.dat'
    print('Using ' + str(selfile) + ' as active site selection for trimming procedure')
    logger.info('Using ' + str(selfile) + ' as active site selection for trimming procedure')
    logger.info('section done\n\n')
        
        
    ### TRIMMING AND CAPPING MODEL, WRITING INPUT FILE
    logger.info('CREATING MODEL AND WRITING INPUT FILE:')
    if 'must_add' in checked_dict.keys():
        seednum,resnum,addnum,totnum = res_atom_count(selfile,checked_dict['must_add'],checked_dict['seedlist'])
    else:
        seednum,resnum,addnum,totnum = res_atom_count(selfile,'',checked_dict['seedlist'])
    if addnum == 0:
        minsize = seednum+1
    else:
        minsize = seednum+addnum

    logger.info('Size of maximal model is: ' + str(totnum))
    print('Size of maximal model is: ' + str(totnum))
    
    option = list(range(minsize,totnum+1))
    x = False
    try:
        model_num = checked_dict['model'].lower()
    except:
        model_num = checked_dict['model']
    logger.info('Valid model sizes ' + str(option))
    logger.info('The user selected the model option in driver_input: '+ str(model_num))
    print('Model selected in driver input: ' + str(model_num))
    
    option.append('all')
    if (model_num.strip().isnumeric() == True and int(model_num) in option) or model_num in ['all','max','maximal']:
        if model_num in ['max','maximal']:
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
    checked_dict['model'] = model_num
    
    # check for protonation exclusions if not specified in input file
    if not {'model_prot_ignore_ids', 'model_prot_ignore_atoms', 'model_prot_ignore_atnames'} & set(checked_dict.keys()):
        logger.info('No "model_prot_ignore..." options specified in input. Command line input requested')
        print(f'No "model_prot_ignore..." options specified in input. Warning: PyMOL may incorrectly change protonation of the ligand/non-canonical residues.')
        freezeinp = input('Exclude seed (S) or manually select exclusions (M) or continue anyway (): \n')
        if freezeinp.lower() == 's':
            checked_dict['model_prot_ignore_ids'] = checked_dict['seed']
        elif freezeinp.lower() == 'm':
            print('Manual input of exclusions (see documentation for pymol_protonate.py for help):')
            checked_dict['model_prot_ignore_ids'] = input('-ignore_ids (whole residues to ignore): ')
            checked_dict['model_prot_ignore_atoms'] = input('-ignore_atoms (specific atoms to ignore): ')
            checked_dict['model_prot_ignore_atnames'] = input('-ignore_atnames (atom types to ignore): ')
        else:
            logger.info('Continuing without anything excluded from protonation')
    if 'model_prot_ignore_ids' in checked_dict.keys():
        logger.info('Residues excluded from protonation: '+ checked_dict['model_prot_ignore_ids'])
    if 'model_prot_ignore_atoms' in checked_dict.keys():
        logger.info('Specific atoms excluded from protonation: '+ checked_dict['model_prot_ignore_atoms'])
    if 'model_prot_ignore_atnames' in checked_dict.keys():
        logger.info('Atom types excluded from protonation: '+ checked_dict['model_prot_ignore_atnames'])

    # check F-SAPT selections
    if 'qm_input_format' in checked_dict.keys() and checked_dict['qm_input_format'] == 'psi4-fsapt':
        if 'fsapt_fa' not in checked_dict.keys():
            fachk=input(f'F-SAPT selected but fA not defined in input file. Please specify fA or press enter to use seed ({checked_dict["seed"]})\n')
            if fachk:
                checked_dict['fsapt_fa'] = fachk
                logger.info(f'F-SAPT fA command line input: {fachk}')
            else:
                checked_dict['fsapt_fa'] = checked_dict['seed']
                logger.info('F-SAPT fA command line input: seed')
        if checked_dict['model']=='all':
            print('Warning: you have selected all models and F-SAPT input files.')
            fsaptcheck = input('Do you really need F-SAPT input files for all models or only the maximal model? (all/max) \n')
            logger.info(f'Checked if F-SAPT input files for all models needed. Command line input: {fsaptcheck}')

    if checked_dict['model']=='all':
        for num in range(minsize,totnum+1):
            logger.info(f'Making model {num}')
            print(f'Making model {num}')
            trim_model(checked_dict,num,selfile,logger)
            protonate_model(checked_dict,num,logger)
            make_temp_pdb(num,checked_dict['path_to_scripts'],logger)
            if 'qm_input_format' in checked_dict.keys():
                if checked_dict['qm_input_format'] != 'psi4-fsapt' or fsaptcheck.lower() == 'all':
                    create_input_file(checked_dict,str(num),logger)
                elif num == totnum:
                    create_input_file(checked_dict,str(num),logger)
    else:
        print(f'Making model {checked_dict["model"]}')
        trim_model(checked_dict,checked_dict['model'],selfile,logger)
        protonate_model(checked_dict,checked_dict['model'],logger)
        make_temp_pdb(checked_dict['model'],checked_dict['path_to_scripts'],logger)
        if 'qm_input_format' in checked_dict.keys():
            create_input_file(checked_dict,checked_dict['model'],logger)       
        
    logger.info('section done\n')

    logger.info('RINRUS run complete\n')
    return


if __name__ == '__main__':
    """ Usage: RINRUS_driver.py -i rinrus.inp """

    ml = max([len(key) for key in opts.keys()])
    opthelp = 'Recognized keywords and values in input file:\n'
    for key in opts.keys():
        opthelp += f'  {key+":":<{ml+1}} {opts[key][1]}\n'

    parser = argparse.ArgumentParser(description='RINRUS driver module', epilog=opthelp, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', dest='driver_input', default='rinrus.inp', help='RINRUS driver input file')

    args = parser.parse_args()
    
    scriptpath = Path(__file__).parent
    run_rinrus_driver(args.driver_input,scriptpath)
