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

def driver_file_reader(file,logger):
    red = 'True'
    pdb = ''
    Seed = []
    seed = ''
    RIN_program = ''
  #  Histidine = ''
    model_num = ''
    charge = ''
    multi = ''
  #  residues_not_protonated = ''
    Computational_program = ''
    template_path = ''
    basis_set_library = ''
    path_to_RIN = ''
    with open(file,'r') as fp:
        data = fp.readlines()
        for line in data:
            #print(line)
            if '#' in line:
                line
            else:
                if 'PDB' in line:
                    pdb += line.replace('PDB:','').replace('\n','').replace(' ','')
                if 'Protonate_initial' in line:
                    red = line.replace('Protonate_initial:','').replace('\n','').replace(' ','')
                if 'Seed:' in line:
                    if ',' in line:
                        line = line.replace('Seed:','').replace('\n','').replace(' ','')
                        seed+= line
                        line = line.split(',')
                        for i in line:
                            if i=='':
                                pass
                            else:
                                Seed.append(i)
                        
                    else:
                        seed+= line.replace('Seed:','').replace('\n','').replace(' ','')
                        Seed.append(line.replace('Seed:','').replace('\n','').replace(' ',''))
                if 'RIN_program' in line:
                    RIN_program += line.replace('RIN_program:','').replace('\n','').replace(' ','')
                if 'Seed_charge' in line:
                    charge += line.replace('Seed_charge:','').replace('\n','').replace(' ','')
                if 'Multiplicity' in line:
                    multi+= line.replace('Multiplicity:','').replace('\n','').replace(' ','')
                if 'Computational_program' in line:
                    Computational_program += line.replace('Computational_program:','').replace('\n','').replace(' ','')
                if 'input_template_path' in line:
                    template_path += line.replace('input_template_path:','').replace('\n','').replace(' ','')
                if 'basisset_library' in line:
                    basis_set_library += line.replace('basisset_library:','').replace('\n','').replace(' ','')
                if 'Model(s)' in line:
                    model_num += line.replace('Model(s):','').replace('\n','').replace(' ','')
                if 'path_to_type_of_RIN:' in line:
                    path_to_RIN += line.replace('path_to_type_of_RIN:', '').replace('\n','').replace(' ','')
    
    
    if red.lower() in ['false', 'f', '0', 'n']:
        red = 'False'
    elif red.lower() in ['true', 't', '1', 'y']:
        red = 'True'

    logger.info('PDB Name: ' + pdb)
    logger.info('Protonate initial PDB: ' + red)
    logger.info('RIN program: ' + RIN_program)
    logger.info('Substrate Charge: ' + str(charge))
    logger.info('Multiplicity: ' + str(multi))
    logger.info('Computational Program: ' + str(Computational_program))
    logger.info('Path to the input template file: ' + str(template_path))
    logger.info('Path to the basis set library: ' + str(basis_set_library))
    logger.info('Seed: ' + str(seed))
    logger.info('Model number selection: ' + str(model_num))
    logger.info('Path to the RINRUS scripts bin directory: ' + str(path_to_RIN))
    
    return pdb,red,Seed,RIN_program,charge,multi,Computational_program,template_path,basis_set_library,seed,model_num,path_to_RIN


def res_atom_count(filename):
    num = 0
    with open(filename,'r') as fp:
        data = fp.readlines()
        num += len(data)
        for i in data:
            if i == '':
                num -= 1
    return num


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
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('Probe run as: '+ str(' '.join(args)))
    logger.info('Return code: '+ str(out.returncode))
    logger.info('Output from Probe: '+ out.stdout)

    print('Analyzing probe RIN')
    path = os.path.expanduser(path_to_RIN+'/probe2rins.py')
    args =  sys.executable + ' ~/git/RINRUS/bin/probe2rins.py -f '+ str(probe)+ ' -s ' + seed
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    logger.info('RIN analysis run as: '+ str(args))
    logger.info('Return code= '+ str(out.returncode))
    #logger.info('Output: \n'+ out.stdout)
    return probe
    

def select_by_distance(calc_type,hydro,pdb,seed,cut,logger,path_to_RIN):
    path = os.path.expanduser(path_to_RIN+'/pdb_dist_rank.py')
    if hydro.lower() == "nohydro":
        arg = [sys.executable, path ,'-pdb',str(pdb),'-s',str(seed),'-cut',cut,'-type',calc_type,'-nohydro']
        result = subprocess.run(arg)
        logger.info('Distance selection run as: '+ str(' '.join(args)))
        logger.info('Return code: '+ str(result.returncode))
        logger.info('Output:\n'+ str(result.stdout))
    else:
        arg = [sys.executable, path ,'-pdb',str(pdb),'-s',str(seed),'-cut',cut,'-type',calc_type]
        result = subprocess.run(arg)
        logger.info('Distance selection run as: '+ str(' '.join(args)))
        logger.info('Return code=: '+ str(result.returncode))
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
    

def trim_model(seed,pdb,selfile,model_num,path_to_RIN,RIN_program,logger):
    path = os.path.expanduser(path_to_RIN+'/rinrus_trim2_pdb.py')
    args = [sys.executable,path, '-s',str(seed).replace('\n',''), '-pdb',str(pdb),'-c',str(selfile), '-model', str(model_num)]
    result = subprocess.run(args)
    logger.info('Model trimming run as: ' + str(' '.join(args)))
    return


def protonate_model(freeze,model_num,path_to_RIN,logger):
    path = os.path.expanduser(path_to_RIN+'/pymol_protonate.py')
    name = 'res_' + str(model_num)+'.pdb'
    arg= [sys.executable,path, '-ignore_ids',str(freeze),'-pdb', name]
    out = subprocess.run(arg)
    logger.info('Model protonation run as: '+ str(' '.join(arg)))
    return


def create_input_file(template,format,basisinfo,charge,model_num,path_to_RIN,logger):
    path = os.path.expanduser(path_to_RIN+'/write_input.py')
    path_2 =os.path.expanduser(basisinfo.replace('\n','').replace(' ','')) 
    noh =  ' res_'+model_num+'.pdb '
    adh = ' res_'+model_num+'_h.pdb'
    arg= [sys.executable, path ,'-intmp',str(template),'-format',str(format),'-basisinfo',path_2,'-c',str(charge),'-type', 'hopt','-noh',str(noh).replace(' ',''),'-adh',str(adh).replace(' ','')]
    result = subprocess.run(arg)
    logger.info('Write_input run as: '+ str(' '.join(arg)))
    return



def run_rinrus_driver(file):
    ### SET UP LOGGING
    dt = datetime.now().strftime("%Y-%m-%d")
    #logging.basicConfig(filename="rinrus_log.out",
    #                format='%(asctime)s %(message)s',
    #                filemode='w')
    lf = f"rinrus_log_{dt}.out"
    logging.basicConfig(level=logging.DEBUG,
                    filename=lf,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')
    logger = logging.getLogger()
    logger.info('RUNNING RINRUS DRIVER:')
    logger.info('Python executable being used: ' + sys.executable + '\n')

    
    ### PARSE INPUTS
    logger.info('PARSING DRIVER INPUT FILE:')
    driver_input_file = file
    logger.info('Inputs from: ' + file)
    pdb,red,Seed,RIN_program,charge,multi,Computational_program,template_path,basis_set_library,seed,model_num,path_to_RIN = driver_file_reader(file,logger)
    RIN_program = RIN_program.lower()
    amountofseed = len(Seed)
    model_num = model_num.strip()
    logger.info('section done\n')

    
    ### INITIAL PROTONATION
    if red == 'True':
        logger.info('PROTONATION OF INITIAL PDB:')
        mod_pdb = run_reduce(pdb,logger,path_to_RIN)
        logger.info('section done\n')
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
        calc_type = input("Do you want distance based calc to use average Cartesian coordinates or center of mass of the seed? (avg or mass): \n")
        logger.info('Selected calc_type: '+ str(calc_type))
        hydro = input("Do you want to include protons in the distance calc? (yes=hydro, no=nohydro) \n")
        logger.info('Selected proton option (nohydro or hydro): '+ str(hydro))
        cut = input("What is the cutoff distance in angstroms? \n")
        logger.info('Selected cut off distance: '+ str(cut))
        select_by_distance(calc_type,hydro,pdb,seed,cut,logger,path_to_RIN)
        cut = float(cut)
        selfile = 'res_atom-%.2f.dat'%cut
    elif RIN_program.lower() == 'manual':
        selfile = 'res_atoms.dat'
    print('Using ' + str(selfile) + 'as active site selection for trimming procedure')
    logger.info('Using ' + str(selfile) + 'as active site selection for trimming procedure')
    logger.info('section done\n')
        
        
    ### TRIMMING AND CAPPING MODEL, WRITING INPUT FILE
    logger.info('CREATING MODEL AND WRITING INPUT FILE:')
    num_lines = res_atom_count(selfile)
    logger.info('Size of maximal model is: ' + str(num_lines))
    print('Size of maximal model is: ' + str(num_lines))
    
    option = list(range(amountofseed+1,num_lines+1))
    option.append('all')
    x = False
    logger.info('Valid options for building models ' + str(option))
    logger.info('The user selected the model option in driver_input: '+ str(model_num))
    print('Model selected in driver input: ' + str(model_num))
    
    if (model_num.strip().isnumeric() == True and int(model_num) in option) or model_num.lower() in ['all','max','maximal']:
        if model_num.lower() in ['max','maximal']:
            model_num = num_lines
        x = True
    else:
        print('Model selected in input file is not a valid option. Valid options are:')
        print(option)
        logger.info("The user did not input a correct model number in driver_input file. Requesting input.")    
    while x != True:
        model_num = input('Which model would you like (max = ' + str(num_lines) + ')?\n')
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
    for i in Seed:
        seed_name+=i + ','
    
    print('Should anything be avoided in the capping protonation step? Typically the seed (' + seed_name[0:-1] + ') is avoided')
    freezeinp = input('Exclusion options: nothing/seed/[give specific list] \n')
    if freezeinp.lower() == 'nothing':
        freeze = ''
    elif freezeinp.lower() == 'seed':
        freeze = seed_name[0:-1]
    else:
        freeze = freezeinp
    logger.info('Residues excluded from protonation: '+ freeze)

    if model_num=='all':
        tot = []
        for num in range(amountofseed+1,num_lines+1):
            tot.append(num)
            trim_model(seed,mod_pdb,selfile,num,path_to_RIN,RIN_program,logger)
            protonate_model(freeze,num,path_to_RIN,logger)
            create_input_file(template_path,Computational_program,basis_set_library,charge,str(num),path_to_RIN,logger)
            shutil.copy('1.inp',str(num)+'.inp')
            shutil.copy('template.pdb','template_'+str(num)+'.pdb')
    else:
        trim_model(seed,mod_pdb,selfile,model_num,path_to_RIN,RIN_program,logger)
        protonate_model(freeze,model_num,path_to_RIN,logger)
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
