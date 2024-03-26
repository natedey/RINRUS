import os
import argparse
import shutil
import shlex,subprocess
from pathlib import Path
import sys
from subprocess import Popen, PIPE,STDOUT
import logging
import io


def driver_file_reader(file,logger):
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
                if 'Substrate(s)_charge' in line:
                    charge += line.replace('Substrate(s)_charge:','').replace('\n','').replace(' ','')
                if 'Multiplicity' in line:
                    multi+= line.replace('Multiplicity:','').replace('\n','').replace(' ','')
                if 'Computational_program' in line:
                    Computational_program += line.replace('Computational_program:','').replace('\n','').replace(' ','')
                if 'input_template_path' in line:
                    template_path += line.replace('input_template_path:','').replace('\n','').replace(' ','')
                if 'basisset_library' in line:
                    basis_set_library += line.replace('basisset_library:','')
                if 'Model(s)' in line:
                    model_num += line.replace('Model(s):','')
                if 'path_to_type_of_RIN:' in line:
                    path_to_RIN += line.replace('path_to_type_of_RIN:', '').replace('\n','').replace(' ','')

    logger.info('PDB Name: ' + pdb)
    logger.info('RIN program: ' + RIN_program)
    logger.info('Substrate Charge: ' + str(charge))
    logger.info('Multiplicity: ' + str(multi))
    logger.info('Computational Program: ' + str(Computational_program))
    logger.info('Path to the input template file: ' + str(template_path))
    logger.info('Path to the basis set library: ' + str(basis_set_library))
    logger.info('Seed: ' + str(seed))
    logger.info('Model number selection: ' + str(model_num))
    logger.info('Path to the RINRUS scripts bin directory: ' + str(path_to_RIN))
    
                    
                    #print(line)
    return pdb,Seed,RIN_program,charge,multi,Computational_program,template_path,basis_set_library,seed,model_num,path_to_RIN
def commands_step1(pdb,logger,path_to_RIN):
    """
    runs reduce with -NOFLIP flag and uses logger to log output
    """
    path = os.path.expanduser(path_to_RIN+'/reduce')
    pdb_2 = pdb.replace('.pdb','')
    args = path_to_RIN + '/reduce -NOFLIP -Quiet  '+ str(pdb)+ ' > '+ str(pdb_2)+'_h.pdb'
    io.StringIO(initial_value='', newline='\r')
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    a = vars(out)
    print(out.stdout)
    logger.info('The reduce command inputted: '+ str(out.args))
    logger.info('return code= '+ str(out.returncode))
    logger.info('Output from Reduce:\n'+ out.stdout)
    shutil.copy(str(pdb_2)+'_h.pdb',str(pdb_2)+'_h_modify.pdb')
    mod_pdb = str(pdb_2)+'_h_modify.pdb'
    return mod_pdb


def commands_step2(pdb,logger,path_to_RIN):
    """_summary_
    Hardcoded Flags:
        -unformated
        -MC
        -self "all"
        -Quiet
    Args:
        pdb (_type_): _h modified pdb
        logger (_type_): shows the probe input, return code and output from the probe program
        path_to_RIN (_type_): path to where file is

    Returns:
        _type_: .probe file
    """
    
    probe = pdb.replace('.pdb','')
    args = [path_to_RIN+'/probe -unformated -MC -self "all" -Quiet '+ pdb +' > '+ probe + '.probe']
    probe = probe + '.probe'
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    a = vars(out)
    print(out.stdout)
    logger.info('The inputted probe: '+ str(out.args))
    logger.info('return code= '+ str(out.returncode))
    logger.info('Output from Probe:\n'+ out.stdout)
    return probe


def commands_step3(probe,seed,logger,path_to_RIN):
    """_summary_

    Args:
        probe (_type_): _description_
        seed (_type_): _description_
        logger (_type_): _description_
        path_to_RIN (_type_): _description_

    Returns:
        nothing
    """
    print(probe)
    path = os.path.expanduser(path_to_RIN+'/probe2rins.py')
    #os.system('python3 ~/git/RINRUS/bin/probe2rins.py -f '+ str(probe)+ ' -s ' + seed)
    #args = ['python3',path, '-f',str(probe),'-s',seed]
    #args = ['python3',path, '-f',str(probe),'-s','A:263']
    #args_2 = ['python3',path, '-f',probe,'-s',seed]
    args =  'python3 ~/git/RINRUS/bin/probe2rins.py -f '+ str(probe)+ ' -s ' + seed
    #print(args)
    #out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    a = vars(out)
    print(out.stdout)
    logger.info('The inputted probe2rin command: '+ str(out.args))
    logger.info('return code= '+ str(out.returncode))
    logger.info('Output :\n'+ out.stdout)
    
    return

def res_atom_count(filename):
    num = 0
    with open(filename,'r') as fp:
        data = fp.readlines()
        num += len(data)
        for i in data:
            if i == '':
                num -= 1
    return num

def commands_step4(seed,pdb,model_num,path_to_RIN,RIN_program,logger):
    #print(pdb)
    if RIN_program.lower() == 'arpeggio':
        path = os.path.expanduser(path_to_RIN+'/rinrus_trim2_pdb.py')
        args = ['python3',path, '-s',str(seed).replace('\n',''), '-pdb',str(pdb),'-c','contact_counts.dat', '-model', str(model_num)]
        logger.info('The model number is:' + str(model_num))
        result = subprocess.run(args)
        logger.info('The inputted rinrus trim path: ' + str(result.args))
    else: 
        path = os.path.expanduser(path_to_RIN+'/rinrus_trim2_pdb.py')
        args = ['python3',path, '-s',str(seed).replace('\n',''), '-pdb',str(pdb), '-model', str(model_num)]
        result = subprocess.run(args)
        logger.info('The model number is:' + str(model_num))
        logger.info('The inputted rinrus trim path: ' + str(result.args))
    #out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    #logger.info('The inputted rinrus_trim2_pdb.py command: '+ str(out.args))
    #logger.info('return code= '+ str(out.returncode))
    #logger.info('Output :\n'+ out.stdout)
    
    return

def commands_step5(freeze,model_num,path_to_RIN,logger):
    
    path = os.path.expanduser(path_to_RIN+'/pymol_scripts.py')
    name = 'res_' + str(model_num)+'.pdb'
    #loc = str(freeze).find(':')
    freeze_2 = ''
    for i in freeze:
        if i.isnumeric()==True:
            freeze_2 += i
        elif i == ',':
            freeze_2 += i
    #print(freeze_2)
    arg= ['python3',path, '-ignore_ids',str(freeze_2),'-pdbfilename', name]
    #out = subprocess.run(args,shell=True,stdout=PIPE,stderr=STDOUT,universal_newlines=True)
    out = subprocess.run(arg)
    logger.info('The inputted pymol script path: '+ str(out.args))
    #logger.info('The inputted pymol_scripts.py command: '+ str(out.args))
    #logger.info('return code= '+ str(out.returncode))
    #logger.info('Output :\n'+ out.stdout)
    #args.append('\n')
    return

def command_step6(template,format,basisinfo,charge,model_num,path_to_RIN,logger):
    
    path = os.path.expanduser(path_to_RIN+'/write_input.py')
    path_2 =os.path.expanduser(basisinfo.replace('\n','').replace(' ','')) 
    noh =  ' res_'+model_num+'.pdb '
    adh = ' res_'+model_num+'_h.pdb'
    arg= ['python3', path ,'-intmp',str(template),'-format',str(format),'-basisinfo',path_2,'-c',str(charge),'-noh',str(noh).replace(' ',''),'-adh',str(adh).replace(' ','')]
    result = subprocess.run(arg)
    logger.info('write inputs output:'+ str(result.args))
   # arg= ['python3', path ,'-intmp',str(template),'-format',str(format),'-basisinfo',path_2,'-c','-2','-noh',str(noh).replace(' ',''),'-adh',str(adh).replace(' ','')]
   # args.append('\n')
    return


def distance(calc_type,hydro,pdb,seed,cut,logger,path_to_RIN):
    path = os.path.expanduser(path_to_RIN+'/pdb_dist_rank.py')
    if hydro.lower() == "nohydro":
        arg = ['python3', path ,'-pdb',str(pdb),'-s',str(seed),'-cut',cut,'-type',calc_type,'-nohydro']
        result = subprocess.run(arg)
        logger.info('The inputted distance command: '+ str(result.args))
        logger.info('return code= '+ str(result.returncode))
        logger.info('Output :\n'+ str(result.stdout))
    else:
        arg = ['python3', path ,'-pdb',str(pdb),'-s',str(seed),'-cut',cut,'-type',calc_type]
        result = subprocess.run(arg)
        logger.info('The inputted distance command: '+ str(result.args))
        logger.info('return code= '+ str(result.returncode))
        logger.info('Output :\n'+ str(result.stdout))
    return

def arpreggio(pdb,seed,path_to_RIN,logger,):
    print(path_to_RIN)
    path = os.path.expanduser(path_to_RIN+'/arpeggio/arpeggio.py')
    arg = ['python3',path,str(pdb)]
    logger.info(str(arg))
    result = subprocess.run(arg)
    #logger.info('arpeggio.py output' + result.args)
    path = os.path.expanduser(path_to_RIN+'/arpeggio2rins.py')
    arg = ['python3',path,'-f',str(pdb).replace('pdb','contacts'),'-s',seed]
    result = subprocess.run(arg)
    logger.info(str(arg))
    #logger.info('arpeggio2rins.py output' + result.args)
    return

def main(file,nor):
    logging.basicConfig(filename="rinrus_log.out",
                    format='%(asctime)s %(message)s',
                    filemode='w')
    logger = logging.getLogger()
    # Creating an object
    logger = logging.getLogger()
 
# Setting the threshold of logger to DEBUG
    logger.setLevel(logging.DEBUG)
    print(nor)
    driver_input_file = file
    pdb,Seed,RIN_program,charge,multi,Computational_program,template_path,basis_set_library,seed,model_num,path_to_RIN= driver_file_reader(file,logger)
    RIN_program = RIN_program.lower()
    amountofseed = len(Seed)
    print(amountofseed)
    mod_pdb = pdb
    if nor == 'False':
        mod_pdb = commands_step1(pdb,logger,path_to_RIN)
    if nor == 'True':
        mod_pdb=pdb
    if RIN_program.lower()== 'probe':
        probe = commands_step2(mod_pdb,logger,path_to_RIN)
        commands_step3(probe,seed,logger,path_to_RIN)
        num_lines = res_atom_count('res_atoms.dat')
        logger.info('The Maximum amount of lines in res_atom.dat file is ' + str(num_lines))
        option = list(range(amountofseed+1,num_lines+1))
        option.append('all')
        logger.info('The options for building models ' + str(option))
        print('Insert '+ str(num_lines)+ ' for the largest model')
        print('Other options are listed below')
        print(option)
        x = False
        a = model_num.strip().isnumeric()
        print(model_num)
        print(a)
        if a == True:
            if int(model_num) in option:
                    logger.info('The user selected the model option: '+ str(model_num))
                    x = True
            else:
                print("The user did not input a correct model mumber in driver_input file")
                logger.info("The user did not input a correct model mumber in driver_input file")
        while x != True:
            model_num = input('What models number would you like? (type "all" if you want all of the models ) ')
            a = model_num.isnumeric()
            if a==True:
                    model_num = int(model_num)
                    if model_num in option:
                        model_num = str(model_num)
                        logger.info('The user selected the model option: '+ str(model_num))
                        x = True
                    else:
                        logger.error('User did not select a valid option. User selection was: '+ str(model_num))
                        logger.info('The options are '+ str(option))
                        print(option)
                        print("try again")
                        x =  False
            elif model_num.lower()== 'all':
                logger.info('The user selected the model option: '+ str(model_num))
                x = True
            else:
                logger.error('User did not select a valid option. User selection was: '+ str(model_num))
                print('4')
                logger.info('The options are '+ str(option))
                print(option)
                print("try again")
                x = False
        seed_name = ''
        for i in Seed:
            seed_name+=i + ','
        print('copy and paste to tell pymol not to protonate the seed: '+  seed_name[0:-1])
        logger.info('copy and paste to tell pymol not to protonate the seed: '+  seed_name[0:-1])
        freeze = input("What residues do you not want PyMol to protonate? (Typically, this is the seed) ")
        logger.info('The user froze: '+  freeze)
        if model_num=='all':
            num_lines = res_atom_count('res_atoms.dat')
            tot = []
            for num in range(amountofseed+1,num_lines+1):
                print(num)
                tot.append(num)
                commands_step4(seed,mod_pdb,num,path_to_RIN,RIN_program,logger)
                commands_step5(freeze,num,path_to_RIN,logger)
                command_step6(template_path,Computational_program,basis_set_library,charge,str(num),path_to_RIN,logger)
                shutil.copy('1.inp',str(num)+'.inp')
                shutil.copy('template.pdb','template_'+str(num)+'.pdb')
        else:
            commands_step4(seed,mod_pdb,model_num,path_to_RIN,RIN_program,logger)
            commands_step5(freeze,model_num,path_to_RIN,logger)
            command_step6(template_path,Computational_program,basis_set_library,charge,str(model_num),path_to_RIN,logger)
        

    if RIN_program.lower() == 'arpeggio':
        seed_name = ''
        for i in Seed:
            seed_name+=i + ','
        print('copy and paste to tell pymol not to protonate the seed: '+  seed_name[0:-1])
        logger.info('copy and paste to tell pymol not to protonate the seed: '+  seed_name[0:-1])
        freeze = input("What residues do you not want PyMol to protonate? (Typically, this is the seed) ")
        logger.info('The user froze: '+  freeze)
        arpreggio(pdb,seed,path_to_RIN,logger)
        model_num = input('What model number would you like? (type "all" if you want all of the models ) ')
        if model_num=='all':
            num_lines = res_atom_count('contact_counts.dat')
            tot = []
            for num in range(amountofseed+1,num_lines+1):
                print(num)
                tot.append(num)
                commands_step4(seed,pdb,num,path_to_RIN,RIN_program,logger)
                commands_step5(freeze,num,path_to_RIN,logger)
                command_step6(template_path,Computational_program,basis_set_library,charge,str(num),path_to_RIN,logger)
                shutil.copy('1.inp',str(num)+'.inp')
                shutil.copy('template.pdb','template_'+str(num)+'.pdb')
        else:
            commands_step4(seed,pdb,model_num,path_to_RIN,RIN_program,logger)
            commands_step5(freeze,model_num,path_to_RIN,logger)
            command_step6(template_path,Computational_program,basis_set_library,charge,model_num,path_to_RIN,logger)
    if RIN_program.lower() == 'manual':
        num_lines = res_atom_count('res_atoms.dat')
        option = list(range(amountofseed+1,num_lines+1))
        option.append('all')
        print('Insert '+ str(num_lines)+ ' for the largest model')
        print('Other options are listed below')
        print(option)
        model_num = input('What model number would you like? (type "all" if you want all of the models ) ')
        seed_name = ''
        for i in Seed:
            seed_name+=i + ','
        print('copy and paste to tell pymol not to protonate the seed: '+  seed_name[0:-1])
        logger.info('copy and paste to tell pymol not to protonate the seed: '+  seed_name[0:-1])
        freeze = input("What residues do you not want PyMol to protonate? (Typically, this is the seed) ")
        logger.info('The user froze: '+  freeze)
        if model_num=='all':
            tot = []
            for num in range(amountofseed+1,num_lines+1):
                print(num)
                tot.append(num)
                commands_step4(seed,mod_pdb,num,path_to_RIN,RIN_program,logger)
                commands_step5(freeze,num,path_to_RIN,logger)
                command_step6(template_path,Computational_program,basis_set_library,charge,str(num),path_to_RIN,logger)
                shutil.copy('1.inp',str(num)+'.inp')
                shutil.copy('template.pdb','template_'+str(num)+'.pdb')
        else:
            commands_step4(seed,mod_pdb,model_num,path_to_RIN,RIN_program,logger)
            commands_step5(freeze,model_num,path_to_RIN,logger)
            command_step6(template_path,Computational_program,basis_set_library,charge,model_num,path_to_RIN,logger)
        print("Program assumes user has already created there own res_atom.dat file")
        
    if RIN_program.lower() == 'distance':
        calc_type = input("Do you want distance based calc to be average Cartesian coordinates or center of mass of the seed? (avg or mass): ")
        logger.info('User inputted calc_type '+ str(calc_type))
        hydro = input("Do you want to include protons in the distance calc? (yes=hydro, no=nohydro) ")
        logger.info('User inputted whether nohydro or hydro '+ str(hydro))
        cut = input("What is the cutoff distance in angstroms? ")
        logger.info('User inputted cut off distance '+ str(cut))
        distance(calc_type,hydro,pdb,seed,cut,logger,path_to_RIN)
        print('For cluster model creation, create your own res_atoms and change RIN program to manual. This program does not run center of mass or avg on center_atoms')
        logger.info('For cluster model creation, create your own res_atoms and change RIN program to manual. This program does not run center of mass or avg on center_atoms')
        
         
    return

if __name__ == '__main__':
    """ Usage: input_processing.py -i driver_input """
    parser = argparse.ArgumentParser(description='generate output containing script lines from driver_input')
    parser.add_argument('-i',
                        dest='driver_input',
                        default='driver_input',
                        help=' This is the driver_input file')
    parser.add_argument('-nor',
                        dest='NOR',
                        default='False',
                        help='If true the the pdb will not be protonated')

args = parser.parse_args()
nor = args.NOR

main(args.driver_input,nor)
