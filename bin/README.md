# How to make QM-cluster models and input files with RINRUS

The RINRUS workflow can be run all in one go with the driver or done step by step with the individual scripts. Both ways of using RINRUS are described in this document.

# Using the driver
The driver runs the entire RINRUS model building procedure from initial PDB structure to input file with one command (and a few interactive inputs along the way). The structure pre-processing described below in part 1 of the step-by-step usage instructions needs to be done before using the driver (except for protonation with reduce).

```bash
# Usage of the RINRUS driver
python $HOME/git/RINRUS/bin/RINRUS_driver.py -i driver_input

# Arguments:
-i FILE    rinrus driver input file (default: driver_input)
``` 

### Input file
The `driver_input` file should contain the following details:
```
PDB: [starting pdb file name]
Protonate_initial: [true/t/y or false/f/n]
Seed: [seed residues]
RIN_program: [probe/arpeggio/distance/manual]
Model(s): [model number or maximal or all] 
Seed_charge: [overall charge of seed]
Multiplicity: [model multiplicity]
Computational_program: [gaussian/gau-xtb/orca/qchem]
input_template_path: [path to input file template]
basisset_library: [path to Gaussian basis set file]
path_to_type_of_RIN: [path to RINRUS bin directory]
```

### Log file
The driver writes a log file called `rinrus_log_[date].out` with the details of the run/commands called.

<details>
<summary>Commands run by the driver:</summary>
```bash
# If reduce protonation selected
$HOME/git/RINRUS/bin/reduce -NOFLIP -Quiet PDB.pdb > PDB_h.pdb 

# If probe RIN selected:
$HOME/git/RINRUS/bin/probe -unformated -MC -self "all" -Quiet PDB_h.pdb > PDB.probe
$HOME/git/RINRUS/bin/probe2rins.py -f PDB.probe -s [seed residues]

# If arpeggio RIN selected
$HOME/git/RINRUS/bin/arpeggio/arpeggio.py PDB_h.pdb
$HOME/git/RINRUS/bin/arpeggio2rins.py -f PDB.contacts -s [seed residues]

# If distance RIN selected
$HOME/git/RINRUS/bin/pdb_dist_rank.py -pdb PDB_h.pdb -s [seed residues] -cut [cutoff] -type [avg/mass]

# Model trimming and capping
$HOME/git/RINRUS/bin/rinrus_trim2_pdb.py -s [seed residues] -pdb PDB_h.pdb -model N -c [res_atoms.dat/contact_counts.dat/res_atom-X.dat depending on RIN program]
$HOME/git/RINRUS/bin/pymol_protonate.py -ignore_ids [excluded ids] -pdb res_N.pdb

# Input file creation
$HOME/git/RINRUS/bin/write_input.py -intmp input_template_path -format computational_program -basisinfo basisset_library -c seed_charge -type hopt -noh res_N.pdb -adh res_N_h.pdb
```
</details>

# Step by step usage of RINRUS

Examples throughout this section are based on the enzyme COMT (PDB:3BWM), with the active site being built around the Mg2+ ion, SAM cofactor and catechol substrate. 

## 1. Preparing structure files
If your structure is not appropriately pre-processed then go through the following steps to protonate and clean up substrate molecules. This is the user's responsibility to verify, since RINRUS will generate models from provided structure and cannot yet run sanity checks.

1. After getting a raw PDB file (`3bwm.pdb`), check for nonstandard fragments and clean up atoms or residues with alternate positions/conformations. Listed below are some possible tasks that need to be done before starting the workflow:

    a) If there are multiple conformations for a residue then keep only one conformation.

    b) If starting from an MD simulation, check that metal coordination is correct. You may also want to delete neutralizing counterions if they have migrated into the active site.

    c) If crystallographic symmetry is embedded into the pdb, one needs to be very careful and manually "unfold" the multimers from the PDB website. This might be more common with older X-ray crystal structures. RINRUS will eventually be improved to recognize this and automate symmetry unpacking. 

    d) We have seen situations where AMBER truncates atom names of two-letter elements when writing to PDB format (turning FE into F for example). If the PDB file comes from a source other than PDB repository or similar database, check all atom names, residue names, and residue IDs.

2. If the protein is not yet protonated, run `reduce` to generate a new protonated PDB file  (`3bwm_h.pdb`) or protonate with other program of your choice. Skip this step if the pdb file is from an MD simulation since it has been protonated already.
    ```bash
    # Example usage of reduce
    $HOME/git/RINRUS/bin/reduce -NOFLIP 3bwm.pdb > 3bwm_h.pdb 
    ``` 

3. Check the new PDB file. Probe doesn't recognize covalent bonds between fragments, which becomes a problem with metalloenzymes and metal-ligand coordination. Our workaround is to replace a metal coordination center with an atom that probe recognizes as a H-bond donor or acceptor to capture interactions with the coordinating ligand atoms. For example in COMT, we replace Mg with O, and save as a new PDB file (`3bwm_h_modify.pdb`).

4. Check all ligands, make sure H atoms were added correctly (may need to delete or add more H based on certain conditions). The user needs to protonate ligand and substrate properly since reduce does not recognize most substrates and may add H improperly. Check the ligand 2D drawing for the PDB entry on the RCSB website for some guidance on substrate protonation if starting from scratch. Performing this step correctly is the user's responsibility! RINRUS will generate models from provided structure and does not yet have comprehensive sanity checks.

## 2. Selecting the active site/generating the active site RIN

***Defining the "seed" is one of the most important steps of the QM-cluster model building process.*** What should you select as the seed? Typically, the seed will be the substrate(s) (or ligand in biochemical terms) participating in the chemical reaction. Any amino acid residues, co-factors, or fragments which participate in the active site catalytic breaking and forming of chemical bonds may also need to be included as part of the seed, but this will generate much larger models compare to only using the substrate.

The seed is specified as a comma-separated list of colon-separated Chain:ResidueID pairs. In the example of 3BWM, we select the seed as A:300 (Mg2+), A:301 (SAM) and A:302 (catechol). 

If the PDB you're using does not have chain identifiers, you will need to specify ":XXX" where XXX is the residue ID number in this step and beyond. Our current defaults are wonky in these cases and need to be improved. If the protein is multimeric, use the chain of your choice for seed fragments. Note that some multimeric x-ray crystal structures may not necessarily have equivalent active sites!

The active site RIN can be automatically determined from the seed fragments based on probe contacts, arpeggio contacts or distance. For all selection metrics the key output is a file, usually called `res_atoms.dat` or similar, containing a list of the identified active site residues ranked by the chosen metric and the atoms identified by the selection procedure. This file will be used as the input for the trimming procedure in part 3.

<details>
  <summary> <b><font size="+1">2A. Using Probe contact count ranking</font></b> </summary>

First run `probe` on the (modified) PDB file to generate a *.probe file of all contacts in the enzyme
``` bash
# Example usage of probe:
$HOME/git/RINRUS/bin/probe -unformated -MC -self "all" 3bwm_h_modify.pdb > 3bwm_h_modify.probe 
```

Then generate the active site RIN from the probe contacts with `probe2rins.py`
``` bash
# Example usage of probe2rins:
python3 $HOME/git/RINRUS/bin/probe2rins.py -f 3bwm_h_modify.probe -s A:300,A:301,A:302

# All arguments for probe2rins:
-f FILE   probe contacts file
-s SEED   seed fragment(s) (e.g. A:300,A:301,A:302)
```

This produces `freq_per_res.dat`, `rin_list.dat`, `res_atoms.dat`, and `*.sif`.

**Use `res_atoms.dat` as the input for the trimming procedure in part 3.**

**NOTE:** Remember to replace the metal atom in the PDB before continuing/use the unmodified PDB for the remaining steps if it was replaced with O in the preprocessing.

</details>

<br>

<details>
  <summary> <b><font size="+1">2B. Using Arpeggio contact count or contact type ranking</font></b> </summary>

*Make sure openbabel libraries are available to properly use RINRUS with arpeggio.*

First, run `arpeggio` to generate the contact file (you need to make sure config.py is in the same directory as arpeggio.py)
````bash
# Example usage of arpeggio
python3 ~/git/RINRUS/bin/arpeggio/arpeggio.py 3bwm_h.pdb
````  

Then generate the active site RIN from arpeggio contacts using `arpeggio2rins.py`
````bash
# Example usage of arpeggio2rins
python3 ~/git/RINRUS/bin/arpeggio2rins.py -f 3bwm_h.contacts -s A:300,A:301,A:302

# All arguments for arpeggio2rins:
-f FILE   arpeggio contacts file
-s SEED   seed fragment(s) (A:300,A:301,A:302)
````
This produces the files `contact_counts.dat`, `contype_counts.dat` and, `node_info.dat`. Both `contact_counts.dat` and `contype_counts.dat` have the same format as `res_atoms.dat`.

**Use `contact_counts.dat` (residues ranked by number of contacts) or `contype_counts.dat` (residues ranked by number of interaction types) as the input for the trimming procedure in part 3.**

</details>

<br>

<details>
  <summary> <b><font size="+1">2C. Using distance ranking</font></b> </summary>

There are two key options that determine how RINRUS calculates the distance between residues and the seed for distance-based selection and ranking. 
* Seed centre: centre of mass or the average of the Cartesian coordinates.
* Hydrogens: distances can be calculated using all atoms or only heavy atoms (no hydrogens)

Use `pdb_dist_rank.py` to select all fragments with any atoms within a cutoff radius of the seed centre. The seed centre can be determined from only a few atoms instead of the full fragments by using the "-satom" flag instead of "-s".
```bash
# Example usage of pdb_dist_rank selecting residues within 5A of the full seed COM
python3 ~/git/RINRUS/bin/pdb_dist_rank.py -pdb 3bwm_h.pdb -s A:300,A:301,A:302 -cut 5 -type mass

# All arguments for pdb_dist_rank
-type CENTRE  how to determine centre of seed ("mass" or "avg")
-nohydro      ignore hydrogen atoms (true if flag present)
-cut COFF     cut off distance (default: 3 Å)
-s SEED       seed fragment(s) (e.g. A:300,A:301,A:302)
-satom ATOMS  seed atoms (e.g. A:301:C8,A:301:N9,A:302:C1,A:302:N1)
```

This produces the files `res_atom-[COFF].dat`, `res_atom-[COFF]_new.dat`, `res_atom-[COFF]_new.csv` and `data_analysis_output_[COFF].csv`. The latter 3 files separate out the side chain atoms and main chain atoms of each residue. This is currently not an ideal functional group-based partitioning and an updated workflow is in development.

**Use `res_atom-[COFF].dat` as the input for the trimming procedure in part 3.**

</details>

<br>

<details>
  <summary> <b><font size="+1">2D. Manual selection and ranking</font></b> </summary>

You can generate your own `res_atoms.dat` file using an existing `res_atoms.dat` file as a template or from scratch.
* The first two columns should list the chain and residue ID of a given residue
* The third column should contain an arbitrary ranking value
* The rest of the line should be the selected atom(s). 

The residues should be listed in the order you want them to be added to the model.

Example format:
```
A    300      554      O   
A    40       491      CB   CE   CG   HB2  HE1  HE3  HG3  O    SD  
A    141      478      CB   CG   HB3  O    OD1  OD2 
A    143      378      CD1  CD2  CE2  CE3  CG   CH2  CZ2  CZ3  H    NE1 
A    91       335      CB   CG2  H    HB   HG21 HG22 N   
A    170      304      CG   HD21 ND2  OD1
```

**Use your new `res_atoms.dat` as the input for the trimming procedure in part 3.**

</details>

<br>

<details>
  <summary> <b><font size="+1">2E/2F. Combinatorial model building from probe and arpeggio</font></b> </summary>

These instructions have not been updated/checked recently. Ignore for now...

<details>
  <summary> <b>Combinatorial model building from arpeggio</b> </summary>

1. Generate arpeggio contact files as in section 2B to compute the combinations

````bash
python3 ~/git/RINRUS/bin/arpeggio/arpeggio.py 2cht_h-TS.pdb -s /A/203/
````  
After generating the arpeggio contacts, the contact file is used for computing the combinations

2. Run combifromcontacts.py script with the defined seed (chain/residue numbers) which takes combinations of the different interactions 
```bash
python3 ~/git/RINRUS/combi_script/combifromcontacts.py 2cht_h.contacts A/203 2cht_h.sift 
```
This generates LongCombi.dat, SimpCombi.dat, and ModSimpCombi.dat

3. Run genmodelfiles.py to remove redundant models ModSimpCombi.dat file
```bash
python3 ~/git/RINRUS/combi_script/genmodelfiles.py ModSimpCombi.dat
```

4. Use res_atoms_*.dat files generated to prepare modified list of models. The res_atoms_*.dat files are then translated   into their corresponding PDB files 
```bash 
ls res_atoms_*.dat > list
```

5. Open the list and remove "res_atoms_" and ".dat" to leave only the model numbers within list, run the command below after you 
```bash
for i in `cat list`; do python3 ~/git/RINRUS/bin/rinrus_trim2_pdb.py -pdb 2cht_h.pdb -s A:203 -ratom res_atoms_${i}.dat; mv res_*_atom_info.dat atom_info_${i}.dat; mv res_*_froz_info.dat froz_info_${i}.dat; mv res_*.pdb model_${i}.pdb; done
```

6. To identify which pdbs are identical, create a new list of the various model pdb names and run the identifiles.py script
```bash
ls model_*.pdb > list
```
```bash
python3 ~/git/RINRUS/combi_script/identifiles.py list
```
The generated file (UniqueModels.dat) lists all the unique models and removes redundant models

7. Complete the valences by adding protons to severed bonds and waters via PyMol
 ```bash 
 for i in `cat list`; do python3 ~/git/RINRUS/bin/pymol_scripts.py ${i} 203; ~qcheng1/bin/pymol -qc log.pml; done
 ```

</details>
 
<details>
  <summary> <b>Combinatorial model building from probe</b> </summary>
1. Refer to section 2A to generate probe files to compute combinations

2. Run gen-probe-combi.py script with the defined seed (chain/residue numbers/atom(s)) which takes combinations of the different interactions 
```bash
python3 ~/git/RINRUS/git/combi_script/gen-probe-combi.py -f 2cht_h_ac_aligned.probe -seed A/203/C1,O1,C2,O2,C3,O3,C4,O4,C5,O5,C6,O7,C8,C9, C10,C11,HO5,H01,H02,H03,H04,H05,H06,H07
```
Note: Multiple seed indices can be indicated by space separation as A/202 A/202 

3. Run the next step which combine multiple step but remember to edit the necessary part

```bash
ls -lrt| grep -v slurm |awk '{print $9}'|grep -E res_atoms_|cut -c 11-12|cut -d. -f1>list; mkdir pdbs; for i in `cat list`; do mkdir model-${i}-01; cd model-${i}-01; mv ../res_atoms_${i}.dat .; python3 ~/git/RINRUS/bin/rinrus_trim2_pdb.py -s A:203 -ratom res_atoms_${i}.dat -pdb ../2cht_h_ac_aligned.pdb; python3 ~/git/RINRUS/bin/pymol_scripts.py -ignore_ids 203 -pdbfilename *.pdb; cp *_h.pdb model-${i}_h.pdb; cp model-${i}_h.pdb ../pdbs/; cp res_atoms_${i}.dat ../pdbs/${i}.dat ; cd ..; done
```
Be aware that while the res_atoms_#.dat model sets generated are technically unique to each other, once RINRUS generates the full trimmed QM-models, a lot of the them will no longer be unique and will be identical to others. So there a need to check to determine which QM-models are still unique (and not redundant)

</details>
</details>


## 3. Trimming and capping the models

### Trimming the models

Use `rinrus_trim2_pdb.py` to generate trimmed cluster models based on the atoms/residues listed in `res_atoms.dat` (or an equivalent file). If any atom from a functional group (amino acid side chain, peptide bond, water molecule, etc) is included in `res_atoms.dat`, the entire functional group will be included; the trimming algorithm only truncates the model at alpha carbons. 

```bash
# Example usage of rinrus_trim2_pdb using a probe "res_atoms.dat" file
python3 ~/git/RINRUS/bin/rinrus_trim2_pdb.py -pdb 3bwm_h.pdb -s A:300,A:301,A:302

# Example usage of rinrus_trim2_pdb using an arpeggio "contact_counts.dat" file and making only the maximal model
python3 ~/git/RINRUS/bin/rinrus_trim2_pdb.py -pdb 3bwm_h.pdb -s A:300,A:301,A:302 -c contact_counts.dat -model max

# All arguments for rinrus_trim2_pdb
-pdb FILE       pre-processed PDB file (e.g. 3bwm_h.pdb)
-s SEED         seed fragment(s)) (e.g. A:300,A:301,A:302)
-c FILE         atom info file (default: res_atoms.dat)
-unfrozen ATS   residues to unfreeze (Chain:resID:<CA/CB/CACB>)
-model N        specify which model to create (number or 'max')
```

If no model is specified, the script will generate the entire "ladder" of possible models by adding residues based on their order in `res_atoms.dat`, otherwise only the model containing N residues will be created (or the maximal model if 'max' is specified). For each model, the files `res_N.pdb`, `res_N_froz_info.dat` and `res_N_atom_info.dat` are created. 

If canonical residues are included in the seed, you can unfreeze CA and/or CB in these residues with the unfrozen flag. Specify the carbon atoms to unfreeze as CA, CB or CACB.


### Capping the truncation sites with hydrogens

Use `pymol_protonate.py` to add capping hydrogens to each `res_N.pdb` file where bonds were broken in the trimming procedure. You must have a local copy of PyMOL installed! The script generates a `log.pml` input file containing commands that perform the hydrogen addition and then runs it in PyMOL. 
```bash
# Example usage of pymol_protonate
python3 $HOME/git/RINRUS/bin/pymol_protonate.py -pdb res_N.pdb -ignore_ids "A:300,A:301,A:302" 

# All arguments for pymol_scripts
-pdb FILE           pdb of model to be protonated
-ignore_ids CH:ID   fragment(s) to exclude from protonation (e.g. all seed fragments)
-ignore_ats ATOMS   atom name(s) to exclude (e.g. "ND1,NE2" to avoid changing histidine protonation. Be very careful with this!)
```

The "ignore_ids" and "ignore_ats" flags are used to specify residue IDs and atom types that should not be protonated. You will most likely want to put the seed fragments in the "ignore_ids" list, otherwise PyMOL might reprotonate your noncanonical amino acids/substrate molecules (and make very poor decisions). By default the script avoids re-protonating the end nitrogens in arginine (NH1 and NH2). If you use "ignore_ats" to avoid any other atom types (for example the nitrogens in histidine side chains), check your models carefully before proceeding!

## 4. Generating input files

Use `write_input.py` to generate input files for quantum chemistry packages. Currently input files can be created for Gaussian, xTB (through Gaussian), ORCA and Q-Chem. Some of the options for this script are specific to our QM-cluster modelling workflow and may not be relevant to other users.

Basic usage of `write_input.py` to create an input file to do a geometry optimisation + frequency calculation with a given cluster model `res_N_h.pdb`:
```bash
# Example usage of write_input: Gaussian
python3 $HOME/git/RINRUS/bin/write_input.py -pdb res_N_h.pdb -c 2 -format gaussian -intmp $HOME/git/RINRUS/bin/gaussian_input_template.txt -basisinfo ~/git/RINRUS/template_files/basisinfo

# Example usage of write_input: xTB (in Gaussian)
python3 $HOME/git/RINRUS/bin/write_input.py -pdb res_N_h.pdb -c 2 -format gau-xtb -intmp $HOME/git/RINRUS/bin/xtb_input_template.txt

# Example usage of write_input: ORCA
python3 $HOME/git/RINRUS/bin/write_input.py -pdb res_N_h.pdb -c 2 -format orca -intmp $HOME/git/RINRUS/bin/orca_input_template.txt 

# Example usage of write_input: Q-Chem
python3 $HOME/git/RINRUS/bin/write_input.py -pdb res_N_h.pdb -c 2 -format qchem -intmp $HOME/git/RINRUS/bin/qchem_input_template.txt 

# Useful arguments for basic usage of write_input:
-pdb FILE        pdb file to use as structure in input file
-m MULT          model multiplicity
-c LIGCHRG       ligand charge
-format PROG     software package (gaussian/gau-xtb/orca/qchem)
-intmp FILE      input template for selected software package
-inpn NAME       name of input file (default: 1.inp)
-basisinfo FILE  basis library file (only needed for Gaussian calculations)
```
<br> 

Full set of arguments for using `write_input.py`:
```bash
# General arguments:
-m MULT          model multiplicity
-c LIGCHRG       ligand charge
-format PROG     software package (gaussian/gau-xtb/orca)
-intmp FILE      input template for selected software package
-inpn NAME       name of input file (default: 1.inp)
-basisinfo FILE  basis library file (only needed for Gaussian calculations)
-wdir PATH       working directory
-type TYPE       type of structure processing for input file:
                 'pdb': use pdb file normally (default)
                 'hopt': freeze all heavy atoms so only hydrogen atoms optimized
                 'gauout': take structure from Gaussian output
                 'replacecoords': update selected atom coords e.g. to create TS guess

# If using -type 'pdb' (default if no type specified), these are required:
-pdb FILE        pdb file

# If using -type 'hopt', these are required:
-noh FILE        cluster model pdb file before H added (e.g. res_N.pdb)
-adh FILE        cluster model pdb file after H added (e.g. res_N_h.pdb)
-tmp FILE        template pdb file

# If using -type 'gauout', these are required:
-tmp FILE        template pdb file
-outf FILE       Gaussian output file
-ckp FILE        Gaussian checkpoint file

# If using -type 'replacecoords', these are required:
-pdb1 FILE       starting pdb file
-pdb2 FILE       section of structure with new coordinates in pdb format
-parts FILE      text file containing species in pdb1 to be replaced with pdb2 coordinates
```
