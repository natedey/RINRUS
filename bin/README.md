# How to make QM-cluster models and input files with RINRUS

The RINRUS workflow can be run all in one go with the driver or done step by step with the individual scripts. Both ways of using RINRUS are described in this document.

> [!IMPORTANT]  
> As of July 2025, all selection schemes (probe/arpeggio/distance) have been updated to group the selected atoms/contact counts from standard amino acid residues by protein functional groups (side chains and peptide bonds) instead of by residue IDs. This matches the partitioning used for F-SAPT analysis/ranking and the structural building blocks used in the trimming procedure. The different parts of each residue are now listed and ranked independently in `res_atoms.dat`, so they can be added at different stages when doing incremental model building. The contents of the final/maximal model will remain exactly the same but the `res_N`/`model_N` size label will usually be higher as this number reflects the number of separate groups in `res_atoms.dat`. [See examples for more information](../examples/2025_NEW_EXAMPLES). 


# Using the driver
The driver runs the entire RINRUS model building procedure from initial PDB structure to input file with one command (and a few interactive inputs along the way). The structure pre-processing described below in part 1 of the step-by-step usage instructions needs to be done before using the driver (except for protonation with reduce).

```bash
# Usage of the RINRUS driver
python $HOME/git/RINRUS/bin/RINRUS_driver.py -i rinrus.inp

# Arguments:
-i FILE    rinrus driver input file (default: rinrus.inp)
``` 

### Input file
The driver recognises the following (case independent) options in `rinrus.inp`.
```
# required
pdb:                       starting PDB filename
seed:                      ch:ID[,ch:ID,...]
rin_program:               probe OR arpeggio OR distance OR manual
model:                     all OR maximal OR max OR number

# optional
path_to_scripts:           path to RINRUS bin directory				defaults to same folder as RINRUS_driver.py
protonate_initial:         true OR false 					defaults to false
res_atoms_file:            filename 						only used with "rin_program: manual"
arpeggio_rank:             contacts OR types 					only used with "rin_program: arpeggio"
dist_type:                 avg OR com OR closest				only used with "rin_program: distance"
dist_satom:                ch:ID:atom[,ch:ID:atom,...]				only used with "rin_program: distance"
dist_max:                  number						only used with "rin_program: distance"
dist_noh:                  true OR false					only used with "rin_program: distance"
must_add:                  ch:ID[:S/:N/:C] etc
unfrozen:		   ch:ID[,ch:ID:CA,cd:ID:CB,...]
nc_res_info:		   filename
model_prot_ignore_ids:     ch:ID[,ch:ID,...]
model_prot_ignore_atoms:   ch:ID:atom[,ch:ID:atom,...]
model_prot_ignore_atnames: atom[,atom,...]				
qm_input_format:           gaussian OR orca OR qchem OR gau-xtb OR psi4-fsapt	if no value, no QM inputs and remaining options ignored
qm_input_template:         filename						defaults to those in RINRUS/template_files/
gaussian_basis_intmp:      [true/false]						defaults to false
qm_calc_hopt:              [true/false]						defaults to false
seed_charge:               integer						defaults to 0
multiplicity:              integer						defaults to 1
fsapt_fa:                  ch:ID[,ch:ID,...]
```
A [template rinrus.inp file](../template_files/rinrus.inp) listing all options is provided in the template_files directory. Unneeded options can be deleted, commented out or just given no value. 


### Log file
The driver writes a log file called `rinrus_log_[date].out` with the details of the run/commands called.

<details>
    <summary>Examples of commands run by the driver:</summary>
    
```bash
# If reduce protonation selected
$HOME/git/RINRUS/bin/reduce -NOFLIP -Quiet PDB.pdb > PDB_h.pdb 

# If probe RIN selected:
$HOME/git/RINRUS/bin/probe -unformated -MC -self "all" -Quiet PDB.pdb > PDB.probe
$HOME/git/RINRUS/bin/probe2rins.py -f PDB.probe -s seed

# If arpeggio RIN selected
$HOME/git/RINRUS/bin/arpeggio/arpeggio.py PDB.pdb
$HOME/git/RINRUS/bin/arpeggio2rins.py -f PDB.contacts -s seed

# If distance RIN selected
$HOME/git/RINRUS/bin/dist_rank.py -pdb PDB.pdb -s seed -type dist_type [-max dist_max] [-satom dist_satom] [-noH]

# Model trimming and capping
$HOME/git/RINRUS/bin/rinrus_trim2_pdb.py -s seed -pdb PDB.pdb -model model -ra (res_atoms file for selected RIN program) [-mustadd most_add]
$HOME/git/RINRUS/bin/pymol_protonate.py -pdb res_N.pdb [-ignore_ids model_prot_ignore_ids] [-ignore_atoms model_prot_ignore_atoms] [-ignore_atnames model_prot_ignore_atnames]
$HOME/git/RINRUS/bin/make_template_pdb.py -model N

# Input file creation (if QM_input_format defined)
$HOME/git/RINRUS/bin/write_input.py -format qm_input_format -c seed_charge -pdb model_N_template.pdb [-intmp qm_input_template] [-basisinfo intmp] [-type hopt] [-fA fsapt_fa]
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
-s SEED   seed fragment(s) (e.g. "A:300,A:301,A:302")
```

This produces `FG_probe_counts.dat` and `res_atoms.dat`.

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
python3 ~/git/RINRUS/bin/arpeggio2rins.py -f 3bwm_h.contacts -s A:300,A:301,A:302 -pdb 3bwm_h.pdb

# All arguments for arpeggio2rins:
-f FILE   arpeggio contacts file
-s SEED   seed fragment(s) (e.g. "A:300,A:301,A:302")
-pdb PDB  pdb file
-prox     include proximal contacts (ignored by default, "-proximal" also works as a flag)
````
This produces the files `FG_arpeggio_counts.dat`, `res_atoms.dat` and `res_atoms_types.dat`.

The usage of `arpeggio2rins.py` is very similar to that of `probe2rins.py`, but the pdb file is additionally required to facilitate the FG-based grouping as the residue names are not printed in the arpeggio output. 
By default, proximal contacts are ignored as these do not represent "meaningful" interactions, just that atoms are within 5 Å of each other. 
These can be included with the optional '-prox' or '-proximal' flag. 
Note that the original version of `arpeggio2rins.py` only changed the counts when ignoring proximal contacts; atoms which only had proximal contacts were still included in `contact_counts.dat` (equivalent to `res_atoms.dat`) so unnecessary functional groups were often added to the models. This has been fixed in the new version (published July 2025).

**Use `res_atoms.dat` (ranked first by number of contacts, then by number of types) or `res_atoms_types.dat` (ranked by types then contacts) as the input for the trimming procedure in part 3.**

</details>

<br>

<details>
  <summary> <b><font size="+1">2C. Using distance ranking</font></b> </summary>

There are two key options that determine how RINRUS calculates the distance between residues and the seed for distance-based selection and ranking. 
* Distance type: distance can be calculated to the seed's centre of mass or average Cartesian coordinates, or the closest seed atom.
* Hydrogens: distances can be calculated using all atoms or only heavy atoms (all hydrogens are ignored). This applies to both the residue-seed distance calculation and calculation of the seed COM/average centre.

Use `dist_rank.py` to select all fragments with any atoms within a cutoff radius of the seed. A limited set of seed atoms can be selected for the distance calculations by using the '-satom' flag instead of '-s'.
```bash
# Example usage of dist_rank selecting residues within 5A of the full seed COM
python3 ~/git/RINRUS/bin/dist_rank.py -pdb 3bwm_h.pdb -s A:300,A:301,A:302 -max 5 -type mass

# All arguments for dist_rank
-s SEED       seed fragment(s) (e.g. A:300,A:301,A:302)
-satom ATOMS  seed atoms (e.g. A:301:C8,A:301:N9,A:302:C1,A:302:N1)
-type TYPE    how to calculate distance from seed ('closest' or 'avg' or 'mass')
-max CUTOFF   cut off distance in Ã… (default: 5)
-noH          ignore hydrogen atoms (true if flag present)
-store        store distances and pdb/seed/noH selections so they can be reanalysed with a new distance type/cutoff value
-recut        use stored distances and pdb/seed/noH selections (overwrites any -s/-satom/-type/-noH flags used)
-byres        group by residue instead of functional group for comparison with old selection schemes
```

This produces a file `all_dists_table.dat` listing all distances of all atoms in the pdb from the given seed. 
The atoms are filtered by distance type and cutoff value and grouped by functional group to give `sorted_FGs_[type]_[max].dat` and a matching `res_atoms.dat` file.

If you want to test different cutoffs/distance types, use the '-store' flag to save the distances/info about how they were calculated to a python dictionary in `all_dists.pkl`. 
Then use the '-recut' flag to load the saved info instead of doing the initial distance calculations again.

**Use `res_atoms.dat` as the input for the trimming procedure in part 3.**

</details>

<br>

<details>
  <summary> <b><font size="+1">2D. Manual selection and ranking</font></b> </summary>

You can generate your own `res_atoms.dat` file using an existing `res_atoms.dat` file as a template or from scratch.
* The first two columns should list the chain and residue ID of a given residue.
* The third column is where the ranking value would go. This isn't actually used by the trimming script so the value doesn't matter as long as something is there.
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
python3 ~/git/RINRUS/bin/rinrus_trim2_pdb.py -pdb 3bwm_h.pdb -s A:300,A:301,A:302 -ra contact_counts.dat -model max

# All arguments for rinrus_trim2_pdb
-pdb FILE       pre-processed PDB file (e.g. 3bwm_h.pdb)
-s SEED         seed fragment(s)) (e.g. A:300,A:301,A:302)
-ra FILE        atom info file (default: res_atoms.dat)
-model N        specify which model to create (number or 'max' or 'all')
-mustadd RES	fragment(s) that must be in model (ch:ID for whole res or ch:ID:<S/N/C/S+N/S+C/N+C/S+N+C>)
-unfrozen RES   residues/atoms to avoid freezing (ch:ID or ch:ID:CA or ch:ID:CB)
-ncres FILE	non-canonical res info (file containing lines of "resname [all atom names]")
```

If no model is specified, the script will generate the entire "ladder" of possible models by adding residues based on their order in `res_atoms.dat` (same as specifying '-model all'), otherwise only the model containing N residues will be created (or the maximal model if 'max' is specified). For each model, the files `res_N.pdb`, `res_N_froz_info.dat` and `res_N_atom_info.dat` are created. 

RINRUS automatically freezes the CA atoms of recognised (standard or specified in the non-canonical res info file) residues, as well as CB of residues with large side chains (Arg, Lys, Glu, Gln, Met, Trp, Tyr and Phe). These can be removed with the unfrozen flag. Specifying just the residue will remove all automatically applied constraints, or CA/CB can be specifically unfrozen. RINRUS will print a warning if any part of the seed has been frozen in case these are part of substrate peptides or free amino acids etc that shouldn't be constrained.

The mustadd flag is used to define groups that must be in the model but are not part of the seed (e.g. metal coordinating residues or groups covalently bonded to the seed). These are included independently of what is in `res_atoms.dat` so the model numbers might double count some groups. 


### Capping the truncation sites with hydrogens

Use `pymol_protonate.py` to add capping hydrogens to each `res_N.pdb` file where bonds were broken in the trimming procedure. You must have a local copy of PyMOL installed! The script generates a `log.pml` input file containing commands that perform the hydrogen addition and then runs it in PyMOL. 
```bash
# Example usage of pymol_protonate
python3 $HOME/git/RINRUS/bin/pymol_protonate.py -pdb res_N.pdb -ignore_ids "A:300,A:301,A:302" 

# All arguments for pymol_scripts
-pdb FILE               pdb of model to be protonated
-ignore_ids CH:ID       fragment(s) to exclude from protonation (e.g. all seed fragments)
-ignore_atoms CH:ID:AT  specific atom(s) to exclude from protonation (e.g. A:25:C+O,A:25:N) 
-ignore_atnames ATOMS   atom type(s) to exclude (e.g. "ND1,NE2" to avoid changing histidine protonation. Be very careful with this!)
```

The "ignore_ids", "ignore_atoms" and "ignore_atnames" flags are used to specify residue IDs, specific atoms and atom types that should not be protonated. You will most likely want to put the seed fragments in the "ignore_ids" list, otherwise PyMOL might reprotonate your noncanonical amino acids/substrate molecules (and make very poor decisions). By default the script avoids re-protonating the end nitrogens in arginine (NH1 and NH2). If you use "ignore_ats" to avoid any other atom types generally (for example the nitrogens in histidine side chains), check your models carefully before proceeding!

### Making template model pdb files

Use `make_template_pdb.py` to reformat the protonated model pdbs. This makes the file `model_N_template.pdb` which encodes the frozen atom info into the last column of each line. Important for correctly applying constraints in the QM input files made in the next step!
```bash
# Usage of make_template_pdb with standard name formats (res_N.pdb and res_N_h.pdb)
python3 $HOME/git/RINRUS/bin/make_template_pdb.py -model N

# Usage of make_template_pdb by specifying both files separately
python3 $HOME/git/RINRUS/bin/make_template_pdb.py -noh res_N.pdb -addh res_N_h.pdb
```


## 4. Generating input files

Use `write_input.py` to generate input files for quantum chemistry packages. Currently input files can be created for Gaussian, xTB (through Gaussian), ORCA and Q-Chem. Some of the options for this script are specific to our QM-cluster modelling workflow and may not be relevant to other users.

**Make sure to use a model pdb file that has the frozen atom info encoded (see template pdb section above) to ensure the atom constraints are applied properly**

Basic usage of `write_input.py` to create an input file to do a basic geometry optimisation + frequency calculation:
```bash
# Example usage of write_input: Gaussian
python3 $HOME/git/RINRUS/bin/write_input.py -pdb model_N_template.pdb -c 2 -format gaussian

# Example usage of write_input: xTB (in Gaussian)
python3 $HOME/git/RINRUS/bin/write_input.py -pdb model_N_template.pdb -c 2 -format gau-xtb

# Example usage of write_input: ORCA
python3 $HOME/git/RINRUS/bin/write_input.py -pdb model_N_template.pdb -c 2 -format orca

# Example usage of write_input: Q-Chem
python3 $HOME/git/RINRUS/bin/write_input.py -pdb model_N_template.pdb -c 2 -format qchem

# Useful arguments for basic usage of write_input:
-pdb FILE           pdb file to use as structure in input file
-m MULT             model multiplicity
-c LIGCHRG          ligand charge
-format PROG        software package (gaussian/gau-xtb/orca/qchem)
-intmp FILE         input template for selected software package (if none specified, the ones in $HOME/git/RINRUS/template_files/ are used)
-inpn NAME          name of input file (default: 1.inp)
-basisinfo 'intmp' (Gaussian only) use basis set info from the input template file instead of the library in $HOME/git/RINRUS/lib3/gaussian_basis_dict.py
```
<br> 
F-SAPT input files can be set up in one of two ways:
```bash
# Specify format as psi4-fsapt
python3 $HOME/git/RINRUS/bin/write_input.py -pdb model_N_template.pdb -c -2 -format psi4-fsapt -fA A:203

# Specify fsapt as type (overwrites any format option to 'psi4-fsapt')
python3 $HOME/git/RINRUS/bin/write_input.py -pdb model_N_template.pdb -c 2 -type fsapt -fA A:203
```


Full set of arguments for using `write_input.py`:
```bash
# General arguments:
-m MULT            model multiplicity
-c LIGCHRG         ligand charge
-format PROG       software package (gaussian/gau-xtb/orca/qchem/psi4-fsapt)
-intmp FILE        input template for selected software package (defaults to the ones in $HOME/git/RINRUS/template_files/ if none specified)
-inpn NAME         name of input file (default: 1.inp or input.dat for psi4-fsapt)
-basisinfo 'intmp' (Gaussian only) use basis set info from the input template file instead of default RINRUS basis set library
-wdir PATH         working directory
-type TYPE         type of structure processing for input file:
                   'pdb': use pdb file normally (default)
                   'hopt': freeze all heavy atoms so only hydrogen atoms optimized
                   'gauout': take structure from Gaussian output
                   'replacecoords': update selected atom coords e.g. to create TS guess
		   'fsapt': write fsapt input with struc separated into fA and fB

# If using -type 'pdb' (default if no type specified), 'hopt' or 'fsapt', these are required:
-pdb FILE          pdb file

# If using -type 'gauout', these are required:
-tmp FILE          template pdb file for creating new pdb from output
-outf FILE         Gaussian output file
-ckp FILE          Gaussian checkpoint file

# If using -type 'replacecoords', these are required:
-pdb1 FILE         starting pdb file
-pdb2 FILE         pdb to take coordinates from
-parts FILE        text file containing species in pdb1 to be replaced with pdb2 coordinates

# If using -type 'fsapt' or -format 'psi4-fsapt', these are required:
-fA SEED           res ID(s) to define as fragment A (e.g. A:203), rest of enzyme will be fragment B 
```

## (Optional) 5. Convert optimized coords back to pdb format

Use `gopt_to_pdb.py` or `xyz_to_pdb.py` to convert coordinates from QM outputs into pdb format. 
```bash
# Example usage of gopt_to_pdb
python3 $HOME/git/RINRUS/bin/gopt_to_pdb.py -gout 1.out -pdb model_N_template.pdb -name model_N_opt

# Example usage of xyz_to_pdb
python3 $HOME/git/RINRUS/bin/xyz_to_pdb.py -xyz orca.xyz -pdb model_N_template.pdb -name model_N_opt

# Arguments for gopt_to_pdb/xyz_to_pdb
-gout OUTFILE   Gaussian output file (gopt_to_pdb only)
-xyz XYZFILE    XYZ file (xyz_to_pdb only)
-frame FRAME    geometry to use: -1 for final (default), integer or "all". NOTE: COUNT STARTS AT 0.
-name SAVENAME  name for new pdb file (defaults: "gopt" or base name of specified XYZ file)
```
If making only final frame, structures are saved as `[name].pdb`. If option for -frame is "all" or an integer, file(s) will be called `[name].f[num].pdb`.

