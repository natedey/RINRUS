## Usage example 1 - generating a single or a few input files with probe interaction count ranking 

If your structure is not appropriately pre-processed then go through steps 1-4 to protonate and clean up substrate molecules. This is the user's responsibility to verify, since RINRUS will generate models from provided structure and cannot yet run sanity checks

1. After getting a raw PDB file (`3bwm.pdb`), check for nonstandard fragments, clean up atoms or residues with alternate positions/conformations. Listed below are some possible tasks that need to be done before starting the workflow:
   a) if there are multiple conformations for a residue then keep only one conformation
   b) If starting from MD simulation, you may want to delete neutralizing counterions if they have migrated into the active site
   c) If crystallographic symmetry is embedded into the pdb, one needs to be very careful and manually "unfold" the multimers from the PDB website. This might be more common with older X-ray crystal structures. RINRUS will eventually be improved to recognize this and automate symmetry unpacking. 
  d) We have seen situations where AMBER truncates atom names of two-letter elements when writing to PDB format (turning FE into F for example). If the PDB file comes from source other than PDB repository or similar database, check all atom names, residue names, and residue ID #s.
  e) If coming from MD, check that metal coordination is correct

2. If the protein is not yet protonated, run `reduce` to generate a new protonated PDB file  (`3bwm_h.pdb`): (or protonate with other program of your choice). Skip this step if the pdb file is from MD simulation since it has been protonated already.
```bash
$HOME/git/RINRUS/bin/reduce -NOFLIP 3bwm.pdb > 3bwm_h.pdb 
``` 

3. copy the protonated PDB to a new filename (`3bwm_h.pdb > 3bwm_h_modify.pdb`)

4. Check the new PDB file. Probe doesn't recognize covalent bonds between fragments, which becomes a problem with metalloenzymes and metal-ligand coordination. Our workaround is to replace a metal coordination center with an atom that probe recognizes as a H-bond donor or acceptor to capture interactions with the coordinating ligand atoms. For example in COMT, we replace Mg with O, and save as a new PDB file.

5. Check all ligands, make sure H atoms were added correctly (may need to delete or add more H based on certain conditions). The user needs to protonate ligand and substrate properly since reduce does not recognize most substrates and may add H improperly. Check the ligand 2D drawing for the PDB entry on rcsb website for some guidance on substrate protonation if starting from scratch. Performing this step correctly is the user's responsibility! RINRUS will generate models from provided structure and does not yet have comprehensive sanity checks.

6. Use the new PDB file (`3bwm_h_modify.pdb`) to run `probe` to generate .probe file	
``` bash
$HOME/git/RINRUS/bin/probe -unformated -MC -self "all" 3bwm_h_modify.pdb > 3bwm_h_modify.probe 
```
###NOTE: When creating QM-cluster models, remember to replace metal atom in PDB if it was replaced with O in step 4

***Step 7 is one of the most important steps of the QM-cluster model building process. The user must now define the "seed". What is the seed? Typically, the seed will be the substrate(s) (or ligand in biochemical terms) participating in the chemical reaction. Any amino acid residues, co-factors, or fragments which participate in the active site catalytic breaking and forming of chemical bonds may also need to be included as part of the seed, but this will generate much larger models compare to only using the substrate***

##Using the example of PDB:3BWM, we will select the PDB residue ID# of three fragments: 300(metal Mg2+), 301 (SAM) and 302 (Catechol - the substrate) as the seed. 3BWM only has one chain, so we must specify chain A throughout. If the PDB does not have chain identifiers, you will need to specify ":XXX" where XXX is residue id number to use them in this step and beyond. If the protein is multimeric, use the chain of your choice for seed fragments. Note that some multimeric x-ray crystal structures may not necessarily have equivalent active sites!
##NOTE: If there is no Chain ID in the PDB file, our defaults are wonky and need to be improved. 

7. Run `probe2rins.py`. The seed is a comma-separated list of colon-separated pairs, the first part being the ID of the PDB subunit, the second part being the residue number in that subunit:(Chain:ResID)
``` bash
python3 $HOME/git/RINRUS/bin/probe2rins.py -f 3bwm_h_modify.probe -s A:300,A:301,A:302
```
```bash
  -f PROBEFILE  probe_file
  -s SEED       seed for select RIN, in the format of A:300,A:301,A:302
```
This produces `freq_per_res.dat`, `rin_list.dat`, `res_atoms.dat`, and `*.sif`.

8. With the res_atoms.dat file generated, use run rinrus_trim2_pdb to generate the trimmed PDB model:
```bash
python3 ~/git/RINRUS/bin/rinrus_trim2_pdb.py -s A:300,A:301,A:302 -pdb 3bwm_h_modify.pdb 
```

This generates automatically generate the entire "ladder" of possible models based on a ranking scheme which contain `res_NNN.pdb`, `res_NNN_froz_info.dat` and `res_NNN_atom_info.dat` for the all models, where `NNN` is the number of residues in that model.

#Note: if you want to generate one model based on a any ranking scheme, you will need to use the flag `-model NNN`:
```bash
python3 ~/git/RINRUS/bin/rinrus_trim2_pdb.py -s A:300,A:301,A:302 -pdb 3bwm_h_modify.pdb -model NNN 
```
```bash
  -pdb R_PDB       protonated pdbfile
  -s SEED          Chain:Resid,Chain:Resid
  -c R_ATOM        atom info for each residue (eg. res_atom.dat)
  -unfrozen        CA/CB/CACB to unfreeze backbone of seed canonical residues
  -model METHOD    generate one or all trimmed models, if "7" is given, then
                   will generate the 7th model
```
#Note: if you include canonical residues in your seed, you may want to unfreeze CA/CB or both backbone atoms in the seed canonical residues. To do so, you will use the flag `-unfrozen Chain:residue_id:<carbon atoms>` where <carbon atoms> is CA for the alpha carbon only, CB for the beta carbon only, and CACB for both. 
   
9. The trimming procedure creates uncapped backbone pieces. Next use pymol_script.py to add capping hydrogens where bonds were broken when the model was trimmed. You must have a local copy of pymol installed! Run `pymol_scripts.py` to add hydrogens to one or more `res_NNN.pdb` files:

which
- generates a `log.pml` PyMOL input file containing commands that perform the hydrogen addition, and then
- runs PyMOL to perform the addition.
If `-ignore_ids` is specified, those residue IDs will not have hydrogens added. NOTE: This is an important part of the process and you will most likely want to put the seed residues in this list. If you don't, pymol might (probably will) reprotonate your noncanonical amino acids/substrate molecules and make very poor decisions. This step will also ignore re-protonating histidine side chain protons. Please check histidine protonation state of your pdb file before proceeding.
```bash
python3 $HOME/git/RINRUS/bin/pymol_scripts.py -ignore_ids "300,301,302" -pdbfilename res_.pdb
```
```bash
-ignore_ids "300,301,302" (if you do not want pymol to re-protonate fragments 300, 301, and 302)
-ignore_ids "300,301,302 and not name NE2" (to ignore fragments 300, 301, 302, and all NE2 atoms in the model. See [pymol boolean algebra](https://www.pymolwiki.org/index.php/Selection_Algebra) page for formatting details. 
```

Note: To properly loop pymol_script.py over all possible models, you will have to write a bash or python script to loop. An example of a bash loop is given here:
```bash
ls -lrt| grep -v slurm |awk '{print $9}'|grep -E _atom_info |cut -c 5-6 |cut -d_ -f1>list; mkdir pdbs; for i in `cat list`; do mkdir model-${i}-01; cd model-${i}-01; mv ../res_${i}.pdb .;mv ../res_${i}_atom_info.dat  .;mv ../res_${i}_froz_info.dat .; python3 ~/git/RINRUS/bin/pymol_scripts.py -ignore_ids "300,301,302" -pdbfilename *.pdb; cp *_h.pdb model-${i}_h.pdb; cp model-${i}_h.pdb ../pdbs/ ; cd ..; done
```

10. Run `write_input.py` for a single model to generate an input file for quantum chemistry packages. You will have to loop this with python or a shell script to iterate over all possible models. Simple DFT/xTB templates are included:
```bash
python3 ~/git/RINRUS/bin/write_input.py -intmp ~/git/RINRUS/bin/gaussian_input_template.txt -format gaussian -basisinfo ~/git/RINRUS/template_files/basisinfo -c -2 -noh res_NN.pdb -adh res_NN_h.pdb 
```
```bash
-noh NO_H_PDB     trimmed_pdb_file
-adh H_ADD_PDB    hadded_pdb_file
-intmp INPUT_TMP  template_for_write_input
-m MULTIPLICITY   multiplicity
-c LIGAND_CHARGE  charge_of_ligand
-format FMAT      input_file_format eg.'gaussian','qchem','gau-xtb'
-basisinfo        basisset library file(ONLY FOR GAUSSIAN- if prvided then it will automatically pull from the library based on atom present in model,other wise it will be from intmp)
"run python3 ~/git/RINRUS/bin/write_input.py --help" to see more information about flags
```

## Usage example 2 - generating a single or a few input files with distance-based ranking: 

1. Follow step 1-5 from example 1 to get protonated pdb. In this case we will work with `3bwm.pdb` and after step 2 of example 1 it will generate `3bwm_h.pdb`. If you have modified anything from `3bwm_h.pdb` then use that for following process.

2. To create a ranking list of residues/fragments/atoms within 5 Angstroms from either the center of mass of the seed fragments/atoms OR the average of Cartesian coordinates of the seed fragments(in this example, A:300,A:301,A:302) or atoms (in this example A:301:C8,A:301:N9,A:302:C1,A:302:N1), run pdb_dist_rank.py. This script can be run with or without hydrogens included in the calculation of distances by adding the "-nohydro" flag to neglect hydrogen.
```bash
python3 ~/git/RINRUS/bin/pdb_dist_rank.py -pdb 3bwm_h.pdb -s A:300,A:301,A:302 -cut 5 -type avg -nohydro
python3 ~/git/RINRUS/bin/pdb_dist_rank.py -pdb 3bwm_h.pdb -s A:300,A:301,A:302 -cut 5 -type avg
python3 ~/git/RINRUS/bin/pdb_dist_rank.py -pdb 3bwm_h.pdb -s A:300,A:301,A:302 -cut 5 -type mass -nohydro
python3 ~/git/RINRUS/bin/pdb_dist_rank.py -pdb 3bwm_h.pdb -s A:300,A:301,A:302 -cut 5 -type mass
python3 ~/git/RINRUS/bin/pdb_dist_rank.py -pdb 3bwm_h.pdb -satom A:301:C8,A:301:N9,A:302:C1,A:302:N1 -cut 5 -type mass -nohydro
python3 ~/git/RINRUS/bin/pdb_dist_rank.py -pdb 3bwm_h.pdb -satom A:301:C8,A:301:N9,A:302:C1,A:302:N1 -cut 5 -type mass
```
```bash
  -type       "mass" or "avg"
  -nohydro     If flag present, ignore hydrogen atoms from distance calculations               
  -cut COFF    cut_off_dist, default = 3 Å
  -s SEED      center_residues, examples: A:300,A:301,A:302
  -satom       center_atoms, examples: A:301:C8,A:301:N9,A:302:C1,A:302:N1
```
This script will generate a file named `res_atom-5.00.dat`, which has same format as `res_atoms.dat`, contains information about all residue atoms within 5 Angstroms in increasing order of distance from the seed, and contains information about important atoms to be included in the identified residues. This script will also generate a `res_atom-5.00_new.dat` that contains the information about the individual side chain and main chain for each residue and a `res_atom-5.00_new.csv` that can be loaded to excel. Lastly, a `data_analysis_output_-5.00.csv` file is created to show the closest distance of the side chain and main chain atom type for each residue from the seed.  

3. Once the `res_atom-5.00.dat` file is generated, generate the trimmed PDB model using RINRUS by following steps 9 and 11 of example 1 to protonate the pdbs and generate input files.

## Usage example 3 - generating a single or a few input files with arpeggio ranking by either interaction counts or number of interaction types: 

#Note: Before starting, ensure that the pdb is cleaned. Use example 1, steps 1 to 7 above to clean your pdb
1. Make sure openbabel libraries are available to properly use RINRUS with arpeggio

2. Run the arpeggio.py script to generate the contact file (you also need to make sure config.py is in the same directory as arpeggio.py)
````bash
python3 ~/git/RINRUS/bin/arpeggio/arpeggio.py 2cht_h-TS.pdb
````  

3. Use the `arpeggio2rins.py` script to generate "contact_counts.dat",  "contype_counts.dat" and, "node_info.dat":
````bash
python3 ~/git/RINRUS/bin/arpeggio2rins.py -f 2cht_h-TS.contacts -s C:202
````
Then you can choose to use either rank by number of arpeggio contact counts between fragments and seed ('-c contact_counts') or by number of arpeggio interaction types ('-c contype_counts') to generate models. Note that `contact_counts.dat` and `contype_counts.dat` have the same format as `freq_per_res.dat`, while `node_info.dat` has the same format as `res_atoms.dat`. These files can be used interchangeably. 

4. Run rinrus_trim2_pdb.py to generate a specific model, for example the model with the seed and seven most important fragments: 
 ```bash
python3 ~/git/RINRUS/bin/rinrus_trim2_pdb.py -pdb 2cht_h-TS.pdb -c contact_counts.dat -s C:202 -model 7
```
#Note: if you want to unfreeze CA or CB or both atoms in noncanonical resides in your models, you will need to use the flag `-free Chain:residue_id:carbon atoms to unfreeze`:
```bash
python3 ~/git/RINRUS_master/bin/rinrus_trim2_pdb.py -s A:203,A:7,C:63 -pdb 2cht_h_ac_aligned.pdb -c res_atoms.dat -free A:7:CACB,C:63:CA -model 26
```
For -free A:7:CACB,C:63:CA in the line above, both the CA and CB carbons are free in Arg7 and in Arg63 only CB atom is free in the model.
5. Run the pymol_script.py to add hydrogens to cap trimmed residues:
 ```bash
python3 ~/git/RINRUS/bin/pymol_scripts.py -ignore_ids 202 -pdbfilename res_7.pdb 
```
#Note: You can write a bash script to loop over all the models. See example below
```bash
ls -lrt| grep -v slurm |awk '{print $9}'|grep -E _atom_info |cut -c 5-6 |cut -d_ -f1>list; mkdir pdbs; for i in `cat list`; do mkdir model-${i}-01; cd model-${i}-01; mv ../res_${i}.pdb .;mv ../res_${i}_atom_info.dat  .;mv ../res_${i}_froz_info.dat .; python3 ~/git/RINRUS/bin/pymol_scripts.py -ignore_ids 202 -pdbfilename *.pdb; cp *_h.pdb model-${i}_h.pdb; cp model-${i}_h.pdb ../pdbs/ ; cd ..; done
```
7. Run `write_input.py` for a single model to generate a template file and input file. You will have to loop this with python or a shell script to iterate over all possible models.
```bash
python3 ~/git/RINRUS/bin/write_input.py -intmp input_temp -c -2 -noh res_7.pdb -adh res_7_h.pdb -format gaussian
```

## Usage example 4 - generating a single or a few input files with manual ranking (from SAPT, ML, or from some scheme that doesn't yet interface with RINRUS automatically)

1. Obtain the a protonated PDB file to use for generating an res_atoms.dat file by hand
 
2. Generate your own res_atoms.dat file from a list of residues, adding side chain and main chain atoms (if any) you want to include in your model.
   Use an existing res_atoms.dat file as a template. More documentation is forthcoming for this part.

3. If you are interested in truncated models (smaller than the "maximal model"), make sure the residues are listed in order that you want to add them,
   with decreasing interaction count in the third column.

## Usage example 5 - GENERATE ALL THE THINGS!!! Combinatorial model building from probe and arpeggio

## Usage example 5a - Combinatorial model building from arpeggio - Dr. D still needs to proofread this part!
1. Refer to usage example 3 steps 1 to 2 to generate arpeggio contact files to compute the combinations

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
## Usage example 5b - Combinatorial model building from probe
1. Refer to usage example 1 steps 1 to 7 to generate probe files to compute combinations

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
