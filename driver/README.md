The driver runs the entire RINRUS model building procedure from initial PDB structure to input file with one command. To run individual steps or for more information on how the procedure works, see the [step-by-step usage instructions](../bin/README.md).

```bash
# Usage of the RINRUS driver
$HOME/git/RINRUS/driver/RINRUS_driver.py -i driver_input -red True

# Arguments:
-i FILE    driver input file (default: driver_input)
-red BOOL  protonate initial structure with reduce true/false (default: True)
``` 

## Input file
The `driver_input` file should contain the following details:
```
PDB: [starting pdb file name]
Seed: [seed residues]
RIN_program: [probe/arpeggio/distance/manual]
Model(s): [model number or maximal or all] 
Substrate(s)_charge: [overall charge of seed]
Multiplicity: [model multiplicity]
Computational_program: [gaussian/gau-xtb/orca/qchem]
input_template_path: [path to input file template]
basisset_library: [path to Gaussian basis set file]
path_to_type_of_RIN: [path to RINRUS bin directory]
```

## Log file
The driver writes a log file called `rinrus_log.out` with the details of the run/commands called.

## Commands run by the driver

### Reduce protonation if selected

```bash
$HOME/git/RINRUS/bin/reduce -NOFLIP -Quiet PDB.pdb > PDB_h.pdb 
```

### RIN generation/atom selection

```bash
# If probe selected:
$HOME/git/RINRUS/bin/probe -unformated -MC -self "all" -Quiet PDB_h.pdb > PDB.probe
$HOME/git/RINRUS/bin/probe2rins.py -f PDB.probe -s [seed residues]

# If arpeggio selected
$HOME/git/RINRUS/bin/arpeggio/arpeggio.py PDB_h.pdb
$HOME/git/RINRUS/bin/arpeggio2rins.py -f PDB.contacts -s [seed residues]

# If distance selected
$HOME/git/RINRUS/bin/pdb_dist_rank.py -pdb PDB_h.pdb -s [seed residues] -cut [cutoff] -type [avg/mass]
```

If using the distance selection metric, the subsequent trimming/capping/input file generation steps need to be run manually. 

### Model trimming and capping
```bash
$HOME/git/RINRUS/bin/rinrus_trim2_pdb.py -s [seed residues] -pdb PDB_h.pdb -model N [-c contact_counts.dat (if arpeggio)]
$HOME/git/RINRUS/bin/pymol_protonate.py -ignore_ids [seed residues] -pdb res_N.pdb
```

### Input file creation
```bash
$HOME/git/RINRUS/bin/write_input.py -intmp input_template_path -format computational_program -basisinfo basisset_library -c substrate_charge -noh res_N.pdb -adh res_N_h.pdb
```
