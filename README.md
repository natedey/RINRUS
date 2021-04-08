# RINRUS
Residue Interaction Network-based ResidUe Selector (RINRUS) is a QM-cluster model building tool which starts from a raw PDB file. After seriers procedures, important residues for certain chemical reactions will be selected and trimmed PDB files as well as quantum mechanics (QM) input files will be generated.

## Installation

Clone this repository, then add the library code under `lib3` to your `PYTHONPATH`. For example, in `~/git`:
``` bash
cd ~/git
git clone git@github.com:MiloCheng17/RINRUS.git
export PYTHONPATH="~/git/RINRUS/lib3:$PYTHONPATH"
```

### Python dependencies

- Python >= 3.x
- NumPy

For certain scripts (optional),
- matplotlib
- BioPython

### External dependencies

- [probe](https://github.com/rlabduke/probe)
- [reduce](https://github.com/rlabduke/reduce)

which both require
- CMake >= 3.10
- Any C/C++ compiler suite with C++11 support

## Usage example

1. After getting a raw PDB file (`3bwm.pdb`), check for ambiguous atoms and do general clean up.

2. Run `reduce` to generate a new H-added PDB file (`3bwm_h.pdb`):
```bash
$HOME/git/RINRUS/bin/reduce -NOFLIP 3bwm.pdb > 3bwm_h.ent
```

3. Check the new PDB file. If there is metal, replace with an atom with same coordination (such as replace Mg with O), and save as a new PDB file.

4. Check all ligands, make sure H atoms were added correctly (may need to delete or add more H based on certain condition)

5. If there is atom "CA" or "CB" in ligands, replace with "CA'" and "CB'"

6. After previous 3 check stpes, save the file to a new file (`3bwm_h_modify.pdb`)

7. Use this new PDB file (`3bwm_h_modify.pdb`) to run `probe` and save the result to `*.probe`:
``` bash
$HOME/git/RINRUS/bin/probe -unformated -MC -self "all" 3bwm_h_modify.pdb > 3bwm_h_modify.probe
```

8. Run `probe2rins.py`. The seed is a comma-separated list of colon-separated pairs, the first part being the ID of the PDB subunit, the second part being the residue number in that subunit:
``` bash
python3 $HOME/git/RINRUS/bin/probe2rins.py -f 3bwm_h.probe -s 'A:300,A:301,A:302'
```
This produces `freq_per_res.dat`, `rin_list.dat`, `res_atoms.dat`, and `*.sif`.

9. Run `freq_per_res.dat` file for trim residues:
``` bash
python3 $HOME/git/RINRUS/bin/probe_freq_2pdb.py 3bwm_h_mg.ent 3bwm_h.probe freq_per_res.dat 'A:300,A:301,A:302'
```
This creates a `model_detail.dat` file that associates model sizes (residue count) with the residue numbers in that model, plus a `res_NNN.pdb` for each model.

10. Run `pymol_scripts.py` to add hydrogens to `res_*.pdb`, which generates a `log.pml` file:
```bash
python3 $HOME/git/RINRUS/bin/pymol_scripts.py res_*.pdb 300,301,302
# no need to run "pymol -qc log.pml" now
```

11. Run `write_input.py` to generate template file and input file
Example:
```bash
python3 bin/write_input.py -noh res_*.pdb -adh res_*_h.pdb -intmp input_templat
```

### Distance based rules

1. Run `python3 ~/git/RINRUS/bin/pdb_2dist_freq.py -pdb 3bwm_h_mg.ent -s A:300,A:301,A:302 -cut 5`
this will generate a file named `dist_per_res-5.00.dat`

2. Then Run `python3 ~/git/RINRUS/bin/pdb_dist_2pdb.py -pdb 3bwm_h_mg.ent -s A:300,A:301,A:302 -cut 5`
This will generate a set of trimmed pdb files such as `dres_3.pdb, dres_4.pdb...`

3. Then run `python3 ~/git/RINRUS/bin/pymol_scripts.py dres_3.pdb dres_3_h.pdb`
This will generate the H added pdb files.
