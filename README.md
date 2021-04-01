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

After get a raw pdb file (`raw.pdb`), check ambiguous atoms, clean up
Run `reduce` to generate a new H-added PDB file (`raw_h.pdb`):
```bash
$HOME/git/RINRUS/bin/reduce -NOFLIP 3bwm.pdb > 3bwm_h.ent
```

Check the new pdb file, if there is metal, replace with an atom with same coordination (such as replace Mg with O, and save as a new pdb file)
Check all ligands, make sure H atoms were added correctly (may need to delete or add more H based on certain condition)
If there is atom "CA" or "CB" in ligands, replace with "CA'" and "CB'"
After previous 3 check stpes, save the file to a new file (`raw_h_modify.pdb`)
Use this new pdb file (`raw_h_modify.pdb`), run probe and save result to `*.probe`
Run probe to generate correct probe file
Example:
``` bash
$HOME/git/RINRUS/bin/probe -unformated -MC -self "all" 3bwm_h.ent > 3bwm_h.probe
```

Run probe2rins.py to generate `freq_per_res.dat`. The seed is comma-separated list of colon-separated pairs, the first part being the ID of the PDB subunit, the second part being the residue number in that subunit:
``` bash
python3 $HOME/git/RINRUS/bin/probe2rins.py -f 3bwm_h.probe -s 'A:300,A:301,A:302'
```

Run `freq_per_res.dat` file for trim residues
Example:
``` bash
python3 $HOME/git/RINRUS/bin/probe_freq_2pdb.py 3bwm_h_mg.ent 3bwm_h.probe freq_per_res.dat 'A:300,A:301,A:302'
```
Now it will also print out a `model_detail.dat` file, tells which residues in various size model
If you like to keep the information printed on the screen, you can select it and copy past to a txt file

```
Seeds are: [('A', 300), ('A', 301), ('A', 302)]
3 {'A': [300, 301, 302]}
4 {'A': [300, 301, 302, 40]}
5 {'A': [300, 301, 302, 40, 141]}
6 {'A': [300, 301, 302, 40, 141, 143]}
7 {'A': [300, 301, 302, 40, 141, 143, 91]}
8 {'A': [300, 301, 302, 40, 141, 143, 91, 170]}
9 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68]}
10 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411]}
11 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72]}
12 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66]}
13 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142]}
14 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90]}
15 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199]}
16 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441]}
17 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119]}
18 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144]}
19 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402]}
20 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42]}
21 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458]}
22 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174]}
23 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67]}
24 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169]}
25 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41]}
26 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41, 89]}
27 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41, 89, 118]}
28 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41, 89, 118, 71]}
29 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41, 89, 118, 71, 198]}
30 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41, 89, 118, 71, 198, 120]}
31 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41, 89, 118, 71, 198, 120, 38]}
32 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41, 89, 118, 71, 198, 120, 38, 139]}
33 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41, 89, 118, 71, 198, 120, 38, 139, 146]}
34 {'A': [300, 301, 302, 40, 141, 143, 91, 170, 68, 411, 72, 66, 142, 90, 199, 441, 119, 144, 402, 42, 458, 174, 67, 169, 41, 89, 118, 71, 198, 120, 38, 139, 146, 117]}
```

Use `raw_h.pdb`, `freq_per_res.dat`, and other needed files run `probe_freq_2pdb.py`
Run `pymol_scripts.py` to add H to the trimmed pdb files for differnt size model

Run `pymol_scripts.py` to add H to `res_*.pdb`, which generate a `log.pml` file
Example:
```bash
python3 bin/pymol_scripts.py res_*.pdb 300,301,302
### not need to run "pymol -qc log.pml" now
```

Run `write_input.py` to generate template file and input file
Example:
```bash
python3 bin/write_input.py -noh res_*.pdb -adh res_*_h.pdb -intmp input_templat
```

### Distance based rules

Run `python3 ~/git/RINRUS/bin/pdb_2dist_freq.py -pdb 3bwm_h_mg.ent -s A:300,A:301,A:302 -cut 5`
this will generate a file named `dist_per_res-5.00.dat`
Then Run `python3 ~/git/RINRUS/bin/pdb_dist_2pdb.py -pdb 3bwm_h_mg.ent -s A:300,A:301,A:302 -cut 5`
This will generate a set of trimmed pdb files such as `dres_3.pdb, dres_4.pdb...`
Then run `python3 ~/git/RINRUS/bin/pymol_scripts.py dres_3.pdb dres_3_h.pdb`
This will generate the H added pdb files.
