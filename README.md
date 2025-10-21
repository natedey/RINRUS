# RINRUS

Residue Interaction Network-based ResidUe Selector (RINRUS) is a QM-cluster model building tool for biomolecular systems. Starting from a raw PDB file, after running a series of preparation tasks, the tool will
- select important residues for chemical reactions, and
- generate trimmed PDB files with the corresponding quantum chemical inputs.

RINRUS is the first tool available that performs automated and algorithmic trimming and capping of enzyme models. Reproducibility is embedded into the model construction workflow, setting new community standards.

A software review paper is in preparation. For now, the best way to acknowledge RINRUS is to cite:
[DOI: 10.1016/j.bpj.2021.07.029](https://doi.org/10.1016/j.bpj.2021.07.029)
and
[DOI: 10.1039/D3CP06100K](https://doi.org/10.1039/D3CP06100K)

The development of RINRUS has been supported by the National Science Foundation Division of Biological Infrastructure
(CAREER BIO-1846408) and the Department of Energy Basic Energy Sciences (SBIR DE-SC0021568).
If you have questions or would like to contribute to the development of RINRUS, please email 
Prof. Nathan DeYonker, Department of Chemistry, University of Memphis at ndyonker@memphis.edu

## Installation

Clone this repository, then add `bin` to your `PATH` and library code under `lib3` to your `PYTHONPATH`. For example, in `~/git`:
``` bash
cd ~/git
git clone https://github.com/natedey/RINRUS.git
export PATH="~/git/RINRUS/bin:$PATH"
export PYTHONPATH="~/git/RINRUS/lib3:$PYTHONPATH"
```

### Python dependencies

- Python >= 3.x
- NumPy
- pymol
  - If installing via conda, it's under `-c conda-forge pymol-open-source`.
- openbabel (required for Arpeggio)
- BioPython (required for Arpeggio)
- pandas
  
For certain scripts (optional),
- matplotlib


### External dependencies

- [probe](https://github.com/rlabduke/probe) - version 2.16.130520 is packaged with RINRUS
- [reduce](https://github.com/rlabduke/reduce) - version 3.23 is packaged with RINRUS
- [arpeggio](https://github.com/harryjubb/arpeggio) - packaged with RINRUS
- obabel version = openbabel/2.4.1

Currently, precompiled binaries of probe, reduce, and arpeggio are present in `bin/`.

which require
- CMake >= 3.10
- Any C/C++ compiler suite with C++11 support

PDB HET groups (some ligands and noncanonical amino acids) can be protonated by reduce. 
The connectivity table file is included with RINRUS: bin/reduce_wwPDB_het_dict.txt
You must set a shell environment variable to allow reduce to use the reduce_wwPDB_het_dict file, 
as the default location is /local:
```bash
setenv REDUCE_HET_DICT /home/$USER/git/RINRUS/bin/reduce_wwPDB_het_dict.txt 
```
or
```bash
export REDUCE_HET_DICT=/home/$USER/git/RINRUS/bin/reduce_wwPDB_het_dict.txt
```
Further preprocessing of ligands likely required!


## How to use

The entire processed-PDB-to-input-file RINRUS workflow can be run at once using the driver, or steps can be run individually.

Full usage instructions are described in [bin/](bin/README.md).


## Contributors
This code was conceptualized in the DeYonker group at the University of Memphis Department of Chemistry.
Prof. Qianyi Cheng was the primary contributor to early versions of RINRUS, who now runs her independent research group at University of Memphis.
Additional code, conceptualization, and documentation have been provided by
Prof. Nathan DeYonker, Dr. Dominique Wappett, Dr. Thomas Summers, Dr. Donatus Agbaglo, Dr. Tejas Suhagia, Dr. Taylor Santaloci, Prof. Jose Fernando Ruggiero Bachega, (Universidade Federal de Ciências da Saúde de Porto Alegre), and Dr. Eric Berquist (Q-Chem, Inc.)
