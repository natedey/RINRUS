# RINRUS

Residue Interaction Network-based ResidUe Selector (RINRUS) is a QM-cluster model building tool.  Starting from a raw PDB file, after running a series of preparation tasks, the tool will
- select important residues for chemical reactions, and
- generate trimmed PDB files with the corresponding quantum chemical inputs.

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
- pymol
  - If installing via conda, it's under `-c conda-forge pymol-open-source`.
- openbabel (required for Arpeggio)

For certain scripts (optional),
- matplotlib
- BioPython

### External dependencies

- [probe](https://github.com/rlabduke/probe) - version 2.16.130520 is packaged with RINRUS
- [reduce](https://github.com/rlabduke/reduce) - version 3.23 is packaged with RINRUS
- [arpeggio](http://biosig.unimelb.edu.au/arpeggioweb)
- obabel version = openbabel/2.4.1

Currently, precompiled binaries of probe, reduce, and arpeggio are present in `bin/`.

which require
- CMake >= 3.10
- Any C/C++ compiler suite with C++11 support

PDB HET groups (some ligands and noncanonical amino acids) can be protonated by reduce. 
The connectivity table file is included with RINRIS: bin/reduce_wwPDB_het_dict.txt
You must set a shell environment variable to allow reduce to use the reduce_wwPDB_het_dict file, 
as the default location is /local:

setenv REDUCE_HET_DICT /home/$USER/git/RINRUS/bin/reduce_wwPDB_het_dict.txt 

or

export REDUCE_HET_DICT=/home/$USER/git/RINRUS/bin/reduce_wwPDB_het_dict.txt

Further preproccesing of ligands likely required!

Current production-level use cases are described in `bin/`.

## Usage example 1 - generating a single or a few input files with probe interaction count ranking 

## Usage example 2 - generating a single or a few input files with distance-based ranking

## Usage example 3 - generating a single or a few input files with arpeggio interaction-type ranking

## Usage example 4 - generating a single or a few input files with manual ranking (from SAPT, ML, or from some scheme that doesn't yet interface with RINRUS automatically)

## GENERATE ALL THE THINGS!!! Combinatorial model building from arpeggio or probe
## Usage example 5a - Combinatorial model building from arpeggio
## Usage example 5b - Combinatorial model building from probe

### Contributors
This code was conceptualized in the DeYonker group at the University of Memphis Department of Chemistry.
Prof. Qianyi Cheng is the primary coder, and now contributes to RINRUS with her independent research group.
Additional conceptualization, documentation and code have been provided by
Prof. Nathan DeYonker, Dr. Thomas Summers (University of Nevada-Reno), Donatus Agbaglo, Tejas Suhagia, Dr. Taylor Santaloci, Prof. Jose Fernando Ruggiero Bachega, (Universidade Federal de Ciências da Saúde de Porto Alegre), and Dr. Eric Berquist (Q-Chem, Inc.)
