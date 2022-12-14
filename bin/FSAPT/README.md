Background: The user has a cluster model PDB file that has already been 
cleaned/trimmed/protonated, preferably the "maximal" model from the 
probe or arpeggio RIN (probe in this example). A Psi4 FSAPT computation has 
been run on this PDB structure wherein the "seed" species of interest 
(e.g. the chorismate of chorismate mutase or glycine of GNMT) is listed 
as the first interacting body and the enzyme is the second interacting 
body in the Psi4 input file. The atom ordering of the first body is in 
the same order as in the PDB file, and the atom ordering of the second 
body is in the same order as in the PDB file. The overall atom indices 
do not need to be the same (e.g. it is okay that the chorismate 
is listed as the first body in the FSAPT calculation, but it does not 
appear first in the PDB file, as long as the order of the chorismate 
atoms is the same and the order of the rest of the enzyme atoms is the 
same). NOTE: If there is more than one seed fragment, the fragments must be in the
same order in the first body as they are in the pdb, even if they are noncontiguous 
fragments. 
The user wishes to identify the side/main chains and waters present in the cluster model 
and compute the FSAPT interaction energies between the first body (seed) and each 
of these protein's active site functional groups.
The provided example uses GNMT with free glycine (A:294) as the seed and is 
contained in RINRUS/examples/FSAPT/.

1) Identify functional group atoms

The gen-FG-analysis.py script is a rudimentary script that can be used 
to identify which atoms of the cluster model PDB correspond to residue main 
chains, residue side chains, and waters. Because the atom indices between 
the PDB and the FSAPT computation may be different (e.g. with glycine 
located at the end of the PDB ordering but listed as the first body in 
the FSAPT calculation) the atom indices printed in the output file will be 
shifted according to the expected atom indices in the FSAPT calculation. 
Example script usage to be run in the directory where the Psi4 FSAPT calculation
was performed:

	input : python3 ~/git/RINRUS/bin/FSAPT/gen-FG-atomIDs.py -p template_27.pdb -s A:294
	output: pdbFG.dat

The output pdbFG.dat file contains info on: how many atoms are in 
the structure (line 1), the atom names of the first body (line 2), and the 
names and corresponding shifted atom indices for the functional groups 
identified by the script (for example, the side chain of residue 7 of 
chain A would correspond to atoms 25-46 in the geom.xyz script generated 
by the FSAPT calculation in the fsapt directory for this computation).

2) Calculate all of the functional group atom interactions with the first body

Now that the atom indices of the functional groups have been 
identified and summarized in the generated pdbFG.dat file, we can 
automate the calculation of their interaction energies with the 
first body using the analyze-FG-SAPT.py script. The script needs 
to be run within the ./fsapt directory output by the FSAPT computation 
in the computation's scratch directory as this directory contains the 
data files required for calculating the energies among the user-specified 
functional groups/partitions. Assuming that the pdbFG.dat file is located 
in the top directory before the fsapt directory (this location can be changed 
with the flag -p ) the script can be run as-is. 

analyze-FG-SAPT.py uses the information within pdbFG.dat to 
generate file fA.dat, which contains the functional group information for 
the first body (by default using the whole first body [i.e. the whole 
glycine substrate], though this can be changed to indicate only a particular 
functional group of the first body [i.e. only one of the carboxylates] 
using the -a flag), and then file fB.dat, which contains the functional 
group information for the second body. The script will begin by setting 
fB.dat to the first functional group (specifically the enzyA within 
fB.dat will correspond to the atom indices of the first functional group 
and enzyB will correspond to the atom indices of the rest of the enzyme 
of the second body). The Psi4 fsapt.py script will then automatically 
calculate the FSAPT interaction energies, and then the 
specific interactions between the two bodies of interest (specifically 
enzyA of fB.dat and seedA of fA.dat) will be saved. The fB.dat file is 
then re-written using the next functional group in the pdbFG.dat file 
as the next enzyA functional group, and the process is iterated until 
all of the functional groups interaction energies has been computed. 

NOTE: analyze-FG-sapt.py requires you to point to the location of the Psi4 fsapt.py script.
You must have Psi4 1.5 or later installed for RINRUS to work with FSAPT-D calculations.

Example usage (within the /fsapt directory):

	input: python3 ~/git/RINRUS/bin/FSAPT/analyze-FG-SAPT.py -path <your local path where share/psi4/fsapt scripts are located> 
	output: ../FG-SAPT.dat 

3) (Optional) Gather probe counts and arpeggio interactions for functional group atoms
The gen-FG-analysis-probe.py script counts the contacts for 
the functional groups interacting with a user-specified seed. Requires 
a probe file to have been run on the PDB. 

Example usage:

	input: python3 gen-FG-analysis-probe.py -p 1nbh.probe -s A:294
	output: FG-probe.dat which lists the functional groups and the number of contacts with the seed

Similarly, the gen-FG-analysis-arpeggio.py script counts the interactions for 
the functional groups interacting with a user-specified seed. 	
Requires an arpeggio contacts file to have been run on the cluster PDB. 

Example usage: 

	input: python3 gen-FG-analysis-arpeggio.py -c contact_counts.dat -p template_27.pdb -s A:294 
	output: FG-arpeggio.dat which lists the functional groups and the number of given interaction types. The order of the interaction 
	types column is the same as in the contacts file, which can be found in the README of Arpeggio

4) Generate QM-cluster model input files for quantum chemistry software packages based on the unsigned magnitude of the fragment FSAPT interaction energy. This
should be equivalent to step 8 of the original workflow. The sapt2rins.py script generates res_atoms-fsapt.dat, which then allows you to create models for quantum chemistry software packages. 
Note, in this case, we want the A:293 SAM cofactor to be included in the seed even though our FSAPT calculation only had A:294 GLY as fragA. 

Example usage: 

	input: python3 sapt2rins.py -p FG-SAPT.dat -c contact_counts.dat -s A:293,A:294 
	output: res_atoms-fsapt.dat 

This brings us to Step 9 on the original workflow, but you will need an additional argument "-c" when running rinrus_trim2_pdb.py to tell RINRUS to use the 
FSAPT-based ranking.

	input: python3 rinrus_trim2_pdb.py -s A:293,A:294 -pdb template_27.pdb -c res_atoms-fsapt.dat 
	output: All of the FSAPT-ranked models in PDB format

Then the workflow should proceed as normal. 
