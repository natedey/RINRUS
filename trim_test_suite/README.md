## RINRUS trimming test suite

This test suite is used to check that the RINRUS trimming procedure is working properly.
Scripts are designed for a slurm queueing system using array jobs to effectively parallelise the tests. 

The atoms in each residue are treated as belonging to four sections: side chain (SC), N-terminus (NT), C-terminus (CT) and the alpha carbon/hydrogen(s) (A).
Each residue can be trimmed to a maximum of 7 possible models depending on which atoms are selected to be included:
* SC: side chain only
* NT: peptide bond containing the residue's N-terminus
* CT: peptide bond containing the residue's C-terminus
* SC_NT: side chain and N-terminus peptide bond
* SC_CT: side chain and C-terminus peptide bond
* NT_CT: N-terminus peptide bond and C-terminus peptide bond
* SC_NT_CT: side chain, N-terminus peptide bond and C-terminus peptide bond
CA/HA will be included in all possible representations of the residue. If CA and/or its hydrogen(s) are the only atoms selected from a residue, it is trimmed to the NT_CT model so that CA is anchored by the main chain instead of represented as methane.


### Testing trimming of individual residues

The script 'test_single_residues.sh' tests whether each set of possible atom combinations creates the model expected (e.g. does RINRUS generate the SC_NT model from any and every selection of at least one side chain atom and at least one N-terminus atom). 

Optional arguments:
* -r "res" (e.g. -r "Pro" or -r "Cys,Met")
* -s "atom combination subset" (e.g. -s "SC,SC_NT,SC_CT,SC_NT_CT" to test all combination subsets containing any side chain atoms)
* -c "keep"/"reuse" (options for the combinations files used in the tests)

WARNING: for a residue with N atoms, there are 2^N - 1 atom combinations. The files containing the lists of combinations can get big! To save space, the default behaviour is to create only the ones needed for the tests and then delete them at the end. They don't take that long to generate compared to how long it takes to actually do the tests but if you want to keep them to look at them after the tests are done or something then use the optional keep flag. If you've already got the combinations files, you can reuse them in the next run. But be careful that you haven't deleted any or modified them accidentally because the script doesn't check that the reused ones are correct.


### Testing trimming of two or more adjacent residues

Coming soon...
