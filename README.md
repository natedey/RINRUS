# RINRUS
residue interaction network-based residue selector

# After get a raw pdb file (raw.pdb), check ambiguous atoms, clean up
# Run reduce, generate a new H added pdb file (raw_h.pdb) 
# Check the new pdb file, if there is metal, replace with an atom with same coordination (such as replace Mg with O)
# Check all ligands, make sure H atoms were added correctly (may need to delete or add more H based on certain condition)
# If there is atom "CA" or "CB" in ligands, replace with "CA'" and "CB'"
# After previous 3 check stpes, save the file to a new file (raw_h_modify.pdb)
# Use this new pdb file (raw_h_modify.pdb), run probe and save result to *.probe
# Run probe2rins.py
# Modify freq_per_res.dat file for trim residues
# Use raw_h.pdb, freq_per_res.dat, and other needed files run probe_freq_2pdb.py
# Run pymol_scripts.py to add H to the trimmed pdb files for differnt size model
