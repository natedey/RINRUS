load res_25.pdb
cmd.select("sel","res_25 and not (chain A and resi 203) and not (name NH1 or name NH2)")
cmd.h_add("sel")
cmd.save("./res_25_h.pdb")
