load res_15.pdb
cmd.select("sel","res_15 and not (chain A and resi 203) and not (name NH1 or name NH2)")
cmd.h_add("sel")
cmd.save("./res_15_h.pdb")
