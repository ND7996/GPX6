from Qpyl.qmakefep import make_fep
fepstring = make_fep(qmap_file="fepsechuman.qmap", 
                        pdb_file="GPX6sec_human.pdb", 
                        forcefield="oplsaa",
                        parm_files=["qoplsaa_all2.prm"],
                        lib_files=["qoplsaa.lib", "GPX.lib"])
open("GPX6_wthumansec.fep", "w").write(fepstring)