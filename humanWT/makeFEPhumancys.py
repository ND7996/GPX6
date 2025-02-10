from Qpyl.qmakefep import make_fep
fepstring = make_fep(qmap_file="fepqmaphumancys.qmap", 
                        pdb_file="GPX6cys_human.pdb", 
                        forcefield="oplsaa",
                        parm_files=["qoplsaa_all.prm"],
                        lib_files=["qoplsaa.lib", "GPX.lib"])
open("GPX6_wthumancys.fep", "w").write(fepstring)