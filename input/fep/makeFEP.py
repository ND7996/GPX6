from Qpyl.qmakefep import make_fep
fepstring = make_fep(qmap_file="fepqmapmousecys.qmap", 
                        pdb_file="GPX6cys_mouse.pdb", 
                        forcefield="oplsaa",
                        parm_files=["qoplsaa_all.prm"],
                        lib_files=["qoplsaa.lib", "GPX.lib"])
open("GPX6_wtmousecys.fep", "w").write(fepstring)