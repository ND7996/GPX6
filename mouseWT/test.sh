for system in GPX6mousecys 47_48_52_54 
do 
mkdir $SCRATCH/$system
rm -rf $SCRATCH/$system/replica*
q_genfeps.py genfeps.proc --pdb $system.pdb relax_012.inp relax --repeats 100 --frames 51 --fromlambda 1.0 --prefix $SCRATCH/$system/replica --rs run_qdyn_5.sh
done

