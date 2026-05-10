#!/bin/bash -l
#SBATCH -p CLUSTER-AMD
#SBATCH -n 64
#SBATCH -t 3-00:00:00
#SBATCH -J Q-md

## Load modules for qdynp
module load Q/5.10

OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

steps=( $(ls -1v *inp | sed 's/.inp//') )

for step in ${steps[@]}
do
 echo "Running step ${step}"
 if srun -n 16 qdyn5p ${step}.inp > ${step}.log
 then echo -e "$OK"
   cp ${step}.re  ${step}.re.rest
 else 
  echo -e "$FAILED"
  echo "Check output (${step}.log) for more info."
  exit 1
 fi
done
