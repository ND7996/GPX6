#!/bin/bash
#SBATCH --job-name=Q     # Job name
#SBATCH --error=messages.err.txt      # Standard error file
#SBATCH --output=messages.out.txt      # Standard output file
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=1        # Number of CPU cores per task (adjust as needed)
#SBATCH --time=24:00:00           # Time limit hrs:min:sec (adjust as needed)
#SBATCH --mem=2G
#SBATCH --partition=normal1,normal2,normal3,normal4,normal5,highmem,gpu

OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"
steps=`ls -1v *inp | sed 's/.inp//'`
echo $steps
#steps=( $(ls -1v *inp | sed 's/.inp//') )

#rs=$((1 + $RANDOM % 1000000))
#sed -i s/987654321/$rs/ equil_000_*.inp

for step in $steps
do
  echo "Running equilibration step ${step}"

  #if mpirun -np 6 Qdyn6p ${step}.inp > ${step}.log
  if  qdyn5_r8 ${step}.inp > ${step}.log
  then
    echo -e "$OK"
  else 
    echo -e "$FAILED"
    echo "Check output (${step}.log) for more info."
    exit 1
  fi
done
