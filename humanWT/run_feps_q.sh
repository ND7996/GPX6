#!/bin/bash 

OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

steps=( $(ls -1v *inp | sed 's/.inp//') )

#rs=$((1 + $RANDOM % 1000000))
#sed -i s/987654321/$rs/ equil_000_*.inp

for step in ${steps[@]}
do
  echo "Running equilibration step ${step}"
  if mpirun -np 6 Qdyn6p ${step}.inp > ${step}.log
  then
    echo -e "$OK"
  else 
    echo -e "$FAILED"
    echo "Check output (${step}.log) for more info."
    exit 1
  fi
done
