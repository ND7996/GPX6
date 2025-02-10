#!/bin/bash
#############################

#SBATCH -J relax
#SBATCH -n 1
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -D .

#############################
echo '##################################################'
echo 'running run_qdyn'
echo '##################################################'

cd /home/ndas/cysM/replicamousecys000
sh run_qdyn.sh

