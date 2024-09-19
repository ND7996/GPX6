#!/bin/bash

# Base directory where all the input files are located
input_dir=~/nayanika/github/PhD_Thesis/EVB/protein_stepwise/GPX6WT/mousecys/1-prep

# List the systems to be processed
for system in GPX6cys_mouserestart.pdb
do
    # Strip the `.pdb` extension to use the system name in directory creation
    system_name=$(basename "$system" .pdb)

    # Create directory for the system in $SCRATCH
    mkdir -p "$SCRATCH/$system_name"

    # Remove any existing replicas for a clean start
    rm -rf "$SCRATCH/$system_name/replica*"

    # Run q_genfeps.py with the system's pdb file and other input files from the specified folder
    q_genfeps.py "$input_dir/genfeps.proc" --pdb "$input_dir/$system" "$input_dir/relax_012.inp" relax --repeats 100 --frames 51 --fromlambda 1.0 --prefix "$SCRATCH/$system_name/replica" --rs "$input_dir/run_qdyn_5.sh"
done
