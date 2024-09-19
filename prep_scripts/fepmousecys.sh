#!/bin/sh

# List of directories to process
directories="/home/hp/nayanika/github/PhD_Thesis/EVB/protein_stepwise/GPX6MUT/individual_mutants/mousecys/1-prep/S47A
/home/hp/nayanika/github/PhD_Thesis/EVB/protein_stepwise/GPX6MUT/individual_mutants/mousecys/1-prep/F48Y
/home/hp/nayanika/github/PhD_Thesis/EVB/protein_stepwise/GPX6MUT/individual_mutants/mousecys/1-prep/T52A
/home/hp/nayanika/github/PhD_Thesis/EVB/protein_stepwise/GPX6MUT/individual_mutants/mousecys/1-prep/T54Q
/home/hp/nayanika/github/PhD_Thesis/EVB/protein_stepwise/GPX6MUT/individual_mutants/mousecys/1-prep/R99C"

# Base directory where you want to create folders
base_scr_dir="/home/hp/results/mousecys"

# Create the mousecys directory if it does not exist
mkdir -p "$base_scr_dir"

# Loop through each directory and process .pdb files
for dir in $directories
do
    # Construct the system name based on the directory name
    system_name=$(basename "$dir")

    # Create directory for the system in $base_scr_dir
    mkdir -p "$base_scr_dir/$system_name"

    # Remove any existing replicas for a clean start
    rm -rf "$base_scr_dir/$system_name/replica*"

    # Find the .pdb file inside the current directory
    for pdb_file in "$dir"/*.pdb
    do
        if [ -f "$pdb_file" ]; then
            # Run q_genfeps.py with the system's pdb file and other input files from the current directory
            q_genfeps.py "$dir/genfeps.proc" --pdb "$pdb_file" "$dir/relax_012.inp" relax --repeats 100 --frames 51 --fromlambda 1.0 --prefix "$base_scr_dir/$system_name/replica" --rs "$dir/run_qdyn_5.sh"
        else
            echo "No .pdb file found in $dir"
        fi
    done
done
