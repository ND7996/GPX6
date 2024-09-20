#!/bin/sh

# List of directories to process
directories="/home/hp/nayanika/github/GPX6/structures/S47A_F48Y.pdb"

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
            q_genrelax.py "$dir/genrelax.proc" --top "" --pdb "$pdb_file" --fep "$dir/GPX6_wtmousecys.fep" --outdir "$base_scr_dir/$system_name/minim" --rs "$dir/run_qdyn_5.sh"
        else
            echo "No .pdb file found in $dir"
        fi
    done
done

