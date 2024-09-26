#!/bin/sh

# List of directories to process
directories="/home/hp/results/humansec/C99R/minim"
directories="$directories /home/hp/nayanika/github/GPX6/input"  # Combine both directories

# Base directory where you want to create folders
base_scr_dir="/home/hp/results/humansec/C99R"

# Create the base directory if it does not exist
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
            # Run q_genfeps.py with the specified path for selgenfeps.proc
            q_genfeps.py "/home/hp/nayanika/github/GPX6/input/selgenfeps.proc" --pdb "$pdb_file" "$dir/relax_012.inp" relax --repeats 100 --frames 51 --fromlambda 1.0 --prefix "$base_scr_dir/$system_name/replica" --rs "$dir/run_qdyn_5.sh"
        else
            echo "No .pdb file found in $dir"
        fi
    done
done

echo "Processing completed for all directories."
