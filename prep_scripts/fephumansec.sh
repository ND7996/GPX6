#!/bin/sh

# List of directories to process
directories="/home/hp/results/humansec/T87K/minim"
directories="$directories /home/hp/nayanika/github/GPX6/input"  # Combine both directories

# Base directory where you want to create folders
base_scr_dir="/home/hp/results/humansec/T87K"

# Create the base directory if it does not exist
mkdir -p "$base_scr_dir"

# Path to cluster scripts directory
cluster_scripts_dir="/home/hp/nayanika/github/GPX6/cluster_scripts"

# Check if run_qdyn_5.sh exists in the cluster scripts directory
if [ -f "$cluster_scripts_dir/run_qdyn_5.sh" ]; then
    run_qdyn_file="$cluster_scripts_dir/run_qdyn_5.sh"
else
    echo "run_qdyn_5.sh not found in $cluster_scripts_dir"
    exit 1
fi

# Check if run_replica.sh exists in the cluster scripts directory
if [ -f "$cluster_scripts_dir/run_replica.sh" ]; then
    run_replica_file="$cluster_scripts_dir/run_replica.sh"
else
    run_replica_file=""  # Leave empty if not found
    echo "run_replica.sh not found, proceeding without it."
fi

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
             q_genfeps.py "/home/hp/nayanika/github/GPX6/input/selgenfeps.proc" --pdb "$pdb_file" "$dir/relax_012.inp" relax --repeats 100 --frames 51 --fromlambda 1.0 --prefix "$base_scr_dir/$system_name/replica" --rs "$run_qdyn_file"

            # Copy the run_qdyn_5.sh script to each replica folder
            for replica_dir in "$base_scr_dir/$system_name/replica"*; do
                if [ -d "$replica_dir" ]; then
                    # Copy run_qdyn_5.sh from the cluster scripts directory to each replica
                    cp "$run_qdyn_file" "$replica_dir/"
                else
                    echo "Replica directory $replica_dir not found"
                fi
            done

            # Optionally, run run_replica.sh if it's present
            if [ -n "$run_replica_file" ]; then
                for replica_dir in "$base_scr_dir/$system_name/replica"*; do
                    if [ -d "$replica_dir" ]; then
                        bash "$run_replica_file" "$replica_dir"
                    else
                        echo "Replica directory $replica_dir not found, skipping."
                    fi
                done
            fi
        else
            echo "No .pdb file found in $dir"
        fi
    done
done

echo "Processing completed for all directories."
