#!/bin/sh

# Path to the .pdb file
pdb_file="/home/hp/nayanika/github/GPX6/pdb_top/S47A_F48Y.pdb"

# Base directory where you want to create folders
base_scr_dir="/home/hp/results/mousecys"

# Directories for the input files
run_dir="/home/hp/nayanika/github/GPX6/cluster_scripts"
genrelax_dir="/home/hp/nayanika/github/GPX6/input"
topology_dir="/home/hp/nayanika/github/GPX6/pdb_top"
fep_dir="/home/hp/nayanika/github/GPX6/input/fep"  # Correct FEP directory path

# Create the mousecys directory if it does not exist
mkdir -p "$base_scr_dir"

# Extract the system name from the PDB file (removing the extension)
system_name=$(basename "$pdb_file" .pdb)

# Create a directory for the system in the base output directory
mkdir -p "$base_scr_dir/$system_name"

# Remove any existing replicas for a clean start
#rm -rf "$base_scr_dir/$system_name/replica*"

# Check if the .pdb file exists
if [ -f "$pdb_file" ]; then
    # Run q_genrelax.py with the specified directories and inputs
    q_genrelax.py "$genrelax_dir/genrelax.proc" \
    --top "$topology_dir/S47A_F48Y.top" \
    --pdb "$pdb_file" \
    --fep "$fep_dir/GPX6_wtmousecys.fep" \
    --outdir "$base_scr_dir/$system_name/minim" \
    --rs "$run_dir/run_qdyn_5.sh"
else
    echo "No .pdb file found at $pdb_file"
fi
