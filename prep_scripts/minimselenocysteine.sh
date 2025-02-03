#!/bin/bash

# Directories setup
pdb_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"
base_scr_dir="/home/hp/results/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"
run_dir="/home/hp/nayanika/github/GPX6/cluster_scripts"
genrelax_dir="/home/hp/nayanika/github/GPX6/input"
topology_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"
fep_dir="/home/hp/nayanika/github/GPX6/input/fep"

# Create the base directory if it does not exist
mkdir -p "$base_scr_dir"

# Iterate through all .pdb files in the pdb_dir
for pdb_file in "$pdb_dir"/*.pdb; do
    [ -e "$pdb_file" ] || continue  # Skip if no .pdb files found
    
    system_name=$(basename "$pdb_file" .pdb)  # Extract system name
    topology_file="$topology_dir/$system_name.top"
    system_dir="$base_scr_dir/$system_name"

    # Check if topology file exists
    if [ ! -f "$topology_file" ]; then
        echo "Warning: Topology file not found for $system_name: $topology_file"
        continue
    fi

    # Create system directory
    mkdir -p "$system_dir"

    # Construct and run q_genrelax.py command
    q_genrelax.py "$genrelax_dir/selgenrelax.proc" \
        --top "$topology_file" \
        --pdb "$pdb_file" \
        --fep "$fep_dir/GPX6_wtmousesec.fep" \
        --outdir "$system_dir/minim" \
        --rs "$run_dir/run_qdyn_5.sh"
    
    # Check if command succeeded
    if [ $? -eq 0 ]; then
        echo "Successfully processed: $system_name"
    else
        echo "Error processing $system_name"
    fi
done
