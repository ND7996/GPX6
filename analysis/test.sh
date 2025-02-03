#!/bin/bash

# Function to create folders using a greedy decreasing approach
generate_folders() {
    local base_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE"
    mkdir -p "$base_dir"  # Ensure base directory exists

    # Given variant numbers in the specified sequence
    local variants=(48 47 52 99 54 144 177 74 178 143 87 142 104 102 139 24 181 4 60 107)

    # Associative array to track created combinations
    declare -A used_combinations
    local combinations_count=0
    local total_variants=${#variants[@]}

    echo "Starting folder creation..."

    # Iterate over each variant as a starting point
    for ((i=0; i<total_variants && combinations_count<210; i++)); do
        local combination=""

        # Generate decreasing combinations from the current variant
        for ((j=i; j<total_variants && combinations_count<210; j++)); do
            if [ -z "$combination" ]; then
                combination="${variants[j]}"
            else
                combination="${combination}_${variants[j]}"  # Add underscore separator between variants
            fi

            # Create folder only if combination hasn't been used
            if [ -z "${used_combinations[$combination]}" ]; then
                folder_name="$base_dir/$combination"
                mkdir -p "$folder_name"  # Create folder
                used_combinations[$combination]=1
                combinations_count=$((combinations_count + 1))
                echo "Created folder: $folder_name"
            fi
        done
    done

    echo "Folder generation complete!"
}

# Run the folder generation
generate_folders

# PREP_STRUCTURES (Generate PDB files)
echo "Running prep_structure_mouse.py..."
python3 /home/hp/nayanika/github/GPX6/prep_structures/prep_structure_mouse.py
echo "PDB files generated!"

# PREP_TOPOLOGY (Generate topology files inside the respective folders)
# Specify the target directory for saving PDB and Topology files
target_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"

# Ensure the directory exists
mkdir -p "$target_dir"

# Running prep.sh inside the folder where it's located
echo "Running prep.sh in $target_dir..."
cd "$target_dir"  # Change to the specific folder

# Run the prep.sh script to generate the topology and PDB files
sh "$target_dir/prep.sh"
echo "Topology and PDB files generated in $target_dir!"

# (CHANGE SOME ATOMS)(EXTRA STEP) 
# Run the atomnames script to generate the atom names in PDB files
sh "$target_dir/atomnames.sh"
echo "Atom names changed in $target_dir!" 

cp -r /home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107 /home/hp/results/MOUSE
