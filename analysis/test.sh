#!/bin/bash

# Function to create folders using a greedy decreasing approach
generate_folders() {
    local base_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE"
    mkdir -p "$base_dir"
    
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
                combination="${combination}_${variants[j]}"
            fi
            
            # Create folder only if combination hasn't been used
            if [ -z "${used_combinations[$combination]}" ]; then
                folder_name="$base_dir/$combination"
                mkdir -p "$folder_name"
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

# Execute the Python notebook for structure preparation
echo "Running prep_structure_mouse.ipynb..."
jupyter nbconvert --to notebook --execute /home/hp/nayanika/github/GPX6/prep_structures/prep_structure_mouse.ipynb --output /home/hp/nayanika/github/GPX6/prep_structures/prep_structure_mouse_out.ipynb

echo "Notebook execution complete!"
