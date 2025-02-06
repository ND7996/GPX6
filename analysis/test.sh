#!/bin/bash

# Function to create folders with structured greedy decreasing approach
generate_folders() {
    local base_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE"
    mkdir -p "$base_dir"  # Ensure base directory exists

    # Given variant numbers
    local variants=(48 47 52 99 54 144 177 74 178 143 87 142 104 102 139 24 181 4 60 107)
    local total_variants=${#variants[@]}

    echo "Starting folder creation..."

    local folder_count=0

    # Iterate over each variant to create exactly 20 folders
    for ((i=0; i<total_variants; i++)); do
        local start_variant="${variants[i]}"
        
        # First, create the single variant folder
        single_folder="$base_dir/$start_variant"
        mkdir -p "$single_folder"
        echo "Created folder: $single_folder"
        folder_count=$((folder_count + 1))

        local combination="$start_variant"  # Start with the single variant

        # Create 19 more unique folders for this variant
        for ((j=i+1, count=1; j<total_variants && count<20; j++, count++)); do
            combination="${combination}_${variants[j]}"  # Append next variant
            folder_name="$base_dir/$combination"
            mkdir -p "$folder_name"  # Create folder
            echo "Created folder: $folder_name"
            folder_count=$((folder_count + 1))
        done
    done

    echo "Folder generation complete! Total folders created: $folder_count"
}

# Run the folder generation
generate_folders
