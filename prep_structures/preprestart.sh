#!/bin/bash

# Base directory containing all system folders
base_dir="/home/hp/results/C49U"

# Paths to required files
lib1="/home/hp/nayanika/github/GPX6/parameters/qoplsaa.lib"
lib2="/home/hp/nayanika/github/GPX6/parameters/GPX.lib"
prm="/home/hp/nayanika/github/GPX6/parameters/qoplsaa_all2.prm"
relax_file="relax_012.re"  # Make sure this file is accessible in the same directory or provide the full path

# Iterate through each subdirectory in the base directory
for system_dir in "$base_dir"/*; do
    if [ -d "$system_dir/minim" ]; then
        minim_dir="$system_dir/minim"
        
        # Find the .top file in the minim directory
        top_file=$(find "$minim_dir" -maxdepth 1 -name "*.top" | head -n 1)
        
        # Check if the .top file exists
        if [ -f "$top_file" ]; then
            # Create a temporary input file for the commands
            input_file="$minim_dir/generate_minim.inp"
            output_pdb="$minim_dir/minim.pdb"

            # Write the commands to the input file
            cat <<EOF > "$input_file"
! Read in library files
readlib $lib1
readlib $lib2
! Read in parameter files
readprm $prm
! Read in structure
readtop $(basename "$top_file")
rx $relax_file
writepdb $(basename "$output_pdb") y
quit
EOF

            # Navigate to the minim directory and run qdyn (or appropriate command)
            cd "$minim_dir" || exit
            qprep5 < generate_minim.inp  
            
            # Log the status
            if [ -f "$output_pdb" ]; then
                echo "Generated $output_pdb successfully."
            else
                echo "Failed to generate $output_pdb in $minim_dir."
            fi
            
            # Cleanup
            rm -f "$input_file"
        else
            echo "No .top file found in $minim_dir."
        fi
    else
        echo "No 'minim' directory in $system_dir."
    fi
done
