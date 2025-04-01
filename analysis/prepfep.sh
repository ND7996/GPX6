#!/bin/bash

# Base directories
base_dir="/home/hp/results/MOUSE/level1"
param_dir="/home/hp/nayanika/github/GPX6/parameters"

# List of mutation folders
folders=("D148E" "F139L" "G102S" "H177Q" "K87T" "P142S" "R99C" "S4R" "T178A" "T54Q" "Y104F" "E143S" "F48Y" "H144Q" "I24L" "N107S" "R181S" "S47A" "T52A" "T60A")

for folder in "${folders[@]}"
do
    # Define paths
    system_dir="$base_dir/$folder"
    minim_dir="$system_dir/minim"
    fep_file="$system_dir/replica000_fep_050_0.000.re"
    output_pdb="$system_dir/replica000_fep_050_0.000.pdb"
    top_file="$minim_dir/${folder}_solvated.top"
    input_file="$minim_dir/generate_fep.inp"

    # Check if required files exist
    if [[ -f "$fep_file" && -f "$top_file" ]]; then
        echo "Processing $folder..."

        # Create input script for qprep5
        cat <<EOF > "$input_file"
readlib $param_dir/qoplsaa.lib
readlib $param_dir/GPX.lib
readprm $param_dir/qoplsaa_all2.prm
! Read in structure
readtop $top_file
rx $fep_file
writepdb $output_pdb y
quit
EOF

        # Run qprep5
        cd "$minim_dir" || exit
        qprep5 < generate_fep.inp

        # Check if pdb was created
        if [[ -f "$output_pdb" ]]; then
            echo "Generated $output_pdb successfully."
        else
            echo "Failed to generate $output_pdb for $folder."
        fi

        # Cleanup input file
        rm -f "$input_file"

    else
        echo "Missing files for $folder: Skipping..."
    fi
done

echo "Processing complete!"

