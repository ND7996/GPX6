#!/bin/bash

# Define directories and parameter files
param_dir="/home/hp/nayanika/github/GPX6/parameters"
pdb_dir="/home/hp/nayanika/github/GPX6/prep_structures/mousecys"
output_dir="$pdb_dir"  # Use the same directory for output files

# Loop over each specified PDB file
for pdb_file in "$pdb_dir"/*.pdb; do
    # Extract the base name (e.g., S4R from S4R.pdb)
    base_name=$(basename "$pdb_file" .pdb)

    # Write commands to a temporary file for each PDB
    cat <<EOF > "$output_dir/${base_name}.in"
! Load libraries and parameters
readlib $param_dir/qoplsaa.lib
readlib $param_dir/GPX.lib
readprm $param_dir/qoplsaa_all.prm

! Read in structure
readpdb $pdb_file

! Set parameters
set solvent_pack 2.7
boundary sphere 49:SG 25.
solvate 49:SG 25. grid HOH

! Write output files
maketop $output_dir/${base_name}.top
writetop $output_dir/${base_name}.top
writepdb $output_dir/${base_name}.pdb y

quit
EOF

    # Run the command with the input file
    qprep5 < "$output_dir/${base_name}.in"  # Replace `your_tool_command` with the actual command
done
