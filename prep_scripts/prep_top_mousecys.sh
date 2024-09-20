#!/bin/sh
#working script

# Base directory containing all the PDB files, libraries, and parameter files
pdb_dir="/home/hp/nayanika/github/PhD_Thesis/EVB/protein_stepwise/GPX6MUT/individual_mutants/mousecys/1-prep"

# Full path to the necessary library and parameter files
qoplsaa_lib="$pdb_dir/qoplsaa.lib"
gpx_lib="$pdb_dir/GPX.lib"
qoplsaa_prm="$pdb_dir/qoplsaa_all.prm"

# Corresponding PDB file names
pdb_files="47_nohyd.pdb 48_nohyd.pdb 52_nohyd.pdb 54_nohyd.pdb 99_nohyd.pdb"

# Base directory where you want to create folders and output files
base_scr_dir="/home/hp/results/topology/mousecys"

# Create the base directory if it does not exist
mkdir -p "$base_scr_dir"

# Loop through each PDB file
for pdb_file in $pdb_files
do
    echo "Processing file: $pdb_file"
    
    # Construct the system name based on the PDB file (e.g., "47" from "47_nohyd.pdb")
    system_name=$(echo $pdb_file | cut -d'_' -f1)

    # Create directory for the system in $base_scr_dir
    mkdir -p "$base_scr_dir/$system_name"

    # Check if the specified PDB file exists in the pdb_dir
    if [ -f "$pdb_dir/$pdb_file" ]; then
        echo "Found PDB file: $pdb_file"

        pdb_filename=$(basename "$pdb_file" .pdb)

        # Generate the correct filenames based on your format (e.g., "47_mouse.top", "47_mouse.pdb")
        output_top="${system_name}_mouse.top"
        output_pdb="${system_name}_mouse.pdb"

        # Generate input script for each system
        inp_file="$base_scr_dir/$system_name/${system_name}_mouse.inp"
        cat > "$inp_file" <<EOF
!new opls parameters and substrate parameters (ffld?)
readlib  $qoplsaa_lib
readlib  $gpx_lib

! parameters
readprm  $qoplsaa_prm

! read in structure
readpdb $pdb_file

set solvent_pack 2.7

! Solvate around C49
boundary   sphere 49:SG 25.
solvate           49:SG 25. grid HOH

! write out topology and pdb with corresponding indices 
maketop $output_top
writetop $output_top
writepdb $output_pdb y

quit
EOF

        # Output information
        echo "Generated input script for $pdb_filename in $base_scr_dir/$system_name"
        echo "Expected output files: $output_top and $output_pdb"

        # Check for atom 49:SG in the PDB file
        if grep -q "49 SG" "$pdb_dir/$pdb_file"; then
            echo "Atom 49:SG found in PDB."
        else
            echo "Warning: Atom 49:SG not found in PDB. Skipping this file."
            continue
        fi

        # Run qprep5 on the generated input file to produce the topology and PDB files
        cd "$base_scr_dir/$system_name"
        qprep5 < "$inp_file"

    else
        echo "PDB file $pdb_file not found in $pdb_dir"
    fi
done
