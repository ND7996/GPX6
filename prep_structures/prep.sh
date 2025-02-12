#!/bin/bash

# Directory containing PDB files
pdb_dir="../../GPX6/prep_structures/MOUSE/level0"

# Loop through each PDB file
for pdb_file in "$pdb_dir"/*.pdb; do
    base_name=$(basename "$pdb_file" .pdb)
    
    # Create input file for qprep5 in the target directory
    cat <<EOF > "$pdb_dir/${base_name}.inp"
readlib ../../GPX6/parameters/qoplsaa.lib
readlib ../../GPX6/parameters/GPX.lib
readprm ../../GPX6/parameters/qoplsaa_all2.prm

readpdb "$pdb_file"

set solvent_pack 2.7
boundary sphere 49:SG 25.
solvate 49:SG 25. grid HOH

maketop "$pdb_dir/${base_name}_solvated.top"
writetop "$pdb_dir/${base_name}_solvated.top"
writepdb "$pdb_dir/${base_name}_solvated.pdb" y

quit
EOF

    # Run qprep5 for each PDB, specifying the input file path within the directory
    qprep5 "$pdb_dir/${base_name}.inp"
done
