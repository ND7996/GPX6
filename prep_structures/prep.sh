#!/bin/bash

# Directory containing PDB files
pdb_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"

# Loop through each PDB file
for pdb_file in "$pdb_dir"/*.pdb; do
    base_name=$(basename "$pdb_file" .pdb)
    
    # Create input file for qprep5
    cat <<EOF > "${base_name}.inp"
readlib /home/hp/nayanika/github/GPX6/parameters/qoplsaa.lib
readlib /home/hp/nayanika/github/GPX6/parameters/GPX.lib
readprm /home/hp/nayanika/github/GPX6/parameters/qoplsaa_all2.prm

readpdb "$pdb_file"

set solvent_pack 2.7
boundary sphere 49:SE 25.
solvate 49:SE 25. grid HOH

maketop "${base_name}_solvated.top"
writetop "${base_name}_solvated.top"
writepdb "${base_name}_solvated.pdb" y

quit
EOF

    # Run qprep5 for each PDB
    qprep5 "${base_name}.inp"
done
