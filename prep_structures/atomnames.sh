#!/bin/bash

# Directory containing PDB files
pdb_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"

# Process each PDB file
for pdb_file in "$pdb_dir"/*.pdb; do
    echo "Processing: $pdb_file"
    
    # Create a temporary file
    temp_file="${pdb_file}.tmp"
    
    # Make the replacements in correct order
    # First change original HB2 to HB1
    sed '/SEC    49/s/ HB2 / HB1 /' "$pdb_file" > "$temp_file"
    
    # Then change HB3 to HB2 (and this won't get converted to HB1)
    sed -i '/SEC    49/s/ HB3 / HB2 /' "$temp_file"
    
    # Replace H to HN
    sed -i '/SEC    49/s/ H   / HN  /' "$temp_file"
    
    # Replace HG to HG1
    sed -i '/SEC    49/s/ HG  / HG1 /' "$temp_file"
    
    # Move temporary file back to original
    mv "$temp_file" "$pdb_file"
    
    echo "Completed processing: $pdb_file"
done

echo "All files processed!"