#!/bin/bash

# Directory containing PDB files
PDB_DIR="/home/hp/nayanika/github/GPX6/prep_structures/humansec"

# Loop through each PDB file in the directory
for pdb_file in "$PDB_DIR"/*.pdb; do
    # Create a temporary file for modifications
    temp_file="${pdb_file%.pdb}_modified.pdb"
    
    # Process each line in the PDB file with proper formatting
    awk '{
        if ($4 == "SEC" && $5 == 49) {
            if ($3 == "HB3") $3 = "HB2";
            else if ($3 == "HB2") $3 = "HB1";
            else if ($3 == "HG") $3 = "HG1";
            else if ($3 == "H") $3 = "HN";
        }
        
        # Fix formatting and spacing
        printf "%-6s%5s %-4s %-3s%4s    %-8s%-8s%-8s\n", $1, $2, $3, $4, $5, $6, $7, $8
    }' "$pdb_file" > "$temp_file"
    
    # Optionally, replace the original file with the modified one
    mv "$temp_file" "$pdb_file"
done

echo "Atom name modifications and spacing fixed for all PDB files in $PDB_DIR."
