#!/bin/bash

# Directory containing PDB files
pdb_directory="/home/hp/nayanika/github/GPX6/prep_structures/humansec"

# Loop through each PDB file in the directory
for pdb_file in "$pdb_directory"/*.pdb; do
    # Temporary file to store modified lines
    temp_file="${pdb_file}.tmp"
    
    # Process each line and apply changes only to residue 49, and remove unwanted lines
    awk '
    {
        # Check for unwanted atoms and skip those lines
        if ($1 == "ATOM" && ($4 == "HOH" && ($5 == 1448 || $5 == 1449 || $5 == 1450))) next;
        if ($1 == "ATOM" && $5 == 49 && $6 == 0.000 && $7 == 0.000 && $8 == 0.000) next;
        
        # Change atom names in residue 49
        if ($1 == "ATOM" && $5 == 49) {
            if ($3 == "H")  {$3 = "HN"}
            else if ($3 == "HG") {$3 = "HG1"}
            else if ($3 == "HB3") {$3 = "HB2"}
            else if ($3 == "HB2") {$3 = "HB1"}
        }
        
        # Print the modified line with only atom name and residue info (no coordinates)
        printf "%-6s%5d %-4s %-3s %4d\n", $1, $2, $3, $4, $5
    }' "$pdb_file" > "$temp_file"

    # Replace the original file with the modified file
    mv "$temp_file" "$pdb_file"
    echo "Processed $pdb_file"
done
