#!/bin/bash

# Define input files
file1="overlaps_allatoms_mousecys"  # Contains data from mouse Cys
file2="overlaps_allatoms_humansec"   # Contains data from human Sec

# Define output file for overlaps
overlap_output="common_overlaps.txt"

# Header information to include in output
header="Allowed overlap: -0.4\nH-bond overlap reduction: 0\nIgnore contacts between atoms separated by 2 bonds or less\nDetect intra-residue contacts: True"

# Extract contact data from both files
# Here we are saving only the residues and their positions for comparison
contacts1=$(grep -E "^[A-Z]" "$file1" | awk '{print $1, $2, $3, $4, $5}')  # Capture residues and their positions from file1
contacts2=$(grep -E "^[A-Z]" "$file2" | awk '{print $1, $2, $3, $4, $5}')  # Capture residues and their positions from file2

# Create temporary files to store contacts
temp1=$(mktemp)
temp2=$(mktemp)

# Save contacts to temporary files
echo "$contacts1" > "$temp1"
echo "$contacts2" > "$temp2"

# Sort and prepare for comparison, focusing on residue numbers (ignoring the types)
awk '{print $2, $3, $4}' "$temp1" | sort -u > /tmp/overlaps1.txt
awk '{print $2, $3, $4}' "$temp2" | sort -u > /tmp/overlaps2.txt

# Find common overlaps based on residue number
common_overlaps=$(comm -12 /tmp/overlaps1.txt /tmp/overlaps2.txt)

# Initialize an array to hold full lines of common overlaps
declare -A common_lines

# Find the full lines in both files corresponding to the common overlaps
while read -r overlap; do
    # Split the overlap into components
    res1=$(echo "$overlap" | awk '{print $1}')
    atom1=$(echo "$overlap" | awk '{print $2}')
    atom2=$(echo "$overlap" | awk '{print $3}')
    
    # Search for lines in both files that match the overlap
    line1=$(grep -E "$res1.*$atom1.*$atom2" "$file1")
    line2=$(grep -E "$res1.*$atom1.*$atom2" "$file2")
    
    # Store the lines in an associative array using residue and atom as key
    common_lines["$overlap"]+="$line1\n$line2\n"
done <<< "$common_overlaps"

# Write the output to the output file
{
    echo -e "$header"  # Write the header to output
    echo "Common overlaps:"  # Add a label for common overlaps
    echo ""  # Add an empty line for formatting
    if [[ -z "$common_overlaps" ]]; then
        echo "No common overlaps found."
    else
        # Output the full lines of common overlaps
        for overlap in "${!common_lines[@]}"; do
            echo -e "${common_lines[$overlap]}"
        done
    fi
} > "$overlap_output"

# Clean up temporary files
rm "$temp1" "$temp2" /tmp/overlaps1.txt /tmp/overlaps2.txt

# Output results
echo "Common overlaps saved to '$overlap_output'"
