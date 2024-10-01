#!/bin/bash

# Define the input files
file1="hbond_mousecys.info"
file2="hbond_humansec.info"

# Define output files for lost and gained H-bonds
lost_output="lost_hbonds.txt"
gained_output="gained_hbonds.txt"

# Extract H-bonds and store them in temporary files, excluding unwanted lines
grep -vE "(Finding|Models used|#|Constraints|no hydrogen|N/A)" "$file1" | awk '{print $1, $2, $3, "->", $4, $5, $6}' | sort > parsed_mousecys.txt
grep -vE "(Finding|Models used|#|Constraints|no hydrogen|N/A)" "$file2" | awk '{print $1, $2, $3, "->", $4, $5, $6}' | sort > parsed_humansec.txt

# Find lost H-bonds (in file1 but not in file2)
comm -13 parsed_humansec.txt parsed_mousecys.txt > "$lost_output"

# Find gained H-bonds (in file2 but not in file1)
comm -23 parsed_humansec.txt parsed_mousecys.txt > "$gained_output"

# Clean up temporary files
rm parsed_mousecys.txt parsed_humansec.txt

# Output results
echo "Lost H-bonds saved to '$lost_output'"
echo "Gained H-bonds saved to '$gained_output'"
