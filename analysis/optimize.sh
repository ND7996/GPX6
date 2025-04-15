#!/bin/bash

# Input file
input_file="/home/hp/results/MOUSE/level5/combined_latex_table.tex"

# Initialize variables
lowest_value=99999
lowest_mutation=""

# Process the file line by line
while IFS= read -r line; do
    # Check if line contains "kcal/mol" (indicating it's a data line)
    if [[ $line == *"kcal/mol"* ]]; then
        # Extract the mutation name
        mutation=$(echo "$line" | grep -o -E '[A-Z][0-9]+[A-Z]|Cys')
        if [[ -z "$mutation" ]]; then
            continue
        fi
        
        # Extract the first energy value
        energy=$(echo "$line" | grep -o -E '[0-9]+\.[0-9]+' | head -n 1)
        if [[ -z "$energy" ]]; then
            continue
        fi
        
        # Convert to a comparable number and check if it's the lowest
        if (( $(echo "$energy < $lowest_value" | bc -l) )); then
            lowest_value=$energy
            lowest_mutation=$mutation
        fi
    fi
done < "$input_file"

# Output the result
if [[ -n "$lowest_mutation" ]]; then
    echo "Mutation with the lowest dG*: $lowest_mutation"
    echo "Value: $lowest_value kcal/mol"
else
    echo "No valid mutation data found."
fi
