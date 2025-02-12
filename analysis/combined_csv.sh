#!/bin/bash

# Define the directory and output file
input_dir="/home/hp/results"
output_file="/home/hp/results/combined_results.csv"

# Find all CSV files in the directory
csv_files=($(find "$input_dir" -maxdepth 1 -type f -name "*.csv"))

# Check if CSV files are found
if [ ${#csv_files[@]} -eq 0 ]; then
    echo "No CSV files found in $input_dir."
    exit 1
fi

# Extract header from the first file
head -n 1 "${csv_files[0]}" > "$output_file"

# Loop through CSV files and append data (excluding header)
for file in "${csv_files[@]}"; do
    tail -n +2 "$file" >> "$output_file"
done

# Notify the user
echo "Combined CSV file created: $output_file"
