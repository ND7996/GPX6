#!/bin/bash

# List of mutation folders
folders=("D148E" "E143S" "F139L" "F48Y" "G102S" "H144Q" "H177Q" "I24L" "K87T" "N107S" "P142S" "R181S" "R99C" "S47A" "S4R" "T178A" "T52A" "T54Q" "T60A" "Y104F")

# Output LaTeX file
output_file="latex_table.tex"

# Start writing the LaTeX table
{
  echo "\\begin{table}[ht]"
  echo "\\centering"
  echo "\\begin{tabular}{|c|c|c|}"
  echo "\\hline"
  echo "Mutation & Mean dG* & Mean dG0 \\\\"
  echo "\\hline"
} > "$output_file"

# Loop through each mutation folder
for folder in "${folders[@]}"; do
  if [ -d "$folder" ]; then
    echo "Processing $folder..."
    cd "$folder" || continue  # Move into the mutation directory

    # Run q_mapper and q_analysefeps.py
    q_mapper.py 13.52 1.30 --bins 50 --skip 100 --min 10 --temp 300.0 --dirs rep{00..15} --qfep_exec qfep5
    q_analysefeps.py rep{00..15} > test.out

    # Extract statistics with grep and output next 5 lines
    grep Statistics test.out -A 5 > stats_output.txt

    # Extract Mean and Std.error for dG* and dG0
    dg_star_mean=$(grep 'dG\*' stats_output.txt | awk '{print $2}')
    dg_star_err=$(grep 'dG\*' stats_output.txt | awk '{print $5}')
    dg0_mean=$(grep 'dG0' stats_output.txt | awk '{print $2}')
    dg0_err=$(grep 'dG0' stats_output.txt | awk '{print $5}')

    # LaTeX formatted values
    dg_star="${dg_star_mean} \$\\pm\$ ${dg_star_err} kcal/mol"
    dg0="${dg0_mean} \$\\pm\$ ${dg0_err} kcal/mol"

    # Append result to LaTeX table
    echo "    $folder & $dg_star & $dg0 \\\\" >> "../$output_file"
    echo "\\hline" >> "../$output_file"

    cd ..  # Move back to the parent directory
  else
    echo "Skipping $folder (not found)."
  fi
done

# Finish LaTeX table
{
  echo "\\end{tabular}"
  echo "\\caption{Free Energy Changes}"
  echo "\\end{table}"
} >> "$output_file"

# Notify the user
echo "LaTeX table written to $output_file"

# Show LaTeX file content
cat "$output_file"

