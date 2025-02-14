#!/bin/bash

# List of mutation folders
folders=("D148E" "E143S" "F139L" "F48Y" "G102S" "H144Q" "H177Q" "I24L" "K87T" "N107S" "P142S" "R181S" "R99C" "S47A" "S4R" "T178A" "T52A" "T54Q" "T60A" "Y104F")

# Output combined LaTeX file
combined_output="combined_latex_table.tex"

# Start writing the combined LaTeX file
{
  echo "\\documentclass{article}"
  echo "\\usepackage{booktabs}"
  echo "\\begin{document}"
  echo "\\begin{table}[ht]"
  echo "\\centering"
  echo "\\begin{tabular}{|c|c|c|}"
  echo "\\hline"
  echo "Mutation & Mean dG* & Mean dG0 \\\\"
  echo "\\hline"
} > "$combined_output"

# Loop through each folder and append its LaTeX table content
for folder in "${folders[@]}"; do
  if [ -f "$folder/latex_table.tex" ]; then
    echo "Adding results from $folder..."
    tail -n +2 "$folder/latex_table.tex" >> "$combined_output"  # Skip the first header line
  else
    echo "Skipping $folder (latex_table.tex not found)."
  fi
done

# Finish LaTeX document
{
  echo "\\end{tabular}"
  echo "\\caption{Combined Free Energy Changes}"
  echo "\\end{table}"
  echo "\\end{document}"
} >> "$combined_output"

# Notify the user
echo "Combined LaTeX table written to $combined_output"

# Show LaTeX file content
cat "$combined_output"
