#!/bin/bash

# Step 1: Run q_mapper and q_analysefeps.py and store output in test.out
q_mapper.py 13.52 1.30 --bins 50 --skip 100 --min 10 --temp 300.0 --dirs rep* --qfep_exec qfep5
q_analysefeps.py rep* > test.out

# Step 2: Extract statistics with grep and output next 5 lines
grep Statistics test.out -A 5 > stats_output.txt

# Output LaTeX table file
output_file="latex_table.tex"

# Step 3: Extract Mean and Std.error for dG* and dG0
dg_star=$(grep 'dG\*' stats_output.txt | awk '{printf "%s +- %s kcal/mol", $2, $5}')
dg0=$(grep 'dG0' stats_output.txt | awk '{printf "%s +- %s kcal/mol", $2, $5}')

# Step 4: Write LaTeX table to file with proper LaTeX formatting
{
  echo "\\begin{table}[ht]"
  echo "\\centering"
  echo "\\begin{tabular}{|c|c|c|}"
  echo "\\hline"
  echo "  & Mean dG* & Mean dG0 \\\\"
  echo "\\hline"
  echo "WTmousecys & $dg_star & $dg0 \\\\"
  echo "\\hline"
  echo "\\end{tabular}"
  echo "\\caption{Free Energy Changes}"
  echo "\\end{table}"
} > "$output_file"

# Notify the user
echo "LaTeX table written to $output_file"
