#!/bin/bash

# Step 1: Run q_mapper and q_analysefeps.py and store output in test.out
q_mapper.py 50 48.5 --bins 50 --skip 100 --min 10 --temp 300.0 --dirs rep* --qfep_exec qfep5
q_analysefeps.py rep* > test.out

# Step 2: Extract statistics with grep and output next 5 lines
grep Statistics test.out -A 5 > stats_output.txt

# Output LaTeX table file
output_file="latex_table.tex"

# Step 3: Extract Mean and Std.error for dG* and dG0
dg_star=$(grep 'dG\*' stats_output.txt | awk '{printf "%s +- %s kcal/mol", $2, $5}')
dg0=$(grep 'dG0' stats_output.txt | awk '{printf "%s +- %s kcal/mol", $2, $5}')

# Step 4: Write LaTeX table to file with proper LaTeX formatting
cat <<EOL > $output_file
\\begin{table}[ht]
\\centering
\\begin{tabular}{|c|c|c|}
\\hline
  & Mean dG* & Mean dG0 \\\\
\\hline
WTmousecys & $dg_star & $dg0 \\\\
\\hline
\\end{tabular}
\\caption{Free Energy Changes}
\\end{table}
EOL

# Notify the user
echo "LaTeX table written to $output_file"

