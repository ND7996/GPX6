#!/bin/bash

# Directory containing the stats file
DATA_DIR="/home/hp/nayanika/github/GPX6/analysis"
# File containing the existing tables
STATS_FILE="$DATA_DIR/stats.tex"
# File where the combined table will be stored
TABLE_FILE="/home/hp/nayanika/github/GPX6/table/Free_Energy.tex"

# Start writing the complete LaTeX document to the output file
{
    echo "\documentclass{article}"
    echo "\usepackage{amsmath}"  # For math symbols
    echo "\begin{document}"

    echo "\begin{table}[ht]"
    echo "    \centering"
    echo "    \begin{tabular}{|c|c|c|}"
    echo "    \hline"
    echo "    Sample & Mean dG* (kcal/mol) & Mean dG0 (kcal/mol) \\\\"
    echo "    \hline"
} > "$TABLE_FILE"

# Read the existing stats file and extract the relevant lines
while IFS= read -r line; do
    # Sanitize the line
    clean_line=$(echo "$line" | tr -d '\r' | sed 's/[^[:print:]]//g')

    # Check if the line contains table data (i.e., if it contains '&' indicating columns)
    if echo "$clean_line" | grep -q '\&'; then
        # Escape underscores and backslashes, and handle \pm
        clean_line=$(echo "$clean_line" | sed -e 's/_/\\_/g' -e 's/\\pm/\$\\pm\$/g' -e 's/\$/\\$/g')

        # Debug output
        echo "Debug: $clean_line"  # Print to terminal for debugging

        # Add line to the table
        echo "    $clean_line \\\\" >> "$TABLE_FILE"
    fi
done < <(grep -A 100 '\begin{tabular}' "$STATS_FILE" | grep -B 100 '\end{table}')

# Close the table and document
{
    echo "    \hline"  # Add final horizontal line
    echo "    \end{tabular}"
    echo "    \caption{Free Energy Changes}"
    echo "\end{table}"
    echo "\end{document}"
} >> "$TABLE_FILE"

# Print the final LaTeX file for inspection
echo "Contents of $TABLE_FILE:"
cat "$TABLE_FILE"

# End of script
