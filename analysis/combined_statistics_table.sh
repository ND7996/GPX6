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
    echo "\usepackage{graphicx}" # If you need additional graphics support
    echo "\begin{document}"
    echo "\begin{table}[ht]"
    echo "    \centering"
    echo "    \begin{tabular}{|c|c|c|}"
    echo "    \hline"
    echo "    Distance & 10 Å &  \\\\"  # Add the row with distance info
    echo "    System & Mean dG* (kcal/mol) & Mean dG0 (kcal/mol) \\\\"
    echo "    \hline"
    echo "    \hline"
} > "$TABLE_FILE"

# Read the stats file and extract relevant data
while IFS= read -r line; do
    # Clean up the line and substitute +- with the proper LaTeX \pm
    clean_line=$(echo "$line" | tr -d '\r' | sed 's/\+-/\\pm/g')

    # Check if the line contains data (includes keywords like WT, S47A, F48Y, T54Q, R99C, Cys, Sec, etc.)
    if echo "$clean_line" | grep -q "WT\|S47A\|F48Y\|T52A\|T54Q\|R99C\|Cys\|Sec\|Mouse\|Human"; then
        # Append the cleaned line to the table, preserving its format, and add \hline for row separation
        echo "$clean_line \\\\" >> "$TABLE_FILE"
        echo "    \hline" >> "$TABLE_FILE"
    fi
done < "$STATS_FILE"

# Close the table and document
{
    echo "    \end{tabular}"
    echo "    \caption{Free Energy Changes}"
    echo "\end{table}"
    echo "\end{document}"
} >> "$TABLE_FILE"

# Print the final LaTeX file for inspection
echo "Contents of $TABLE_FILE:"
cat "$TABLE_FILE"
