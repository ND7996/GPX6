#!/bin/bash

# Directory containing the stats file
DATA_DIR="/home/hp/nayanika/github/GPX6/analysis"
# File containing the existing tables
STATS_FILE="$DATA_DIR/stats.tex"
# File where the combined table will be stored
TABLE_FILE="/home/hp/nayanika/github/GPX6/table/Free_Energy.tex"

# Start writing the LaTeX document to the output file
{
    echo "\documentclass{article}"
    echo "\usepackage{amsmath}" # for \pm
    echo "\usepackage{booktabs}" # optional for better table rules
    echo "\usepackage{caption}"  # for table captions
    echo "\begin{document}"
    echo ""
    
    echo "\begin{table}[ht]"
    echo "    \centering"
    echo "    \begin{tabular}{|c|c|c|}"
    echo "    \hline"
    echo "    System & Mean dG* (kcal/mol) & Mean dG0 (kcal/mol) \\\\"
    echo "    \hline"
} > "$TABLE_FILE"

# Read the stats file and extract relevant data
while IFS= read -r line; do
    # Clean up the line and substitute +- with the proper LaTeX \pm
    clean_line=$(echo "$line" | tr -d '\r' | sed 's/+/\\pm/g')

    # Check if the line contains data (includes keywords like WT, S47A, F48Y, T54Q, R99C, Cys, Sec, etc.)
    if echo "$clean_line" | grep -qE "WT|S47A|F48Y|T54Q|R99C|Cys|Sec|Human"; then
        # Append the cleaned line to the table, preserving its format, and add \hline for row separation
        echo "$clean_line \\\\" >> "$TABLE_FILE"
        echo "    \hline" >> "$TABLE_FILE"
    fi
done < "$STATS_FILE"

# Close the table and the document
{
    echo "    \end{tabular}"
    echo "    \caption{Free Energy Changes}"
    echo "\end{table}"
    echo ""
    echo "\end{document}"
} >> "$TABLE_FILE"

# Print the final LaTeX file for inspection
echo "Contents of $TABLE_FILE:"
cat "$TABLE_FILE"
