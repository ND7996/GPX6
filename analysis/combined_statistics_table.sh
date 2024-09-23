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

# Read the stats file and extract relevant data
while IFS= read -r line; do
    # Sanitize the line to remove unwanted characters and add \pm where needed
    clean_line=$(echo "$line" | tr -d '\r' | sed 's/\+/\$\\pm\$/g' | sed 's/\\//g')

    # Check if the line contains the sample data
    if echo "$clean_line" | grep -q "WTmousecys"; then
        # Extract Mean dG* and Mean dG0 values using consistent field separation
        mean_dg_star=$(echo "$clean_line" | awk -F'&' '{print $2}' | sed 's/kcal\/mol//g' | sed 's/^[ \t]*//;s/[ \t]*$//')
        mean_dg0=$(echo "$clean_line" | awk -F'&' '{print $3}' | sed 's/kcal\/mol//g' | sed 's/^[ \t]*//;s/[ \t]*$//')

        # Add the extracted values to the table
        echo "    WTMOUSECYS & \$$mean_dg_star\$ kcal/mol & \$$mean_dg0\$ kcal/mol \\\\" >> "$TABLE_FILE"
    fi
done < "$STATS_FILE"

# Close the table and document
{
    echo "    \hline"
    echo "    \end{tabular}"
    echo "    \caption{Free Energy Changes}"
    echo "\end{table}"
    echo "\end{document}"
} >> "$TABLE_FILE"

# Print the final LaTeX file for inspection
echo "Contents of $TABLE_FILE:"
cat "$TABLE_FILE"
