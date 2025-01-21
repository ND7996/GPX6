#!/bin/bash

# Directory containing the stats file
DATA_DIR="/home/hp/nayanika/github/GPX6/analysis"
# File containing the existing tables
STATS_FILE="$DATA_DIR/statsmouse1.tex"
# File where the combined table will be stored
TABLE_FILE="/home/hp/nayanika/github/GPX6/table/Free_Energy1.tex"

# Start writing the complete LaTeX document to the output file
{
    echo "\documentclass{article}"
    echo "\usepackage{amsmath}"  # For math symbols
    echo "\usepackage{graphicx}" # If you need additional graphics support
    echo "\begin{document}"
} > "$TABLE_FILE"

# Read the stats file and extract relevant data
while IFS= read -r line; do
    # Clean up the line and substitute +- with the proper LaTeX \pm
    clean_line=$(echo "$line" | tr -d '\r' | sed 's/\+-/\\pm/g')

    # Check if the line contains one of the specific mutations (e.g., N107S, T60A)
    if echo "$clean_line" | grep -q -E "N107S|T60A|S4R|R181S|I24L|G102S|Y104F|P142S|K87T|E143S|T178A|G74A|H177Q|H144Q|Q54T|C99R|T52A|S47A|F48Y|Sec C49U|Cys"; then
        # Add the cleaned line to the table for specific mutations
        echo "\\begin{table}[ht]" >> "$TABLE_FILE"
        echo "    \\centering" >> "$TABLE_FILE"
        echo "    \\begin{tabular}{|c|c|c|}" >> "$TABLE_FILE"
        echo "    \\hline" >> "$TABLE_FILE"
        echo "    System & Mean dG* (kcal/mol) & Mean dG0 (kcal/mol) \\\\" >> "$TABLE_FILE"
        echo "    \\hline" >> "$TABLE_FILE"

        # Add the cleaned line to the table
        echo "    $clean_line \\\\" >> "$TABLE_FILE"
        echo "    \\hline" >> "$TABLE_FILE"

        # Close the table
        echo "    \\end{tabular}" >> "$TABLE_FILE"
        echo "    \\caption{Free energy changes for selected mutations}" >> "$TABLE_FILE"
        echo "\\end{table}" >> "$TABLE_FILE"
    fi
done < "$STATS_FILE"

# Close the document
{
    echo "\end{document}"
} >> "$TABLE_FILE"

# Print the final LaTeX file for inspection
echo "Contents of $TABLE_FILE:"
cat "$TABLE_FILE"
