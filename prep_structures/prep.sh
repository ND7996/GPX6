#!/bin/bash
set -euo pipefail

# Directory containing PDB files
pdb_dir="${PDB_DIR:-./prep_structures/HUMAN/level0}"
param_dir="${PARAM_DIR:-./parameters}"

# Create a log file
log_file="$pdb_dir/solvation_log.txt"
echo "Starting solvation process at $(date)" > "$log_file"

# Loop through each PDB file
for pdb_file in "$pdb_dir"/*.pdb; do
    # Skip already solvated PDB files
    if [[ "$pdb_file" == *"_solvated.pdb" ]]; then
        continue
    fi
    
    base_name=$(basename "$pdb_file" .pdb)
    echo "Processing $base_name..." | tee -a "$log_file"
    
    # Create input file for qprep5 in the target directory
    cat <<EOF > "$pdb_dir/${base_name}.inp"
readlib $param_dir/qoplsaa.lib
readlib $param_dir/GPX.lib
readprm $param_dir/qoplsaa_all2.prm
readpdb "$pdb_file"
set solvent_pack 2.7
boundary sphere 49:SE 25.
solvate 49:SE 25. grid HOH
maketop "$pdb_dir/${base_name}_solvated.top"
writetop "$pdb_dir/${base_name}_solvated.top"
writepdb "$pdb_dir/${base_name}_solvated.pdb" y
quit
EOF
    # Run qprep5 for each PDB with error checking
    qprep5 "$pdb_dir/${base_name}.inp" > "$pdb_dir/${base_name}_qprep.log" 2>&1
    
    # Check if the solvated files were created
    if [ -f "$pdb_dir/${base_name}_solvated.pdb" ] && [ -f "$pdb_dir/${base_name}_solvated.top" ]; then
        echo "Successfully created solvated structure for $base_name" | tee -a "$log_file"
    else
        echo "ERROR: Failed to create solvated structure for $base_name" | tee -a "$log_file"
        echo "Check $pdb_dir/${base_name}_qprep.log for details" | tee -a "$log_file"
    fi
done

# Report summary
echo "Solvation process completed at $(date)" | tee -a "$log_file"
echo "Summary:" | tee -a "$log_file"
total_pdb=$(find "$pdb_dir" -maxdepth 1 -name "*.pdb" ! -name "*_solvated.pdb" | wc -l)
solvated_pdb=$(find "$pdb_dir" -maxdepth 1 -name "*_solvated.pdb" | wc -l)
echo "Total PDB files: $total_pdb" | tee -a "$log_file"
echo "Successfully solvated: $solvated_pdb" | tee -a "$log_file"
echo "Failed: $((total_pdb - solvated_pdb))" | tee -a "$log_file"
