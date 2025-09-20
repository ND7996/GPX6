#!/bin/bash
# Directory containing PDB files
pdb_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/level16"
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
readlib /home/hp/nayanika/github/GPX6/parameters/qoplsaa.lib
readlib /home/hp/nayanika/github/GPX6/parameters/GPX.lib
readprm /home/hp/nayanika/github/GPX6/parameters/qoplsaa_all2.prm
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
echo "Total PDB files: $(ls "$pdb_dir"/*.pdb | grep -v "_solvated.pdb" | wc -l)" | tee -a "$log_file"
echo "Successfully solvated: $(ls "$pdb_dir"/*_solvated.pdb | wc -l)" | tee -a "$log_file"
echo "Failed: $(($(ls "$pdb_dir"/*.pdb | grep -v "_solvated.pdb" | wc -l) - $(ls "$pdb_dir"/*_solvated.pdb | wc -l)))" | tee -a "$log_file"
