#!/bin/bash

# Define partition mappings
declare -A partition_map=(
    [A178T]="normal1"
    [A52T]="normal1"
    [A60T]="normal2"
    [A74G]="normal2"
    [F104Y]="normal3"
    [H173R]="normal3"
    [L139F]="normal4"
    [L24I]="normal4"
    [N3K]="normal5"
    [Q144H]="normal5"
    [Q177H]="normal5"
    [Q54T]="gpu"
    [S102G]="gpu"
)

# Base directory
base_dir="/home/nsekhar/stepwise/MUT/step1/HUMAN/level3/"

# Loop through all relevant subdirectories
for sub_dir in "${!partition_map[@]}"; do
  full_path="${base_dir}${sub_dir}"
  
  if [ -d "$full_path" ]; then
    for replica_dir in "$full_path"/replica*/; do
      # Specify the slurm file name
      slurm_file="${replica_dir}run_qdyn_5.sh"

      # Check if the slurm job file exists
      if [ -f "$slurm_file" ]; then
        # Update the partition based on the mapping
        sed -i "s/^#SBATCH --partition=.*$/#SBATCH --partition=${partition_map[$sub_dir]}/" "$slurm_file"

        # Add the memory line if it doesn't already exist
        if ! grep -q "^#SBATCH --mem=" "$slurm_file"; then
          sed -i '/^#SBATCH --partition=/a #SBATCH --mem=2G              # Memory per node (adjust as needed)' "$slurm_file"
          echo "Added memory specification to $slurm_file"
        fi

        echo "Updated partition in $slurm_file to ${partition_map[$sub_dir]}"
      else
        echo "run_qdyn_5.sh not found in $replica_dir"
      fi
    done
  else
    echo "Directory $full_path does not exist."
  fi
done

