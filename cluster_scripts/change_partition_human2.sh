#!/bin/bash

# Define high-memory mutation directories
highmem_dirs=("S142P" "S143E" "S181R" "Y48F")

# Base directory
base_dir="/home/nsekhar/stepwise/MUT/step1/HUMAN/level3"

# Loop through specified high-memory directories
for sub_dir in "${highmem_dirs[@]}"; do
  full_path="${base_dir}/${sub_dir}"

  if [ -d "$full_path" ]; then
    for replica_dir in "$full_path"/replica*/; do
      # Specify the slurm file name
      slurm_file="${replica_dir}run_qdyn_5.sh"

      # Check if the slurm job file exists
      if [ -f "$slurm_file" ]; then
        # Update the partition to highmem
        sed -i 's/^#SBATCH --partition=.*$/#SBATCH --partition=highmem/' "$slurm_file"

        # Add the memory line if it doesn't already exist
        if ! grep -q "^#SBATCH --mem=" "$slurm_file"; then
          sed -i '/^#SBATCH --partition=/a #SBATCH --mem=2G              # Memory per node (adjust as needed)' "$slurm_file"
          echo "Added memory specification to $slurm_file"
        fi

        # Remove the time limit line
        sed -i '/^#SBATCH --time=24:00:00/d' "$slurm_file"

        echo "Updated partition and removed time limit in $slurm_file"
      else
        echo "run_qdyn_5.sh not found in $replica_dir"
      fi
    done
  else
    echo "Directory $full_path does not exist."
  fi
done

