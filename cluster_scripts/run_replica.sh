#!/bin/bash
#SBATCH --job-name=Q     # Job name
#SBATCH --error=messages.err.txt      # Standard error file
#SBATCH --output=messages.out.txt      # Standard output file
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=1         # Number of CPU cores per task (adjust as needed)
#SBATCH --mem=16G              # Memory per node (adjust as needed)
#SBATCH --time=24:00:00           # Time limit hrs:min:sec (adjust as needed)
#SBATCH --partition=normal3

#Loop through each directory from replicamousecys000 to replicamousecys009
for i in {000..099}; do
  dir="replica$i"
  echo "Processing directory: $dir"
  # Change to the replica directory
  cd "$dir" || { echo "Failed to enter directory $dir"; continue; }
  # Submit the run_qdyn_5.sh script to SLURM
  sbatch run_qdyn_5.sh
  # Change back to the base directory
  cd ..
done

