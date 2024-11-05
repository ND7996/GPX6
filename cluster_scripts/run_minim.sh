#!/bin/bash

# Base directory where the folders are located
base_dir="/home/nsekhar/stepwise/MUT/step1/mousecys"

# List of directories to process
directories=(
    "D119E" "E137D" "I182T" "K3N" "N107S" "N192K" "N31Q" "P67N" 
    "Q33K" "R181S" "S195N" "T194F" "T71I" "V30I" "F29Y" "I24L" 
    "M188L" "N120K" "N22L" "N69G" "Q184K" "R173H" "S126T" "S4R" 
    "T60A" "V16I"
)

# Loop through each specified directory
for dir_name in "${directories[@]}"; do
    dir_path="$base_dir/$dir_name"
    
    # Check if the directory exists
    if [[ -d "$dir_path" ]]; then
        minim_dir="$dir_path/minim"
        run_qdyn_script="$minim_dir/run_qdyn_5.sh"

        # Check for the minim subdirectory and script
        if [[ -d "$minim_dir" && -f "$run_qdyn_script" ]]; then
            # Submit a job for this directory
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=$dir_name     # Job name for the directory
#SBATCH --error=${dir_name}_err.txt      # Standard error file
#SBATCH --output=${dir_name}_out.txt      # Standard output file
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=1          # Number of CPU cores per task (adjust as needed)
#SBATCH --mem=16G              # Memory per node (adjust as needed)
#SBATCH --time=24:00:00           # Time limit hrs:min:sec (adjust as needed)
#SBATCH --partition=normal1,normal2,normal3,normal4,normal5,highmem,gpu

# Change to the minim directory
cd "$minim_dir"

# Add execute permission and run the script
chmod +x "$run_qdyn_script"
./run_qdyn_5.sh

EOF
        else
            echo "minim directory or run_qdyn_5.sh not found in $dir_path"
        fi
    else
        echo "Directory $dir_path does not exist"
    fi
done

echo "Completed submitting jobs for all directories."

