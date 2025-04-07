#!/bin/bash

# Define the base directory
BASE_DIR="/home/nsekhar/stepwise/MUT/step1/HUMAN/level3"

# Define the mutation folders
MUTATIONS=(
    "A178T" "A52T" "A60T" "A74G" "F104Y" "H173R" "L139F" 
    "L24I" "N3K" "Q144H" "Q177H" "Q54T" "S102G" "S142P" "S143E" 
    "S181R" "Y48F"
)

# Iterate through each mutation folder
for mutation in "${MUTATIONS[@]}"; do
    mutation_path="${BASE_DIR}/${mutation}"
    
    # Check if mutation directory exists
    if [[ ! -d "$mutation_path" ]]; then
        echo "Mutation directory $mutation_path not found, skipping..."
        continue
    fi
    
    echo "Processing mutation: $mutation"
    
    # Change to the mutation directory
    cd "$mutation_path"
    
    # Iterate over replicas from replica000 to replica015 only
    for i in $(seq -f "replica%03g" 0 15); do
        if [[ -d "$i" && -f "$i/run_qdyn_5.sh" ]]; then
            echo "Submitting job in $mutation/$i..."
            (cd "$i" && sbatch run_qdyn_5.sh)
        else
            echo "Skipping $mutation/$i (not found or missing run_qdyn_5.sh)"
        fi
    done
done

echo "All jobs have been submitted!"

