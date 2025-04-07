#!/bin/bash

# Base directory
BASE_DIR="/home/nsekhar/stepwise/MUT/step1/HUMAN/level3"

# Updated mutation folders
MUTATIONS=("A178T" "A47S" "A52T" "A60T" "C99R" "A74G" "F104Y" "L139F"
           "L24I" "Q144H" "Q177H" "Q54T" "N3K" "S102G" "H173R" "S142P"
           "S143E" "S181R" "Y48F")

# Loop through each mutation folder
for MUT in "${MUTATIONS[@]}"; do
    MINIM_DIR="$BASE_DIR/$MUT/minim"
    RUN_SCRIPT="$MINIM_DIR/run_qdyn_5.sh"

    echo "Checking directory: $MINIM_DIR"

    # Check if the minim directory exists
    if [ -d "$MINIM_DIR" ]; then
        echo "Entering $MINIM_DIR"

        # Check if the script exists and is executable
        if [ -f "$RUN_SCRIPT" ]; then
            chmod +x "$RUN_SCRIPT"  # Ensure script is executable
            echo "Submitting job for $MUT"
            (cd "$MINIM_DIR" && sbatch run_qdyn_5.sh)  # Submit job via sbatch
        else
            echo "Error: run_qdyn_5.sh not found in $MINIM_DIR"
        fi
    else
        echo "Skipping $MUT: minim directory not found."
    fi
done

echo "All jobs submitted."

