#!/bin/bash
# Define the main directory containing mutation folders
BASE_DIR="/home/nsekhar/stepwise/MUT/step1/HUMAN/level2"
# List of mutation folders
MUTATIONS=(A178T A47S A52T A60T C99R E148D F104Y L139F L24I Q144H Q177H Q54T R4S S102G S107N S142P S143E S181R Y48F)
# Log file for missing .en files
LOG_FILE="missing_en_files.log"
# Output directory for saving .en files
OUTPUT_DIR="/home/nsekhar/stepwise/MUT/collected_en_files"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Clear the log file before starting
> "$LOG_FILE"

echo "Checking for missing .en files..." > "$LOG_FILE"

# Loop through each mutation folder
for MUT in "${MUTATIONS[@]}"; do
    MUT_DIR="$BASE_DIR/$MUT"
    
    # Check if the mutation directory exists
    if [ -d "$MUT_DIR" ]; then
        # Create mutation-specific output directory
        mkdir -p "$OUTPUT_DIR/$MUT"
        
        # Loop through all repXX and replicaXXX folders
        for replica_folder in "$MUT_DIR"/rep[0-9][0-9] "$MUT_DIR"/replica[0-9][0-9][0-9]; do
            if [ -d "$replica_folder" ]; then
                replica_name=$(basename "$replica_folder")
                en_files=("$replica_folder"/*.en)
                
                # Check if any .en files exist in this replica folder
                if [ -e "${en_files[0]}" ]; then
                    # Create replica-specific directory in output
                    mkdir -p "$OUTPUT_DIR/$MUT/$replica_name"
                    
                    # Copy all .en files to the output directory
                    cp "$replica_folder"/*.en "$OUTPUT_DIR/$MUT/$replica_name/"
                else
                    # Log only missing .en files for this replica
                    echo "$MUT/$replica_name" >> "$LOG_FILE"
                fi
            fi
        done
    else
        # Log if mutation folder does not exist
        echo "$MUT (folder not found)" >> "$LOG_FILE"
    fi
done

echo "Process completed. Folders with missing .en files are logged in $LOG_FILE."

