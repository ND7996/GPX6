#!/bin/bash

# Define the main directory containing mutation folders
BASE_DIR="/home/nsekhar/stepwise/MUT/step1/HUMAN/level0"

# List of mutation folders
MUTATIONS=(A178T A47S A52T A60T C99R E148D F104Y L139F L24I Q144H Q177H Q54T R4S S102G S107N S142P S143E S181R Y48F)

# Loop through each mutation folder
for MUT in "${MUTATIONS[@]}"; do
    MUT_DIR="$BASE_DIR/$MUT"
    
    # Check if the mutation directory exists
    if [ -d "$MUT_DIR" ]; then
        cd "$MUT_DIR" || { echo "Failed to enter $MUT_DIR"; continue; }
        
        # Loop through all replicaXXX folders
        for replica_folder in replica[0-2][0-9][0-9]; do
            if [ -d "$replica_folder" ]; then
                replica_number=${replica_folder#replica}
                replica_number=$((10#$replica_number))  # Convert to decimal to avoid octal issues
                destination_folder=$(printf "traj%02d" "$replica_number")
                
                mkdir -p "$destination_folder"
                
                # Copy only .dcd files
                cp "$replica_folder"/*.dcd "$destination_folder"/ 2>/dev/null
            else
                echo "Folder $replica_folder does not exist in $MUT_DIR."
            fi
        done
    else
        echo "Mutation folder $MUT_DIR does not exist."
    fi

done

echo "All .dcd files copied successfully!"

