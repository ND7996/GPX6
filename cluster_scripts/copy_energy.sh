#!/bin/bash

# Define the main directory containing replica folders
BASE_DIR="/home/nsekhar/stepwise/MUT/step1/MOUSE_OLD/level0/R181S"

# Navigate to the base directory
cd "$BASE_DIR" || { echo "Failed to enter $BASE_DIR"; exit 1; }

# Loop through all replicaXXX folders
for replica_folder in replica[0-9][0-9][0-9]; do
    # Check if the replica folder exists before processing
    if [ -d "$replica_folder" ]; then
        # Extract numeric part and format as repXX
        replica_number=${replica_folder#replica}
        destination_folder=$(printf "rep%02d" "$replica_number")

        # Create destination folder if it doesn't exist
        mkdir -p "$destination_folder"

        # Copy only .en files to the destination folder
        cp "$replica_folder"/*.en "$destination_folder"/ 2>/dev/null
    else
        echo "Folder $replica_folder does not exist in $BASE_DIR."
    fi

done

echo "All .en files copied successfully!"
~                                                                                                                                                                                      
~                                                        
