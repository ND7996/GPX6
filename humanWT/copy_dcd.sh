#!/bin/bash

# Loop through 0 to 99
for i in $(seq 0 99); do
    # Format the folder names with leading zeros
    replica_folder=$(printf "replica%03d" $i)
    destination_folder=$(printf "traj%02d" $i)

    # Create destination folder if it doesn't exist
    mkdir -p $destination_folder

    # Check if the replica folder exists before copying files
    if [ -d "$replica_folder" ]; then
        echo "Copying from $replica_folder to $destination_folder"
        cp $replica_folder/*.dcd $destination_folder/ 2>/dev/null
    else
        echo "Folder $replica_folder does not exist."
    fi
done

echo "All files copied successfully."

