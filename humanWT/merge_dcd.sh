#!/bin/bash

# Destination folder where you want to copy the files
DESTINATION="combined_dcd_files"

# Create destination folder if it doesn't exist
mkdir -p $DESTINATION

# Loop through traj00 to traj99
for i in $(seq -w 0 99); do
  # Format the trajXX folder name
  traj_folder="traj$i"
  
  # Check if the directory exists
  if [ -d "$traj_folder" ]; then
    echo "Processing $traj_folder"

    # Loop through the .dcd files in the folder
    for file in "$traj_folder"/*.dcd; do
      if [ -f "$file" ]; then
        # Get the base filename
        base_name=$(basename "$file")
        
        # Copy the file with the new name prefixed with trajXX_
        cp "$file" "$DESTINATION/${traj_folder}_${base_name}"
      fi
    done
  else
    echo "$traj_folder does not exist"
  fi
done

echo "All files have been copied with the trajXX prefix."

