#!/usr/bin/env python
"""
PyMOL Trajectory Distance Analysis Script
----------------------------------------
This script measures specific atom-pair distances across all frames
of a loaded trajectory and saves them to a CSV file.

Usage:
1. Load your structure and trajectory in PyMOL
2. Run this script from within PyMOL:
   run /path/to/trajectory_distances.py
3. Call the function:
   analyze_trajectory_distances("output.csv")
"""

from pymol import cmd
import os
import csv

def analyze_trajectory_distances(output_file="trajectory_distances.csv"):
    """
    Analyze distances between specific atom pairs across all frames of a trajectory
    
    Parameters:
    - output_file (str): Path to the output CSV file
    
    Returns:
    - bool: True if successful
    """
    # Define the atom pairs to measure
    atom_pairs = [
        ("49.SG", "196.O1"),
        ("83.OE1", "196.O1"),
        ("49.SG", "49.HG1"),
        ("196.O1", "196.H1"),
        ("83.OE1", "196.H1"),
        ("196.O1", "49.HG1")
    ]
    
    # Get the number of states (frames) in the trajectory
    n_states = cmd.count_states()
    if n_states <= 0:
        print("Error: No trajectory loaded or no states found.")
        return False
    
    print(f"Analyzing distances across {n_states} frames...")
    
    # Open CSV file for writing
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        # Write header row
        header = ["Frame"] + [f"{pair[0]}-{pair[1]}" for pair in atom_pairs]
        csvwriter.writerow(header)
        
        # Process each frame
        for state in range(1, n_states + 1):
            # Create a row for this frame
            row = [state]
            
            # Measure each distance in this frame
            for atom1, atom2 in atom_pairs:
                # Parse the atom info
                resi1, name1 = atom1.split('.')
                resi2, name2 = atom2.split('.')
                
                # Measure the distance in the current state
                try:
                    dist = cmd.get_distance(
                        f"resi {resi1} and name {name1}", 
                        f"resi {resi2} and name {name2}",
                        state=state
                    )
                    row.append(f"{dist:.3f}")
                except Exception as e:
                    print(f"Error measuring {atom1}-{atom2} in frame {state}: {str(e)}")
                    row.append("N/A")
            
            # Write the row for this frame
            csvwriter.writerow(row)
            
            # Display progress for long trajectories
            if state % 100 == 0 or state == n_states:
                print(f"Processed {state}/{n_states} frames")
    
    print(f"Analysis complete! Distances saved to {os.path.abspath(output_file)}")
    return True

# Print usage information when script is loaded
print("\nTrajectory Distance Analysis Script loaded.")
print("To analyze distances across all frames, use:")
print("analyze_trajectory_distances(\"output_filename.csv\")")
print("\nExample:")
print("analyze_trajectory_distances(\"my_trajectory_distances.csv\")")