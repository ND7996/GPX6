import time
from pymol import cmd
import os

# Define file paths
pdb_path = "/home/hp/results/C49U/S47A/replica000/equil.pdb"
traj_path = "/home/hp/results/C49U/S47A/replica000/combined_trajectory.dcd"

# Check if the files exist
print(f"Checking files at: {pdb_path} and {traj_path}")
print(f"File exists (PDB): {os.path.exists(pdb_path)}")
print(f"File exists (Trajectory): {os.path.exists(traj_path)}")

# Load the static structure first
if os.path.exists(pdb_path):
    cmd.load(pdb_path)
else:
    print(f"Error: PDB file not found at {pdb_path}")

# Load the combined trajectory file
if os.path.exists(traj_path):
    cmd.load_traj(traj_path)
else:
    print(f"Error: Trajectory file not found at {traj_path}")

# Define the mutant residue (replace with your mutant residue)
mutant_residue = "X47"  # Example: "A100" or "X47" (adjust to your mutant residue)
residue_49 = "49"  # Residue 49

# Define a distance threshold for hydrogen bonds (e.g., 3.5 Å)
distance_threshold = 3.5

# Select the mutant residue (adjust chain and residue as needed)
cmd.select("mutant", f"resi {mutant_residue}")  # Replace with your mutant residue ID
cmd.select("residue49", f"resi {residue_49}")  # Select residue 49

# Iterate through the trajectory frames
for frame in range(1, cmd.count_frames() + 1):
    # Load the current frame from the trajectory
    cmd.frame(frame)

    # Select mutant and neighboring residues within a 3.5 Å radius
    cmd.select("hbond_mutant_neighbors", "(byres (mutant around 3.5)) and (hydrogen)")

    # Visualize the hydrogen bonds
    cmd.distance("hbonds", "mutant", "hbond_mutant_neighbors")
    cmd.color("red", "hbonds")  # Color the hydrogen bonds red

    # Show the selected atoms (mutant and residue 49)
    cmd.show("spheres", "mutant")
    cmd.show("sticks", "residue49")

    # Label the mutant and residue 49
    cmd.label("mutant", "\"Mutant\"")
    cmd.label("residue49", "\"Residue 49\"")

    # Color the mutant and residue 49 for distinction
    cmd.color("yellow", "mutant")
    cmd.color("blue", "residue49")

    # Optional: Adjust the view to zoom in on the relevant residues
    cmd.zoom("all", 5)  # Adjust the zoom level
    
    # Pause for 100 milliseconds (0.1 seconds)
    time.sleep(0.1)

# Optional: Save the resulting frames as an image or movie
# cmd.save_movie("/home/hp/results/hbonds_movie.mp4", frame_range=(1, cmd.count_frames()))
