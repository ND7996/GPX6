import pymol
from pymol import cmd
import glob
import os

# Set the paths
pdb_file = "/home/hp/results/MOUSE_OLD/MOUSE/level0/R181S/minim/minim.pdb"
dcd_path = "/home/hp/results/MOUSE_OLD/MOUSE/level0/R181S/traj00"


dcd_files = sorted(glob.glob(os.path.join(dcd_path, "fep*.dcd")))

# Initialize PyMOL
pymol.finish_launching(['pymol', '-q'])  # Start in quiet mode

# Load reference structure
cmd.load(pdb_file, "structure")

# Ensure the structure object exists before loading trajectories
if not dcd_files:
    print("Error: No DCD files found!")
else:
    print(f"Found {len(dcd_files)} DCD files.")

# Load DCD files and append trajectories
for dcd in dcd_files:
    try:
        cmd.load_traj(dcd, "structure", state=0)  # Removed 'quiet'
        print(f"Loaded trajectory: {dcd}")
    except Exception as e:
        print(f"Error loading {dcd}: {e}")

# Check if trajectories loaded
num_frames = cmd.count_states("structure")
if num_frames <= 1:
    print("Warning: Only one frame detected! Check if trajectories are loading correctly.")


