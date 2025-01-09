import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter
import glob
import os

# Define the path to the directory containing relax files and minim.pdb
data_path = "/home/hp/results/C49U/T60A/minim"

# Find all relax.dcd files in the specified directory
dcd_files = sorted(glob.glob(os.path.join(data_path, "relax_*.dcd")))  # Sort to maintain numerical order

# Reference PDB file
pdb_file = os.path.join(data_path, "minim.pdb")  # The PDB file in the same folder

# Check if files are found
if not dcd_files:
    raise FileNotFoundError(f"No relax_*.dcd files found in the directory: {data_path}")

# Load the reference structure using the PDB file
u = mda.Universe(pdb_file, dcd_files[0])  # Initialize with the first .dcd file

# Output file name for the combined trajectory
output_file = os.path.join(data_path, "combined_trajectory.dcd")

# Write the combined trajectory
with DCDWriter(output_file, u.atoms.n_atoms) as w:
    for dcd in dcd_files:
        print(f"Processing {dcd}...")
        u.load_new(dcd)  # Load the next .dcd file
        for ts in u.trajectory:  # Iterate over all frames in the trajectory
            w.write(u.atoms)

print(f"Combined trajectory saved as {output_file}.")
