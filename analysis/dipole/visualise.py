import os
from pymol import cmd
from pymol.cgo import CYLINDER

# Function to load dipole data from file
def extract_dipole_coordinates(dipole_file):
    """Extract dipole vector from the file."""
    if not os.path.exists(dipole_file):
        print(f"Warning: File {dipole_file} does not exist.")
        return None
    
    with open(dipole_file, 'r') as f:
        lines = f.readlines()
        
    # Check for a valid dipole vector (assuming it's the first line for simplicity)
    for line in lines:
        if line.startswith("dipole"):
            try:
                # Extract x, y, z components of the dipole vector
                dipole_vector = list(map(float, line.split()[1:]))
                if len(dipole_vector) == 6:  # Ensure it's a valid 3D dipole vector
                    return dipole_vector
            except ValueError:
                print(f"Error reading dipole data in {dipole_file}.")
                return None
    return None

# Function to visualize dipole
def visualize_dipole(dipole_vector, name="dipole"):
    """Visualize the dipole vector in PyMOL."""
    if dipole_vector is None:
        return

    # Create the dipole as a cylinder (start and end points)
    dipole_obj = [
        CYLINDER,
        dipole_vector[0], dipole_vector[1], dipole_vector[2],  # Start coordinates
        dipole_vector[3], dipole_vector[4], dipole_vector[5],  # End coordinates
        0.3,  # Radius
        1.0, 0.0, 0.0,  # Color start (red)
        1.0, 0.0, 0.0   # Color end (red)
    ]

    # Load the dipole into PyMOL
    cmd.load_cgo(dipole_obj, name)

# Path to PDB file
pdb_file = '/home/hp/results/MOUSE/level3/F48Y/minim/minim.pdb'
cmd.load(pdb_file, 'minim')

# Loop through replicas and load their dipole data
for i in range(16):  # Assuming replica000 to replica015
    reactant_dipole_file = f"/home/hp/results/MOUSE/level3/F48Y/replica{i:03d}/reactant_dipole.py"
    product_dipole_file = f"/home/hp/results/MOUSE/level3/F48Y/replica{i:03d}/product_dipole.py"

    print(f"Processing replica {i:03d}...")

    # Load and visualize reactant dipole
    reactant_dipole = extract_dipole_coordinates(reactant_dipole_file)
    if reactant_dipole:
        visualize_dipole(reactant_dipole, name=f"reactant_dipole_{i:03d}")
    else:
        print(f"Warning: No dipole vector found in {reactant_dipole_file}")

    # Load and visualize product dipole
    product_dipole = extract_dipole_coordinates(product_dipole_file)
    if product_dipole:
        visualize_dipole(product_dipole, name=f"product_dipole_{i:03d}")
    else:
        print(f"Warning: No dipole vector found in {product_dipole_file}")

# Display the protein as cartoon and recenter
cmd.show("cartoon", "minim")
cmd.zoom()
cmd.recenter()

print("Visualization completed.")
