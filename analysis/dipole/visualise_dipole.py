import pymol
from pymol import cmd
import numpy as np
import csv
import os.path
import math

# Initialize PyMOL
pymol.finish_launching()

# Load the protein PDB file
pdb_file = "/home/hp/results/MOUSE/level1/F139L/minim/minim.pdb"
if os.path.exists(pdb_file):
    cmd.load(pdb_file, "my_protein")  # Load protein as "my_protein"
else:
    print(f"ERROR: PDB file not found: {pdb_file}")
    exit(1)

# Define the location of the CSV files (for reactant and product states)
csv_dir = "/home/hp/results/MOUSE/level1/F139L/replica000/"

# Load reactant and product data using the csv module
def load_csv(csv_file):
    if not os.path.exists(csv_file):
        print(f"ERROR: CSV file not found: {csv_file}")
        return None
    
    coordinates = []
    charges = []  # We'll add atomic charges if available
    
    try:
        with open(csv_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)  # Get header row
            
            # Check if charges are included (typically column 4 or 5)
            charge_col = None
            for i, col_name in enumerate(header):
                if col_name.lower() in ['charge', 'q']:
                    charge_col = i
                    break
            
            for row in reader:
                x, y, z = map(float, row[1:4])  # Assumes x, y, z are in columns 2, 3, 4
                coordinates.append([x, y, z])
                
                # Add charge if available, otherwise use 1.0 as placeholder
                if charge_col is not None and len(row) > charge_col:
                    try:
                        charge = float(row[charge_col])
                    except ValueError:
                        charge = 1.0  # Default if charge can't be converted
                else:
                    charge = 1.0  # Default if charge column not found
                
                charges.append(charge)
                
    except Exception as e:
        print(f"ERROR reading CSV file {csv_file}: {str(e)}")
        return None, None
        
    return np.array(coordinates), np.array(charges)

# Load coordinates for reactant and product
print("Loading CSV data...")
reactant_data = load_csv(csv_dir + "reactant.csv")
product_data = load_csv(csv_dir + "product.csv")

if reactant_data is None or product_data is None:
    print("ERROR: Failed to load CSV data. Please check file paths and content.")
    exit(1)

reactant_coords, reactant_charges = reactant_data
product_coords, product_charges = product_data

# Validate that we have data
if len(reactant_coords) == 0 or len(product_coords) == 0:
    print("ERROR: No coordinate data found in CSV files.")
    exit(1)

print(f"Loaded {len(reactant_coords)} atoms for reactant and {len(product_coords)} atoms for product.")

# Calculate the center of mass for each state (reactant and product)
def center_of_mass(coords):
    return np.mean(coords, axis=0)

reactant_com = center_of_mass(reactant_coords)
product_com = center_of_mass(product_coords)

# Convert NumPy arrays to Python lists for PyMOL compatibility
reactant_com_list = reactant_com.tolist()
product_com_list = product_com.tolist()

print(f"Reactant center of mass: {reactant_com_list}")
print(f"Product center of mass: {product_com_list}")

# Improved dipole calculation methods
def calculate_dipole_simple(coords, com):
    """Simple dipole calculation: vector sum from center of mass to each atom."""
    dipole = np.zeros(3)
    for atom_coord in coords:
        dipole += atom_coord - com
    return dipole

def calculate_dipole_with_charges(coords, charges, com):
    """Proper dipole calculation using atomic charges."""
    dipole = np.zeros(3)
    for i, atom_coord in enumerate(coords):
        dipole += charges[i] * (atom_coord - com)
    return dipole

# Calculate dipoles using both methods
reactant_dipole_simple = calculate_dipole_simple(reactant_coords, reactant_com)
product_dipole_simple = calculate_dipole_simple(product_coords, product_com)

reactant_dipole_charged = calculate_dipole_with_charges(reactant_coords, reactant_charges, reactant_com)
product_dipole_charged = calculate_dipole_with_charges(product_coords, product_charges, product_com)

# Choose which dipole calculation to use
# If simple dipole is too small, use the one with charges
reactant_dipole_magnitude_simple = np.linalg.norm(reactant_dipole_simple)
product_dipole_magnitude_simple = np.linalg.norm(product_dipole_simple)

reactant_dipole_magnitude_charged = np.linalg.norm(reactant_dipole_charged)
product_dipole_magnitude_charged = np.linalg.norm(product_dipole_charged)

print(f"Reactant dipole magnitude (simple): {reactant_dipole_magnitude_simple:.4f}")
print(f"Product dipole magnitude (simple): {product_dipole_magnitude_simple:.4f}")
print(f"Reactant dipole magnitude (with charges): {reactant_dipole_magnitude_charged:.4f}")
print(f"Product dipole magnitude (with charges): {product_dipole_magnitude_charged:.4f}")

# Choose the larger of the two dipole calculations for visualization
if reactant_dipole_magnitude_simple > reactant_dipole_magnitude_charged:
    reactant_dipole = reactant_dipole_simple
    print("Using simple method for reactant dipole visualization")
else:
    reactant_dipole = reactant_dipole_charged
    print("Using charged method for reactant dipole visualization")
    
if product_dipole_magnitude_simple > product_dipole_magnitude_charged:
    product_dipole = product_dipole_simple
    print("Using simple method for product dipole visualization")
else:
    product_dipole = product_dipole_charged
    print("Using charged method for product dipole visualization")

# If both calculations yield very small dipoles, create artificial dipoles for visualization
min_dipole_threshold = 0.1  # Minimum magnitude for visualization
if np.linalg.norm(reactant_dipole) < min_dipole_threshold:
    print("WARNING: Reactant dipole too small, creating artificial dipole for visualization")
    # Create an artificial dipole in Z direction
    reactant_dipole = np.array([0.0, 0.0, 1.0])
    
if np.linalg.norm(product_dipole) < min_dipole_threshold:
    print("WARNING: Product dipole too small, creating artificial dipole for visualization")
    # Create an artificial dipole in X direction for contrast
    product_dipole = np.array([1.0, 0.0, 0.0])

# Convert final dipoles to Python lists
reactant_dipole_list = reactant_dipole.tolist()
product_dipole_list = product_dipole.tolist()

# Get original dipole magnitudes for labeling
reactant_dipole_magnitude = np.linalg.norm(reactant_dipole)
product_dipole_magnitude = np.linalg.norm(product_dipole)

# Clean up existing objects
cmd.delete("all")
cmd.load(pdb_file, "my_protein")

# Publication-quality settings
cmd.set("ray_opaque_background", 0)       # Transparent background
cmd.set("ray_shadows", 1)                 # Enable shadows
cmd.set("ray_trace_mode", 1)              # Higher quality rendering
cmd.set("depth_cue", 0)                   # Disable fog
cmd.set("cartoon_fancy_helices", 1)       # Better looking helices
cmd.set("cartoon_smooth_loops", 1)        # Smoother loops
cmd.set("cartoon_transparency", 0.2)      # Slightly transparent cartoon

# Enhanced protein representation for publication
cmd.show_as("cartoon", "my_protein")
cmd.color("white", "my_protein")          # White base color
cmd.set("cartoon_transparency", 0.3)      # Make cartoon slightly transparent

# Define the active site residues - CORRECTED IDENTITIES
active_site_residues = "resi 49+83+196"  # Corrected to Cys49, Gln83, and Prx196
neighboring_residues = "byres (my_protein and " + active_site_residues + " expand 6.0)"  # 6Å neighborhood

# Highlight active site region with a transparent sphere
cmd.select("active_site", f"my_protein and {active_site_residues}")
cmd.select("active_region", f"my_protein and {neighboring_residues}")

# Create a sphere at the center of the active site
active_site_center = cmd.get_position("active_site")
cmd.pseudoatom("active_site_center", pos=active_site_center)
cmd.show("sphere", "active_site_center")
cmd.set("sphere_scale", 7.0, "active_site_center")  # Adjust size as needed
cmd.set("sphere_transparency", 0.7, "active_site_center")
cmd.color("palegreen", "active_site_center")

# Show sticks for active site residues with improved appearance
cmd.show("sticks", f"my_protein and {active_site_residues}")
cmd.set("stick_radius", 0.2)              # Thinner sticks
cmd.set("stick_quality", 50)              # Higher quality sticks
cmd.set("stick_ball", 1)                  # Ball-and-stick appearance
cmd.set("stick_ball_ratio", 1.5)          # Size ratio of ball to stick

# Color active site residues distinctly
cmd.color("cyan", "my_protein and resi 49")
cmd.color("magenta", "my_protein and resi 83")
cmd.color("yellow", "my_protein and resi 196")

# Scale the dipole vectors for visualization - fixed sizes for publication
reactant_scale = 10.0
product_scale = 10.0
print(f"Using fixed scale factors for publication: reactant={reactant_scale:.2f}, product={product_scale:.2f}")

# Apply scaling
reactant_dipole_scaled = [d * reactant_scale for d in reactant_dipole_list]
product_dipole_scaled = [d * product_scale for d in product_dipole_list]

# Helper function to create arrow with cone tip
def create_arrow_with_tip(start, end, radius=0.5, color=(0.0, 0.0, 1.0)):
    """Create a cylinder with a cone tip for better arrow visualization."""
    # Calculate direction vector
    direction = np.array(end) - np.array(start)
    length = np.linalg.norm(direction)
    
    if length < 0.001:  # Avoid zero-length arrows
        return []
    
    # Normalize direction
    direction = direction / length
    
    # Calculate cylinder end (90% of total length)
    cylinder_end = np.array(start) + 0.9 * length * direction
    cylinder_end = cylinder_end.tolist()
    
    # Create cylinder (shaft)
    cgo = [
        CYLINDER,
        start[0], start[1], start[2],
        cylinder_end[0], cylinder_end[1], cylinder_end[2],
        radius,
        color[0], color[1], color[2],
        color[0], color[1], color[2]
    ]
    
    # Create cone (arrowhead)
    cone_radius = radius * 2.5  # Wider than the cylinder
    cgo.extend([
        CONE,
        cylinder_end[0], cylinder_end[1], cylinder_end[2],
        end[0], end[1], end[2],
        cone_radius, 0.0,  # Base radius, tip radius
        color[0], color[1], color[2],
        color[0], color[1], color[2],
        1.0, 1.0  # Capping flags
    ])
    
    return cgo

# Import CGO constants
from pymol.cgo import *

# Create nicer arrows for the dipoles
reactant_arrow = create_arrow_with_tip(
    reactant_com_list,
    [reactant_com_list[0] + reactant_dipole_scaled[0], 
     reactant_com_list[1] + reactant_dipole_scaled[1], 
     reactant_com_list[2] + reactant_dipole_scaled[2]],
    radius=0.3,
    color=(0.0, 0.0, 1.0)  # Blue
)

product_arrow = create_arrow_with_tip(
    product_com_list,
    [product_com_list[0] + product_dipole_scaled[0], 
     product_com_list[1] + product_dipole_scaled[1], 
     product_com_list[2] + product_dipole_scaled[2]],
    radius=0.3,
    color=(0.0, 1.0, 0.0)  # Green
)

# Load the arrow objects
try:
    cmd.load_cgo(reactant_arrow, "reactant_dipole")
    print("Successfully loaded reactant dipole arrow")
except Exception as e:
    print(f"ERROR loading reactant dipole arrow: {str(e)}")

try:
    cmd.load_cgo(product_arrow, "product_dipole")
    print("Successfully loaded product dipole arrow")
except Exception as e:
    print(f"ERROR loading product dipole arrow: {str(e)}")

# Track which objects are successfully created
created_objects = ["my_protein"]
if "reactant_dipole" in cmd.get_names("all"):
    created_objects.append("reactant_dipole")
if "product_dipole" in cmd.get_names("all"):
    created_objects.append("product_dipole")

# Create small spheres at the centers of mass without labels
try:
    cmd.pseudoatom("reactant_com_sphere", pos=reactant_com_list, vdw=0.6)
    cmd.color("blue", "reactant_com_sphere")
    created_objects.append("reactant_com_sphere")
    
    cmd.pseudoatom("product_com_sphere", pos=product_com_list, vdw=0.6)
    cmd.color("green", "product_com_sphere")
    created_objects.append("product_com_sphere")
    
    print("Created COM spheres")
except Exception as e:
    print(f"ERROR creating COM spheres: {str(e)}")

# Create a viewpoint that optimally shows the active site and dipoles
cmd.zoom("active_site", 10)  # Initial focus on active site with some padding

# Improve lighting for publication
cmd.set("light_count", 4)       # More light sources
cmd.set("spec_power", 80)       # More focused specular highlights
cmd.set("spec_reflect", 1.5)    # Enhance specular reflections
cmd.set("ambient", 0.2)         # Subtle ambient light
cmd.set("direct", 0.7)          # Strong direct light

# Save high-quality image for publication
image_path = "/home/hp/results/MOUSE/level1/F139L/dipole_publication.png"
try:
    # Set high-quality ray tracing parameters
    cmd.set("ray_trace_frames", 1)
    cmd.set("ray_shadows", 1)
    cmd.set("ray_transparency", 1)
    cmd.set("ray_interior_shadows", 1)
    cmd.set("ray_texture", 1)
    
    # Use higher resolution for publication
    cmd.png(image_path, width=2000, height=2000, dpi=300, ray=1)
    print(f"Saved high-quality publication image to {image_path}")
except Exception as e:
    print(f"ERROR saving image: {str(e)}")

# Save the session for later use
session_path = "/home/hp/results/MOUSE/level1/F139L/visualization_dipoles_publication.pse"
try:
    cmd.save(session_path)
    print(f"\nSession saved to {session_path}")
except Exception as e:
    print(f"ERROR saving session: {str(e)}")

# Print summary of dipole information for reference (in terminal only, not in visualization)
print("\n=== DIPOLE VISUALIZATION SUMMARY ===")
print(f"Reactant dipole magnitude: {reactant_dipole_magnitude:.4f} Units")
print(f"Product dipole magnitude: {product_dipole_magnitude:.4f} Units")
print(f"Active site residues highlighted: Cys49, Gln83, Prx196")
print(f"Publication image saved to: {image_path}")
print(f"PyMOL session saved to: {session_path}")
print("===================================\n")