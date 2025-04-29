#!/usr/bin/env python3
import csv
import os
import subprocess
import tempfile
import argparse

# Create a temporary PyMOL script
def create_pymol_script(csv_file, pdb_file):
    script_content = f"""
import csv
from pymol import cmd, cgo

def draw_arrow(start, end, color=[1, 0, 0], name='dipole_arrow', radius=0.2, head_size=0.7):
    \"\"\"
    Draw an arrow from start to end with the given color.
    The arrow points FROM start TO end (direction: start -> end)
    \"\"\"
    arrow = [
        cgo.CYLINDER, *start, *end, radius, *color, *color,
        cgo.CONE, *end, *start, head_size, 0.0, *color, *color, 1.0, 0.0
    ]
    cmd.load_cgo(arrow, name)

# Clear previous objects
cmd.delete("all")

# Load the protein structure
cmd.load("{pdb_file}", "my_protein")
cmd.show("cartoon", "my_protein")
cmd.color("gray", "my_protein")

# Create a list to store arrow centers for focusing
center_points = []

# Read the CSV file and create arrows
with open("{csv_file}", 'r') as f:
    reader = csv.reader(f)
    header = next(reader, None)  # Skip header
    
    for i, row in enumerate(reader):
        try:
            # Get atom coordinates
            x1, y1, z1 = float(row[0]), float(row[1]), float(row[2])  # Atom_49
            x2, y2, z2 = float(row[3]), float(row[4]), float(row[5])  # Atom_83
            x3, y3, z3 = float(row[6]), float(row[7]), float(row[8])  # Atom_196
            
            # Store the actual atom coordinates for later use
            atom_49 = [x1, y1, z1]
            atom_83 = [x2, y2, z2]
            atom_196 = [x3, y3, z3]
            
            # Calculate vector midpoints directly between the atoms (but don't offset them)
            # This places arrows exactly on the atoms instead of offset
            center_points.append([(x1+x2+x3)/3, (y1+y2+y3)/3, (z1+z2+z3)/3])
            
            # Create arrows directly between atoms
            # FROM Atom 49 TO Atom 83 (red)
            draw_arrow(atom_49, atom_83, color=[1, 0, 0], name=f"dipole_arrow_49_83_{{i}}", radius=0.15)
            
            # FROM Atom 49 TO Atom 196 (blue)
            draw_arrow(atom_49, atom_196, color=[0, 0, 1], name=f"dipole_arrow_49_196_{{i}}", radius=0.15)
            
            # FROM Atom 83 TO Atom 196 (green)
            draw_arrow(atom_83, atom_196, color=[0, 1, 0], name=f"dipole_arrow_83_196_{{i}}", radius=0.15)
            
            # Create visible spheres at the exact atom positions
            cmd.pseudoatom("atom_49", pos=atom_49, label="49", color=[1, 1, 0])  # Yellow
            cmd.pseudoatom("atom_83", pos=atom_83, label="83", color=[1, 0.5, 0])  # Orange
            cmd.pseudoatom("atom_196", pos=atom_196, label="196", color=[0.5, 0, 0.5])  # Purple
            
            # Set label options
            cmd.set("label_size", 14)
            cmd.set("label_position", [0, 0, 0])  # Place label at atom center
            
        except (ValueError, IndexError) as e:
            print(f"Error with row {{i}}: {{e}}")

# Select the actual residues by ID instead of just using atom numbers
# First, try to find the atoms directly
try:
    cmd.select("atom_49", "id 49")
    cmd.select("atom_83", "id 83")
    cmd.select("atom_196", "id 196")
except:
    pass

# If direct atom IDs don't work, try residue numbers
try:
    cmd.select("res_49", "resi 49")
    cmd.select("res_83", "resi 83")
    cmd.select("res_196", "resi 196")
except:
    pass

# MODIFICATION: Enhance active site representation with sticks
cmd.show("sticks", "res_49 or res_83 or res_196")
cmd.color("yellow", "res_49")
cmd.color("orange", "res_83")
cmd.color("purple", "res_196")

# MODIFICATION: Highlight side chains to better show orientation
cmd.select("active_site", "res_49 or res_83 or res_196")
cmd.show("sticks", "active_site")
cmd.set("stick_radius", 0.25, "active_site")  # Make sticks thicker for better visibility

# MODIFICATION: Show orientation of side chains in active site
cmd.set("stick_color", "0")  # Force sticks to use assigned colors

# Enhance visibility of active site
cmd.set("stick_ball", "on", "active_site")
cmd.set("stick_ball_ratio", 2.0, "active_site")
cmd.set("stick_ball_color", "atomic", "active_site")

# Add surface representation around active site to show spatial orientation
cmd.select("around_active", "byres (active_site around 6.0)")
cmd.show("surface", "around_active")
cmd.set("transparency", 0.85, "around_active")
cmd.color("white", "around_active and not active_site")

# Enhance visibility
cmd.set("cgo_line_width", 2)
cmd.set("cgo_sphere_quality", 4)
cmd.set("sphere_scale", 0.4)

# Focus view on the center of all arrows and active site
if center_points:
    avg_x = sum(p[0] for p in center_points) / len(center_points)
    avg_y = sum(p[1] for p in center_points) / len(center_points)
    avg_z = sum(p[2] for p in center_points) / len(center_points)
    cmd.pseudoatom("center_focus", pos=[avg_x, avg_y, avg_z], visible=0)
    cmd.zoom("center_focus", 15)  # Tighter zoom
    cmd.delete("center_focus")

# Apply nice visualization settings
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_transparency", 0.8)  # More transparent cartoon for better active site visibility
cmd.set("sphere_transparency", 0.2)
cmd.set("sphere_quality", 3)
cmd.set("depth_cue", 1)
cmd.set("ray_shadows", 0)
cmd.set("ray_trace_mode", 1)

# MODIFICATION: Add dihedral angle labels to show orientation
cmd.dihedral("dihe1", "atom_49", "atom_83", "atom_196", "atom_49")
cmd.set("label_dihedral_digits", 1)
cmd.set("label_size", 14)

# Optional: Add distance labels between atoms
cmd.distance("dist_49_83", "atom_49", "atom_83")
cmd.distance("dist_49_196", "atom_49", "atom_196")
cmd.distance("dist_83_196", "atom_83", "atom_196")
cmd.set("dash_radius", 0.1)  # Thinner distance lines
cmd.set("dash_color", "gray50")

# MODIFICATION: Create alternate views to better show active site orientation
cmd.scene("front_view", "store")
cmd.turn("x", 90)
cmd.scene("top_view", "store")
cmd.turn("x", -90)
cmd.turn("y", 90)
cmd.scene("side_view", "store")
cmd.scene("front_view", "recall")

# Optional: Save images from different angles to show orientation
cmd.bg_color("white")
cmd.ray(1200, 1200)
cmd.png("{os.path.expanduser('~/pymol_dipole_front.png')}")
cmd.scene("top_view", "recall")
cmd.ray(1200, 1200)
cmd.png("{os.path.expanduser('~/pymol_dipole_top.png')}")
cmd.scene("side_view", "recall")
cmd.ray(1200, 1200)
cmd.png("{os.path.expanduser('~/pymol_dipole_side.png')}")
cmd.scene("front_view", "recall")

# Print legend information
print("\\nActive Site Orientation Visualization:")
print("- Red arrows: Atom 49->83")
print("- Blue arrows: Atom 49->196")
print("- Green arrows: Atom 83->196")
print("- Yellow sticks/sphere: Residue/Atom 49")
print("- Orange sticks/sphere: Residue/Atom 83")
print("- Purple sticks/sphere: Residue/Atom 196")
print("\\nMultiple view images saved:")
print(f"- Front view: {os.path.expanduser('~/pymol_dipole_front.png')}")
print(f"- Top view: {os.path.expanduser('~/pymol_dipole_top.png')}")
print(f"- Side view: {os.path.expanduser('~/pymol_dipole_side.png')}")
print("\\nUse keyboard shortcuts F1-F3 to switch between saved views")
"""
    
    # Create a temporary file
    fd, path = tempfile.mkstemp(suffix='.py', prefix='pymol_script_')
    with os.fdopen(fd, 'w') as f:
        f.write(script_content)
    
    return path

def launch_pymol_with_script(script_path):
    """Launch PyMOL with the specified script"""
    try:
        # Try different possible PyMOL command names
        pymol_commands = ["pymol", "pymol.exe", "/usr/local/bin/pymol", "/opt/homebrew/bin/pymol"]
        
        for cmd in pymol_commands:
            try:
                # Try to launch PyMOL with the script
                subprocess.run([cmd, "-r", script_path], check=True)
                print(f"Successfully launched PyMOL with {cmd}")
                return True
            except (subprocess.SubprocessError, FileNotFoundError):
                continue
        
        # If we get here, none of the commands worked
        print("Could not find PyMOL. Make sure it's installed and in your PATH.")
        return False
            
    except Exception as e:
        print(f"Error launching PyMOL: {e}")
        return False

def check_csv_format(csv_file):
    """Check if the CSV file has the expected format"""
    try:
        with open(csv_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            first_row = next(reader, None)
            
            if not first_row or len(first_row) < 9:
                print("WARNING: CSV file might not have enough columns (need at least 9 for x,y,z coordinates of 3 atoms)")
                print(f"Found {len(first_row)} columns in the first data row")
                return False
            
            # Try parsing the numbers
            for i in range(9):
                float(first_row[i])
                
            return True
    except Exception as e:
        print(f"Error checking CSV file format: {e}")
        return False

if __name__ == "__main__":
    # Default file paths from original script
    default_csv_file = "/home/hp/results/MOUSE/level1/F139L/replica000/product_atom_positions_with_charges.csv"
    default_pdb_file = "/home/hp/results/MOUSE/level1/F139L/minim/minim.pdb"
    
    # Add command-line argument support while keeping backward compatibility
    parser = argparse.ArgumentParser(description="Visualize dipoles in PyMOL based on atom coordinates")
    parser.add_argument("--csv", help="Path to the CSV file with atom coordinates", default=default_csv_file)
    parser.add_argument("--pdb", help="Path to the PDB structure file", default=default_pdb_file)
    
    args = parser.parse_args()
    
    csv_file = args.csv
    pdb_file = args.pdb
    
    # Check if files exist
    if not os.path.exists(csv_file):
        print(f"CSV file not found: {csv_file}")
        exit(1)
    
    if not os.path.exists(pdb_file):
        print(f"PDB file not found: {pdb_file}")
        exit(1)
    
    # Verify CSV format
    if not check_csv_format(csv_file):
        proceed = input("CSV format might be incorrect. Proceed anyway? (y/n): ")
        if proceed.lower() != 'y':
            exit(1)
    
    print(f"Creating PyMOL script to visualize dipoles...")
    script_path = create_pymol_script(csv_file, pdb_file)
    print(f"Temporary script created at: {script_path}")
    
    print("Launching PyMOL...")
    success = launch_pymol_with_script(script_path)
    
    if success:
        print("PyMOL launched successfully!")
    else:
        print("\nCouldn't automatically launch PyMOL.")
        print(f"Open PyMOL manually and run: run {script_path}")
    
    # Keep temporary file for manual use if needed
    print(f"\nYou can also run this script manually in PyMOL with:")
    print(f"run {script_path}")