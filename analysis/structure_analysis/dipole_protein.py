#!/usr/bin/env python3
import csv
import os
import subprocess
import tempfile

# Create a temporary PyMOL script
def create_pymol_script(csv_file, pdb_file):
    script_content = f"""
import csv
from pymol import cmd, cgo

def draw_arrow(start, end, color=[1, 0, 0], name='dipole_arrow'):
    \"\"\"
    Draw an arrow from start to end with the given color.
    The arrow points FROM start TO end (direction: start -> end)
    \"\"\"
    radius = 0.2
    head_size = 0.7
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
            
            # Calculate vector midpoints to place arrows properly on residues
            mid_49_83 = [(x1+x2)/2, (y1+y2)/2, (z1+z2)/2]
            mid_49_196 = [(x1+x3)/2, (y1+y3)/2, (z1+z3)/2]
            mid_83_196 = [(x2+x3)/2, (y2+y3)/2, (z2+z3)/2]
            
            # Store midpoint for camera focus
            center_points.append(mid_49_83)
            
            # Create arrows pointing FROM Atom 49 to Atom 83 (red)
            draw_arrow([x1, y1, z1], [x2, y2, z2], color=[1, 0, 0], name=f"dipole_arrow_49_83_{{i}}")
            
            # Create arrows pointing FROM Atom 49 to Atom 196 (blue)
            draw_arrow([x1, y1, z1], [x3, y3, z3], color=[0, 0, 1], name=f"dipole_arrow_49_196_{{i}}")
            
            # Create arrows pointing FROM Atom 83 to Atom 196 (green)
            draw_arrow([x2, y2, z2], [x3, y3, z3], color=[0, 1, 0], name=f"dipole_arrow_83_196_{{i}}")
            
            # Optional: Label the atoms for better understanding
            cmd.pseudoatom("label_49", pos=[x1, y1, z1], label="49")
            cmd.pseudoatom("label_83", pos=[x2, y2, z2], label="83")
            cmd.pseudoatom("label_196", pos=[x3, y3, z3], label="196")
            cmd.set("label_size", 14)
            cmd.set("label_position", [0, 0, 2])
            
        except (ValueError, IndexError) as e:
            print(f"Error with row {{i}}: {{e}}")

# Enhance visibility
cmd.set("cgo_line_width", 3)
cmd.set("cgo_sphere_quality", 4)

# Focus view on center of all arrows
if center_points:
    avg_x = sum(p[0] for p in center_points) / len(center_points)
    avg_y = sum(p[1] for p in center_points) / len(center_points)
    avg_z = sum(p[2] for p in center_points) / len(center_points)
    cmd.pseudoatom("center_atom", pos=[avg_x, avg_y, avg_z])
    cmd.zoom("center_atom", 20)  # Adjust zoom level as needed
    cmd.delete("center_atom")

# Ensure arrows are visible
cmd.show("cgo", "*arrow*")

# Show relevant atoms as spheres for better context
# Try both methods of atom selection to ensure compatibility with different PDB formats
try:
    # Try first by atom ID
    cmd.select("dipole_atoms", "id 49+83+196")
    if cmd.count_atoms("dipole_atoms") == 0:
        # If that fails, try by residue number
        cmd.select("dipole_atoms", "resi 49+83+196")
    if cmd.count_atoms("dipole_atoms") == 0:
        # If that also fails, don't highlight any atoms
        print("Could not select atoms 49, 83, and 196 - check your PDB file structure")
    else:
        cmd.show("spheres", "dipole_atoms")
        cmd.color("yellow", "dipole_atoms")
        cmd.set("sphere_scale", 0.5, "dipole_atoms")
except:
    print("Could not select atoms - check your PDB file structure")

# Apply nice visualization settings
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_transparency", 0.5)
cmd.set("sphere_transparency", 0.5)
cmd.set("stick_radius", 0.2)
cmd.set("sphere_quality", 3)
cmd.set("depth_cue", 1)
cmd.set("ray_shadows", 0)

# Optional: Save image
cmd.bg_color("white")
cmd.ray(1200, 1200)
cmd.png("{os.path.expanduser('~/pymol_dipole_visualization.png')}")

# Replace cmd.comment() with print() for PyMOL message output
print("Red: 49->83, Blue: 49->196, Green: 83->196")
print("\\nVisualization complete!")
print("- Red arrows: Atom 49->83")
print("- Blue arrows: Atom 49->196")
print("- Green arrows: Atom 83->196")
print(f"Image saved to: {os.path.expanduser('~/pymol_dipole_visualization.png')}")
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
    # Replace these with your actual paths
    csv_file = "/home/hp/results/MOUSE/level1/D148E/replica000/product_atom_positions_with_charges.csv"
    pdb_file = "/home/hp/results/MOUSE/level1/D148E/minim/minim.pdb"
    
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
