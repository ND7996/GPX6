from pymol import cmd

# Path to the PDB file
pdb_path = "/home/hp/results/C49U/T60A/minim/minim.pdb"  # Replace with actual filename

# Residues to select (include chain if needed, e.g., ("A", 3))
residues = [
    (None, 3), (None, 4), (None, 16), (None, 22), (None, 24), (None, 25),
    (None, 27), (None, 29), (None, 30), (None, 31), (None, 33), (None, 35),
    (None, 40), (None, 47), (None, 48), (None, 52), (None, 54), (None, 60),
    (None, 67), (None, 69), (None, 71), (None, 74), (None, 87), (None, 99),
    (None, 102), (None, 104), (None, 107), (None, 119), (None, 120), (None, 126),
    (None, 137), (None, 139), (None, 142), (None, 143), (None, 144), (None, 148),
    (None, 173), (None, 177), (None, 178), (None, 181), (None, 182), (None, 184),
    (None, 188), (None, 192), (None, 194), (None, 195)
]

# Helper function to build the selection string
def build_selection(residues):
    selections = []
    for chain, resi in residues:
        if chain:
            selections.append(f"(chain {chain} and resi {resi})")
        else:
            selections.append(f"resi {resi}")
    return " or ".join(selections)

# Build the selection string
selection = build_selection(residues)

# Load the PDB file
cmd.load(pdb_path, "mousecys")

# Print feedback if file fails to load
if not cmd.get_object_list("mousecys"):
    print("Error: PDB file failed to load. Check the file path.")

# Select residues
cmd.select("mutated_residues", selection)

# Ensure residues are selected
if cmd.count_atoms("mutated_residues") == 0:
    print("Error: No atoms selected. Check residue numbers and chain identifiers.")

# Show selected residues
cmd.show("sticks", "mutated_residues")

# Optional: Color the selected residues
cmd.color("brown", "mutated_residues")

# Set the background to white for better visibility
cmd.bg_color("white")

# Save the visualization as a PNG image
image_path = "/home/hp/nayanika/github/GPX6/figures/visualise_all_mutants.png"
cmd.png(image_path, width=1920, height=1080, dpi=300)
print(f"Visualization saved as image: {image_path}")
