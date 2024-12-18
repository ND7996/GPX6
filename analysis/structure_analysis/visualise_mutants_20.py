from pymol import cmd

# Path to the PDB file
pdb_path = "/home/hp/results/mousecys/N22L/minim/minim.pdb"  # Replace with actual filename

# Residues to select (only the 20 specified mutant positions)
residues = [
    (None, 4), (None, 24), (None, 47), (None, 48), (None, 52), (None, 54), 
    (None, 60), (None, 74), (None, 87), (None, 99), (None, 102), (None, 104),
    (None, 107), (None, 139), (None, 142), (None, 143), (None, 144), (None, 177),
    (None, 178), (None, 181)
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
cmd.color("red", "mutated_residues")

# Set the background to white for better visibility
cmd.bg_color("white")

# Save the visualization as a PNG image
image_path = "/home/hp/nayanika/github/GPX6/figures/visualise_20_mutants.png"
cmd.png(image_path, width=1920, height=1080, dpi=300)
print(f"Visualization saved as image: {image_path}")
