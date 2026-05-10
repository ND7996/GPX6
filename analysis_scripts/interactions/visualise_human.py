from pymol import cmd
import os

# Create the results directory if it doesn't exist
results_dir = "/home/hp/nayanika/github/GPX6/results/"
os.makedirs(results_dir, exist_ok=True)

# Load the structure
cmd.load("/home/hp/nayanika/github/GPX6/prep_structures/original_humansec.pdb", "humansec")

# List of variant residue numbers
variant_residues = [3, 24, 47, 48, 49, 52, 54, 60, 74, 87, 99, 102, 104, 173, 139, 142, 143, 144, 177, 178, 181]

# Create selection for all variant residues
variant_sel = "resi " + "+".join(str(r) for r in variant_residues)
cmd.select("variants", variant_sel)

# Show structure and highlight variants
cmd.show("cartoon", "humansec")
cmd.color("gray", "humansec")
cmd.color("orange", "variants")
cmd.show("sticks", "variants")

# Set up hydrogen bond detection
cmd.set("h_bond_cutoff_center", 3.6)
cmd.set("h_bond_cutoff_edge", 3.2)

# Show hydrogen bonds
cmd.delete("hbonds")
cmd.distance("hbonds", "variants", "all", mode=2, cutoff=3.6)
cmd.set("dash_color", "blue", "hbonds")
cmd.set("dash_width", 2.0)
cmd.set("dash_gap", 0.2)
cmd.hide("labels", "hbonds")

# Save the PyMOL session and image with updated filenames for human secretory structure
cmd.save(os.path.join(results_dir, "humansec_variants_hbonds.pse"))
cmd.bg_color("white")
cmd.zoom("variants", buffer=6)
cmd.ray(1200, 1200)
cmd.png(os.path.join(results_dir, "humansec_variants_hbonds.png"), dpi=300)

# Get hydrogen bond information using a different approach
hbond_file = os.path.join(results_dir, "humansec_hbond_info.txt")
with open(hbond_file, "w") as f:
    f.write("Hydrogen Bonds Involving Variant Residues in Human Secretory GPX6:\n")
    f.write("==================================================================\n")
    
    # Use PyMOL's built-in hbond function
    cmd.select("variants_hbonds", f"byres ({variant_sel}) within 5 of {variant_sel}")
    hbonds = cmd.find_pairs(f"{variant_sel} and (elem N or elem O)", 
                           f"variants_hbonds and (elem N or elem O)", 
                           mode=1, cutoff=3.6, angle=55)
    
    if hbonds:
        for pair in hbonds:
            atom1 = cmd.get_model(pair[0]).atom[0]
            atom2 = cmd.get_model(pair[1]).atom[0]
            distance = cmd.get_distance(pair[0], pair[1])
            f.write(f"H-Bond: {atom1.resn} {atom1.resi} {atom1.name} (Chain {atom1.chain}) → ")
            f.write(f"{atom2.resn} {atom2.resi} {atom2.name} (Chain {atom2.chain}) ")
            f.write(f"(Distance: {distance:.2f} Å)\n")
    else:
        f.write("No hydrogen bonds detected using criteria: cutoff=3.6Å, angle=55°\n")
        f.write("Hydrogen bonds are visualized in the saved PyMOL session.\n")

print("Analysis complete. Files saved to", results_dir)
print("Filenames: humansec_variants_hbonds.pse, humansec_variants_hbonds.png, humansec_hbond_info.txt")