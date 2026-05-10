from pymol import cmd
import os
import math
import json

# Analyze the specific residue
ref_resi = "83"
ref_chain = "X"
cmd.select("ref_residue", f"resi {ref_resi} and chain {ref_chain}")
ref_atoms = cmd.get_model("ref_residue")

# Calculate center of mass for reference residue
ref_com = cmd.centerofmass(f"resi {ref_resi} and chain {ref_chain}")

# Get reference CA atom if it exists
cmd.select("ref_ca", f"resi {ref_resi} and chain {ref_chain} and name CA")
ref_ca_coords = None
if cmd.count_atoms("ref_ca") > 0:
    ref_ca_coords = cmd.get_model("ref_ca").atom[0].coord

# Visual selection - show the reference residue
cmd.show("sticks", f"resi {ref_resi} and chain {ref_chain}")
cmd.color("red", f"resi {ref_resi} and chain {ref_chain}")

# List of mutations to analyze
mutations = [
    {"mutation": "mousecys", "dg_star": 17.68, "dg_star_error": 0.91, "dg0": -56.96, "dg0_error": 1.03},
    {"mutation": "C49U", "dg_star": 14.84, "dg_star_error": 0.52, "dg0": -62.60, "dg0_error": 0.59},
    {"mutation": "F48Y", "dg_star": 17.38, "dg_star_error": 2.54, "dg0": -55.26, "dg0_error": 1.54},
    {"mutation": "S47A", "dg_star": 21.57, "dg_star_error": 2.46, "dg0": -56.18, "dg0_error": 3.49},
    {"mutation": "T52A", "dg_star": 16.16, "dg_star_error": 0.89, "dg0": -61.77, "dg0_error": 1.79},
    {"mutation": "R99C", "dg_star": 24.29, "dg_star_error": 2.79, "dg0": -55.41, "dg0_error": 3.13},
    {"mutation": "T54Q", "dg_star": 14.62, "dg_star_error": 0.68, "dg0": -55.68, "dg0_error": 3.78},
    {"mutation": "H144Q", "dg_star": 17.66, "dg_star_error": 2.18, "dg0": -53.40, "dg0_error": 3.89},
    {"mutation": "H177Q", "dg_star": 18.66, "dg_star_error": 2.09, "dg0": -55.40, "dg0_error": 3.89},
    {"mutation": "T178A", "dg_star": 14.77, "dg_star_error": 0.51, "dg0": -61.57, "dg0_error": 0.57},
    {"mutation": "E143S", "dg_star": 19.04, "dg_star_error": 1.40, "dg0": -57.10, "dg0_error": 1.56},
    {"mutation": "P142S", "dg_star": 16.30, "dg_star_error": 0.87, "dg0": -59.56, "dg0_error": 1.10},
    {"mutation": "F139L", "dg_star": 15.72, "dg_star_error": 0.48, "dg0": -54.00, "dg0_error": 1.27},
    {"mutation": "Y104F", "dg_star": 16.94, "dg_star_error": 1.28, "dg0": -58.41, "dg0_error": 1.58},
    {"mutation": "K87T", "dg_star": 18.21, "dg_star_error": 2.52, "dg0": -53.05, "dg0_error": 3.25},
    {"mutation": "G102S", "dg_star": 15.83, "dg_star_error": 1.30, "dg0": -59.61, "dg0_error": 1.52},
    {"mutation": "T60A", "dg_star": 18.09, "dg_star_error": 1.04, "dg0": -58.09, "dg0_error": 2.08},
    {"mutation": "R181S", "dg_star": 16.38, "dg_star_error": 0.97, "dg0": -60.55, "dg0_error": 1.26},
    {"mutation": "S4R", "dg_star": 17.37, "dg_star_error": 2.94, "dg0": -56.97, "dg0_error": 3.69},
    {"mutation": "N107S", "dg_star": 20.05, "dg_star_error": 2.12, "dg0": -56.09, "dg0_error": 3.33},
    {"mutation": "D148E", "dg_star": 15.68, "dg_star_error": 0.94, "dg0": -60.30, "dg0_error": 1.09},
]

# Special case for "mousecys" mutation - you may need to adjust this based on what it actually represents
mousecys_position = None  # You'll need to define this if mousecys refers to a specific residue

# Extract the residue numbers from mutation strings
mutation_residues = []
for mutation in mutations:
    mutation_name = mutation["mutation"]
    
    # Skip mousecys for now
    if mutation_name == "mousecys":
        continue
    
    # Extract the residue number - each mutation is in format like "F48Y"
    # This extracts all digits from the string
    residue_number = ''.join(filter(str.isdigit, mutation_name))
    if residue_number:
        mutation_residues.append({
            "residue": residue_number,
            "mutation": mutation_name,
            "dg_star": mutation["dg_star"],
            "dg0": mutation["dg0"]
        })

# Calculate and output distances
current_dir = os.getcwd()
distance_file_path = os.path.join(current_dir, "mutation_distances.txt")
with open(distance_file_path, "w") as distance_file:
    distance_file.write(f"Distances from GLN{ref_resi} (Chain {ref_chain}) to mutation residues:\n")
    distance_file.write("-" * 90 + "\n")
    distance_file.write(f"{'Mutation':>10} | {'Residue':>8} | {'Chain':>6} | {'Residue Name':>12} | {'COM Distance (Å)':>15} | {'CA Distance (Å)':>15} | {'ΔG*':>8} | {'ΔG0':>8}\n")
    distance_file.write("-" * 90 + "\n")
    
    # Store results for JSON output
    json_results = []
    
    # Set for highlighting mutated residues in visualization
    cmd.select("mutated_residues", "none")
    
    for mutation_info in mutation_residues:
        residue = mutation_info["residue"]
        mutation_name = mutation_info["mutation"]
        dg_star = mutation_info["dg_star"]
        dg0 = mutation_info["dg0"]
        
        # Check all chains for this residue
        for chain in cmd.get_chains():
            cmd.select("check_res", f"resi {residue} and chain {chain}")
            if cmd.count_atoms("check_res") > 0:
                # Get residue name
                res_model = cmd.get_model("check_res")
                res_name = res_model.atom[0].resn if res_model.atom else "Unknown"
                
                # Calculate center of mass
                res_com = cmd.centerofmass(f"resi {residue} and chain {chain}")
                
                # Calculate distance between centers of mass
                com_dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(ref_com, res_com)))
                
                # Calculate CA-CA distance if possible
                ca_dist = "N/A"
                cmd.select("res_ca", f"resi {residue} and chain {chain} and name CA")
                if cmd.count_atoms("res_ca") > 0 and ref_ca_coords is not None:
                    res_ca = cmd.get_model("res_ca").atom[0].coord
                    ca_dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(ref_ca_coords, res_ca)))
                    ca_dist = f"{ca_dist:.2f}"
                
                # Add to visualization selection
                cmd.select("mutated_residues", f"mutated_residues or (resi {residue} and chain {chain})")
                
                # Write to file
                distance_file.write(f"{mutation_name:>10} | {residue:>8} | {chain:>6} | {res_name:>12} | {com_dist:>15.2f} | {ca_dist:>15} | {dg_star:>8.2f} | {dg0:>8.2f}\n")
                
                # Add to JSON results
                json_results.append({
                    "mutation": mutation_name,
                    "residue_number": residue,
                    "chain": chain,
                    "residue_name": res_name,
                    "com_distance": round(com_dist, 2),
                    "ca_distance": ca_dist if ca_dist != "N/A" else None,
                    "dg_star": dg_star,
                    "dg0": dg0
                })

# Show mutated residues in visualization
cmd.show("sticks", "mutated_residues")
cmd.color("yellow", "mutated_residues")

# Add a line from reference residue to each mutated residue
for i, mutation_info in enumerate(mutation_residues):
    residue = mutation_info["residue"]
    for chain in cmd.get_chains():
        cmd.select("check_res", f"resi {residue} and chain {chain}")
        if cmd.count_atoms("check_res") > 0:
            # Draw distance measurement
            cmd.distance(f"dist_{i}", 
                         f"resi {ref_resi} and chain {ref_chain} and name CA", 
                         f"resi {residue} and chain {chain} and name CA")

# Save as JSON for further analysis
json_file_path = os.path.join(current_dir, "mutation_distances.json")
with open(json_file_path, "w") as json_file:
    json.dump(json_results, json_file, indent=2)

print(f"\nDistances from GLN{ref_resi} to mutation residues saved to {distance_file_path}")
print(f"Distance data also saved as JSON to {json_file_path}")
print("\nVisualization created:")
print("- Reference residue GLN83 shown in red")
print("- Mutated residues shown in yellow")
print("- Distance measurements displayed as dashed lines")