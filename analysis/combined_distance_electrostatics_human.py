from pymol import cmd, stored
import os
import csv

# ==============================================================================
# SETTINGS
# ==============================================================================
# Choose dataset: 'mouse' or 'human'
DATASET = 'mouse'  # Change to 'human' to use human data

# Active site residues with their chains
ACTIVE_SITE = [
    ('49', 'X'),   # SEC49 (or CYS49) in chain X
    ('83', 'X'),   # GLN83 in chain X
    ('198', 'Y'),  # Residue 198 in chain Y (corrected from 196)
]

# For backward compatibility, also keep a list of just residue numbers
ACTIVE_SITE_RESI = [res[0] for res in ACTIVE_SITE]

# Specific atoms for distance measurements
DISTANCE_ATOM_SELECTIONS = {
    '198': 'O2',  # Corrected from 196
    '49': 'CA',
    '83': 'CA',
}

# Mutation chain (where your mutations are)
MUTATION_CHAIN = 'X'

# H-bond parameters
HBOND_CUTOFF = 3.5  # Angstroms

# Output settings
OUTPUT_DIR = "D:/PhD_Thesis/GPX6/analysis"
OUTPUT_CSV = "mutation_distances_human.csv"
OUTPUT_HBOND_CSV = "mutation_hbonds_human.csv"

# ==============================================================================
# HEATMAP DATA - ΔG* BASED CLASSIFICATION
# ==============================================================================
# Averaged ΔG* values calculated from ALL mutation levels
# WT baseline: Mouse = 17.68 kcal/mol, Human = 16.40 kcal/mol

# MOUSE averaged ΔG* (from actual data)
MOUSE_AVG_DG = {
    '3': 17.4,    # S4R (corrected - this is residue 3 in mouse, not 4)
    '24': 15.7,   # L24I - BLUE
    '47': 21.6,   # S47A - RED (worst!)
    '48': 17.4,   # F48Y
    '49': 14.8,   # C49U - BLUE (active site mutation)
    '52': 16.2,   # T52A
    '54': 14.6,   # T54Q - BLUE (best!)
    '60': 18.1,   # T60A
    '74': 15.1,   # G74A - BLUE
    '87': 18.2,   # K87T
    '99': 24.3,   # R99C - RED (second worst!)
    '102': 15.8,  # G102S - BLUE
    '104': 16.9,  # Y104F
    '173': 20.1,  # N107S->H173R - RED (corrected residue number)
    '139': 15.7,  # F139L - BLUE
    '142': 16.3,  # P142S
    '143': 19.0,  # E143S
    '144': 17.7,  # H144Q
    '177': 18.7,  # H177Q
    '178': 14.8,  # T178A - BLUE
    '181': 16.4,  # R181S
}

# HUMAN averaged ΔG* (calculated from actual human data)
HUMAN_AVG_DG = {
    '3': 16.5,    # N3K
    '24': 18.6,   # L24I (NOTE: highly variable 13.7-23.6, avg ~18.6)
    '47': 15.9,   # A47S - BLUE
    '48': 17.6,   # Y48F
    '52': 17.3,   # A52T
    '54': 18.4,   # Q54T (NOTE: highly variable 12.1-28.9)
    '60': 19.9,   # A60T
    '74': 16.5,   # A74G
    '87': 15.6,   # T87K - BLUE (best!)
    '99': 16.7,   # C99R
    '102': 17.0,  # S102G
    '104': 17.7,  # F104Y (NOTE: variable 14.9-20.3)
    '139': 17.5,  # L139F
    '142': 17.6,  # S142P
    '143': 19.2,  # S143E
    '144': 18.0,  # Q144H
    '173': 15.0,  # H173R - BLUE (best in this position!)
    '177': 19.9,  # Q177H - RED
    '178': 16.6,  # A178T
    '181': 17.8,  # S181R
}

# Select which dataset to use
if DATASET == 'human':
    MUTATION_AVG_DG = HUMAN_AVG_DG
    WT_BASELINE = 16.40
    print("\nUsing HUMAN dataset (WT = 16.40 kcal/mol)")
else:
    MUTATION_AVG_DG = MOUSE_AVG_DG
    WT_BASELINE = 17.68
    print("\nUsing MOUSE dataset (WT = 17.68 kcal/mol)")


# Function to convert ΔG* to RGB color matching your heatmap
# Uses matplotlib's RdBu_r colormap (Red-Blue reversed)
def dg_to_color(dg_value, vmin=14.0, vmax=20.0):
    """Convert ΔG* value to RGB color matching heatmap"""
    # Normalize to 0-1 range
    norm_value = (dg_value - vmin) / (vmax - vmin)
    norm_value = max(0, min(1, norm_value))  # Clamp to [0, 1]
    
    # RdBu_r colormap (red=high, blue=low)
    # Approximate RGB values
    if norm_value > 0.75:  # Very high (RED)
        r, g, b = 0.8 + 0.2*((norm_value-0.75)/0.25), 0.1, 0.1
    elif norm_value > 0.6:  # High (light red/pink)
        t = (norm_value - 0.6) / 0.15
        r, g, b = 0.9, 0.5 - 0.4*t, 0.5 - 0.4*t
    elif norm_value > 0.4:  # Middle (white/light grey)
        t = (norm_value - 0.4) / 0.2
        r, g, b = 0.9 + 0.1*t, 0.9 - 0.4*t, 0.9 - 0.4*t
    elif norm_value > 0.25:  # Low (light blue)
        t = (norm_value - 0.25) / 0.15
        r, g, b = 0.5 - 0.4*t, 0.5 + 0.4*t, 0.9
    else:  # Very low (BLUE)
        t = norm_value / 0.25
        r, g, b = 0.1, 0.1 + 0.4*t, 0.8 + 0.2*(1-t)
    
    return [r, g, b]

# Assign colors to each mutation based on ΔG*
MUTATION_COLORS = {}
for resi, avg_dg in MUTATION_AVG_DG.items():
    rgb = dg_to_color(avg_dg)
    MUTATION_COLORS[resi] = rgb

# For backward compatibility, create the old dictionaries
HIGH_BARRIER_MOUSE = {}
LOW_BARRIER_MOUSE = {}
MIXED_NEUTRAL_MOUSE = {}

# But now use RGB colors instead
for resi in MUTATION_AVG_DG.keys():
    MIXED_NEUTRAL_MOUSE[resi] = MUTATION_COLORS[resi]  # All use continuous colors now

print("\nColor assignment based on averaged ΔG* (continuous colormap):")
for resi, dg in sorted(MUTATION_AVG_DG.items(), key=lambda x: x[1], reverse=True):
    color_type = "RED" if dg > 18.0 else ("BLUE" if dg < 16.0 else "WHITE/GREY")
    print("  Residue %s: ΔG* = %.1f kcal/mol → %s" % (resi, dg, color_type))

# Collect all mutation sites
ALL_MUTATIONS = (list(HIGH_BARRIER_MOUSE.keys()) + 
                 list(LOW_BARRIER_MOUSE.keys()) + 
                 list(MIXED_NEUTRAL_MOUSE.keys()))

print("\n" + "="*70)
print("COMPREHENSIVE MUTATION ANALYSIS WITH ELECTROSTATICS")
print("="*70)
print(f"Active site: {ACTIVE_SITE}")
print(f"Analyzing {len(ALL_MUTATIONS)} mutation sites")

# ==============================================================================
# GET RESIDUE INFORMATION
# ==============================================================================
stored.residues = []
cmd.iterate_state(1, "name CA", "stored.residues.append((resi, resn, chain, x, y, z))")

residue_names = {}
for resi, resn, chain, x, y, z in stored.residues:
    key = (resi, chain)
    residue_names[key] = resn

# Check if active site residues exist
print("\nChecking active site residues:")
for resi, chain in ACTIVE_SITE:
    key = (resi, chain)
    if key in residue_names:
        print("  Found residue %s chain %s (%s)" % (resi, chain, residue_names[key]))
    else:
        print("  WARNING: Residue %s chain %s has no CA atom (might be ligand/cofactor)" % (resi, chain))

# Calculate distances from active site center (using CA atoms that exist)
active_coords = []
for resi, chain in ACTIVE_SITE:
    coords = cmd.get_coords("resi %s and chain %s and name CA" % (resi, chain))
    if coords is not None and len(coords) > 0:
        active_coords.append(coords[0])

print("\nActive site residues for calculations: %s" % str(ACTIVE_SITE))

distances = {}
if len(active_coords) > 0:
    cx = sum([c[0] for c in active_coords]) / len(active_coords)
    cy = sum([c[1] for c in active_coords]) / len(active_coords)
    cz = sum([c[2] for c in active_coords]) / len(active_coords)
    
    for resi, resn, chain, x, y, z in stored.residues:
        dx = x - cx
        dy = y - cy
        dz = z - cz
        dist = (dx*dx + dy*dy + dz*dz) ** 0.5
        key = (resi, chain)
        distances[key] = dist

# ==============================================================================
# VISUALIZATION SETUP
# ==============================================================================
print("\nSetting up visualization...")

cmd.hide("everything", "all")
cmd.delete("dist_*")
cmd.delete("hbond_*")
cmd.bg_color("white")

# Protein backbone - transparent
cmd.show("cartoon", "all")
cmd.color("grey95", "all")
cmd.set("cartoon_transparency", 0.85)
cmd.set("cartoon_fancy_helices", 1)

# ==============================================================================
# SHOW ACTIVE SITE - YELLOW
# ==============================================================================
for resi, chain in ACTIVE_SITE:
    cmd.show("sticks", "resi %s and chain %s" % (resi, chain))
    cmd.show("spheres", "resi %s and chain %s and name CA" % (resi, chain))
    cmd.color("yellow", "resi %s and chain %s" % (resi, chain))
    cmd.set("sphere_scale", 1.3, "resi %s and chain %s and name CA" % (resi, chain))
    cmd.set("stick_radius", 0.25, "resi %s and chain %s" % (resi, chain))
    cmd.set("cartoon_transparency", 0.0, "resi %s and chain %s" % (resi, chain))

# ==============================================================================
# SHOW MUTATION SITES - COLOR BY ΔG* FROM MOUSE HEATMAP (CONTINUOUS)
# ==============================================================================
print("\nColoring mutations by ΔG* (matching heatmap colorscale)...")

for resi in ALL_MUTATIONS:
    if resi in MUTATION_COLORS:
        rgb_color = MUTATION_COLORS[resi]
        
        cmd.show("sticks", "resi %s and chain %s" % (resi, MUTATION_CHAIN))
        cmd.show("spheres", "resi %s and chain %s and name CA" % (resi, MUTATION_CHAIN))
        
        # Set RGB color
        cmd.set_color("color_%s" % resi, rgb_color)
        cmd.color("color_%s" % resi, "resi %s and chain %s" % (resi, MUTATION_CHAIN))
        
        # Size based on ΔG* (larger sphere = more extreme effect)
        dg_value = MUTATION_AVG_DG[resi]
        if dg_value > 18.0 or dg_value < 16.0:
            sphere_size = 0.9
            stick_size = 0.22
        else:
            sphere_size = 0.7
            stick_size = 0.18
            
        cmd.set("sphere_scale", sphere_size, "resi %s and chain %s and name CA" % (resi, MUTATION_CHAIN))
        cmd.set("stick_radius", stick_size, "resi %s and chain %s" % (resi, MUTATION_CHAIN))
        cmd.set("cartoon_transparency", 0.0, "resi %s and chain %s" % (resi, MUTATION_CHAIN))

# ==============================================================================
# DISTANCE MEASUREMENTS TO ACTIVE SITE
# ==============================================================================
print("\nCalculating distances to active site...")

distance_results = []
closest_distances = {}  # Track closest active site residue for each mutation

for mut_resi in ALL_MUTATIONS:
    key = (mut_resi, MUTATION_CHAIN)
    mut_resn = residue_names.get(key, "UNK")
    min_dist = 999
    closest_active = None
    
    # Calculate distances to ALL active site residues
    for active_resi, active_chain in ACTIVE_SITE:
        key_active = (active_resi, active_chain)
        active_resn = residue_names.get(key_active, "UNK")
        active_atom = DISTANCE_ATOM_SELECTIONS.get(active_resi, "CA")
        
        try:
            # First check if the residue exists at all
            residue_atom_count = cmd.count_atoms("resi %s and chain %s" % (active_resi, active_chain))
            if residue_atom_count == 0:
                print("  WARNING: Residue %s chain %s not found in structure - skipping" % (active_resi, active_chain))
                continue
            
            # Check if specific atom exists, fall back to CA if not
            atom_exists = cmd.count_atoms("resi %s and chain %s and name %s" % (active_resi, active_chain, active_atom))
            if atom_exists == 0:
                print("  Warning: Atom %s not found in residue %s chain %s, trying CA" % (active_atom, active_resi, active_chain))
                active_atom = "CA"
                atom_exists = cmd.count_atoms("resi %s and chain %s and name %s" % (active_resi, active_chain, active_atom))
                
                if atom_exists == 0:
                    # Try any heavy atom
                    print("  Warning: CA not found in residue %s chain %s, using first heavy atom" % (active_resi, active_chain))
                    model = cmd.get_model("resi %s and chain %s and not elem H" % (active_resi, active_chain))
                    if model.atom:
                        active_atom = model.atom[0].name
                        print("  Using atom %s from residue %s chain %s" % (active_atom, active_resi, active_chain))
                    else:
                        print("  ERROR: No atoms found in residue %s chain %s - skipping" % (active_resi, active_chain))
                        continue
            
            # Distance from mutation CA to active site specific atom
            dist_value = cmd.get_distance(
                "resi %s and chain %s and name CA" % (mut_resi, MUTATION_CHAIN),
                "resi %s and chain %s and name %s" % (active_resi, active_chain, active_atom))
            
            # Also get CA-to-CA distance if CA exists
            ca_atom_count = cmd.count_atoms("resi %s and chain %s and name CA" % (active_resi, active_chain))
            if ca_atom_count > 0 and active_atom != "CA":
                ca_dist = cmd.get_distance(
                    "resi %s and chain %s and name CA" % (mut_resi, MUTATION_CHAIN),
                    "resi %s and chain %s and name CA" % (active_resi, active_chain))
            else:
                ca_dist = dist_value
            
            # Store ALL distances for CSV (complete dataset)
            distance_results.append({
                'Mutation_Residue': mut_resi,
                'Mutation_ResName': mut_resn,
                'Mutation_Chain': MUTATION_CHAIN,
                'Active_Site_Residue': active_resi,
                'Active_Site_ResName': active_resn,
                'Active_Site_Chain': active_chain,
                'Active_Site_Atom': active_atom,
                'Distance_to_Atom_Angstrom': round(dist_value, 2),
                'CA_to_CA_Distance_Angstrom': round(ca_dist, 2)
            })
            
            # Track closest for visualization
            if dist_value < min_dist:
                min_dist = dist_value
                closest_active = (active_resi, active_chain, active_atom)
                
        except Exception as e:
            print("  ERROR: Could not measure %s to %s chain %s: %s" % (mut_resi, active_resi, active_chain, str(e)))
    
    # Draw line only to closest active site residue (for clean visualization)
    if closest_active:
        closest_distances[mut_resi] = (closest_active, min_dist)
        dist_name = "dist_%s_closest" % mut_resi
        active_resi_close, active_chain_close, active_atom_close = closest_active
        
        try:
            cmd.distance(dist_name,
                        "resi %s and chain %s and name CA" % (mut_resi, MUTATION_CHAIN),
                        "resi %s and chain %s and name %s" % (active_resi_close, active_chain_close, active_atom_close),
                        mode=0)
        except:
            pass

# Style distance lines - MORE VISIBLE
cmd.hide("labels", "dist_*")  # Hide distance values
cmd.set("dash_color", "grey40", "dist_*")
cmd.set("dash_width", 1.5, "dist_*")
cmd.set("dash_gap", 0.2, "dist_*")
cmd.set("dash_radius", 0.05, "dist_*")

print("Measured %d total distances (all mutations to all active site residues)" % len(distance_results))
print("Showing %d lines to closest active site residue for clean visualization" % len(closest_distances))

# ==============================================================================
# H-BOND DETECTION
# ==============================================================================
print("\nDetecting hydrogen bonds...")

# Select all measured residues
all_residues_list = [(resi, MUTATION_CHAIN) for resi in ALL_MUTATIONS] + list(ACTIVE_SITE)
all_res_sel = " or ".join(["(resi %s and chain %s)" % (r, c) for r, c in all_residues_list])
cmd.select("measured_residues", all_res_sel)

# Add hydrogens
cmd.h_add("measured_residues")

# Get all atoms from each residue
residue_atoms = {}
for res, chain in all_residues_list:
    model = cmd.get_model("resi %s and chain %s" % (res, chain))
    residue_atoms[(res, chain)] = {
        'resn': model.atom[0].resn if model.atom else 'UNK',
        'atoms': model.atom
    }

# Check for H-bonds
hbond_results = []
DONOR_PREFIXES = ['N', 'S']
ACCEPTOR_PREFIXES = ['O']

for i, res_chain1 in enumerate(all_residues_list):
    for res_chain2 in all_residues_list[i+1:]:
        atoms1 = residue_atoms[res_chain1]['atoms']
        atoms2 = residue_atoms[res_chain2]['atoms']
        
        for atom1 in atoms1:
            for atom2 in atoms2:
                # Check if polar atoms
                is_atom1_polar = (atom1.name[0] in DONOR_PREFIXES + ACCEPTOR_PREFIXES)
                is_atom2_polar = (atom2.name[0] in DONOR_PREFIXES + ACCEPTOR_PREFIXES)
                
                if not (is_atom1_polar and is_atom2_polar):
                    continue
                
                if atom1.name.startswith('H') or atom2.name.startswith('H'):
                    continue
                
                try:
                    sel1 = "resi %s and name %s" % (atom1.resi, atom1.name)
                    sel2 = "resi %s and name %s" % (atom2.resi, atom2.name)
                    dist = cmd.get_distance(sel1, sel2)
                    
                    if dist <= HBOND_CUTOFF:
                        # Determine donor vs acceptor
                        if atom1.name[0] in DONOR_PREFIXES:
                            donor, acceptor = atom1, atom2
                        elif atom2.name[0] in DONOR_PREFIXES:
                            donor, acceptor = atom2, atom1
                        else:
                            donor, acceptor = atom1, atom2
                        
                        hbond_results.append({
                            'Donor_Residue': donor.resi,
                            'Donor_ResName': donor.resn,
                            'Donor_Atom': donor.name,
                            'Acceptor_Residue': acceptor.resi,
                            'Acceptor_ResName': acceptor.resn,
                            'Acceptor_Atom': acceptor.name,
                            'Distance_Angstrom': round(dist, 2)
                        })
                except:
                    pass

# Remove duplicates
if hbond_results:
    seen = set()
    unique = []
    for hb in hbond_results:
        key = tuple(sorted([
            (hb['Donor_Residue'], hb['Donor_Atom']),
            (hb['Acceptor_Residue'], hb['Acceptor_Atom'])
        ]))
        if key not in seen:
            seen.add(key)
            unique.append(hb)
    hbond_results = unique

# Display H-bonds
if hbond_results:
    print(f"Found {len(hbond_results)} hydrogen bonds")
    
    for idx, hb in enumerate(hbond_results):
        hbond_name = "hbond_%d" % (idx+1)
        donor_sel = "resi %s and name %s" % (hb['Donor_Residue'], hb['Donor_Atom'])
        acceptor_sel = "resi %s and name %s" % (hb['Acceptor_Residue'], hb['Acceptor_Atom'])
        
        try:
            cmd.distance(hbond_name, donor_sel, acceptor_sel)
            cmd.set("dash_color", "red", hbond_name)
            cmd.set("dash_width", 4, hbond_name)
            cmd.set("dash_gap", 0.2, hbond_name)
            cmd.set("dash_length", 0.4, hbond_name)
        except:
            pass
else:
    print("No hydrogen bonds detected")

# Hide hydrogens from view
cmd.hide("everything", "measured_residues and elem H")

# ==============================================================================
# LABELS - BLACK SQUARES
# ==============================================================================
print("\nAdding labels...")

# Label active site
for resi, chain in ACTIVE_SITE:
    key = (resi, chain)
    if key in residue_names:
        resn = residue_names[key]
        cmd.label("resi %s and chain %s and name CA" % (resi, chain), "'%s%s'" % (resn, resi))

# Label mutations
for resi in ALL_MUTATIONS:
    key = (resi, MUTATION_CHAIN)
    if key in residue_names:
        resn = residue_names[key]
        cmd.label("resi %s and chain %s and name CA" % (resi, MUTATION_CHAIN), "'%s%s'" % (resn, resi))

cmd.set("label_size", 10)
cmd.set("label_font_id", 7)
cmd.set("label_color", "white")
cmd.set("label_bg_color", "black")
cmd.set("label_bg_transparency", 0.0)
cmd.set("label_outline_color", "white")
cmd.set("label_bg_outline", 1)
cmd.set("label_position", (2, 2, 2))  # Position labels away from structure

# ==============================================================================
# VIEW SETTINGS
# ==============================================================================
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("orthoscopic", 1)
cmd.set("depth_cue", 0)

cmd.select("all_sites_sel", 
           "(%s) or (%s)" % (
               " or ".join(["(resi %s and chain %s)" % (r, c) for r, c in ACTIVE_SITE]),
               " or ".join(["(resi %s and chain %s)" % (r, MUTATION_CHAIN) for r in ALL_MUTATIONS])
           ))
cmd.center("all_sites_sel")
cmd.zoom("all_sites_sel", buffer=8)
cmd.orient("all_sites_sel")
cmd.deselect()

# ==============================================================================
# SAVE DATA
# ==============================================================================
print("\nSaving data...")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Save distances
csv_file = os.path.join(OUTPUT_DIR, OUTPUT_CSV)
try:
    with open(csv_file, 'w', newline='') as f:
        if distance_results:
            writer = csv.DictWriter(f, fieldnames=distance_results[0].keys())
            writer.writeheader()
            writer.writerows(distance_results)
    print("Saved: %s" % csv_file)
except Exception as e:
    print("Error saving distances: %s" % e)

# Save H-bonds
if hbond_results:
    hbond_csv = os.path.join(OUTPUT_DIR, OUTPUT_HBOND_CSV)
    try:
        with open(hbond_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=hbond_results[0].keys())
            writer.writeheader()
            writer.writerows(hbond_results)
        print("Saved: %s" % hbond_csv)
    except Exception as e:
        print("Error saving H-bonds: %s" % e)

# ==============================================================================
# SUMMARY
# ==============================================================================
print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print("\nDataset: %s (WT baseline = %.2f kcal/mol)" % (DATASET.upper(), WT_BASELINE))
print("\nColor scheme (CONTINUOUS COLORMAP):")
print("  🟡 YELLOW = Active site (49, 83, 198)")
print("  🔴 RED = High ΔG* (>%.1f kcal/mol) - worse than WT" % (WT_BASELINE + 0.5))
print("  ⚪ WHITE = Around WT (~%.1f kcal/mol)" % WT_BASELINE)
print("  🔵 BLUE = Low ΔG* (<%.1f kcal/mol) - better than WT" % (WT_BASELINE - 1.5))
print("  Colors match your heatmap using RdBu_r colormap")

print("\nVisualization:")
print("  • Thin grey lines = distances to active site")
print("  • Thick red lines = hydrogen bonds (< 3.5 Å)")
print("  • Black labels = residue identifiers")
print("  • Color intensity = deviation from WT")

print("\nData files:")
print("  • %s - all distance measurements" % OUTPUT_CSV)
if hbond_results:
    print("  • %s - %d hydrogen bonds detected" % (OUTPUT_HBOND_CSV, len(hbond_results)))
else:
    print("  • No hydrogen bonds detected")

print("\nΔG* values and colors:")
for resi in sorted(MUTATION_AVG_DG.keys(), key=lambda x: MUTATION_AVG_DG[x], reverse=True):
    key = (resi, MUTATION_CHAIN)
    dg_value = MUTATION_AVG_DG[resi]
    if dg_value > 18.0:
        status = "RED - HIGH barrier"
    elif dg_value < 16.0:
        status = "BLUE - LOW barrier"
    else:
        status = "WHITE - near WT"
    
    if key in distances and key in residue_names:
        print("  %s%s: ΔG*=%.1f kcal/mol at %.1f Å (%s)" % 
              (residue_names[key], resi, dg_value, distances[key], status))

print("="*70 + "\n")