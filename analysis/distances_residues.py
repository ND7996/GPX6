"""
PyMOL Distance & H-Bond Measurement - ALL Possible H-Bonding Atoms
Comprehensive search of all donor/acceptor combinations
"""

from pymol import cmd
import os
import csv

# ============== SETTINGS ==============
TARGET_RESIDUE = 102
REFERENCE_RESIDUES = [49, 83, 196]
CHAIN = ""  # Leave empty for all chains, or specify like "A"

# Specific atom selections for DISTANCE measurements
DISTANCE_ATOM_SELECTIONS = {
    196: "O2",
    49: "CA",
    83: "CA",
    102: "CA"
}

# H-bond parameters
HBOND_CUTOFF = 3.5  # Angstroms

# Output settings
OUTPUT_DIR = "D:/PhD_Thesis/GPX6/analysis"
OUTPUT_CSV = "distance_measurements.csv"
OUTPUT_HBOND_CSV = "hydrogen_bonds.csv"
# ======================================

print("\n" + "="*60)
print("Starting comprehensive H-bond analysis...")
print("="*60)

# Build selections
if CHAIN:
    target_sel = f"resi {TARGET_RESIDUE} and chain {CHAIN}"
else:
    target_sel = f"resi {TARGET_RESIDUE}"

target_atom = DISTANCE_ATOM_SELECTIONS.get(TARGET_RESIDUE, "CA")

cmd.select("target_residue", target_sel)
cmd.show("sticks", "target_residue")
cmd.color("yellow", "target_residue")

ref_colors = ["cyan", "magenta", "orange"]
distance_results = []

# ============== DISTANCE MEASUREMENTS ==============

for idx, ref_res in enumerate(REFERENCE_RESIDUES):
    if CHAIN:
        ref_sel = f"resi {ref_res} and chain {CHAIN}"
    else:
        ref_sel = f"resi {ref_res}"
    
    cmd.show("sticks", ref_sel)
    cmd.color(ref_colors[idx % len(ref_colors)], ref_sel)
    
    ref_atom = DISTANCE_ATOM_SELECTIONS.get(ref_res, "CA")
    dist_name = f"dist_{TARGET_RESIDUE}_{target_atom}_to_{ref_res}_{ref_atom}"
    
    cmd.distance(dist_name,
                f"({target_sel}) and name {target_atom}",
                f"({ref_sel}) and name {ref_atom}")
    
    try:
        dist_value = cmd.get_distance(f"({target_sel}) and name {target_atom}",
                                     f"({ref_sel}) and name {ref_atom}")
        
        target_resn = "UNK"
        ref_resn = "UNK"
        
        if cmd.count_atoms(f"{target_sel} and name {target_atom}") > 0:
            target_resn = cmd.get_model(f"{target_sel} and name {target_atom}").atom[0].resn
        
        if cmd.count_atoms(f"{ref_sel} and name {ref_atom}") > 0:
            ref_resn = cmd.get_model(f"{ref_sel} and name {ref_atom}").atom[0].resn
        
        distance_results.append({
            'From_Residue': TARGET_RESIDUE,
            'From_Residue_Name': target_resn,
            'From_Atom': target_atom,
            'To_Residue': ref_res,
            'To_Residue_Name': ref_resn,
            'To_Atom': ref_atom,
            'Distance_Angstrom': round(dist_value, 2)
        })
        
        cmd.set("dash_color", ref_colors[idx % len(ref_colors)], dist_name)
        
    except Exception as e:
        print(f"Warning: {e}")

print("\n" + "="*60)
print("Distance Measurements")
print("="*60)
for result in distance_results:
    print(f"{result['From_Residue_Name']}{result['From_Residue']}:{result['From_Atom']} → " + 
          f"{result['To_Residue_Name']}{result['To_Residue']}:{result['To_Atom']}: " +
          f"{result['Distance_Angstrom']:.2f} Å")
print("="*60)

# ============== COMPREHENSIVE H-BOND DETECTION ==============

print("\n" + "="*60)
print("Comprehensive H-bond search on ALL polar atoms...")
print("="*60)

all_residues = [TARGET_RESIDUE] + REFERENCE_RESIDUES
if CHAIN:
    all_res_sel = " or ".join([f"resi {r} and chain {CHAIN}" for r in all_residues])
else:
    all_res_sel = " or ".join([f"resi {r}" for r in all_residues])

cmd.select("measured_residues", all_res_sel)

# Add hydrogens for accurate detection
cmd.h_add("measured_residues")

# Get all atoms from each residue
residue_atoms = {}
for res in all_residues:
    if CHAIN:
        sel = f"resi {res} and chain {CHAIN}"
    else:
        sel = f"resi {res}"
    
    model = cmd.get_model(sel)
    residue_atoms[res] = {
        'resn': model.atom[0].resn if model.atom else 'UNK',
        'atoms': model.atom
    }
    
    print(f"Residue {res} ({residue_atoms[res]['resn']}): {len(model.atom)} atoms")

print(f"\nSearching for H-bonds (cutoff: {HBOND_CUTOFF} Å)...")

hbond_results = []
hbond_count = 0

# ALL possible H-bond donor/acceptor atoms
# Donors: N, NH, NH2, NH3, NE, NE1, NE2, ND1, ND2, NZ, SG, OH
# Acceptors: O, OE1, OE2, OD1, OD2, OG, OG1, OH, O1, O2, O3, OXT

DONOR_PREFIXES = ['N', 'S']  # N and S atoms can be donors
ACCEPTOR_PREFIXES = ['O']     # O atoms are acceptors

# Check all pairs of residues
for i, res1 in enumerate(all_residues):
    for res2 in all_residues[i+1:]:
        
        print(f"\nChecking {residue_atoms[res1]['resn']}{res1} ↔ {residue_atoms[res2]['resn']}{res2}...")
        
        atoms1 = residue_atoms[res1]['atoms']
        atoms2 = residue_atoms[res2]['atoms']
        
        # Check all atom pairs
        for atom1 in atoms1:
            for atom2 in atoms2:
                
                # Skip if neither is a potential H-bond atom
                is_atom1_polar = (atom1.name[0] in DONOR_PREFIXES + ACCEPTOR_PREFIXES)
                is_atom2_polar = (atom2.name[0] in DONOR_PREFIXES + ACCEPTOR_PREFIXES)
                
                if not (is_atom1_polar and is_atom2_polar):
                    continue
                
                # Skip hydrogen atoms themselves
                if atom1.name.startswith('H') or atom2.name.startswith('H'):
                    continue
                
                # Build selections
                if CHAIN:
                    sel1 = f"resi {atom1.resi} and chain {CHAIN} and name {atom1.name}"
                    sel2 = f"resi {atom2.resi} and chain {CHAIN} and name {atom2.name}"
                else:
                    sel1 = f"resi {atom1.resi} and name {atom1.name}"
                    sel2 = f"resi {atom2.resi} and name {atom2.name}"
                
                # Measure distance
                try:
                    dist = cmd.get_distance(sel1, sel2)
                    
                    if dist <= HBOND_CUTOFF:
                        # Determine donor vs acceptor
                        if atom1.name[0] == 'N' or atom1.name[0] == 'S':
                            donor, acceptor = atom1, atom2
                        elif atom2.name[0] == 'N' or atom2.name[0] == 'S':
                            donor, acceptor = atom2, atom1
                        elif atom1.name[0] == 'O' and atom2.name[0] == 'O':
                            # Both oxygens - could be either way
                            donor, acceptor = atom1, atom2
                        else:
                            donor, acceptor = atom1, atom2
                        
                        hbond_results.append({
                            'Donor_Residue': donor.resi,
                            'Donor_Residue_Name': donor.resn,
                            'Donor_Atom': donor.name,
                            'Donor_Chain': donor.chain if donor.chain else '',
                            'Acceptor_Residue': acceptor.resi,
                            'Acceptor_Residue_Name': acceptor.resn,
                            'Acceptor_Atom': acceptor.name,
                            'Acceptor_Chain': acceptor.chain if acceptor.chain else '',
                            'Distance_Angstrom': round(dist, 2),
                            'Bond_Type': 'H-Bond/Polar'
                        })
                        
                        print(f"  FOUND: {donor.resn}{donor.resi}:{donor.name} - {acceptor.resn}{acceptor.resi}:{acceptor.name} = {dist:.2f} Å")
                        
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

# Display H-bonds visually
if hbond_results:
    print(f"\n{'='*60}")
    print(f"FOUND {len(hbond_results)} POTENTIAL H-BONDS!")
    print(f"{'='*60}")
    
    # Create visual H-bond objects
    for idx, hb in enumerate(hbond_results):
        hbond_name = f"hbond_{idx+1}"
        
        donor_sel = f"resi {hb['Donor_Residue']} and name {hb['Donor_Atom']}"
        acceptor_sel = f"resi {hb['Acceptor_Residue']} and name {hb['Acceptor_Atom']}"
        
        if CHAIN:
            donor_sel += f" and chain {CHAIN}"
            acceptor_sel += f" and chain {CHAIN}"
        
        try:
            cmd.distance(hbond_name, donor_sel, acceptor_sel)
            cmd.set("dash_color", "red", hbond_name)
            cmd.set("dash_width", 6, hbond_name)
            cmd.set("dash_gap", 0.2, hbond_name)
            cmd.set("dash_length", 0.4, hbond_name)
        except:
            pass
    
    print("H-bonds displayed as THICK RED dashed lines!")
else:
    print("\n" + "="*60)
    print("NO H-BONDS FOUND")
    print("This means no N, O, or S atoms are within 3.5 Å between these residues")
    print("="*60)

print("\n" + "="*60)
print("Hydrogen Bonds Detected")
print("="*60)
if hbond_results:
    for hb in hbond_results:
        print(f"{hb['Donor_Residue_Name']}{hb['Donor_Residue']}:{hb['Donor_Atom']} ↔ " +
              f"{hb['Acceptor_Residue_Name']}{hb['Acceptor_Residue']}:{hb['Acceptor_Atom']}: " +
              f"{hb['Distance_Angstrom']:.2f} Å")
else:
    print("None detected")
print("="*60)

# ============== FORMATTING ==============

cmd.set("label_size", 16)
cmd.set("label_color", "black")
cmd.set("label_bg_color", "white")
cmd.set("label_bg_transparency", 0.3)
cmd.set("dash_width", 3)
cmd.set("dash_gap", 0)

# Labels
if cmd.count_atoms(f"{target_sel} and name CA") > 0:
    target_resn = cmd.get_model(f"{target_sel} and name CA").atom[0].resn
    cmd.label(f"{target_sel} and name CA", f'"{target_resn} {TARGET_RESIDUE}"')

for ref_res in REFERENCE_RESIDUES:
    if CHAIN:
        ref_sel = f"resi {ref_res} and chain {CHAIN}"
    else:
        ref_sel = f"resi {ref_res}"
    
    label_atom = DISTANCE_ATOM_SELECTIONS.get(ref_res, "CA")
    label_sel = f"{ref_sel} and name {label_atom}"
    
    if cmd.count_atoms(label_sel) > 0:
        try:
            resn = cmd.get_model(label_sel).atom[0].resn
            cmd.label(label_sel, f'"{resn} {ref_res}"')
        except:
            pass

cmd.set("label_outline_color", "black")
cmd.set("label_font_id", 7)
cmd.set("label_position", (0, 0, 2))

cmd.bg_color("white")
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("orthoscopic", 1)

if CHAIN:
    zoom_sel = " or ".join([f"resi {r} and chain {CHAIN}" for r in all_residues])
else:
    zoom_sel = " or ".join([f"resi {r}" for r in all_residues])

cmd.zoom(zoom_sel, buffer=8)
cmd.orient(zoom_sel)

# Hide hydrogens from view
cmd.hide("everything", "measured_residues and elem H")

# ============== SAVE CSV ==============

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

csv_file = os.path.join(OUTPUT_DIR, OUTPUT_CSV)
try:
    with open(csv_file, 'w', newline='') as f:
        if distance_results:
            writer = csv.DictWriter(f, fieldnames=distance_results[0].keys())
            writer.writeheader()
            writer.writerows(distance_results)
    print(f"\nSaved: {csv_file}")
except Exception as e:
    print(f"Error: {e}")

if hbond_results:
    hbond_csv_file = os.path.join(OUTPUT_DIR, OUTPUT_HBOND_CSV)
    try:
        with open(hbond_csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=hbond_results[0].keys())
            writer.writeheader()
            writer.writerows(hbond_results)
        print(f"Saved: {hbond_csv_file}")
    except Exception as e:
        print(f"Error: {e}")

print("\n" + "="*60)
print("COMPLETED!")
print("="*60)
print(f"Files saved to: {OUTPUT_DIR}")
print(f"  - {OUTPUT_CSV} (Distances)")
if hbond_results:
    print(f"  - {OUTPUT_HBOND_CSV} ({len(hbond_results)} H-bonds)")
    print(f"\nH-bonds shown as THICK RED dashed lines (width=6)")
else:
    print("  - No H-bonds detected between these residues")
print("="*60 + "\n")