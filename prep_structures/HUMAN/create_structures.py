#!/usr/bin/env python3
"""
PyMOL script to generate mutant PDB files for GPX6 human protein
Creates progressive mutations from level 0 to level 19
Uses direct residue replacement method
"""

import pymol
from pymol import cmd
import os

# Initialize PyMOL
pymol.finish_launching(['pymol', '-c'])

# Define the base structure path
base_pdb = "/home/hp/nayanika/github/GPX6/humanWT/humanWT.pdb"

# Define mutations for each level (progressive accumulation)
mutations_by_level = {
     # humansec - base structure
    1: ["T87K"],
    2: ["T87K", "A47S"],
    3: ["T87K", "A47S", "S143E"],
    4: ["T87K", "A47S", "S143E", "A60T"],
    5: ["T87K", "A47S", "S143E", "A60T", "F104Y"],
    6: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P"],
    7: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F"],
    8: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R"],
    9: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T"],
    10: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F"],
    11: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F", "Q144H"],
    12: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F", "Q144H", "Q54T"],
    13: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F", "Q144H", "Q54T", "Q177H"],
    14: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F", "Q144H", "Q54T", "Q177H", "S102G"],
    15: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F", "Q144H", "Q54T", "Q177H", "S102G", "L24I"],
    16: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F", "Q144H", "Q54T", "Q177H", "S102G", "L24I", "N3K"],
    17: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F", "Q144H", "Q54T", "Q177H", "S102G", "L24I", "N3K", "H173R"],
    18: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F", "Q144H", "Q54T", "Q177H", "S102G", "L24I", "N3K", "H173R", "A178T"],
    19: ["T87K", "A47S", "S143E", "A60T", "F104Y", "S142P", "L139F", "S181R", "A52T", "Y48F", "Q144H", "Q54T", "Q177H", "S102G", "L24I", "N3K", "H173R", "A178T", "A74G"]
}

# Amino acid mapping
aa_mapping = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
    'U': 'CYS'  # TEMPORARY: Use CYS first, then manually change to SEC
}

def apply_mutation_simple(mutation_str, object_name="protein"):
    """
    Apply mutation using simple approach
    """
    if len(mutation_str) < 3:
        print(f"Invalid mutation format: {mutation_str}")
        return False
        
    original_aa = mutation_str[0]
    new_aa = mutation_str[-1]
    position = mutation_str[1:-1]
    
    try:
        # Check if residue exists
        selection = f"resi {position} and {object_name}"
        if cmd.count_atoms(selection) == 0:
            print(f"Warning: Residue {position} not found")
            return False
        
        # Get 3-letter code for new amino acid
        new_resn = aa_mapping.get(new_aa, new_aa)
        
        # Use mutagenesis wizard
        cmd.wizard("mutagenesis")
        cmd.refresh_wizard()
        cmd.get_wizard().do_select(f"({object_name} and resi {position})")
        cmd.get_wizard().set_mode(new_resn)
        cmd.get_wizard().apply()
        cmd.set_wizard()
        
        print(f"Applied mutation: {original_aa}{position}{new_aa}")
        return True
        
    except Exception as e:
        print(f"Error applying mutation {mutation_str}: {str(e)}")
        return False

def create_mutant_pdb_simple(level, mutations, output_dir="mutant_pdbs"):
    """
    Create mutant PDB using simple approach
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Load fresh structure
    cmd.reinitialize()
    cmd.load(base_pdb, "protein")
    
    print(f"\n=== Creating Level {level} mutant ===")
    print(f"Mutations to apply: {', '.join(mutations)}")
    
    successful_mutations = []
    failed_mutations = []
    
    # Apply mutations
    for mutation in mutations:
        success = apply_mutation_simple(mutation, "protein")
        if success:
            successful_mutations.append(mutation)
        else:
            failed_mutations.append(mutation)
    
    # Basic cleanup
    cmd.remove("solvent")
    
    # Save file
    output_filename = f"GPX6_level{level:02d}.pdb"
    output_path = os.path.join(output_dir, output_filename)
    
    try:
        cmd.save(output_path, "protein")
        print(f"✓ Saved: {output_path}")
        if successful_mutations:
            print(f"  Successful: {', '.join(successful_mutations)}")
        if failed_mutations:
            print(f"  Failed: {', '.join(failed_mutations)}")
    except Exception as e:
        print(f"✗ Error saving: {str(e)}")
    
    # Clean up
    cmd.delete("protein")

def create_pml_script(output_dir="mutant_pdbs"):
    """
    Create a PyMOL script file (.pml) for manual execution
    """
    os.makedirs(output_dir, exist_ok=True)
    pml_content = f"""# PyMOL script for GPX6 human mutations
# Load this script in PyMOL: File -> Run Script -> select this file

load {base_pdb}, original

"""
    
    for level in range(0, 20):  # Levels 0 to 19
        if level in mutations_by_level:
            mutations = mutations_by_level[level]
            pml_content += f"""# Level {level}: {', '.join(mutations)}
copy level{level}, original
"""
            
            for mutation in mutations:
                if len(mutation) >= 3:
                    position = mutation[1:-1]
                    new_aa = mutation[-1]
                    new_resn = aa_mapping.get(new_aa, new_aa)
                    pml_content += f"alter resi {position} and level{level}, resn='{new_resn}'\n"
            
            pml_content += f"remove solvent\n"
            pml_content += f"save {output_dir}/GPX6_level{level:02d}.pdb, level{level}\n"
            pml_content += f"delete level{level}\n\n"
    
    pml_path = os.path.join(output_dir, "generate_mutations.pml")
    with open(pml_path, 'w') as f:
        f.write(pml_content)
    
    print(f"Created PyMOL script: {pml_path}")
    print("To use: Open PyMOL -> File -> Run Script -> select this .pml file")

def main():
    """
    Main function with multiple approaches
    """
    print("GPX6 Human Mutation Generator")
    print("=" * 50)
    print(f"Base structure: {base_pdb}")
    
    if not os.path.exists(base_pdb):
        print(f"ERROR: Base PDB file not found!")
        return
    
    # Method 1: Create PML script for manual execution
    print("\n1. Creating PyMOL script file...")
    create_pml_script()
    
    # Method 2: Try direct Python execution
    print("\n2. Attempting direct mutation...")
    try:
        for level in range(0, 20):  # Levels 0 to 19
            if level in mutations_by_level:
                create_mutant_pdb_simple(level, mutations_by_level[level])
    except Exception as e:
        print(f"Direct method failed: {str(e)}")
    
    print("\n" + "=" * 50)
    print("COMPLETE!")
    print("If Python method failed, use the .pml script in PyMOL")

if __name__ == "__main__":
    main()