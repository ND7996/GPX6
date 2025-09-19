#!/usr/bin/env python3
"""
PyMOL script to generate mutant PDB files for GPX6 mouse protein
Creates progressive mutations from level 1 to level 20
Uses direct residue replacement method
"""

import pymol
from pymol import cmd
import os

# Initialize PyMOL
pymol.finish_launching(['pymol', '-c'])

# Define the base structure path
base_pdb = "/home/hp/nayanika/github/GPX6/mouseWT/GPX6mousecys.pdb"

# Define mutations for each level (progressive accumulation)
mutations_by_level = {
    1: ["T54Q"],
    2: ["T54Q", "C49U"],
    3: ["T54Q", "C49U", "I24L"],
    4: ["T54Q", "I24L", "F137Y"],
    5: ["T54Q", "I24L", "F137Y", "S47A"],
    6: ["T54Q", "I24L", "F137Y", "S47A", "T50A"],
    7: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A"],
    8: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q"],
    9: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C"],
    10: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q"],
    11: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A"],
    12: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A", "Y104F"],
    13: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A", "Y104F", "S5R"],
    14: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A", "Y104F", "S5R", "T52A"],
    15: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A", "Y104F", "S5R", "T52A", "R87T"],
    16: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A", "Y104F", "S5R", "T52A", "R87T", "Q102S"],
    17: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A", "Y104F", "S5R", "T52A", "R87T", "Q102S", "N107S"],
    18: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A", "Y104F", "S5R", "T52A", "R87T", "Q102S", "N107S", "P142S"],
    19: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A", "Y104F", "S5R", "T52A", "R87T", "Q102S", "N107S", "P142S", "R181S"],
    20: ["T54Q", "I24L", "F137Y", "S47A", "T50A", "G74A", "H144Q", "R39C", "H177Q", "T179A", "Y104F", "S5R", "T52A", "R87T", "Q102S", "N107S", "P142S", "R181S", "F80Y"]
}

# Amino acid mapping
aa_mapping = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
    'U': 'SEC'  # Selenocysteine
}

def apply_mutation_direct(mutation_str, object_name="protein"):
    """
    Apply mutation using direct residue name change
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
        
        # Method 1: Direct residue name alteration
        cmd.alter(f"resi {position} and {object_name}", f"resn='{new_resn}'")
        
        print(f"Applied mutation: {original_aa}{position}{new_aa} ({new_resn})")
        return True
        
    except Exception as e:
        print(f"Error applying mutation {mutation_str}: {str(e)}")
        return False

def apply_mutation_wizard(mutation_str, object_name="protein"):
    """
    Apply mutation using PyMOL wizard method (alternative approach)
    """
    if len(mutation_str) < 3:
        return False
        
    original_aa = mutation_str[0]
    new_aa = mutation_str[-1]
    position = mutation_str[1:-1]
    
    try:
        # Check if residue exists
        selection = f"resi {position} and {object_name}"
        if cmd.count_atoms(selection) == 0:
            return False
        
        # Start mutagenesis wizard
        cmd.wizard("mutagenesis")
        cmd.refresh_wizard()
        
        # Set the mutation
        cmd.get_wizard().set_mode(new_aa)
        cmd.get_wizard().do_select(selection)
        cmd.get_wizard().apply()
        
        # End wizard
        cmd.set_wizard()
        
        print(f"Applied mutation (wizard): {original_aa}{position}{new_aa}")
        return True
        
    except Exception as e:
        print(f"Error with wizard method: {str(e)}")
        return False

def create_mutant_pdb_simple(level, mutations, output_dir="mutant_pdbs"):
    """
    Create mutant PDB using simple residue name replacement
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
        success = apply_mutation_direct(mutation, "protein")
        if success:
            successful_mutations.append(mutation)
        else:
            # Try wizard method as backup
            if apply_mutation_wizard(mutation, "protein"):
                successful_mutations.append(mutation)
            else:
                failed_mutations.append(mutation)
    
    # Rebuild and optimize
    cmd.h_add("protein")
    cmd.remove("solvent")
    
    # Save file
    clean_mutations = [m.replace('U', 'Sec') for m in mutations]  # Make filename safe
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

def create_pml_script(output_dir="mutant_pdbs"):
    """
    Create a PyMOL script file (.pml) for manual execution
    """
    os.makedirs(output_dir, exist_ok=True)
    pml_content = f"""# PyMOL script for GPX6 mutations
# Load this script in PyMOL: File -> Run Script -> select this file

load {base_pdb}, original

"""
    
    for level in range(1, 21):
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
            
            pml_content += f"h_add level{level}\n"
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
    print("GPX6 Mutation Generator")
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
        for level in range(1, 21):
            if level in mutations_by_level:
                create_mutant_pdb_simple(level, mutations_by_level[level])
    except Exception as e:
        print(f"Direct method failed: {str(e)}")
    
    print("\n" + "=" * 50)
    print("COMPLETE!")
    print("If Python method failed, use the .pml script in PyMOL")

if __name__ == "__main__":
    main()