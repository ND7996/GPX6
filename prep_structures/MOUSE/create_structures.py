#!/usr/bin/env python3
"""
Simple PyMOL script to generate mutant PDB files
ONLY applies the specified mutations, no extra changes
"""

import pymol
from pymol import cmd
import os

# Initialize PyMOL
pymol.finish_launching(['pymol', '-c'])

# Define the base structure path
#base_pdb = "/home/hp/nayanika/github/GPX6/mouseWT/GPX6mousecys.pdb"
base_pdb = "/home/hp/nayanika/github/PhD_Thesis/EVB/protein_stepwise/GPX6WT/mousesec/1-prep/GPX6sec_mouse.pdb"

# Define mutations for each level
mutations_by_level = {
    2: ["T54Q"],
    3: ["T54Q", "I24L"],
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

def apply_simple_mutation(mutation_str, object_name="mutant"):
    """
    Apply mutation using PyMOL wizard - SIMPLE version
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
        
        # Map to PyMOL's 3-letter codes
        aa_map = {
            'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
            'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
            'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
            'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
            'U': 'CYS'  # TEMPORARY: Use CYS first, then manually change to SEC
        }
        
        new_resn = aa_map.get(new_aa, new_aa)
        
        print(f"Applying: {original_aa}{position}{new_aa}")
        
        # Use mutagenesis wizard
        cmd.wizard("mutagenesis")
        cmd.refresh_wizard()
        cmd.get_wizard().do_select(f"({object_name} and resi {position})")
        cmd.get_wizard().set_mode(new_resn)
        cmd.get_wizard().apply()
        cmd.set_wizard()
        
        # MANUALLY change CYS to SEC for selenocysteine
        if new_aa == 'U':
            cmd.alter(f"{object_name} and resi {position}", "resn='SEC'")
            cmd.alter(f"{object_name} and resi {position} and name SG", "name='SE'")
            print(f"  Manually converted CYS {position} to SEC {position}")
        
        print(f"Successfully applied: {original_aa}{position}{new_aa}")
        return True
        
    except Exception as e:
        print(f"Error applying mutation {mutation_str}: {str(e)}")
        return False

def create_mutant_level(level, mutations, output_dir="mutant_pdbs"):
    """
    Create mutant PDB for a specific level
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Start fresh
    cmd.reinitialize()
    
    # Load the base wild-type structure
    cmd.load(base_pdb, "mutant")
    
    print(f"\n=== Creating Level {level} mutant ===")
    print(f"Mutations to apply: {', '.join(mutations)}")
    
    successful_mutations = []
    failed_mutations = []
    
    # Apply ALL mutations for this level
    for mutation in mutations:
        success = apply_simple_mutation(mutation, "mutant")
        if success:
            successful_mutations.append(mutation)
        else:
            failed_mutations.append(mutation)
    
    # Save file
    output_filename = f"GPX6_level{level:02d}.pdb"
    output_path = os.path.join(output_dir, output_filename)
    
    try:
        cmd.save(output_path, "mutant")
        print(f"✓ Saved: {output_path}")
        if successful_mutations:
            print(f"  Applied: {', '.join(successful_mutations)}")
        if failed_mutations:
            print(f"  Failed: {', '.join(failed_mutations)}")
            
    except Exception as e:
        print(f"✗ Error saving: {str(e)}")
    
    # Clean up
    cmd.delete("mutant")

def main():
    """
    Main function to create all mutant levels
    """
    print("GPX6 Mouse Mutation Generator")
    print("=" * 50)
    print(f"Base structure: {base_pdb}")
    print("Simple mutagenesis only - no extra changes")
    
    if not os.path.exists(base_pdb):
        print(f"ERROR: Base PDB file not found!")
        return
    
    # Create mutants for each level
    try:
        for level in range(1, 21):
            if level in mutations_by_level:
                create_mutant_level(level, mutations_by_level[level])
            else:
                print(f"No mutations defined for level {level}")
                
    except Exception as e:
        print(f"Error during mutation process: {str(e)}")
    
    print("\n" + "=" * 50)
    print("COMPLETE!")
    print("Mutant PDBs created with simple mutagenesis")

if __name__ == "__main__":
    main()