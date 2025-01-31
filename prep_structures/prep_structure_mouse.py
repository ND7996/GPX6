import os
from pymol import cmd

# Mapping of one-letter to three-letter amino acid codes
amino_acid_map = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU', 
    'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 
    'Y': 'TYR', 'V': 'VAL', 'U': 'SEC'  # Adding SEC for selenocysteine
}

# Mutation list (Mouse Cys -> Human Sec)
mutations_list = [
    (4, 'S', 'R'), (47, 'S', 'A'), (48, 'F', 'Y'), (60, 'T', 'A'), (24, 'I', 'L'),
    (52, 'T', 'A'), (54, 'T', 'Q'), (87, 'K', 'T'), (99, 'R', 'C'),
    (102, 'G', 'S'), (104, 'Y', 'F'), (107, 'N', 'S'), (142, 'P', 'S'), (143, 'E', 'S'),
    (144, 'H', 'Q'), (148, 'D', 'E'), (177, 'H', 'Q'), (178, 'T', 'A'), (139, 'F', 'L'),
    (181, 'R', 'S')
]

# Directory to save the modified structures
output_dir = "/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"

# Ensure the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Function to mutate residues in a given PDB file
def mutate_residues(pdb_file, mutations):
    """
    Mutate residues in a PDB file using PyMOL and save each mutation as a separate file.
    
    Parameters:
    - pdb_file: Path to the input PDB file.
    - mutations: List of mutation tuples, e.g., [(49, "C", "U")].
    """
    # Load the PDB file
    cmd.load(pdb_file, "structure")
    
    # Loop through each mutation request
    for resi, mouse_residue, human_residue in mutations:
        selection = f"resi {resi}"  # Apply mutation to the specified residue, without chain specifier
        
        # Check if the current residue matches the intended target
        cmd.wizard("mutagenesis")
        cmd.do(f"select {selection}")
        cmd.get_wizard().do_select(selection)
        
        # Special handling for SEC at position 49
        if resi == 49 and human_residue == 'U':  # Explicitly set residue 49 to SEC
            cmd.alter(selection, "resn='SEC'")
            cmd.alter(f"{selection} and name HB3", "name='HB2'")
            cmd.alter(f"{selection} and name HB2", "name='HB1'")
            cmd.alter(f"{selection} and name HG", "name='HG1'")
            cmd.alter(f"{selection} and name H", "name='HN'")
            print(f"Residue 49 set to SEC with modified atom names.")
        else:
            three_letter_residue = amino_acid_map.get(human_residue, "UNK")
            cmd.get_wizard().set_mode(three_letter_residue)
            cmd.get_wizard().apply()
        
        # Close the mutagenesis wizard
        cmd.set_wizard()
        
        # Remove hydrogens
        cmd.remove("hydro")
        
        # Output file name using residue name and number (e.g., "K3N.pdb" for K->N at residue 3)
        output_file = f"{output_dir}/{mouse_residue}{resi}{human_residue}.pdb"
        
        # Save the mutated structure
        cmd.save(output_file, "structure", format="pdb")
        print(f"Mutation applied at residue {resi}: {mouse_residue} -> {human_residue}. Saved to {output_file}.")
    
    cmd.delete("all")

# Example usage 
pdb_file_path = "/home/hp/nayanika/github/GPX6/prep_structures/C49U/C49U.pdb"
 # Replace with your actual PDB file path

# Run the mutation function
mutate_residues(pdb_file_path, mutations_list)
