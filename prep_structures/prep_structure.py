import os
from pymol import cmd

# Mapping of one-letter to three-letter amino acid codes
amino_acid_map = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU', 
    'Q': 'GLN', 'G': 'GLY', 'H': 'HID', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 
    'Y': 'TYR', 'V': 'VAL', 'U': 'SEC'  # Adding SEC for selenocysteine
}

def read_mutations_from_file(file_path):
    """
    Reads mutation data from a text file.
    :param file_path: Path to the mutations.txt file
    :return: List of mutations in the format [(ResidueNumber, OriginalAA, MutatedAA), ...]
    """
    mutations = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                res_num, original_aa, mutated_aa = parts
                mutations.append((int(res_num), original_aa, mutated_aa))
    return mutations

def fix_histidine_naming(pdb_file, output_file):
    """
    Fixes histidine naming in PDB files to ensure HIS residues are renamed to HID
    """
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    modified_lines = []
    for line in lines:
        if line.startswith('ATOM') and 'HIS' in line[17:20]:
            line = line[:17] + 'HID' + line[20:]
        modified_lines.append(line)
    
    with open(output_file, 'w') as f:
        f.writelines(modified_lines)

def mutate_residues(pdb_file, mutations, output_dir):
    """
    Mutate residues in a PDB file using PyMOL and save each mutation as a separate file.
    :param pdb_file: Path to the input PDB file.
    :param mutations: List of mutation tuples, e.g., [(49, "C", "U")].
    :param output_dir: Directory to save mutated PDB files.
    """
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist. Creating it.")
        os.makedirs(output_dir)
    else:
        print(f"Directory {output_dir} already exists.")
    
    # Load the PDB file
    print(f"Loading PDB file from: {pdb_file}")
    cmd.load(pdb_file, "structure")
    print(f"Structure names: {cmd.get_names()}")  # Check if structure loaded correctly
    
    for resi, mouse_residue, human_residue in mutations:
        selection = f"resi {resi}"
        print(f"Applying mutation at residue {resi}: {mouse_residue} -> {human_residue}")
        
        cmd.wizard("mutagenesis")
        cmd.do(f"select {selection}")
        cmd.get_wizard().do_select(selection)
        
        three_letter_residue = amino_acid_map.get(human_residue, "UNK")
        cmd.get_wizard().set_mode(three_letter_residue)
        cmd.get_wizard().apply()
        
        # Close wizard and remove hydrogens
        cmd.set_wizard()
        cmd.remove("hydro")
        
        # Remove water molecules (solvent)
        cmd.remove("solvent")
        
        # Construct output file path
        temp_output_file = os.path.join(output_dir, f"temp_{mouse_residue}{resi}{human_residue}.pdb")
        final_output_file = os.path.join(output_dir, f"{mouse_residue}{resi}{human_residue}.pdb")
        
        # Save the PDB file
        try:
            cmd.save(temp_output_file, "structure")
            print(f"Temporarily saved to: {temp_output_file}")
            
            # Fix histidine naming if needed
            if human_residue == 'H':
                print(f"Fixing histidine naming for mutation {mouse_residue}{resi}{human_residue}")
                fix_histidine_naming(temp_output_file, final_output_file)
                os.remove(temp_output_file)  # Remove temporary file
                print(f"Successfully saved with fixed histidine: {final_output_file}")
            else:
                # If not histidine, just rename the file
                os.rename(temp_output_file, final_output_file)
                print(f"Successfully saved: {final_output_file}")
        except Exception as e:
            print(f"Error processing PDB file: {e}")
    
    cmd.delete("all")

if __name__ == "__main__":
    # Define input files using absolute paths
    pdb_file_path = "/home/hp/results/MOUSE/level7/H144Q/minim/minim.pdb"  # Absolute path to PDB file
    mutations_file = "/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/level8.txt"  # Absolute path to mutation list file
    
    # Ensure output directory exists using an absolute path
    output_dir = "/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/level8"
    
    # Read mutations from file
    mutations_list = read_mutations_from_file(mutations_file)
    
    # Run mutation function
    mutate_residues(pdb_file_path, mutations_list, output_dir)