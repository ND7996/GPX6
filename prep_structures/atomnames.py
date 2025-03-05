import os

# Define the base directory containing all PDB subdirectories with minim.pdb files
base_directory = '/home/hp/nayanika/github/GPX6/prep_structures/HUMAN/level0'

# Define atom name replacements for residue SEC 49
atom_replacements = {
    'HB3': 'HB2',
    'HB2': 'HB1',  
    'HG':  'HG1',
    'H':  'HN',
}

def replace_atom_names(pdb_file):
    print(f"Processing file: {pdb_file}")
    output_lines = []

    with open(pdb_file, 'r') as file:
        for line in file:
            # Only consider lines with ATOM and SEC 49
            if line.startswith("ATOM") and " SEC " in line and " 49 " in line:
                # If the atom is HB1, replace it directly with HB2
                if 'HB1' in line:
                    print(f"Original line (HB1): {line.strip()}")
                    line = line.replace('HB1', 'HB2')
                    print(f"Replaced 'HB1' with 'HB2' in line: {line.strip()}")
                # Then proceed with other replacements
                for old_name, new_name in atom_replacements.items():
                    if f" {old_name} " in line:
                        line = line.replace(f" {old_name} ", f" {new_name} ")
                        print(f"Replaced '{old_name}' with '{new_name}' in line: {line.strip()}")

            # Append the line (modified or not) to output
            output_lines.append(line)

    # Write back to the same file to overwrite it
    with open(pdb_file, 'w') as out_file:
        out_file.writelines(output_lines)
    
    print(f"Modified file saved as: {pdb_file}")

# Iterate through all files in the base directory and its subdirectories
for root, dirs, files in os.walk(base_directory):
    for file in files:
        # Check if the file has a '.pdb' extension
        if file.endswith('_solvated.pdb'):
            pdb_filepath = os.path.join(root, file)
            replace_atom_names(pdb_filepath)
