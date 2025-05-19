import os
import subprocess
import tempfile
import shutil

# Mapping of one-letter to three-letter amino acid codes
amino_acid_map = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU', 
    'Q': 'GLN', 'G': 'GLY', 'H': 'HID', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',  # Changed 'H': 'HIS' to 'H': 'HID'
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
            line = line.strip()
            if line:  # Skip empty lines
                parts = line.split()
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

def run_pymol_mutation(pdb_file, resi, original_aa, mutated_aa, output_dir):
    """
    Run PyMOL directly to perform a mutation and save the result
    """
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # First, create a temporary PDB with HIS→HID conversion
    temp_input_file = tempfile.mktemp(suffix='.pdb')
    fix_histidine_naming(pdb_file, temp_input_file)
    
    output_file = os.path.join(output_dir, f"{original_aa}{resi}{mutated_aa}.pdb")
    three_letter_residue = amino_acid_map.get(mutated_aa, "UNK")
    
    # Create a one-line PyMOL command to execute
    cmd_line = (
        f'reinitialize; load {temp_input_file}, struct; '
        f'wizard mutagenesis; '
        f'select target, resi {resi}; '
        f'cmd.get_wizard().do_select("target"); '
        f'cmd.get_wizard().set_mode("{three_letter_residue}"); '
        f'cmd.get_wizard().apply(); '
        f'cmd.set_wizard(); '
        f'remove hydro; remove solvent; '
        f'save {output_file}, struct'
    )
    
    try:
        print(f"Running PyMOL for mutation {original_aa}{resi}{mutated_aa}")
        # Run PyMOL in command mode with the command line passed via -x
        result = subprocess.run(
            ['pymol', '-cq', '-d', 'python', '-x', cmd_line],
            capture_output=True, 
            text=True,
            check=False
        )
        
        # Print output for debugging
        print(f"PyMOL stdout: {result.stdout}")
        print(f"PyMOL stderr: {result.stderr}")
        
        # Check if output file was created
        if os.path.exists(output_file):
            print(f"Success: Created {output_file}")
            
            # Always fix histidine naming
            print(f"Fixing histidine naming for mutation {original_aa}{resi}{mutated_aa}")
            temp_file = output_file + ".temp"
            shutil.move(output_file, temp_file)
            fix_histidine_naming(temp_file, output_file)
            os.remove(temp_file)
            
            os.remove(temp_input_file)  # Cleanup temp input file
            return True
        else:
            print(f"Error: Failed to create output file {output_file}")
            os.remove(temp_input_file)  # Cleanup temp input file
            return False
            
    except Exception as e:
        print(f"Error running PyMOL: {e}")
        if os.path.exists(temp_input_file):
            os.remove(temp_input_file)  # Cleanup temp input file
        return False

def create_direct_pymol_script(pdb_file, mutations, output_dir):
    """
    Create a direct PyMOL script approach as alternative
    """
    # Create the output directory if needed
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # First, create a temporary PDB with HIS→HID conversion
    temp_input_file = tempfile.mktemp(suffix='.pdb')
    fix_histidine_naming(pdb_file, temp_input_file)
        
    # Create a temporary script file
    with tempfile.NamedTemporaryFile(suffix='.pml', delete=False, mode='w') as f:
        script_path = f.name
        
        # Write header
        f.write("# PyMOL script to perform mutations\n\n")
        
        # Process each mutation
        for resi, original_aa, mutated_aa in mutations:
            output_file = os.path.join(output_dir, f"{original_aa}{resi}{mutated_aa}.pdb")
            three_letter_residue = amino_acid_map.get(mutated_aa, "UNK")
            
            f.write(f"# Processing mutation {original_aa}{resi}{mutated_aa}\n")
            f.write("reinitialize\n")
            f.write(f"load {temp_input_file}, structure\n")
            f.write("wizard mutagenesis\n")
            f.write(f"select target, resi {resi}\n")
            f.write("cmd.get_wizard().do_select('target')\n")
            f.write(f"cmd.get_wizard().set_mode('{three_letter_residue}')\n")
            f.write("cmd.get_wizard().apply()\n")
            f.write("cmd.set_wizard()\n")
            f.write("remove hydro\n")
            f.write("remove solvent\n")
            f.write(f"save {output_file}, structure\n\n")
    
    # Run the script
    try:
        print(f"Running PyMOL script for all mutations")
        result = subprocess.run(
            ['pymol', '-cq', '-d', 'python', script_path],
            capture_output=True, 
            text=True,
            check=False
        )
        
        print(f"PyMOL stdout: {result.stdout}")
        print(f"PyMOL stderr: {result.stderr}")
        
        # Check for output files
        success = True
        for resi, original_aa, mutated_aa in mutations:
            output_file = os.path.join(output_dir, f"{original_aa}{resi}{mutated_aa}.pdb")
            if not os.path.exists(output_file):
                print(f"Failed to create: {output_file}")
                success = False
            else:
                print(f"Successfully created: {output_file}")
                
                # Always fix histidine naming for all output files
                print(f"Fixing histidine naming for mutation {original_aa}{resi}{mutated_aa}")
                temp_file = output_file + ".temp"
                shutil.move(output_file, temp_file)
                fix_histidine_naming(temp_file, output_file)
                os.remove(temp_file)
        
        # Clean up script file and temp input file
        os.remove(script_path)
        os.remove(temp_input_file)
        return success
        
    except Exception as e:
        print(f"Error running PyMOL script: {e}")
        # Clean up script file and temp input file
        if os.path.exists(script_path):
            os.remove(script_path)
        if os.path.exists(temp_input_file):
            os.remove(temp_input_file)
        return False

def create_pymol_batch_file(pdb_file, mutations, output_dir):
    """
    Create a batch file approach for running each mutation individually
    """
    # Create temporary directory for batch processing
    temp_dir = tempfile.mkdtemp()
    
    # First, create a temporary PDB with HIS→HID conversion
    temp_input_file = tempfile.mktemp(suffix='.pdb')
    fix_histidine_naming(pdb_file, temp_input_file)
    
    # Create individual script files for each mutation
    script_files = []
    for resi, original_aa, mutated_aa in mutations:
        output_file = os.path.join(output_dir, f"{original_aa}{resi}{mutated_aa}.pdb")
        three_letter_residue = amino_acid_map.get(mutated_aa, "UNK")
        
        script_file = os.path.join(temp_dir, f"mutation_{original_aa}{resi}{mutated_aa}.pml")
        with open(script_file, 'w') as f:
            f.write("reinitialize\n")
            f.write(f"load {temp_input_file}, structure\n")
            f.write("wizard mutagenesis\n")
            f.write(f"select target, resi {resi}\n")
            f.write("cmd.get_wizard().do_select('target')\n")
            f.write(f"cmd.get_wizard().set_mode('{three_letter_residue}')\n")
            f.write("cmd.get_wizard().apply()\n")
            f.write("cmd.set_wizard()\n")
            f.write("remove hydro\n")
            f.write("remove solvent\n")
            f.write(f"save {output_file}, structure\n")
            f.write("quit\n")
        
        script_files.append((script_file, original_aa, resi, mutated_aa))
    
    # Create the output directory if needed
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Process each script file
    for script_file, original_aa, resi, mutated_aa in script_files:
        try:
            print(f"Processing mutation {original_aa}{resi}{mutated_aa}")
            subprocess.run(['pymol', '-qc', script_file], check=False)
            
            output_file = os.path.join(output_dir, f"{original_aa}{resi}{mutated_aa}.pdb")
            if os.path.exists(output_file):
                print(f"Successfully created: {output_file}")
                
                # Always fix histidine naming for all output files
                print(f"Fixing histidine naming for output file")
                temp_file = output_file + ".temp"
                shutil.move(output_file, temp_file)
                fix_histidine_naming(temp_file, output_file)
                os.remove(temp_file)
            else:
                print(f"Failed to create: {output_file}")
                
        except Exception as e:
            print(f"Error processing {script_file}: {e}")
    
    # Clean up temporary directory and input file
    shutil.rmtree(temp_dir)
    os.remove(temp_input_file)

if __name__ == "__main__":
    # Define input files using absolute paths
    pdb_file_path = "/home/hp/results/MOUSE/level13/T52A/minim/minim.pdb"  # Absolute path to PDB file
    mutations_file = "/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/level14.txt"   # Absolute path to mutation list file
    
    # Ensure output directory exists using an absolute path
    output_dir = "/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/level14"
    
    # First, create a version of the input PDB with HIS fixed to HID
    print("Converting input PDB file histidines from HIS to HID")
    fixed_input_pdb = pdb_file_path + ".fixed.pdb"
    fix_histidine_naming(pdb_file_path, fixed_input_pdb)
    
    # Read mutations from file
    mutations_list = read_mutations_from_file(mutations_file)
    
    # Method selection - try different approaches
    print("Using batch file approach for mutations")
    create_pymol_batch_file(fixed_input_pdb, mutations_list, output_dir)
    
    # Clean up the fixed input PDB
    os.remove(fixed_input_pdb)