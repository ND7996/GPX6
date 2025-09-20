import os
import subprocess
import tempfile
import shutil

# ================================
# HARDCODED PATHS - CHANGE THESE!
# ================================
BASE_PDB_FILE = "/home/hp/nayanika/github/GPX6/prep_structures/HUMAN/mutant_pdbs/GPX6_level20.pdb"
MUTATIONS_FILE = "/home/hp/nayanika/github/GPX6/prep_structures/HUMAN/level20.txt"
OUTPUT_DIR = "/home/hp/nayanika/github/GPX6/prep_structures/HUMAN/level20"
LEVEL_NUMBER = 3
# ================================

# Mapping of one-letter to three-letter amino acid codes
amino_acid_map = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU', 
    'Q': 'GLN', 'G': 'GLY', 'H': 'HID', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 
    'Y': 'TYR', 'V': 'VAL', 'U': 'SEC'
}

def read_mutations_from_file(file_path):
    """
    Reads mutation data from a text file.
    """
    mutations = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
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

def create_isolated_mutation_script(base_pdb, resi, original_aa, mutated_aa, output_file):
    """
    Create a completely isolated PyMOL script for a single mutation
    """
    three_letter_mutated = amino_acid_map.get(mutated_aa, "UNK")
    
    script_content = f"""# Isolated mutation script for {original_aa}{resi}{mutated_aa}
# Load the original structure
load {base_pdb}, original

# Create a fresh copy for mutation (ensures no accumulation)
create mutant, original

# Select only the mutant structure and perform single mutation
wizard mutagenesis
select target, mutant and resi {resi}
cmd.get_wizard().do_select("target")
cmd.get_wizard().set_mode("{three_letter_mutated}")
cmd.get_wizard().apply()
cmd.set_wizard()

# Clean up
remove hydro
remove solvent

# Save only the mutant structure
save {output_file}, mutant

# Quit completely to ensure no state is preserved
quit
"""
    return script_content

def process_single_mutation_isolated(base_pdb, resi, original_aa, mutated_aa, output_file):
    """
    Process a single mutation in complete isolation
    """
    # Create temporary script file
    script_file = tempfile.mktemp(suffix='.pml')
    script_content = create_isolated_mutation_script(base_pdb, resi, original_aa, mutated_aa, output_file)
    
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    try:
        print(f"  Processing {original_aa}{resi}{mutated_aa}...")
        
        # Run PyMOL with the script - this starts a completely fresh instance
        result = subprocess.run([
            'pymol', '-cq', script_file
        ], capture_output=True, text=True, timeout=120)  # 2 minute timeout
        
        # Check if output was created
        if os.path.exists(output_file):
            # Fix histidine naming
            temp_file = output_file + ".temp"
            shutil.move(output_file, temp_file)
            fix_histidine_naming(temp_file, output_file)
            os.remove(temp_file)
            print(f"  ✅ Success: {original_aa}{resi}{mutated_aa}.pdb")
            return True
        else:
            print(f"  ❌ Failed: {original_aa}{resi}{mutated_aa} - No output file created")
            if result.stderr:
                print(f"    Error: {result.stderr[:200]}...")  # Show first 200 chars of error
            return False
            
    except subprocess.TimeoutExpired:
        print(f"  ❌ Failed: {original_aa}{resi}{mutated_aa} - Timeout")
        return False
    except Exception as e:
        print(f"  ❌ Failed: {original_aa}{resi}{mutated_aa} - {e}")
        return False
    finally:
        # Clean up script file
        if os.path.exists(script_file):
            os.remove(script_file)

def process_single_level():
    """
    Process a single level with hardcoded paths
    """
    # Check if files exist
    if not os.path.exists(BASE_PDB_FILE):
        print(f"Error: Base PDB file not found: {BASE_PDB_FILE}")
        return False
    
    if not os.path.exists(MUTATIONS_FILE):
        print(f"Error: Mutations file not found: {MUTATIONS_FILE}")
        return False
    
    print(f"\n{'='*60}")
    print(f"PROCESSING LEVEL {LEVEL_NUMBER}")
    print(f"{'='*60}")
    print(f"Base PDB: {BASE_PDB_FILE}")
    print(f"Mutations file: {MUTATIONS_FILE}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Read mutations
    mutations_list = read_mutations_from_file(MUTATIONS_FILE)
    print(f"Found {len(mutations_list)} mutations")
    
    # Create output directory
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created output directory: {OUTPUT_DIR}")
    
    # Create fixed version of input PDB (HIS -> HID)
    print("Creating fixed base PDB (HIS -> HID)...")
    fixed_base_pdb = tempfile.mktemp(suffix='.pdb')
    fix_histidine_naming(BASE_PDB_FILE, fixed_base_pdb)
    
    # Process each mutation in complete isolation
    success_count = 0
    total_mutations = len(mutations_list)
    
    for i, (resi, original_aa, mutated_aa) in enumerate(mutations_list, 1):
        output_file = os.path.join(OUTPUT_DIR, f"{original_aa}{resi}{mutated_aa}.pdb")
        
        print(f"\n[{i}/{total_mutations}] Creating {original_aa}{resi}{mutated_aa}.pdb")
        
        if process_single_mutation_isolated(fixed_base_pdb, resi, original_aa, mutated_aa, output_file):
            success_count += 1
    
    # Clean up
    if os.path.exists(fixed_base_pdb):
        os.remove(fixed_base_pdb)
    
    print(f"\nLevel {LEVEL_NUMBER} Summary:")
    print(f"Successful mutations: {success_count}/{total_mutations}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"{'='*60}")
    
    return success_count == total_mutations

def main():
    """
    Main function
    """
    success = process_single_level()
    
    if success:
        print("🎉 All mutations processed successfully!")
    else:
        print("⚠️  Some mutations failed. Check the output above.")

if __name__ == "__main__":
    main()