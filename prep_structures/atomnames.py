#!/usr/bin/env python3

import os
import re

# Define the base directory containing all PDB files
base_directory = "/home/hp/nayanika/github/GPX6/prep_structures/HUMAN/level20"

def fix_pdb_file(pdb_file):
    print(f"Processing file: {pdb_file}")
    
    # Read all lines
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    # Process SEC 49 atoms
    new_lines = []
    hb_index = 0
    
    for line in lines:
        # Only process lines that are ATOM/HETATM and specifically for SEC residue 49
        if line.startswith(('ATOM', 'HETATM')) and 'SEC' in line and '49' in line:
            # Extract residue number more precisely (column 23-26 in PDB format)
            # PDB format: columns 23-26 contain residue number
            if len(line) >= 26:
                residue_num = line[22:26].strip()
                if residue_num == '49':
                    # Fix atom names only for residue 49
                    atom_name = line[12:16].strip()
                    
                    # Fix H → HN
                    if atom_name == 'H':
                        line = line[:12] + ' HN ' + line[16:]
                    
                    # Fix HG → HG1
                    elif atom_name == 'HG':
                        line = line[:12] + ' HG1' + line[16:]
                    
                    # Fix SG → SE (only for residue 49)
                    elif atom_name == 'SG':
                        line = line[:12] + ' SE ' + line[16:]
                    
                    # Fix HB atoms in order
                    elif atom_name.startswith('HB'):
                        hb_index += 1
                        if hb_index == 1:
                            line = line[:12] + ' HB1' + line[16:]
                        elif hb_index == 2:
                            line = line[:12] + ' HB2' + line[16:]
                        elif hb_index == 3:
                            line = line[:12] + ' HB3' + line[16:]
        
        new_lines.append(line)
    
    # Write the corrected file
    with open(pdb_file, 'w') as f:
        f.writelines(new_lines)
    
    print(f"Fixed file saved: {pdb_file}")

def verify_pdb_file(pdb_file):
    issues = 0
    print(f"Verifying: {pdb_file}")
    
    with open(pdb_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # Only check lines that are ATOM/HETATM and specifically for SEC residue 49
            if line.startswith(('ATOM', 'HETATM')) and 'SEC' in line:
                # Extract residue number (columns 23-26 in PDB format)
                if len(line) >= 26:
                    residue_num = line[22:26].strip()
                    if residue_num == '49':
                        # Extract atom name (positions 13-16 in PDB format)
                        atom_name = line[12:16].strip()
                        
                        if atom_name == "H":
                            print(f"  ERROR: H not converted to HN at line {line_num}")
                            issues += 1
                        if atom_name == "HG":
                            print(f"  ERROR: HG not converted to HG1 at line {line_num}")
                            issues += 1
                        if atom_name == "SG":
                            print(f"  ERROR: SG not converted to SE at line {line_num}")
                            issues += 1
                        if atom_name.startswith("HB") and atom_name not in ["HB1", "HB2", "HB3"]:
                            print(f"  ERROR: Incorrect HB atom name '{atom_name}' at line {line_num}")
                            issues += 1
    
    if issues == 0:
        print("  ✓ All checks passed")
        return True
    else:
        print(f"  ✗ Found {issues} issues")
        return False

# Main processing
if __name__ == "__main__":
    print("Starting SEC residue 49 atom name fix...")
    print("========================================")
    
    # Process all PDB files
    all_verified = True
    pdb_files_processed = 0
    
    for root, dirs, files in os.walk(base_directory):
        for file in files:
            if file.endswith('.pdb'):
                pdb_filepath = os.path.join(root, file)
                fix_pdb_file(pdb_filepath)
                pdb_files_processed += 1
    
    print(f"\nProcessed {pdb_files_processed} PDB files")
    
    print("\nVerification...")
    print("========================================")
    
    # Verify all files
    for root, dirs, files in os.walk(base_directory):
        for file in files:
            if file.endswith('.pdb'):
                pdb_filepath = os.path.join(root, file)
                if not verify_pdb_file(pdb_filepath):
                    all_verified = False
    
    print("\n========================================")
    if all_verified:
        print("✓ All files have been successfully fixed!")
    else:
        print("✗ Some issues remain. Please check the verification output.")