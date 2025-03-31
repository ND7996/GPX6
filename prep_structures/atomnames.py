import os
import re

# Define the base directory containing all PDB files
base_directory = '/home/hp/nayanika/github/GPX6/prep_structures/HUMAN/level2'

def fix_pdb_file(pdb_file):
    print(f"Processing file: {pdb_file}")
    
    with open(pdb_file, 'r') as file:
        content = file.read()

    # First, collect all the lines with SEC 49
    sec_lines = []
    pattern = re.compile(r'^(ATOM\s+\d+\s+\S+\s+SEC\s+49.*)$', re.MULTILINE)
    for match in pattern.finditer(content):
        sec_lines.append(match.group(1))
    
    # Create a dictionary to track atom names
    atom_counter = {}
    
    # Step 1: Replace H with HN with exact spacing
    content = re.sub(
        r'^(ATOM\s+\d+\s+)H(\s+SEC\s+49.*)$',
        r'\1HN \2',
        content,
        flags=re.MULTILINE
    )
    
    # Step 2: Fix the HG1 spacing - ensures exactly one space between HG1 and SEC
    content = re.sub(
        r'^(ATOM\s+\d+\s+)HG1\s+(\s*SEC\s+49.*)$',
        r'\1HG1 \2',
        content,
        flags=re.MULTILINE
    )
    
    # Step 3: Fix the HB names in proper sequence
    hb_lines = []
    hb_pattern = re.compile(r'^(ATOM\s+\d+\s+)(HB\d)(\s+SEC\s+49.*)$', re.MULTILINE)
    for match in hb_pattern.finditer(content):
        hb_lines.append((match.group(0), match.group(1), match.group(2), match.group(3)))
    
    # Sort the HB lines by atom number to ensure correct ordering
    hb_lines.sort(key=lambda x: int(re.search(r'ATOM\s+(\d+)', x[0]).group(1)))
    
    # Replace HB atoms in original order
    for i, (original, prefix, atom, suffix) in enumerate(hb_lines):
        if i == 0:  # First HB should be HB1
            new_atom = "HB1"
        elif i == 1:  # Second HB should be HB2
            new_atom = "HB2"
        else:  # Additional HB atoms (if any)
            new_atom = f"HB{i+1}"
            
        # Replace with corrected spacing
        content = content.replace(original, f"{prefix}{new_atom}{suffix}")
    
    # Fix HG atom name
    content = re.sub(
        r'^(ATOM\s+\d+\s+)HG(\s+SEC\s+49.*)$',
        r'\1HG1\2',
        content,
        flags=re.MULTILINE
    )
    
    # Write the corrected file
    with open(pdb_file, 'w') as file:
        file.write(content)
    
    print(f"Fixed file saved: {pdb_file}")

# Process all solvated PDB files
for root, dirs, files in os.walk(base_directory):
    for file in files:
        if file.endswith('_solvated.pdb'):
            pdb_filepath = os.path.join(root, file)
            fix_pdb_file(pdb_filepath)

# Verify the fixes
def verify_pdb_file(pdb_file):
    issues = []
    with open(pdb_file, 'r') as file:
        for line_num, line in enumerate(file, 1):
            # Check SEC 49 lines
            if 'SEC' in line and '49' in line and line.startswith('ATOM'):
                # Check HB atom name spacing
                hb_match = re.match(r'ATOM\s+\d+\s+(HB\d)\s+SEC', line)
                if hb_match and len(hb_match.group(1)) != 3:
                    issues.append(f"Line {line_num}: Incorrect HB spacing - {line.strip()}")
                
                # Check HG1 spacing
                if 'HG1' in line:
                    if not re.match(r'ATOM\s+\d+\s+HG1\s+SEC', line):
                        issues.append(f"Line {line_num}: Incorrect HG1 spacing - {line.strip()}")
                
                # Check for H that should be HN
                if re.match(r'ATOM\s+\d+\s+H\s+SEC', line):
                    issues.append(f"Line {line_num}: H not converted to HN - {line.strip()}")
    
    if issues:
        print(f"\nIssues in {pdb_file}:")
        for issue in issues:
            print(f"  {issue}")
        return False
    return True

print("\nVerifying all fixed files...")
all_verified = True
for root, dirs, files in os.walk(base_directory):
    for file in files:
        if file.endswith('_solvated.pdb'):
            pdb_filepath = os.path.join(root, file)
            if not verify_pdb_file(pdb_filepath):
                all_verified = False

if all_verified:
    print("\nAll files have been successfully fixed!")
else:
    print("\nSome issues remain. Please check the verification output.")
