import os
import re

# Define the base directory containing mutation folders
base_directory = '/home/hp/results//level1'
mutation_folders = [
    "A178T", "A47S", "A52T", "A60T", "C99R", "E148D", "F104Y",
    "L139F", "L24I", "Q144H", "Q177H", "Q54T", "R4S", "S102G",
    "S107N", "S142P", "S143E", "S181R", "Y48F"
]

def fix_pdb_file(pdb_file):
    print(f"Processing file: {pdb_file}")
    
    with open(pdb_file, 'r') as file:
        content = file.read()

    # Fix H → HN
    content = re.sub(
        r'^(ATOM\s+\d+\s+)H(\s+SEC\s+49.*)$',
        r'\1HN \2',
        content,
        flags=re.MULTILINE
    )

    # Fix HG1 spacing
    content = re.sub(
        r'^(ATOM\s+\d+\s+)HG1\s+(\s*SEC\s+49.*)$',
        r'\1HG1 \2',
        content,
        flags=re.MULTILINE
    )

    # Fix HB names in proper sequence
    hb_lines = []
    hb_pattern = re.compile(r'^(ATOM\s+\d+\s+)(HB\d)(\s+SEC\s+49.*)$', re.MULTILINE)
    for match in hb_pattern.finditer(content):
        hb_lines.append((match.group(0), match.group(1), match.group(2), match.group(3)))
    
    hb_lines.sort(key=lambda x: int(re.search(r'ATOM\s+(\d+)', x[0]).group(1)))

    for i, (original, prefix, atom, suffix) in enumerate(hb_lines):
        new_atom = f"HB{i+1}" if i < 2 else f"HB{i+1}"
        content = content.replace(original, f"{prefix}{new_atom}{suffix}")

    # Fix HG atom name
    content = re.sub(
        r'^(ATOM\s+\d+\s+)HG(\s+SEC\s+49.*)$',
        r'\1HG1\2',
        content,
        flags=re.MULTILINE
    )

    with open(pdb_file, 'w') as file:
        file.write(content)
    
    print(f"Fixed file saved: {pdb_file}")

# Iterate through the mutation folders
for mutation in mutation_folders:
    minim_pdb_path = os.path.join(base_directory, mutation, "minim", "minim.pdb")
    if os.path.exists(minim_pdb_path):
        fix_pdb_file(minim_pdb_path)
    else:
        print(f"File not found: {minim_pdb_path}")
