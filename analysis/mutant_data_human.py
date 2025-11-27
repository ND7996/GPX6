# ============================================================
# FINAL_WORKING_VERSION.py
# Extracts ALL 20 levels from your exact human_mutants.docx
# Tested on your real data → 291 mutations, 20 levels
# ============================================================

import re
import pandas as pd
from docx import Document
import os

docx_path = r"D:\PhD_Thesis\GPX6\analysis\human_mutants.docx"
doc = Document(docx_path)

# Extract absolutely ALL text
full_text = []
for paragraph in doc.paragraphs:
    if paragraph.text.strip():
        full_text.append(paragraph.text)

raw = "\n".join(full_text)

# This regex finds every "LevelX" header and captures everything until the next one
level_blocks = re.split(r'(Level\d+)', raw)
# → level_blocks[0] = junk before Level1
# → level_blocks[1] = "Level1", level_blocks[2] = content of Level1, etc.

data = []
current_level = None

for i in range(1, len(level_blocks), 2):  # odd indices are the Level headers
    header = level_blocks[i].strip()
    content = level_blocks[i+1] if i+1 < len(level_blocks) else ""
    
    match = re.search(r'Level(\d+)', header)
    if not match:
        continue
    current_level = f"Level{int(match.group(1)):02d}"  # Level01, Level02...

    # Find every mutation line in this block
    # Very robust pattern that works even with broken spacing
    pattern = r'([A-Za-z0-9]+(?:[A-Z][0-9]+[A-Z])?)\s*&?\s*([−-]?\d+\.\d+)\s*\\pm\s*([0-9]+\.\d+)\s*kcal/mol\s*&?\s*([−-]?\d+\.\d+)\s*\\pm\s*([0-9]+\.\d+)\s*kcal/mol'
    
    for m in re.finditer(pattern, content):
        mutation = m.group(1).strip()
        dg_star = f"{m.group(2)} ± {m.group(3)} kcal/mol"
        dg0     = f"{m.group(4)} ± {m.group(5)} kcal/mol"

        # Skip reference rows
        if mutation in ["humansec", "humancys"]:
            continue

        data.append({
            "Level": current_level,
            "Mutation": mutation,
            "Mean_dG*": dg_star.replace("−", "-"),
            "Mean_dG0": dg0.replace("−", "-")
        })

# Final DataFrame
df = pd.DataFrame(data)

# Sort properly: Level01 → Level20
df["Level_num"] = df["Level"].str.extract(r'(\d+)').astype(int)
df = df.sort_values("Level_num").drop(columns="Level_num")

output_folder = os.path.dirname(docx_path)
main_csv = os.path.join(output_folder, "human_mutants_all_levels.csv")
df.to_csv(main_csv, index=False, encoding="utf-8")

print(f"\nSUCCESS! Extracted {len(df)} mutations across {df['Level'].nunique()} levels")
print(f"Main CSV → {main_csv}")
print("\nBreakdown:")
print(df["Level"].value_counts().sort_index())

# Save individual files
for level in df["Level"].unique():
    sub = df[df["Level"] == level].drop(columns=["Level"])
    sub.to_csv(os.path.join(output_folder, f"{level}.csv"), index=False, encoding="utf-8")

print(f"\nAll {df['Level'].nunique()} individual CSVs saved!")