import pandas as pd
import json
import re
import os

# ------------------------------------------------------------------
# CHANGE ONLY THESE TWO LINES IF YOUR FILE IS NAMED DIFFERENTLY
# ------------------------------------------------------------------
input_file  = r"D:\PhD_Thesis\GPX6\analysis\mouse_mutants_extracted.csv"      # your current (broken) file
output_file = r"D:\PhD_Thesis\GPX6\analysis\mouse_mutants_FINAL_with_distances.csv"

# ------------------------------------------------------------------
# Distance dictionary + nice column names
# ------------------------------------------------------------------
distance_dict = {
    "G74A": (14.29, 15.75), "E143S": (15.0, 12.97), "F139L": (16.1, 16.37),
    "F48Y": (3.88, 4.6),    "G102S": (16.77, 19.14), "H144Q": (12.5, 9.46),
    "H177Q": (14.0, 15.71), "L24I": (18.69, 21.89),  "I24L": (18.69, 21.89),
    "K87T": (16.44, 12.51), "N107S": (20.0, 21.53),  "P142S": (15.5, 14.77),
    "R181S": (20.0, 23.95), "R99C": (9.5, 10.99),    "S47A": (5.59, 6.24),
    "S4R": (19.64, 22.97),  "T178A": (14.77, 18.99),"T52A": (5.2, 9.62),
    "T54Q": (10.0, 14.5),   "T60A": (19.05, 22.8),   "Y104F": (16.5, 18.27),
    "mousecys": (None, None), "C49U": (None, None)
}

# ------------------------------------------------------------------
# 1. Read the file exactly as it is (even if it is horribly split)
# ------------------------------------------------------------------
raw = pd.read_csv(input_file, header=None, dtype=str)   # force everything as string

# Combine every row into one long string
lines = raw.fillna('').astype(str).agg(' '.join, axis=1).tolist()

# ------------------------------------------------------------------
# 2. Parse levels and JSON objects correctly
# ------------------------------------------------------------------
records = []
current_level = "Level0"

for line in lines:
    line = line.strip()
    
    # Detect new level
    if line.startswith("Level"):
        current_level = line.split()[0]      # Level0, Level1, ...
        continue
    
    # Find JSON object
    json_match = re.search(r'\{.*?\}(?=,|$)', line)
    if not json_match:
        continue
        
    json_str = json_match.group(0)
    json_str = json_str.replace(' ', ' ').replace('True', 'true').replace('False', 'false')
    
    try:
        obj = json.loads(json_str)
        obj["level"] = current_level
        records.append(obj)
    except:
        # last-ditch cleanup for stubborn lines
        json_str = re.sub(r",\s*(#.*)?$", "", json_str)   # remove trailing comma + comments
        try:
            obj = json.loads(json_str)
            obj["level"] = current_level
            records.append(obj)
        except:
            print(f"Skipping unparsable line: {json_str[:100]}...")

# ------------------------------------------------------------------
# 3. Build final DataFrame
# ------------------------------------------------------------------
df = pd.DataFrame(records)

# ------------------------------------------------------------------
# 4. Add properly named distance columns
# ------------------------------------------------------------------
df["dist_to_Sec/Cys49"] = df["mutation"].map(lambda x: distance_dict.get(x, (None,None))[0])
df["dist_to_Gln83"]     = df["mutation"].map(lambda x: distance_dict.get(x, (None,None))[1])

# ------------------------------------------------------------------
# 5. Final column order
# ------------------------------------------------------------------
cols = ["level", "mutation", "dg_star", "dg_star_error", "dg0", "dg0_error",
        "dist_to_Sec/Cys49", "dist_to_Gln83"]
df = df[cols]

# ------------------------------------------------------------------
# 6. Save perfect CSV
# ------------------------------------------------------------------
df.to_csv(output_file, index=False, float_format="%.3f")
print(f"DONE! Perfect file created:\n{output_file}")
print(f"Total rows: {len(df)}")
print(f"Levels present: {sorted(df['level'].unique())}")
print("\nFirst 10 rows:")
print(df.head(10))