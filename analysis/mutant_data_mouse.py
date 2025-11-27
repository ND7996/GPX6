import re
import pandas as pd
from docx import Document

def parse_docx_to_csv(docx_path, output_csv):
    print(f"Reading file: {docx_path}")

    doc = Document(docx_path)

    # Extract all text
    full_text = "\n".join([p.text for p in doc.paragraphs])

    # Split by Level0, Level1, etc.
    level_blocks = re.split(r"(Level\d+)", full_text)

    data_rows = []

    # Process every (level, block) pair
    for i in range(1, len(level_blocks), 2):
        level_name = level_blocks[i]
        block_text = level_blocks[i + 1]

        # Extract dictionaries {...}
        dicts = re.findall(r"\{(.*?)\}", block_text, re.DOTALL)

        for entry in dicts:
            cleaned = entry.replace("\n", "")

            kv = {}
            for part in cleaned.split(","):
                if ":" in part:
                    key, val = part.split(":", 1)
                    kv[key.strip(' "')] = val.strip(' "')

            data_rows.append({
                "level": level_name,
                "mutation": kv.get("mutation"),
                "dg_star": float(kv.get("dg_star", "nan")),
                "dg_star_error": float(kv.get("dg_star_error", "nan")),
                "dg0": float(kv.get("dg0", "nan")),
                "dg0_error": float(kv.get("dg0_error", "nan")),
            })

    df = pd.DataFrame(data_rows)
    df.to_csv(output_csv, index=False)

    print(f"Saved CSV: {output_csv}")
    return df


# ---- RUN SCRIPT ----
docx_path = r"D:\PhD_Thesis\GPX6\analysis\mouse_mutants.docx"
output_csv = r"D:\PhD_Thesis\GPX6\analysis\mouse_mutants_merged.csv"

parse_docx_to_csv(docx_path, output_csv)
