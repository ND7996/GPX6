import os

# ── CONFIG ────────────────────────────────────────────────────────────────────
FASTA_PATH = r"D:\PhD_Thesis\analysis\alignment\all_sequences_for_selection.fasta"

# Positions to extract (1-based, using RAW sequence indexing - do NOT strip X)
MOUSE_POSITIONS = [54, 49, 24, 139, 47, 60, 74, 144, 99, 177, 178, 104, 4, 52, 87, 102, 107, 142, 181, 48, 143]
HUMAN_POSITIONS = [87, 99, 47, 143, 60, 104, 142, 139, 181, 52, 48, 144, 54, 177, 102, 24, 3, 173, 178, 74]

# Expected AAs from the reference table: {position: (human_AA, mouse_AA)}
# X in expected means the WT has an unknown residue there - shown as blank
EXPECTED_AA = {
    3:   ('N', 'K'),
    4:   ('R', 'S'),
    16:  ('I', 'V'),
    22:  ('L', 'N'),
    24:  ('L', 'I'),
    25:  ('N', 'D'),
    27:  ('E', 'G'),
    29:  ('Y', 'F'),
    30:  ('I', 'V'),
    31:  ('Q', 'N'),
    33:  ('K', 'Q'),
    35:  ('F', 'Y'),
    40:  ('V', 'I'),
    47:  ('A', 'S'),
    48:  ('Y', 'F'),
    49:  ('U', 'C'),
    52:  ('A', 'T'),
    54:  ('Q', 'T'),
    60:  ('A', 'T'),
    67:  ('N', 'P'),
    69:  ('G', 'N'),
    71:  ('I', 'T'),
    74:  ('A', 'G'),
    87:  ('T', 'K'),
    99:  ('C', 'R'),
    102: ('S', 'G'),
    104: ('F', 'Y'),
    107: ('S', 'N'),
    119: ('E', 'D'),
    120: ('K', 'N'),
    126: ('T', 'S'),
    137: ('D', 'E'),
    139: ('L', 'F'),
    142: ('S', 'P'),
    143: ('S', 'E'),
    144: ('Q', 'H'),
    148: ('E', 'D'),
    173: ('H', 'R'),
    177: ('Q', 'H'),
    178: ('A', 'T'),
    181: ('S', 'R'),
    182: ('T', 'I'),
    184: ('K', 'Q'),
    188: ('L', 'M'),
    192: ('K', 'N'),
    194: ('F', 'T'),
    195: ('N', 'S'),
}

# ── FASTA PARSER ──────────────────────────────────────────────────────────────
def parse_fasta(filepath):
    sequences = {}
    current_id = None
    current_seq = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)
    return sequences

# ── RESIDUE EXTRACTOR ─────────────────────────────────────────────────────────
def extract_residues(sequence, positions, species_key):
    """
    Index the RAW sequence directly (no X removal).
    The X characters in WT sequences act as alignment coordinate holders
    and are essential for correct 1-based position mapping.

    Rules:
      - If residue matches expected AA  → show it
      - If residue is X (unknown)       → blank (genuine ambiguity in WT)
      - If residue doesn't match        → blank (wrong/not this species' AA)
      - If position out of range        → blank
    """
    exp_col = 0 if species_key == 'human' else 1  # tuple index in EXPECTED_AA

    result = []
    for pos in positions:
        idx = pos - 1
        aa = sequence[idx] if idx < len(sequence) else ''

        if pos in EXPECTED_AA:
            expected = EXPECTED_AA[pos][exp_col]
            if aa == '' or aa == 'X':
                result.append('')          # unknown / out of range → blank
            elif aa == expected:
                result.append(aa)          # correct AA → show it
            else:
                result.append('')          # wrong AA → blank
        else:
            # Position not in expected table: show whatever is there
            # (but still blank X)
            result.append('' if aa == 'X' else aa)

    return result

# ── TABLE BUILDER ─────────────────────────────────────────────────────────────
def build_table(seqs, positions, species_key):
    header = ['Sequence'] + [str(p) for p in positions]
    rows = []
    for seq_id, seq in seqs.items():
        residues = extract_residues(seq, positions, species_key)
        rows.append([seq_id] + residues)
    return header, rows

# ── PRETTY PRINTER ────────────────────────────────────────────────────────────
def print_table(header, rows, title):
    col_w = [
        max(len(str(header[i])), max(len(str(r[i])) for r in rows))
        for i in range(len(header))
    ]
    sep = '+' + '+'.join('-' * (w + 2) for w in col_w) + '+'
    fmt = '|' + '|'.join(f' {{:<{w}}} ' for w in col_w) + '|'

    print(f'\n{"=" * len(sep)}')
    print(f' {title}')
    print(sep)
    print(fmt.format(*[str(h) for h in header]))
    print(sep)
    for row in rows:
        print(fmt.format(*[str(c) for c in row]))
    print(sep)

# ── FASTA EXPORT ─────────────────────────────────────────────────────────────
def save_fasta(header, rows, filepath):
    """
    Each sequence in the output FASTA is the extracted residues joined together.
    Blank positions (no match / unknown) are written as '-'.
    The header line is the original sequence ID.
    The first line of the file is a comment listing the positions in order.
    """
    positions = header[1:]  # position labels from header (skip 'Sequence')
    with open(filepath, 'w') as f:
        f.write(f'# Positions: {" ".join(positions)}\n')
        for row in rows:
            seq_id   = row[0]
            residues = row[1:]
            # blank → '-' so alignment length is constant
            seq_str  = ''.join(r if r else '-' for r in residues)
            f.write(f'>{seq_id}\n{seq_str}\n')
    print(f'  Saved → {filepath}')

# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    print(f'Reading: {FASTA_PATH}')
    sequences = parse_fasta(FASTA_PATH)
    print(f'Found {len(sequences)} sequences\n')

    # Sanity check: print lengths
    print(f'{"ID":<45} {"Raw len":>8}')
    print('-' * 55)
    for seq_id, seq in sequences.items():
        print(f'{seq_id:<45} {len(seq):>8}')

    # Classify sequences
    mouse_seqs    = {k: v for k, v in sequences.items() if 'MOUSE'    in k.upper()}
    human_seqs    = {k: v for k, v in sequences.items() if 'HUMAN'    in k.upper()}
    ancestor_seqs = {k: v for k, v in sequences.items() if 'ANCESTOR' in k.upper()}

    print(f'\nMouse sequences:    {len(mouse_seqs)}')
    print(f'Human sequences:    {len(human_seqs)}')
    print(f'Ancestor sequences: {len(ancestor_seqs)}')

    # Build tables
    mouse_header,    mouse_rows    = build_table(mouse_seqs,    MOUSE_POSITIONS, 'mouse')
    human_header,    human_rows    = build_table(human_seqs,    HUMAN_POSITIONS, 'human')
    ancestor_header, ancestor_rows = build_table(ancestor_seqs, HUMAN_POSITIONS, 'human')

    # Print
    print_table(mouse_header,    mouse_rows,    'MOUSE — selected positions')
    print_table(human_header,    human_rows,    'HUMAN — selected positions')
    print_table(ancestor_header, ancestor_rows, 'ANCESTOR — selected positions (human layout)')

    # Save FASTAs next to the input FASTA
    # Ancestor uses human position layout so it is appended to the human file
    out_dir = os.path.dirname(FASTA_PATH)
    save_fasta(mouse_header, mouse_rows,                 os.path.join(out_dir, 'mouse_residues.fasta'))
    save_fasta(human_header, human_rows + ancestor_rows, os.path.join(out_dir, 'human_residues.fasta'))
    print('\nDone.')

if __name__ == '__main__':
    main()