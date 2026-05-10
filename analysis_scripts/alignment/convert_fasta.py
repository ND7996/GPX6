#!/usr/bin/env python3
"""
Convert nucleotide sequences to protein sequences (translation).
"""

from pathlib import Path

# Standard genetic code
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    '---': '-',  # Gap
}

def translate_sequence(nucleotide_seq):
    """
    Translate a nucleotide sequence to protein sequence.
    
    Args:
        nucleotide_seq: String of nucleotide sequence (with spaces removed)
        
    Returns:
        Protein sequence string
    """
    # Remove all spaces
    nucleotide_seq = nucleotide_seq.replace(' ', '')
    
    protein_seq = []
    
    # Translate in codons (3 nucleotides at a time)
    for i in range(0, len(nucleotide_seq) - 2, 3):
        codon = nucleotide_seq[i:i+3]
        
        # Handle gaps
        if '---' in codon or codon == '---':
            protein_seq.append('-')
        elif codon in CODON_TABLE:
            protein_seq.append(CODON_TABLE[codon])
        else:
            # Unknown codon, mark as X
            protein_seq.append('X')
    
    return ''.join(protein_seq)

def parse_alignment_file(input_file):
    """
    Parse the alignment file and extract species names and sequences.
    
    Args:
        input_file: Path to the input alignment file
        
    Returns:
        Dictionary mapping species names to nucleotide sequences
    """
    sequences = {}
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Split on whitespace
            parts = line.split(None, 1)
            if len(parts) == 2:
                species_name = parts[0]
                # Keep sequence with spaces for now
                sequence = parts[1]
                sequences[species_name] = sequence
    
    return sequences

def write_fasta(sequences, output_file, is_protein=False):
    """
    Write sequences to a FASTA format file.
    
    Args:
        sequences: Dictionary mapping species names to sequences
        output_file: Path to the output FASTA file
        is_protein: Boolean indicating if sequences are protein
    """
    with open(output_file, 'w') as f:
        for species_name, sequence in sequences.items():
            # Write header line with '>' prefix
            seq_type = "protein" if is_protein else "nucleotide"
            f.write(f">{species_name}\n")
            
            # Write sequence in lines of 60 characters (standard FASTA format)
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")

def write_individual_fasta_files(sequences, output_dir, is_protein=False):
    """
    Write each sequence to its own FASTA file.
    
    Args:
        sequences: Dictionary mapping species names to sequences
        output_dir: Directory to save individual FASTA files
        is_protein: Boolean indicating if sequences are protein
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    for species_name, sequence in sequences.items():
        # Create safe filename
        safe_name = species_name.replace(' ', '_').replace('#', 'node')
        suffix = "_protein.fasta" if is_protein else "_nucleotide.fasta"
        output_file = output_path / f"{safe_name}{suffix}"
        
        with open(output_file, 'w') as f:
            f.write(f">{species_name}\n")
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")
        
        print(f"Created: {output_file.name}")

def main():
    # Input file path
    input_file = r"D:\PhD_Thesis\GPX6\analysis\alignment\jordi\GPX6_inferred_nodes"
    
    # Parse the alignment file
    print("Parsing alignment file...")
    nucleotide_sequences = parse_alignment_file(input_file)
    print(f"Found {len(nucleotide_sequences)} sequences\n")
    
    # Translate all sequences to protein
    print("Translating nucleotide sequences to protein...")
    protein_sequences = {}
    for species_name, nuc_seq in nucleotide_sequences.items():
        protein_seq = translate_sequence(nuc_seq)
        protein_sequences[species_name] = protein_seq
        print(f"  {species_name}: {len(nuc_seq.replace(' ', ''))} nt -> {len(protein_seq)} aa")
    
    print()
    
    # Create output directory in the same location as input file
    input_path = Path(input_file)
    output_dir = input_path.parent / "protein_sequences"
    output_dir.mkdir(exist_ok=True)
    
    # Write all protein sequences to a single FASTA file
    protein_fasta = output_dir / "all_protein_sequences.fasta"
    write_fasta(protein_sequences, str(protein_fasta), is_protein=True)
    print(f"\nCreated combined protein FASTA file: {protein_fasta}")
    print(f"Total sequences: {len(protein_sequences)}")
    
    # Write individual protein FASTA files
    individual_dir = output_dir / "individual_proteins"
    print("\nCreating individual protein FASTA files...")
    write_individual_fasta_files(protein_sequences, str(individual_dir), is_protein=True)
    
    # Also save original nucleotide sequences
    print("\nSaving original nucleotide sequences...")
    nuc_dir = output_dir / "nucleotide_sequences"
    nuc_fasta = nuc_dir / "all_nucleotide_sequences.fasta"
    nuc_dir.mkdir(exist_ok=True)
    
    # Remove spaces from nucleotide sequences
    clean_nuc_sequences = {name: seq.replace(' ', '') for name, seq in nucleotide_sequences.items()}
    write_fasta(clean_nuc_sequences, str(nuc_fasta), is_protein=False)
    print(f"Created nucleotide FASTA file: {nuc_fasta}")
    
    print(f"\n{'='*60}")
    print(f"Done! All files saved to: {output_dir}")
    print(f"{'='*60}")
    print(f"Protein sequences:")
    print(f"  - Combined: all_protein_sequences.fasta")
    print(f"  - Individual: individual_proteins/")
    print(f"Nucleotide sequences:")
    print(f"  - Combined: nucleotide_sequences/all_nucleotide_sequences.fasta")

if __name__ == "__main__":
    main()