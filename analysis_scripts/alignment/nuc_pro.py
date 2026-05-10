# Standard genetic code
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def codon_alignment_to_protein_fasta(input_txt, output_fasta):
    with open(input_txt, "r") as infile, open(output_fasta, "w") as outfile:
        for line in infile:
            line = line.rstrip()
            if not line:
                continue

            parts = line.split()
            name_parts = []
            codons = []

            for p in parts:
                if set(p) <= set("ACGT-"):
                    codons.append(p)
                else:
                    name_parts.append(p)

            if not codons:
                continue

            name = "_".join(name_parts).replace("#", "")
            protein = []

            for codon in codons:
                if codon == "---":
                    protein.append("-")
                else:
                    protein.append(CODON_TABLE.get(codon, "X"))

            outfile.write(f">{name}\n")
            for i in range(0, len(protein), 60):
                outfile.write("".join(protein[i:i+60]) + "\n")


codon_alignment_to_protein_fasta(
    r"D:\PhD_Thesis\GPX6\analysis\alignment\jordi\GPX6_inferred_nodes",
    r"D:\PhD_Thesis\GPX6\analysis\alignment\jordi\GPX6_inferred_nodes_protein.fasta"
)
