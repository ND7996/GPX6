from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

# ── Sequences ─────────────────────────────────────────────
ANC_SEQ = "PQKMKMDCNKGVTGTIYEYGALTLNGEEYIQFKQYAGKHVLFVNVATYGLTAQYPELNALQEELKHFGVIVLGFPCNQFGKQEPGKNSEILSGLKYVRPGGGFVPNFQLFEKGDVNGEKEQKVFTFLKNSCPPTSDLLGSSSQLFWEPMKVHDIRWNFEKFLVGPDGVPVMRWFHRAPVSTVKSDILEYLKQF"
HUM_SEQ = "PQNRKVDNKGVTGTIYEYGALTLNGEEYIQFKQFAGKVLFVNVAAYLAAQYPELNALQEELKNFGVIVLAFPCNQFGKQEPGTNSEILLGLKYVCPGSGFVPSFQLFEKGDVNGEKEQKVFTFLKNSPPTSDLLGSSSQLFWEPMKVDIRWNFEKFLVGPDGVPVMWFQAPVSTVKSDILEYLKQFNT"
MOU_SEQ = "PQKSKVDNKGVTGTVYEYGANTIDGGEFVNFQQYAGKILFVNVASFCGLTATYPELNTLQEELKPFNVTVLGFPCNQFGKQEPGKNSEILLGLKYVRPGGGYVPNFQLFEKGDVNGDNEQKVFSFLKNSPPTSELFGSPELFWDPMKVDIRWNFEKFLVGPDGVPVMRWFTPVRIVQSDIMEYLNQTS"

# ── Create SeqRecords ─────────────────────────────────────
records = [
    SeqRecord(Seq(ANC_SEQ), id="ANC"),
    SeqRecord(Seq(HUM_SEQ), id="HUM"),
    SeqRecord(Seq(MOU_SEQ), id="MOU"),
]

# ── Simple progressive alignment (pairwise-based) ─────────
aligner = PairwiseAligner()
aligner.mode = 'global'

# Align ANC vs HUM
aln1 = aligner.align(records[0].seq, records[1].seq)[0]

# Extract aligned sequences
anc_aligned = aln1.aligned[0]
hum_aligned = aln1.aligned[1]

# Build strings with gaps
def build_alignment(seq, aligned_coords):
    aligned_seq = ""
    prev_end = 0
    for start, end in aligned_coords:
        aligned_seq += "-" * (start - prev_end)
        aligned_seq += str(seq[start:end])
        prev_end = end
    return aligned_seq

anc_str = build_alignment(records[0].seq, anc_aligned)
hum_str = build_alignment(records[1].seq, hum_aligned)

# Align result with MOU
aln2 = aligner.align(anc_str.replace("-", ""), records[2].seq)[0]

# ── Coloring scheme ───────────────────────────────────────
def color_residue(res):
    if res in "AVILMFWY":
        return "#8dd3c7"  # hydrophobic
    elif res in "STNQC":
        return "#80b1d3"  # polar
    elif res in "KRH":
        return "#fb8072"  # positive
    elif res in "DE":
        return "#fdb462"  # negative
    elif res in "GP":
        return "#ffffb3"  # gly/pro
    else:
        return "#ffffff"

# ── Write HTML output ─────────────────────────────────────
def write_html(seqs, names, filename="msa_colored.html"):
    with open(filename, "w") as f:
        f.write("<html><body><pre style='font-family: monospace;'>\n")
        for name, seq in zip(names, seqs):
            f.write(f"{name} ")
            for aa in seq:
                color = color_residue(aa)
                f.write(f"<span style='background-color:{color}'>{aa}</span>")
            f.write("\n")
        f.write("</pre></body></html>")

write_html([anc_str, hum_str, str(records[2].seq)], ["ANC", "HUM", "MOU"])

print("MSA saved as msa_colored.html")