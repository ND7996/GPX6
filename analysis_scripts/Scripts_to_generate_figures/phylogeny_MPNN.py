"""
GPX6 ProteinMPNN sequence analysis - Heatmap Only
JCIM compliant: Arial 8pt, black text, frameless legend
"""

import sys
import os

ACS_PATH = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_jcim_style
set_jcim_style()

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black"
})

from pathlib import Path
from itertools import combinations
import math
import warnings

import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════
# PARAMETERS
# ═══════════════════════════════════════════════

TOP_PER_SPECIES     = 1671
HEATMAP_PER_SPECIES = 100

# ═══════════════════════════════════════════════
# COLORS
# ═══════════════════════════════════════════════

C_HUMAN    = "#E41A1C"
C_MOUSE    = "#4DAF4A"
C_ANCESTOR = "#377EB8"
COLOR_OF   = {"human": C_HUMAN, "mouse": C_MOUSE, "ancestor": C_ANCESTOR}

heatmap_cmap = LinearSegmentedColormap.from_list(
    "species_cmap", [C_ANCESTOR, C_MOUSE, C_HUMAN], N=256
)

# ═══════════════════════════════════════════════
# PATHS  ←  print these at startup so you can verify
# ═══════════════════════════════════════════════

# Adjust these three paths to match your actual file locations.
# Run the script once; the startup block will print whether each file exists.

FASTA_ANCESTOR = Path(r"D:\PhD_Thesis\analysis\alignment\ancestor_25.fasta")
FASTA_HUMAN    = Path(r"D:\PhD_Thesis\MPNN\Human\GPX6MPNN_generated\all_generated_sequences.fasta")
FASTA_MOUSE    = Path(r"D:\PhD_Thesis\MPNN\Mouse\GPX6MPNN_generated\all_generated_sequences.fasta")

OUTPUT_DIR = Path(r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES\Phylogeny")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def check_paths():
    ok = True
    for label, p in [("ancestor", FASTA_ANCESTOR),
                     ("human",    FASTA_HUMAN),
                     ("mouse",    FASTA_MOUSE)]:
        exists = p.exists()
        status = "OK " if exists else "MISSING"
        print(f"  [{status}] {label:8s} → {p}")
        if not exists:
            ok = False
    return ok

# ═══════════════════════════════════════════════
# ALIGNER
# ═══════════════════════════════════════════════

ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score   = -10
ALIGNER.extend_gap_score = -0.5

def normsim(a, b):
    """
    Normalised BLOSUM62 similarity in [0, 1].
    Guards against:
      - negative self-scores  (unusual AAs or very short seqs)
      - zero denominator      (empty sequences)
      - domain error in sqrt  (product of negatives)
    """
    try:
        ab = float(ALIGNER.score(a, b))
        aa = float(ALIGNER.score(a, a))
        bb = float(ALIGNER.score(b, b))
        product = aa * bb
        if product <= 0:          # one or both self-scores negative/zero
            return 0.0
        return max(0.0, min(1.0, ab / math.sqrt(product)))
    except Exception:
        return 0.0

# ═══════════════════════════════════════════════
# LOAD SEQUENCES
# ═══════════════════════════════════════════════

def load_sequences():
    """
    Returns parallel lists (names, seqs, species).
    Labels are tracked independently of SeqIO record fields.
    Only sequences with length > 0 are kept.
    """
    names, seqs, species = [], [], []

    for fasta_path, label in [
        (FASTA_ANCESTOR, "ancestor"),
        (FASTA_HUMAN,    "human"),
        (FASTA_MOUSE,    "mouse"),
    ]:
        if not fasta_path.exists():
            print(f"  [SKIP] {label}: file not found")
            continue

        loaded = 0
        skipped = 0
        for record in SeqIO.parse(fasta_path, "fasta"):
            if loaded >= TOP_PER_SPECIES:
                break
            seq = str(record.seq).replace("U", "C").strip()
            if len(seq) == 0:
                skipped += 1
                continue
            names.append(record.id)
            seqs.append(seq)
            species.append(label)   # ← stored separately, never overwritten
            loaded += 1

        note = f" ({skipped} empty skipped)" if skipped else ""
        print(f"  {label:8s}: {loaded} sequences{note}")

    return names, seqs, species

# ═══════════════════════════════════════════════
# SIMILARITY MATRIX  (with progress and error safety)
# ═══════════════════════════════════════════════

def compute_similarity(seqs):
    n     = len(seqs)
    sim   = np.zeros((n, n))
    total = n * (n - 1) // 2
    done  = 0
    errs  = 0

    for i, j in combinations(range(n), 2):
        try:
            s = normsim(seqs[i], seqs[j])
        except Exception:
            s = 0.0
            errs += 1
        sim[i, j] = sim[j, i] = s
        done += 1
        if done % 10_000 == 0:
            pct = 100 * done / total
            print(f"  pairs: {done:,}/{total:,}  ({pct:.1f}%)", end="\r", flush=True)

    np.fill_diagonal(sim, 1.0)
    print(f"\n  done — {errs} pairs had errors (set to 0)")
    return sim

# ═══════════════════════════════════════════════
# DIVERSITY SAMPLING
# ═══════════════════════════════════════════════

def diversity_sample(group_indices, sim, k):
    """
    Greedy max-min diversity sampling.
    Iteratively picks the candidate most dissimilar to
    everything already chosen, so k representatives span
    the full sequence space.
    """
    idx = list(group_indices)
    if len(idx) <= k:
        return idx

    sub = sim[np.ix_(idx, idx)]

    selected = [0]
    min_sim  = sub[0].copy()   # similarity of each candidate to nearest selected

    for _ in range(k - 1):
        candidates = [i for i in range(len(idx)) if i not in selected]
        best = candidates[int(np.argmin(min_sim[candidates]))]
        selected.append(best)
        min_sim = np.minimum(min_sim, sub[best])

    return [idx[s] for s in selected]


def get_heatmap_indices(species, sim):
    all_idx = list(range(len(species)))
    groups  = {
        "human":    [i for i in all_idx if species[i] == "human"],
        "mouse":    [i for i in all_idx if species[i] == "mouse"],
        "ancestor": [i for i in all_idx if species[i] == "ancestor"],
    }

    selected = []
    for label, group in groups.items():
        n = len(group)
        if n == 0:
            print(f"  [WARN] '{label}': 0 sequences — skipped")
            continue
        if n == 1:
            # single sequence (e.g. reconstructed ancestor): use it as-is
            print(f"  '{label}': 1 sequence — used directly (no sampling needed)")
            selected.extend(group)
            continue
        reps = diversity_sample(group, sim, HEATMAP_PER_SPECIES)
        print(f"  '{label}': {n} → {len(reps)} representative sequences")
        selected.extend(reps)

    if len(selected) < 2:
        raise ValueError(
            "Fewer than 2 sequences selected for the heatmap.\n"
            "Check that FASTA files exist and paths are correct (run check_paths())."
        )

    return selected

# ═══════════════════════════════════════════════
# HEATMAP
# ═══════════════════════════════════════════════

def plot_heatmap(sim, species, idx):
    sub  = sim[np.ix_(idx, idx)]
    sp   = [species[i] for i in idx]
    dist = 1.0 - sub

    order = leaves_list(linkage(pdist(dist), method="ward"))
    dist  = dist[order][:, order]
    sp    = [sp[i] for i in order]

    colors   = [COLOR_OF[s] for s in sp]
    fig_size = min(7.5, max(3.5, len(idx) / 6))

    cg = sns.clustermap(
        pd.DataFrame(dist),
        row_cluster=True,
        col_cluster=True,
        cmap=heatmap_cmap,
        figsize=(fig_size, fig_size),
        row_colors=colors,
        col_colors=colors,
        xticklabels=False,
        yticklabels=False,
        linewidths=0.3,
        cbar_kws={
            "label":  "Genetic Distance",
            "ticks":  [0, 0.25, 0.5, 0.75, 1],
            "aspect": 30,
            "shrink": 0.4,
        },
        cbar_pos=(1.02, 0.25, 0.02, 0.5),
    )

    cbar = cg.cax
    cbar.yaxis.label.set_color("black")
    cbar.yaxis.label.set_fontsize(8)
    cbar.tick_params(colors="black", labelsize=7, width=0.5, length=3)
    for s in cbar.spines.values():
        s.set_linewidth(0.3)

    for ax_name in ("ax_row_dendrogram", "ax_col_dendrogram"):
        ax = getattr(cg, ax_name, None)
        if ax:
            for line in ax.lines:
                line.set_linewidth(0.5)
                line.set_color("black")

    legend_elements = [
        mpatches.Patch(color=C_HUMAN,    label="Human"),
        mpatches.Patch(color=C_MOUSE,    label="Mouse"),
        mpatches.Patch(color=C_ANCESTOR, label="Ancestor"),
    ]
    legend = cg.fig.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(0.02, 0.98),
        frameon=False,
        fontsize=8,
        handlelength=1.5,
        handleheight=1.0,
        borderaxespad=0.3,
    )
    for text in legend.get_texts():
        text.set_color("black")

    plt.tight_layout()

    for ext in ("pdf", "png", "tiff"):
        out = OUTPUT_DIR / f"mpnn_heatmap.{ext}"
        cg.savefig(out, dpi=600, bbox_inches="tight")
        print(f"  saved: {out}")

    plt.close()

# ═══════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════

def main():
    print("=" * 60)
    print("GPX6 ProteinMPNN Analysis - JCIM Heatmap")
    print("=" * 60)

    print("\n[0/3] Checking file paths ...")
    all_ok = check_paths()
    if not all_ok:
        print("\n  Fix the MISSING paths above before continuing.")
        print("  The script will run with whatever files are present,")
        print("  but missing species will be absent from the heatmap.\n")

    print("\n[1/3] Loading sequences ...")
    names, seqs, species = load_sequences()
    print(f"  total: {len(seqs)} sequences across {len(set(species))} species groups")

    if len(seqs) < 2:
        raise ValueError("Need at least 2 sequences total. Check your FASTA paths.")

    print("\n[2/3] Computing pairwise similarity matrix ...")
    print(f"  ({len(seqs)} sequences → {len(seqs)*(len(seqs)-1)//2:,} pairs)")
    sim = compute_similarity(seqs)

    print("\n[3/3] Sampling representatives and plotting ...")
    idx = get_heatmap_indices(species, sim)
    print(f"  heatmap: {len(idx)} total sequences")
    plot_heatmap(sim, species, idx)

    print("\n✅ Done. Output saved to:")
    print(f"   {OUTPUT_DIR}")


if __name__ == "__main__":
    main()