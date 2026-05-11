from pathlib import Path
import csv, math, re
import matplotlib.pyplot as plt
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────
HUMAN_FASTA   = Path(r"D:\PhD_Thesis\MPNN\results\proteinmpnn_all\proteinmpnn_outputs\combined_HUMAN_labeled.fasta")
MOUSE_FASTA   = Path(r"D:\PhD_Thesis\MPNN\results\proteinmpnn_all\proteinmpnn_outputs\combined_MOUSE_labeled.fasta")
OUTPUT_PREFIX = Path("results/gpx6_mpnn_ternary")

# ── Reference sequences (original GPX6 WTs) ───────────────────────────────────
ANC_SEQ = "PQKMKMDCNKGVTGTIYEYGALTLNGEEYIQFKQYAGKHVLFVNVATYGLTAQYPELNALQEELKHFGVIVLGFPCNQFGKQEPGKNSEILSGLKYVRPGGGFVPNFQLFEKGDVNGEKEQKVFTFLKNSCPPTSDLLGSSSQLFWEPMKVHDIRWNFEKFLVGPDGVPVMRWFHRAPVSTVKSDILEYLKQF"
HUM_SEQ = "PQNRKVDNKGVTGTIYEYGALTLNGEEYIQFKQFAGKVLFVNVAAYLAAQYPELNALQEELKNFGVIVLAFPCNQFGKQEPGTNSEILLGLKYVCPGSGFVPSFQLFEKGDVNGEKEQKVFTFLKNSPPTSDLLGSSSQLFWEPMKVDIRWNFEKFLVGPDGVPVMWFQAPVSTVKSDILEYLKQFNT"
MOU_SEQ = "PQKSKVDNKGVTGTVYEYGANTIDGGEFVNFQQYAGKILFVNVASFCGLTATYPELNTLQEELKPFNVTVLGFPCNQFGKQEPGKNSEILLGLKYVRPGGGYVPNFQLFEKGDVNGDNEQKVFSFLKNSPPTSELFGSPELFWDPMKVDIRWNFEKFLVGPDGVPVMRWFTPVRIVQSDIMEYLNQTS"
REF_SEQS   = [ANC_SEQ, HUM_SEQ, MOU_SEQ]
REF_LABELS = ["Ancestor\n(node_25)", "Human\nGPX6 WT", "Mouse\nGPX6 WT"]
VERTICES   = np.array([[0.5, math.sqrt(3)/2], [0.0, 0.0], [1.0, 0.0]])

# ── BLOSUM62 affine NW ────────────────────────────────────────────────────────
_B62 = {
    ('A','A'):4,('A','R'):-1,('A','N'):-2,('A','D'):-2,('A','C'):0,('A','Q'):-1,('A','E'):-1,('A','G'):0,('A','H'):-2,('A','I'):-1,('A','L'):-1,('A','K'):-1,('A','M'):-1,('A','F'):-2,('A','P'):-1,('A','S'):1,('A','T'):0,('A','W'):-3,('A','Y'):-2,('A','V'):0,
    ('R','R'):5,('R','N'):-1,('R','D'):-2,('R','C'):-3,('R','Q'):1,('R','E'):0,('R','G'):-2,('R','H'):0,('R','I'):-3,('R','L'):-2,('R','K'):2,('R','M'):-1,('R','F'):-3,('R','P'):-2,('R','S'):-1,('R','T'):-1,('R','W'):-3,('R','Y'):-2,('R','V'):-3,
    ('N','N'):6,('N','D'):1,('N','C'):-3,('N','Q'):0,('N','E'):0,('N','G'):0,('N','H'):1,('N','I'):-3,('N','L'):-3,('N','K'):0,('N','M'):-2,('N','F'):-3,('N','P'):-2,('N','S'):1,('N','T'):0,('N','W'):-4,('N','Y'):-2,('N','V'):-3,
    ('D','D'):6,('D','C'):-3,('D','Q'):0,('D','E'):2,('D','G'):-1,('D','H'):-1,('D','I'):-3,('D','L'):-4,('D','K'):-1,('D','M'):-3,('D','F'):-3,('D','P'):-1,('D','S'):0,('D','T'):-1,('D','W'):-4,('D','Y'):-3,('D','V'):-3,
    ('C','C'):9,('C','Q'):-3,('C','E'):-4,('C','G'):-3,('C','H'):-3,('C','I'):-1,('C','L'):-1,('C','K'):-3,('C','M'):-1,('C','F'):-2,('C','P'):-3,('C','S'):-1,('C','T'):-1,('C','W'):-2,('C','Y'):-2,('C','V'):-1,
    ('Q','Q'):5,('Q','E'):2,('Q','G'):-2,('Q','H'):0,('Q','I'):-3,('Q','L'):-2,('Q','K'):1,('Q','M'):0,('Q','F'):-3,('Q','P'):-1,('Q','S'):0,('Q','T'):-1,('Q','W'):-2,('Q','Y'):-1,('Q','V'):-2,
    ('E','E'):5,('E','G'):-2,('E','H'):0,('E','I'):-3,('E','L'):-3,('E','K'):1,('E','M'):-2,('E','F'):-3,('E','P'):-1,('E','S'):0,('E','T'):-1,('E','W'):-3,('E','Y'):-2,('E','V'):-2,
    ('G','G'):6,('G','H'):-2,('G','I'):-4,('G','L'):-4,('G','K'):-2,('G','M'):-3,('G','F'):-3,('G','P'):-2,('G','S'):0,('G','T'):-2,('G','W'):-2,('G','Y'):-3,('G','V'):-3,
    ('H','H'):8,('H','I'):-3,('H','L'):-3,('H','K'):-1,('H','M'):-2,('H','F'):-1,('H','P'):-2,('H','S'):-1,('H','T'):-2,('H','W'):-2,('H','Y'):2,('H','V'):-3,
    ('I','I'):4,('I','L'):2,('I','K'):-1,('I','M'):1,('I','F'):0,('I','P'):-3,('I','S'):-2,('I','T'):-1,('I','W'):-3,('I','Y'):-1,('I','V'):3,
    ('L','L'):4,('L','K'):-2,('L','M'):2,('L','F'):0,('L','P'):-3,('L','S'):-2,('L','T'):-1,('L','W'):-2,('L','Y'):-1,('L','V'):1,
    ('K','K'):5,('K','M'):-1,('K','F'):-3,('K','P'):-1,('K','S'):0,('K','T'):-1,('K','W'):-3,('K','Y'):-2,('K','V'):-2,
    ('M','M'):5,('M','F'):0,('M','P'):-2,('M','S'):-1,('M','T'):-1,('M','W'):-1,('M','Y'):-1,('M','V'):1,
    ('F','F'):6,('F','P'):-4,('F','S'):-2,('F','T'):-2,('F','W'):1,('F','Y'):3,('F','V'):-1,
    ('P','P'):7,('P','S'):-1,('P','T'):-1,('P','W'):-4,('P','Y'):-3,('P','V'):-2,
    ('S','S'):4,('S','T'):1,('S','W'):-3,('S','Y'):-2,('S','V'):-2,
    ('T','T'):5,('T','W'):-2,('T','Y'):-2,('T','V'):0,
    ('W','W'):11,('W','Y'):2,('W','V'):-3,
    ('Y','Y'):7,('Y','V'):-1,('V','V'):4,
}
def _b(a,b):
    k=(a,b) if (a,b) in _B62 else (b,a)
    return _B62.get(k,-1)

def normsim(seqa, seqb):
    seqa=seqa.replace('U','C').replace('X','A')
    seqb=seqb.replace('U','C').replace('X','A')
    GO=-10; GE=-0.5; NEG=float('-inf')
    m,n=len(seqa),len(seqb)
    M=[[NEG]*(n+1) for _ in range(m+1)]
    Ix=[[NEG]*(n+1) for _ in range(m+1)]
    Iy=[[NEG]*(n+1) for _ in range(m+1)]
    M[0][0]=0
    for i in range(1,m+1): Ix[i][0]=GO+(i-1)*GE
    for j in range(1,n+1): Iy[0][j]=GO+(j-1)*GE
    for i in range(1,m+1):
        for j in range(1,n+1):
            M[i][j]=max(M[i-1][j-1],Ix[i-1][j-1],Iy[i-1][j-1])+_b(seqa[i-1],seqb[j-1])
            Ix[i][j]=max(M[i-1][j]+GO,Ix[i-1][j]+GE)
            Iy[i][j]=max(M[i][j-1]+GO,Iy[i][j-1]+GE)
    score=max(M[m][n],Ix[m][n],Iy[m][n])
    aa=sum(_b(c,c) for c in seqa)
    bb=sum(_b(c,c) for c in seqb)
    d=(aa*bb)**0.5
    return max(0.0,min(1.0,score/d)) if d else 0.0

def to_barycentric(seq):
    sims=np.array([normsim(seq,r) for r in REF_SEQS])
    t=sims.sum()
    return sims/t if t>0 else None

def parse_fasta(path):
    records=[]; hf=None; chunks=[]
    with path.open() as f:
        for raw in f:
            line=raw.strip()
            if not line: continue
            if line.startswith(">"):
                if hf:
                    sid=hf.split("|")[0].strip()
                    m=re.search(r"\bscore=([\d.]+)",hf)
                    records.append((sid,float(m.group(1)) if m else None,"".join(chunks).upper()))
                hf=line[1:]; chunks=[]
            else: chunks.append(line)
    if hf:
        sid=hf.split("|")[0].strip()
        m=re.search(r"\bscore=([\d.]+)",hf)
        records.append((sid,float(m.group(1)) if m else None,"".join(chunks).upper()))
    return records

# ── Load — HUM and MOU only, no templates ────────────────────────────────────
all_records=[]
for fp in [HUMAN_FASTA, MOUSE_FASTA]:
    all_records.extend(
        (sid,sc,seq) for sid,sc,seq in parse_fasta(fp)
        if "_template" not in sid.lower()
    )
print(f"Loaded {len(all_records)} designed sequences")

# ── Compute ───────────────────────────────────────────────────────────────────
rows=[]; skipped=0
for i,(sid,score,seq) in enumerate(all_records):
    if i%100==0: print(f"  {i}/{len(all_records)}...")
    lam=to_barycentric(seq)
    if lam is None: skipped+=1; continue
    xy=lam@VERTICES
    closest=["ANC","HUM","MOU"][int(np.argmax(lam))]
    sp="HUM" if sid.upper().startswith("HUM") else "MOU"
    rows.append({
        "sequence_id":sid, "mpnn_score":score or "",
        "species":sp, "closest_ref":closest,
        "lambda_ANC":round(float(lam[0]),5),
        "lambda_HUM":round(float(lam[1]),5),
        "lambda_MOU":round(float(lam[2]),5),
        "x":round(float(xy[0]),5),
        "y":round(float(xy[1]),5),
    })

print(f"Done. {len(rows)} sequences ({skipped} skipped)")

from collections import Counter
for sp in ["HUM","MOU"]:
    c=Counter(r["closest_ref"] for r in rows if r["species"]==sp)
    print(f"  {sp}: closest_ref counts = {dict(c)}")

OUTPUT_PREFIX.parent.mkdir(parents=True,exist_ok=True)
with open(f"{OUTPUT_PREFIX}_coords.csv","w",newline="") as f:
    import csv as _csv
    w=_csv.DictWriter(f,fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

# ── Plot ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({"font.family":"Arial","font.size":11,"figure.dpi":300,
                     "legend.frameon":True,"legend.edgecolor":"#cccccc",
                     "legend.framealpha":0.9})
C_HUM="#C0392B"; C_MOU="#2471A3"
CORNER_COLORS = [C_HUM, C_MOU, C_HUM]  # ANC=red, MOU=blue, HUM=red
CORNER_NAMES  = ["ANCESTOR_node_25", "MOUSE_GPX6_WT", "HUMAN_GPX6_WT"]

fig,ax=plt.subplots(figsize=(8,8))
ax.set_aspect("equal"); ax.axis("off")

# Triangle — white fill, black border
ax.add_patch(plt.Polygon(VERTICES, fill=True, facecolor="white",
             edgecolor="black", linewidth=2.0, zorder=0))

# Large triangle markers at each vertex
for (vx,vy),vc in zip(VERTICES, CORNER_COLORS):
    ax.scatter([vx],[vy], marker="^", s=300, color=vc,
               edgecolors="none", zorder=5)

# Corner labels — bold, coloured, outside triangle
CORNER_OFFSETS = [(0,0.055,"center","bottom"),
                  (-0.07,-0.055,"right","top"),
                  ( 0.07,-0.055,"left","top")]
for (vx,vy),(dx,dy,ha,va),lbl,vc in zip(VERTICES,CORNER_OFFSETS,CORNER_NAMES,CORNER_COLORS):
    ax.text(vx+dx, vy+dy, lbl, ha=ha, va=va,
            fontsize=11, fontweight="bold", color=vc)

# Scatter — circles only, colour by species
hum_r=[r for r in rows if r["species"]=="HUM"]
mou_r=[r for r in rows if r["species"]=="MOU"]
ax.scatter([r["x"] for r in mou_r],[r["y"] for r in mou_r],
           color=C_MOU, s=30, edgecolors="none", alpha=0.70, zorder=3,
           label="Mouse GPX6 variants")
ax.scatter([r["x"] for r in hum_r],[r["y"] for r in hum_r],
           color=C_HUM, s=30, edgecolors="none", alpha=0.70, zorder=3,
           label="Human GPX6 variants")

# Legend — top right, clean
ax.legend(loc="upper right", framealpha=0.9, fontsize=9,
          title="Sequence group", title_fontsize=9,
          markerscale=1.4, handletextpad=0.5, borderpad=0.7)

ax.set_xlim(-0.15,1.15); ax.set_ylim(-0.13,1.02)
fig.tight_layout()

out=Path(f"{OUTPUT_PREFIX}.png")
out.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(out,dpi=300,bbox_inches="tight")
plt.show()
print(f"Saved: {out}")
print(f"CSV:   {OUTPUT_PREFIX}_coords.csv")