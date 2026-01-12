import glob, subprocess, numpy as np, pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
from matplotlib import colors

###########################################
# 1. Merge FASTA files
###########################################

out = open("all_sequences.fasta","w")

for f in glob.glob("/home/hp/nayanika/github/GPX6/analysis/alignment/MOUSE/*.fasta"):
    out.write(open(f).read())

for f in glob.glob("/home/hp/nayanika/github/GPX6/analysis/alignment/HUMAN/*.fasta"):
    out.write(open(f).read())

out.write(open("/home/hp/nayanika/github/GPX6/analysis/alignment/node_25.fasta").read())
out.close()

###########################################
# 2. Build ML phylogenetic tree
###########################################

subprocess.run([
    "iqtree2",
    "-s","all_sequences.fasta",
    "-m","LG+G",
    "-nt","AUTO",
    "-redo"
])

###########################################
# 3. Load and root tree on node_25
###########################################

t = Tree("all_sequences.fasta.treefile", format=1)

ancestor = [x for x in t.get_leaf_names() if "25" in x][0]
t.set_outgroup(t & ancestor)

###########################################
# 4. Compute distances from ancestor
###########################################

A = t & ancestor
dist = {leaf.name: A.get_distance(leaf) for leaf in t}

pd.DataFrame.from_dict(dist,orient="index",columns=["distance"]).to_csv("distances_from_node25.csv")

###########################################
# 5. Color scaling
###########################################

human = [x for x in dist if "HUMAN" in x.upper()]
mouse = [x for x in dist if "MOUSE" in x.upper()]

hvals = np.array([dist[x] for x in human])
mvals = np.array([dist[x] for x in mouse])

hmin,hmax = hvals.min(), hvals.max()
mmin,mmax = mvals.min(), mvals.max()

def scale(x,a,b):
    return (x-a)/(b-a) if b>a else 0.5

colors_map = {}

for k in dist:
    if k in human:
        s = scale(dist[k],hmin,hmax)
        colors_map[k] = colors.to_hex((0.6+0.4*s, 0.0, 0.0))   # dark→light red
    elif k in mouse:
        s = scale(dist[k],mmin,mmax)
        colors_map[k] = colors.to_hex((0.0, 0.4+0.6*(1-s), 0.0))  # dark→light green
    else:
        colors_map[k] = "black"

###########################################
# 6. Draw final tree
###########################################

ts = TreeStyle()
ts.show_leaf_name = False
ts.scale = 120

for leaf in t:
    ns = NodeStyle()
    ns["fgcolor"] = colors_map[leaf.name]
    ns["size"] = 12
    leaf.set_style(ns)
    leaf.add_face(TextFace(leaf.name, fsize=10, fgcolor=colors_map[leaf.name]), column=0)

t.render("GPX6_EVOLUTIONARY_TREE.png", tree_style=ts, dpi=300)

print("Done. Output files:")
print("  GPX6_EVOLUTIONARY_TREE.png")
print("  distances_from_node25.csv")
