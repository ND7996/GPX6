#!/usr/bin/env python3
"""
Generate ROOTED Phylogenetic Tree with Ancestor_node25 as ROOT
Mouse_WT and Human_WT emerge FROM the ancestor
"""

import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
from pathlib import Path
import re


# Core sequences - The evolutionary reference points
CORE_SEQUENCES = {
    'Ancestor_node25': 'IRQFWASCLFPLFLVGFAQQTLKPQKMKMDCNKGVTGTIYEYGALTLNGEEYIQFKQYAGKHVLFVNVATYGLTAQYPELNALQEELKHFGVIVLGFPCNQFGKQEPGKNSEILSGLKYVRPGGGFVPNFQLFEKGDVNGEKEQKVFTFLKNSCPPTSDLLGSSSQLFWEPMKVHDIRWNFEKFLVGPDGVPVMRWFHRAPVSTVKSDILEYLKQF',
    'Mouse_WT': 'AQKLWGSCLFSLFMAALAQETLNPQKSKVDCNKGVTGTVYEYGANTIDGGEFVNFQQYAGKHILFVNVASFCGLTATYPELNTLQEELKPFNVTVLGFPCNQFGKQEPGKNSEILLGLKYVRPGGGYVPNFQLFEKGDVNGDNEQKVFSFLKNSCPPTSELFGSPEHLFWDPMKVHDIRWNFEKFLVGPDGVPVMRWFHHTPVRIVQSDIMEYLNQT',
    'Human_WT': 'FQQFQASCLVLFFLVGFAQQTLKPQNRKVDCNKGVTGTIYEYGALTLNGEEYIQFKQFAGKHVLFVNVAAYGLAAQYPELNALQEELKNFGVIVLAFPCNQFGKQEPGTNSEILLGLKYVCPGSGFVPSFQLFEKGDVNGEKEQKVFTFLKNSCPPTSDLLGSSSQLFWEPMKVHDIRWNFEKFLVGPDGVPVMHWFHQAPVSTVKSDILEYLKQF'
}


def read_fasta_file(filepath):
    """Read a FASTA file and return sequences dictionary."""
    sequences = {}
    current_name = None
    current_seq = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    if current_name:
                        sequences[current_name] = ''.join(current_seq)
                    current_name = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line.upper())
            
            if current_name:
                sequences[current_name] = ''.join(current_seq)
                
    except Exception as e:
        print(f"  ✗ Error reading {filepath}: {e}")
        return {}
    
    return sequences


def read_all_level_files(mouse_dir, human_dir):
    """Read all level files from both directories."""
    all_sequences = {}
    
    mouse_path = Path(mouse_dir)
    if mouse_path.exists():
        fasta_files = list(mouse_path.glob('*.fasta')) + list(mouse_path.glob('*.fa'))
        for fasta_file in sorted(fasta_files):
            seqs = read_fasta_file(fasta_file)
            level_match = re.search(r'level(\d+)', fasta_file.name.lower())
            
            for name, seq in seqs.items():
                if level_match:
                    level_num = int(level_match.group(1))
                    key = f"M2H_Level{level_num:02d}"
                else:
                    key = f"Mouse_{name}"
                all_sequences[key] = seq
            
            if seqs:
                print(f"  ✓ Mouse: {fasta_file.name} ({len(seqs)} sequences)")
    
    human_path = Path(human_dir)
    if human_path.exists():
        fasta_files = list(human_path.glob('*.fasta')) + list(human_path.glob('*.fa'))
        for fasta_file in sorted(fasta_files):
            seqs = read_fasta_file(fasta_file)
            level_match = re.search(r'level(\d+)', fasta_file.name.lower())
            
            for name, seq in seqs.items():
                if level_match:
                    level_num = int(level_match.group(1))
                    key = f"H2M_Level{level_num:02d}"
                else:
                    key = f"Human_{name}"
                all_sequences[key] = seq
            
            if seqs:
                print(f"  ✓ Human: {fasta_file.name} ({len(seqs)} sequences)")
    
    return all_sequences


def calculate_identity(seq1, seq2):
    """Calculate pairwise sequence identity (%)."""
    matches = 0
    total = 0
    
    for aa1, aa2 in zip(seq1, seq2):
        if aa1 in '-XU*' or aa2 in '-XU*':
            continue
        total += 1
        if aa1 == aa2:
            matches += 1
    
    return (matches / total * 100) if total > 0 else 0.0


def calculate_distance_matrix(sequences_dict):
    """Calculate pairwise distance matrix from sequences."""
    seq_names = list(sequences_dict.keys())
    n = len(seq_names)
    
    print(f"  Calculating {n}x{n} distance matrix...")
    distance_matrix = np.zeros((n, n))
    
    for i, name1 in enumerate(seq_names):
        for j, name2 in enumerate(seq_names):
            if i == j:
                distance_matrix[i, j] = 0.0
            else:
                identity = calculate_identity(sequences_dict[name1], sequences_dict[name2])
                distance_matrix[i, j] = 100.0 - identity
    
    return distance_matrix, seq_names


def create_rooted_tree_manual(sequences_dict, start_node_num=23):
    """
    Create a properly rooted tree with Ancestor_node25 as ROOT.
    
    Tree structure:
        Ancestor_node25 (root)
        ├─ Mouse_WT + M2H levels
        └─ Human_WT + H2M levels
    """
    
    print("\n[Building ROOTED Tree with Ancestor as ROOT]")
    
    # Separate sequences into groups
    ancestor_seq = sequences_dict.get('Ancestor_node25')
    mouse_wt_seq = sequences_dict.get('Mouse_WT')
    human_wt_seq = sequences_dict.get('Human_WT')
    
    # Get M2H and H2M sequences
    m2h_seqs = {k: v for k, v in sequences_dict.items() if 'M2H' in k}
    h2m_seqs = {k: v for k, v in sequences_dict.items() if 'H2M' in k}
    
    print(f"  • Ancestor: Ancestor_node25")
    print(f"  • Mouse lineage: Mouse_WT + {len(m2h_seqs)} M2H levels")
    print(f"  • Human lineage: Human_WT + {len(h2m_seqs)} H2M levels")
    
    # Build Mouse lineage subtree (Mouse_WT + M2H)
    mouse_lineage = {'Mouse_WT': mouse_wt_seq}
    mouse_lineage.update(m2h_seqs)
    
    print(f"\n  Building Mouse lineage subtree...")
    dist_matrix_mouse, names_mouse = calculate_distance_matrix(mouse_lineage)
    condensed_mouse = squareform(dist_matrix_mouse, checks=False)
    linkage_mouse = linkage(condensed_mouse, method='average')
    tree_mouse = to_tree(linkage_mouse, rd=False)
    
    # Build Human lineage subtree (Human_WT + H2M)
    human_lineage = {'Human_WT': human_wt_seq}
    human_lineage.update(h2m_seqs)
    
    print(f"  Building Human lineage subtree...")
    dist_matrix_human, names_human = calculate_distance_matrix(human_lineage)
    condensed_human = squareform(dist_matrix_human, checks=False)
    linkage_human = linkage(condensed_human, method='average')
    tree_human = to_tree(linkage_human, rd=False)
    
    # Convert to Newick strings
    node_counter = [start_node_num]
    
    def to_newick(node, names):
        if node.is_leaf():
            return names[node.id]
        left = to_newick(node.left, names)
        right = to_newick(node.right, names)
        label = node_counter[0]
        node_counter[0] += 1
        return f"({left},{right}){label}"
    
    mouse_subtree = to_newick(tree_mouse, names_mouse)
    human_subtree = to_newick(tree_human, names_human)
    
    # Combine into rooted tree with Ancestor at root
    root_node = node_counter[0]
    rooted_tree = f"(Ancestor_node25,({mouse_subtree},{human_subtree}){root_node}){root_node + 1};"
    
    print(f"\n  ✓ Tree rooted at: Ancestor_node25")
    print(f"  ✓ Root node: {root_node + 1}")
    print(f"  ✓ Node range: {start_node_num} to {node_counter[0] + 1}")
    
    return rooted_tree


def main():
    """Main function."""
    
    print("="*90)
    print("ROOTED PHYLOGENETIC TREE GENERATOR")
    print("Ancestor_node25 → Mouse_WT & Human_WT → Mutation Levels")
    print("="*90)
    
    # Configuration
    mouse_dir = r"D:\PhD_Thesis\GPX6\analysis\alignment\MOUSE"
    human_dir = r"D:\PhD_Thesis\GPX6\analysis\alignment\HUMAN"
    
    # Output file
    output_newick = "GPX6_ROOTED_tree.nwk"
    
    print("\n[1] Reading sequences...")
    level_sequences = read_all_level_files(mouse_dir, human_dir)
    
    # Combine all sequences
    all_sequences = CORE_SEQUENCES.copy()
    all_sequences.update(level_sequences)
    
    print(f"\n  Total sequences: {len(all_sequences)}")
    print(f"    • Ancestor: 1")
    print(f"    • Mouse_WT + Human_WT: 2")
    print(f"    • Mutation levels: {len(level_sequences)}")
    
    # Create rooted tree
    print("\n[2] Building rooted phylogenetic tree...")
    rooted_tree = create_rooted_tree_manual(all_sequences, start_node_num=23)
    
    # Save tree
    with open(output_newick, 'w') as f:
        f.write(rooted_tree + "\n")
    
    print(f"\n[3] Tree saved to: {output_newick}")
    
    # Create readable version
    readable = rooted_tree.replace(',', ',\n  ')
    with open("GPX6_ROOTED_tree_readable.txt", 'w') as f:
        f.write(readable)
    
    print("\n" + "="*90)
    print("✅ ROOTED TREE COMPLETE!")
    print("="*90)
    
    print("\n📊 Tree Structure:")
    print("  Ancestor_node25 (ROOT)")
    print("  ├─ Mouse_WT")
    print("  │  └─ M2H_Level01-20 (Mouse→Human mutations)")
    print("  └─ Human_WT")
    print("     └─ H2M_Level01-19 (Human→Mouse mutations)")
    
    print(f"\n📄 Output Files:")
    print(f"  • {output_newick} - Use this in FigTree")
    print(f"  • GPX6_ROOTED_tree_readable.txt - Human-readable version")
    
    print("\n💡 In FigTree:")
    print("  1. File → Open → GPX6_ROOTED_tree.nwk")
    print("  2. ☑ Tip Labels (to see sequence names)")
    print("  3. ☑ Node Labels → Display: 'label'")
    print("  4. You'll see Ancestor_node25 as the root!")
    
    print(f"\n🌳 Tree Preview (first 200 chars):")
    print("  " + rooted_tree[:200] + "...")
    print()


if __name__ == "__main__":
    main()