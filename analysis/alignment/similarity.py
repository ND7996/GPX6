#!/usr/bin/env python3
"""
PHYLOGENY WITH RELIABLE PAIRWISE ALIGNMENT
Uses Bio.Align.PairwiseAligner (more stable than pairwise2)
"""

import os
import glob
import warnings
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

# Suppress warnings
warnings.filterwarnings("ignore", category=Warning)

def clean_sequence(seq_str):
    """Clean sequence to contain ONLY standard amino acids (no X, no gaps)."""
    standard_aa = set('ACDEFGHIKLMNPQRSTVWY')
    
    cleaned = []
    for aa in seq_str.upper():
        if aa in standard_aa:
            cleaned.append(aa)
        # Skip everything else (including X, -, and any other characters)
    
    return ''.join(cleaned)

def find_and_combine_sequences():
    """Find all FASTA files and combine sequences with unique IDs."""
    
    print("Finding FASTA files...")
    
    # Your directories
    human_dir = r"D:\PhD_Thesis\GPX6\analysis\alignment\HUMAN"
    mouse_dir = r"D:\PhD_Thesis\GPX6\analysis\alignment\MOUSE"
    jordi_file = r"D:\PhD_Thesis\GPX6\analysis\alignment\jordi\GPX6_inferred_nodes_protein.fasta"
    
    all_sequences = []
    id_counter = {}
    
    def process_directory(directory, label):
        """Process all FASTA files in a directory."""
        if not os.path.exists(directory):
            print(f"  WARNING: Directory not found: {directory}")
            return
        
        files = glob.glob(os.path.join(directory, "*.fasta")) + \
                glob.glob(os.path.join(directory, "*.fas")) + \
                glob.glob(os.path.join(directory, "*.fa"))
        
        print(f"  Found {len(files)} files in {label}")
        
        for file_path in files:
            try:
                for record in SeqIO.parse(file_path, "fasta"):
                    base_id = record.id.replace(">", "").replace("|", "_").replace(" ", "_")
                    base_id = base_id.replace("(", "_").replace(")", "_").replace(":", "_")
                    base_id = base_id.replace("[", "_").replace("]", "_").replace(",", "_")
                    
                    if base_id.startswith("JORDI_"):
                        base_id = base_id[6:]
                    
                    # Make unique
                    if base_id in id_counter:
                        id_counter[base_id] += 1
                        unique_id = f"{base_id}_{id_counter[base_id]}"
                    else:
                        id_counter[base_id] = 1
                        unique_id = base_id
                    
                    # Clean sequence
                    cleaned_seq = clean_sequence(str(record.seq))
                    
                    record.id = unique_id
                    record.seq = Seq(cleaned_seq)
                    record.description = ""
                    all_sequences.append(record)
            except Exception as e:
                print(f"  WARNING: Could not read {file_path}: {e}")
    
    # Process all directories
    process_directory(human_dir, "HUMAN")
    process_directory(mouse_dir, "MOUSE")
    
    # Process Jordi file
    if os.path.exists(jordi_file):
        try:
            for record in SeqIO.parse(jordi_file, "fasta"):
                base_id = record.id.replace(">", "").replace("|", "_").replace(" ", "_")
                base_id = base_id.replace("(", "_").replace(")", "_").replace(":", "_")
                base_id = base_id.replace("[", "_").replace("]", "_").replace(",", "_")
                
                if base_id.startswith("JORDI_"):
                    base_id = base_id[6:]
                
                if base_id in id_counter:
                    id_counter[base_id] += 1
                    unique_id = f"{base_id}_{id_counter[base_id]}"
                else:
                    id_counter[base_id] = 1
                    unique_id = base_id
                
                # Clean sequence
                cleaned_seq = clean_sequence(str(record.seq))
                
                record.id = unique_id
                record.seq = Seq(cleaned_seq)
                record.description = ""
                all_sequences.append(record)
        except Exception as e:
            print(f"  WARNING: Could not read {jordi_file}: {e}")
    
    if not all_sequences:
        print("ERROR: No sequences found!")
        return None
    
    print(f"SUCCESS: Total sequences: {len(all_sequences)}")
    
    # Save combined FASTA
    SeqIO.write(all_sequences, "COMBINED_SEQUENCES.fasta", "fasta")
    print("Saved: COMBINED_SEQUENCES.fasta")
    
    return all_sequences

def calculate_similarity_with_reliable_aligner(sequences):
    """
    Calculate similarity using Bio.Align.PairwiseAligner.
    This is more stable than pairwise2.
    """
    
    print("\nCalculating similarity using pairwise alignment...")
    print("  Using Bio.Align.PairwiseAligner (stable and reliable)")
    
    n = len(sequences)
    
    if n > 100:
        print(f"\n  WARNING: You have {n} sequences!")
        print(f"  This will require {n*(n-1)//2} pairwise alignments")
        response = input("\n  Continue anyway? (yes/no): ")
        if response.lower() not in ['yes', 'y']:
            return None
    
    similarity_matrix = np.zeros((n, n))
    
    # Set up aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    
    # Try to use BLOSUM62
    try:
        blosum62 = substitution_matrices.load("BLOSUM62")
        aligner.substitution_matrix = blosum62
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        use_blosum = True
        print("  Using BLOSUM62 scoring matrix")
    except:
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        use_blosum = False
        print("  Using simple match/mismatch scoring")
    
    print(f"  Performing {n * (n-1) // 2} optimal pairwise alignments...")
    print("  This may take 10-20 minutes...")
    
    for i in range(n):
        seq1_str = str(sequences[i].seq).upper()
        
        for j in range(i, n):
            if i == j:
                similarity_matrix[i][j] = 1.0
            else:
                seq2_str = str(sequences[j].seq).upper()
                
                if len(seq1_str) == 0 or len(seq2_str) == 0:
                    similarity_matrix[i][j] = 0.0
                    similarity_matrix[j][i] = 0.0
                    continue
                
                try:
                    # Perform alignment
                    alignments = aligner.align(seq1_str, seq2_str)
                    
                    if alignments:
                        best_alignment = alignments[0]
                        aligned_seq1 = str(best_alignment[0])
                        aligned_seq2 = str(best_alignment[1])
                        
                        # Calculate identity
                        matches = 0
                        total_comparable = 0
                        blosum_score = 0
                        max_blosum = 0
                        
                        for aa1, aa2 in zip(aligned_seq1, aligned_seq2):
                            if aa1 != '-' and aa2 != '-':
                                total_comparable += 1
                                
                                if aa1 == aa2:
                                    matches += 1
                                
                                if use_blosum:
                                    try:
                                        blosum_score += blosum62[aa1, aa2]
                                        max_blosum += max(blosum62[aa1, aa1], blosum62[aa2, aa2])
                                    except:
                                        pass
                        
                        if total_comparable > 0:
                            identity = matches / total_comparable
                            
                            if use_blosum and max_blosum > 0:
                                blosum_sim = max(0, blosum_score / max_blosum)
                                combined = 0.7 * identity + 0.3 * blosum_sim
                            else:
                                combined = identity
                            
                            similarity_matrix[i][j] = combined
                            similarity_matrix[j][i] = combined
                        else:
                            similarity_matrix[i][j] = 0.0
                            similarity_matrix[j][i] = 0.0
                    else:
                        similarity_matrix[i][j] = 0.0
                        similarity_matrix[j][i] = 0.0
                        
                except Exception as e:
                    print(f"    Warning: Failed to align {sequences[i].id} vs {sequences[j].id}: {e}")
                    similarity_matrix[i][j] = 0.0
                    similarity_matrix[j][i] = 0.0
        
        if (i + 1) % 5 == 0 or i == n - 1:
            print(f"    Processed {i + 1}/{n} sequences")
    
    return similarity_matrix

def save_similarity_matrix_and_report(similarity_matrix, sequences):
    """Save similarity matrix and create detailed report."""
    
    print("\n  Saving similarity matrix and creating report...")
    
    n = len(sequences)
    seq_names = [seq.id for seq in sequences]
    
    # Save matrix
    np.savetxt("COMBINED_SIMILARITY_MATRIX.csv", similarity_matrix, 
               delimiter=",", fmt="%.6f",
               header=",".join(seq_names), comments='')
    print("  Saved: COMBINED_SIMILARITY_MATRIX.csv")
    
    # Create report
    with open("SIMILARITY_REPORT.txt", "w", encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("SEQUENCE SIMILARITY ANALYSIS REPORT\n")
        f.write("="*80 + "\n\n")
        
        f.write("METHODOLOGY:\n")
        f.write("-" * 80 + "\n")
        f.write("1. Optimal pairwise global alignment (Bio.Align.PairwiseAligner)\n")
        f.write("2. Identity = matches / comparable positions\n")
        f.write("3. BLOSUM62 similarity scoring\n")
        f.write("4. Combined metric: 70% identity + 30% BLOSUM62\n")
        f.write("\n")
        
        # Statistics
        upper_triangle = []
        for i in range(n):
            for j in range(i+1, n):
                upper_triangle.append(similarity_matrix[i][j])
        
        if upper_triangle:
            f.write("SIMILARITY STATISTICS:\n")
            f.write("-" * 80 + "\n")
            f.write(f"Total comparisons: {len(upper_triangle)}\n")
            f.write(f"Mean similarity: {np.mean(upper_triangle):.4f} ({np.mean(upper_triangle)*100:.2f}%)\n")
            f.write(f"Median similarity: {np.median(upper_triangle):.4f} ({np.median(upper_triangle)*100:.2f}%)\n")
            f.write(f"Min similarity: {np.min(upper_triangle):.4f} ({np.min(upper_triangle)*100:.2f}%)\n")
            f.write(f"Max similarity: {np.max(upper_triangle):.4f} ({np.max(upper_triangle)*100:.2f}%)\n")
            f.write(f"Std deviation: {np.std(upper_triangle):.4f}\n")
            f.write("\n")
        
        # Top similar pairs
        f.write("TOP 20 MOST SIMILAR PAIRS:\n")
        f.write("-" * 80 + "\n")
        
        pairs = []
        for i in range(n):
            for j in range(i+1, n):
                pairs.append((similarity_matrix[i][j], i, j))
        
        pairs.sort(reverse=True)
        
        for rank, (sim, i, j) in enumerate(pairs[:20], 1):
            f.write(f"{rank}. {seq_names[i]} <-> {seq_names[j]}\n")
            f.write(f"   Similarity: {sim:.4f} ({sim*100:.2f}%)\n\n")
        
        # Distribution
        f.write("\nSIMILARITY DISTRIBUTION:\n")
        f.write("-" * 80 + "\n")
        bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        bin_labels = ['0-20%', '20-40%', '40-60%', '60-80%', '80-100%']
        for i in range(len(bins)-1):
            count = sum(1 for s in upper_triangle if bins[i] <= s < bins[i+1])
            if i == len(bins)-2:
                count = sum(1 for s in upper_triangle if bins[i] <= s <= bins[i+1])
            pct = 100 * count / len(upper_triangle) if len(upper_triangle) > 0 else 0
            f.write(f"{bin_labels[i]:10s}: {count:4d} pairs ({pct:5.1f}%)\n")
    
    print("  Saved: SIMILARITY_REPORT.txt")

def build_phylogeny_from_similarity(similarity_matrix, sequences):
    """Build phylogenetic tree from similarity matrix using UPGMA."""
    
    print("\nBuilding phylogeny from similarity matrix...")
    
    n = len(sequences)
    
    # Convert similarity to distance
    print("  Converting similarity to phylogenetic distance...")
    distance_matrix = 1.0 - similarity_matrix
    distance_matrix = np.maximum(distance_matrix, 0.0)
    np.fill_diagonal(distance_matrix, 0.0)
    distance_matrix = distance_matrix + 1e-10
    
    # Save distance matrix
    seq_names = [seq.id for seq in sequences]
    np.savetxt("DISTANCE_MATRIX.csv", distance_matrix, delimiter=",", fmt="%.6f",
               header=",".join(seq_names), comments='')
    print("  Saved: DISTANCE_MATRIX.csv")
    
    # Build tree using UPGMA
    condensed_dist = squareform(distance_matrix, checks=False)
    print("  Building tree using UPGMA...")
    Z = linkage(condensed_dist, method='average')
    
    # Build Newick tree
    print("  Converting to Newick format...")
    newick_tree = build_newick_tree(Z, seq_names)
    
    total_tree_length = np.sum(Z[:, 2])
    print(f"  SUCCESS: Tree built!")
    print(f"  Total tree length: {total_tree_length:.6f}")
    
    return newick_tree, Z

def build_newick_tree(linkage_matrix, leaf_names):
    """Build Newick format tree from linkage matrix."""
    
    n = len(leaf_names)
    node_heights = {}
    
    for i in range(n):
        node_heights[i] = 0.0
    
    def get_newick_subtree(node_id):
        if node_id < n:
            return leaf_names[node_id]
        else:
            cluster_idx = node_id - n
            left_id = int(linkage_matrix[cluster_idx, 0])
            right_id = int(linkage_matrix[cluster_idx, 1])
            height = linkage_matrix[cluster_idx, 2]
            
            left_branch = height - node_heights.get(left_id, 0.0)
            right_branch = height - node_heights.get(right_id, 0.0)
            
            node_heights[node_id] = height
            
            left_subtree = get_newick_subtree(left_id)
            right_subtree = get_newick_subtree(right_id)
            
            return f"({left_subtree}:{left_branch:.6f},{right_subtree}:{right_branch:.6f})"
    
    root_id = 2 * n - 2
    tree_str = get_newick_subtree(root_id) + ";"
    
    return tree_str

def visualize_tree(linkage_matrix, sequences):
    """Create visual representation of phylogenetic tree."""
    
    print("\nCreating tree visualization...")
    
    try:
        plt.figure(figsize=(12, max(8, len(sequences) * 0.3)))
        
        dendrogram(
            linkage_matrix,
            labels=[seq.id for seq in sequences],
            orientation='right',
            distance_sort='ascending',
            show_leaf_counts=False
        )
        
        plt.xlabel('Phylogenetic Distance', fontsize=12)
        plt.title('Similarity-Based Phylogenetic Tree (UPGMA)', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        plt.savefig('PHYLOGENETIC_TREE.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("  Saved: PHYLOGENETIC_TREE.png")
        return True
    except Exception as e:
        print(f"  WARNING: Could not create visualization: {e}")
        return False

def create_figtree_instructions():
    """Create instructions for FigTree."""
    
    with open("FIGTREE_INSTRUCTIONS.txt", "w", encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("HOW TO VISUALIZE IN FIGTREE\n")
        f.write("="*80 + "\n\n")
        
        f.write("1. Download FigTree: http://tree.bio.ed.ac.uk/software/figtree/\n")
        f.write("2. Open FigTree\n")
        f.write("3. File -> Open -> Select: SIMILARITY_BASED_PHYLOGENY.nwk\n")
        f.write("4. Customize:\n")
        f.write("   - Check [Branch Labels] to show distances\n")
        f.write("   - Check [Scale Bar] to show scale\n")
        f.write("   - Adjust fonts in Appearance panel\n")
        f.write("5. Export: File -> Export Trees (PDF/PNG/SVG)\n")

def main():
    """Main workflow."""
    
    print("\n" + "="*80)
    print(" RELIABLE SIMILARITY-BASED PHYLOGENY")
    print("="*80)
    print("\nUses stable Bio.Align.PairwiseAligner for accurate results\n")
    
    # Step 1: Load sequences
    print("="*80)
    print("STEP 1: Loading Sequences")
    print("="*80)
    
    sequences = find_and_combine_sequences()
    if not sequences:
        return
    
    # Step 2: Calculate similarity
    print("\n" + "="*80)
    print("STEP 2: Calculating Similarity Matrix")
    print("="*80)
    
    similarity_matrix = calculate_similarity_with_reliable_aligner(sequences)
    if similarity_matrix is None:
        print("\nAlignment cancelled. Exiting.")
        return
    
    # Save matrix and report
    save_similarity_matrix_and_report(similarity_matrix, sequences)
    
    # Step 3: Build tree
    print("\n" + "="*80)
    print("STEP 3: Building Phylogenetic Tree")
    print("="*80)
    
    newick_tree, linkage_matrix = build_phylogeny_from_similarity(similarity_matrix, sequences)
    
    # Step 4: Visualize
    print("\n" + "="*80)
    print("STEP 4: Creating Visualizations")
    print("="*80)
    
    visualize_tree(linkage_matrix, sequences)
    create_figtree_instructions()
    
    # Save tree
    output_file = "SIMILARITY_BASED_PHYLOGENY.nwk"
    with open(output_file, "w", encoding='utf-8') as f:
        f.write(newick_tree + "\n")
    
    print(f"\nSaved: {output_file}")
    
    # Final summary
    print("\n" + "="*80)
    print(" SUCCESS - PHYLOGENY COMPLETE!")
    print("="*80)
    
    print(f"""
OUTPUT FILES:
   SIMILARITY_BASED_PHYLOGENY.nwk - Tree (OPEN IN FIGTREE!)
   PHYLOGENETIC_TREE.png - Preview
   COMBINED_SIMILARITY_MATRIX.csv - All similarities
   DISTANCE_MATRIX.csv - All distances
   SIMILARITY_REPORT.txt - Detailed statistics
   FIGTREE_INSTRUCTIONS.txt - How to use FigTree

This phylogeny is based on proper pairwise sequence alignment!
Open the .nwk file in FigTree to visualize.
""")

if __name__ == "__main__":
    try:
        from Bio import SeqIO, Align
        from Bio.Align import substitution_matrices
        from Bio.Seq import Seq
        import numpy as np
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"\nERROR: Missing package: {e}")
        print("\nInstall: pip install biopython numpy scipy matplotlib")
        exit(1)
    
    main()