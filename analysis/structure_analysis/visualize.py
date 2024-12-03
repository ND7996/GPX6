import pymol
from pymol import cmd

# Initialize PyMOL in script mode
pymol.finish_launching()

# Load the structure
cmd.load("/home/hp/nayanika/github/GPX6/prep_structures/original_mousecys.pdb")  # Replace with the correct path to your PDB file

# Set the background to white
cmd.bg_color('white')

# Color the entire structure grey
cmd.color('grey', 'all')

# Updated cluster assignments (based on your provided clusters)
cluster_residues = {
    0: [142, 143, 144],  # Cluster 0 (P, E, H)
    1: [16, 33, 35],     # Cluster 1 (V, Q, Y)
    2: [181, 182, 184, 188],  # Cluster 2 (R, I, Q, M)
    3: [47, 48, 99],     # Cluster 3 (S, F, R)
    4: [192, 194, 195],  # Cluster 4 (N, T, S)
    5: [119, 120, 126, 148],  # Cluster 5 (D, N, S, D)
    6: [177, 178],       # Cluster 6 (H, T)
    7: [24, 25, 27],     # Cluster 7 (I, D, G)
    8: [40, 71],         # Cluster 8 (I, T)
    9: [3, 4],           # Cluster 9 (K, S)
    10: [52, 54],        # Cluster 10 (T, T)
    11: [137, 139],      # Cluster 11 (E, F)
    12: [102],           # Cluster 12 (G)
    13: [67, 69],        # Cluster 13 (P, N)
    14: [60, 107],       # Cluster 14 (T, N)
    15: [22, 29, 30, 31],# Cluster 15 (N, F, V, N)
    16: [74],            # Cluster 16 (G)
    17: [173],           # Cluster 17 (R)
    18: [104],           # Cluster 18 (Y)
    19: [87],            # Cluster 19 (K)
}

# Assigning colors for each cluster
cluster_colors = {
    0: 'red',
    1: 'blue',
    2: 'green',
    3: 'yellow',
    4: 'cyan',
    5: 'magenta',
    6: 'orange',
    7: 'purple',
    8: 'lime',
    9: 'pink',
    10: 'brown',
    11: 'violet',
    12: 'grey',
    13: 'black',
    14: 'lightgreen',
    15: 'lightblue',
    16: 'lightyellow',
    17: 'lightpink',
    18: 'darkgreen',
    19: 'darkblue',
}

# Apply color and representation to each cluster
for cluster_id, residues in cluster_residues.items():
    print(f"Processing Cluster {cluster_id} with residues {residues}")  # Debug print statement
    
    # Create a selection string for the current cluster
    selection = ' or '.join([f'resid {res}' for res in residues])
    
    # Select and color the cluster
    cmd.select(f'cluster_{cluster_id}', selection)
    cmd.color(cluster_colors[cluster_id], f'cluster_{cluster_id}')
    
    # Show the sticks representation for the selected residues
    cmd.show('sticks', f'cluster_{cluster_id}')

# Optionally, zoom into the structure to focus on the colored residues
cmd.zoom()

# Export image with transparent background
output_path = "/home/hp/nayanika/github/GPX6/figures/visualize_clusters.png"
cmd.ray()  # Render the image
cmd.png(output_path, dpi=300, ray=1, quiet=1, transparent=True)

# Print confirmation message with the file path
print(f"Image saved to {output_path}")
