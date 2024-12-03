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
    0: [22, 24, 25, 27, 29, 87],   # Cluster 0 (hydrophobic and polar residues)
    1: [3, 4, 173, 177, 178],      # Cluster 1 (charged and polar residues)
    2: [16, 30, 31, 33, 35],       # Cluster 2 (hydrophobic and polar residues)
    3: [137, 139, 142, 143, 144],  # Cluster 3 (acidic, aromatic, and imino residues)
    4: [40, 69, 71, 74, 107],      # Cluster 4 (hydrophobic, polar, and flexible residues)
    5: [67, 192, 194, 195],        # Cluster 5 (polar residues and proline)
    6: [102, 104],                 # Cluster 6 (flexible and aromatic residues)
    7: [119, 120, 126, 148],       # Cluster 7 (polar and acidic residues)
    8: [60, 181, 182, 184, 188],   # Cluster 8 (hydrophobic and charged residues)
    9: [47, 48, 52, 54, 99],       # Cluster 9 (polar and hydrophobic residues)
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
}

# Apply color and representation to each cluster
for cluster_id, residues in cluster_residues.items():
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
cmd.ray()  # Render the image
cmd.png("output_image.png", dpi=300, ray=1, quiet=1, transparent=True)
