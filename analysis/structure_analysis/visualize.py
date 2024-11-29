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

# Cluster assignments (based on your provided clusters)
cluster_residues = {
    1: [3, 4, 173, 177, 178],
    2: [16, 30, 31, 33, 35],
    0: [22, 24, 25, 27, 29, 87],
    4: [40, 69, 71, 74, 107],
    9: [47, 48, 52, 54, 99],
    8: [60, 181, 182, 184, 188],
    5: [67, 192, 194, 195],
    6: [102, 104],
    7: [119, 120, 126, 148],
    3: [137, 139, 142, 143, 144]
}

# Assigning colors for each cluster
cluster_colors = {
    1: 'red',
    2: 'blue',
    0: 'green',
    4: 'yellow',
    9: 'cyan',
    8: 'magenta',
    5: 'orange',
    6: 'purple',
    7: 'lime',
    3: 'pink'
}

# Apply color and selection to each cluster
for cluster_id, residues in cluster_residues.items():
    # Create a selection string for the current cluster
    selection = ' or '.join([f'resid {res}' for res in residues])
    
    # Select and color the cluster
    cmd.select(f'cluster_{cluster_id}', selection)
    cmd.color(cluster_colors[cluster_id], f'cluster_{cluster_id}')

# Optionally, zoom into the structure to focus on the colored residues
cmd.zoom()

# Show the representation in a more readable format (sticks or spheres)
#cmd.show('sticks', 'all')
