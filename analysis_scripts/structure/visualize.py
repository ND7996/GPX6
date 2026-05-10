import pymol
from pymol import cmd

# Initialize PyMOL in script mode
pymol.finish_launching()

# Function to calculate the centroid of a list of residues
def calculate_centroid(residues):
    coords = []
    for res in residues:
        selection = f'resid {res}'
        # Debug: Check selection
        print(f"Processing residue: {res}, Selection: {selection}")
        try:
            model = cmd.get_model(selection)
            res_coords = [atom.coord for atom in model.atom]
            print(f"Coordinates for residue {res}: {res_coords}")
            coords.extend(res_coords)
        except Exception as e:
            print(f"Error fetching coordinates for residue {res}: {e}")

    if not coords:
        print("No coordinates were found! Check residue numbers or the loaded structure.")
        return None

    # Calculate centroid
    x, y, z = zip(*coords)
    return sum(x) / len(x), sum(y) / len(y), sum(z) / len(z)

# Load the structure
cmd.load("/home/hp/nayanika/github/GPX6/prep_structures/original_mousecys.pdb")  # Replace with the correct path to your PDB file

# Define residues for clusters
cluster_1_residues = [3, 4, 16, 22, 25, 27, 29, 30, 31, 33, 35, 40, 60, 67, 69, 71, 107, 119, 120, 126, 137, 148, 173, 181, 182, 184, 188, 192, 194, 195]
cluster_2_residues = [24, 47, 48, 52, 54, 74, 87, 99, 102, 104, 139, 142, 143, 144, 177, 178]

# Calculate centroids
centroid_1 = calculate_centroid(cluster_1_residues)
centroid_2 = calculate_centroid(cluster_2_residues)

print(f"Centroid for Cluster 1: {centroid_1}")
print(f"Centroid for Cluster 2: {centroid_2}")

# Set background and initial colors
cmd.bg_color('white')
cmd.color('grey', 'all')

# Color and highlight residues for each cluster
def highlight_cluster(residues, color_name, cluster_name):
    selection = ' or '.join([f'resid {res}' for res in residues])
    cmd.select(cluster_name, selection)
    cmd.color(color_name, cluster_name)
    cmd.show('sticks', cluster_name)

highlight_cluster(cluster_1_residues, 'blue', 'cluster_1')
highlight_cluster(cluster_2_residues, 'red', 'cluster_2')

# Add centroids to visualization
if centroid_1:
    cmd.pseudoatom(object="centroid_1", pos=centroid_1, color="blue", label="Cluster 1 Centroid")
    cmd.show("spheres", "centroid_1")
if centroid_2:
    cmd.pseudoatom(object="centroid_2", pos=centroid_2, color="red", label="Cluster 2 Centroid")
    cmd.show("spheres", "centroid_2")

# Zoom into the structure
cmd.zoom()

# Export image with transparent background
output_path = "/path/to/output/visualize_clusters.png"  # Replace with your desired output path
cmd.ray()  # Render the image
cmd.png(output_path, dpi=300, ray=1, quiet=1, transparent=True)

# Print confirmation
print(f"Image saved to {output_path}")
