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

# Cluster 1: Residue 3-K
cmd.select('cluster_1', 'resid 3')
cmd.color('red', 'cluster_1')

# Cluster 2: Residues 4-S, 177-H, 178-T, 182-I
cmd.select('cluster_2', 'resid 4 or resid 177 or resid 178 or resid 182')
cmd.color('blue', 'cluster_2')

# Cluster 3: Residue 16-V
cmd.select('cluster_3', 'resid 16')
cmd.color('green', 'cluster_3')

# Cluster 4: Residue 22-N
cmd.select('cluster_4', 'resid 22')
cmd.color('yellow', 'cluster_4')

# Cluster 5: Residue 24-I
cmd.select('cluster_5', 'resid 24')
cmd.color('cyan', 'cluster_5')

# Cluster 6: Residue 25-D
cmd.select('cluster_6', 'resid 25')
cmd.color('magenta', 'cluster_6')

# Cluster 7: Residue 27-G
cmd.select('cluster_7', 'resid 27')
cmd.color('orange', 'cluster_7')

# Cluster 8: Residue 29-F
cmd.select('cluster_8', 'resid 29')
cmd.color('purple', 'cluster_8')

# Cluster 9: Residues 30-V, 31-N
cmd.select('cluster_9', 'resid 30 or resid 31')
cmd.color('lime', 'cluster_9')

# Cluster 10: Residue 33-Q
cmd.select('cluster_10', 'resid 33')
cmd.color('brown', 'cluster_10')

# Cluster 11: Residue 35-Y
cmd.select('cluster_11', 'resid 35')
cmd.color('pink', 'cluster_11')

# Cluster 12: Residue 40-I
cmd.select('cluster_12', 'resid 40')
cmd.color('violet', 'cluster_12')

# Cluster 13: Residues 47-S, 99-R
cmd.select('cluster_13', 'resid 47 or resid 99')
cmd.color('teal', 'cluster_13')

# Cluster 14: Residues 48-F, 52-T
cmd.select('cluster_14', 'resid 48 or resid 52')
cmd.color('indigo', 'cluster_14')

# Cluster 15: Residue 54-T
cmd.select('cluster_15', 'resid 54')
cmd.color('olive', 'cluster_15')

# Cluster 16: Residues 60-T, 181-R
cmd.select('cluster_16', 'resid 60 or resid 181')
cmd.color('coral', 'cluster_16')

# Cluster 17: Residues 67-P, 69-N
cmd.select('cluster_17', 'resid 67 or resid 69')
cmd.color('gold', 'cluster_17')

# Cluster 18: Residue 71-T
cmd.select('cluster_18', 'resid 71')
cmd.color('aqua', 'cluster_18')

# Cluster 19: Residue 74-G
cmd.select('cluster_19', 'resid 74')
cmd.color('salmon', 'cluster_19')

# Cluster 20: Residue 87-K
cmd.select('cluster_20', 'resid 87')
cmd.color('lime', 'cluster_20')

# Cluster 21: Residue 102-G
cmd.select('cluster_21', 'resid 102')
cmd.color('forest', 'cluster_21')

# Cluster 22: Residue 104-Y
cmd.select('cluster_22', 'resid 104')
cmd.color('chartreuse', 'cluster_22')

# Cluster 23: Residue 107-N
cmd.select('cluster_23', 'resid 107')
cmd.color('peach', 'cluster_23')

# Cluster 24: Residues 119-D, 120-N, 126-S
cmd.select('cluster_24', 'resid 119 or resid 120 or resid 126')
cmd.color('periwinkle', 'cluster_24')

# Cluster 25: Residue 137-E
cmd.select('cluster_25', 'resid 137')
cmd.color('lavender', 'cluster_25')

# Cluster 26: Residue 139-F
cmd.select('cluster_26', 'resid 139')
cmd.color('plum', 'cluster_26')

# Cluster 27: Residue 142-P
cmd.select('cluster_27', 'resid 142')
cmd.color('tan', 'cluster_27')

# Cluster 28: Residue 143-E
cmd.select('cluster_28', 'resid 143')
cmd.color('lightblue', 'cluster_28')

# Cluster 29: Residue 144-H
cmd.select('cluster_29', 'resid 144')
cmd.color('peachpuff', 'cluster_29')

# Cluster 30: Residue 148-D
cmd.select('cluster_30', 'resid 148')
cmd.color('snow', 'cluster_30')

# Cluster 31: Residue 173-R
cmd.select('cluster_31', 'resid 173')
cmd.color('crimson', 'cluster_31')

# Cluster 32: Residue 184-Q
cmd.select('cluster_32', 'resid 184')
cmd.color('darkorange', 'cluster_32')

# Cluster 33: Residues 188-M, 192-N
cmd.select('cluster_33', 'resid 188 or resid 192')
cmd.color('lightgreen', 'cluster_33')

# Cluster 34: Residues 194-T, 195-S
cmd.select('cluster_34', 'resid 194 or resid 195')
cmd.color('plum', 'cluster_34')

# Optionally, zoom to the structure
cmd.zoom()
