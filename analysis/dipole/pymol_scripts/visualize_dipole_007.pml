# PyMOL script for visualizing dipole in replica 007
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[37.65592, 30.438929, 48.89918], name=R0, color=blue
pseudoatom reactant, pos=[37.691124, 30.398695, 48.814846], name=R1, color=blue
pseudoatom reactant, pos=[37.299477, 31.097425, 48.207783], name=R2, color=blue
pseudoatom reactant, pos=[38.18175, 30.398893, 49.46429], name=R3, color=blue
pseudoatom reactant, pos=[38.18149, 30.703558, 49.102818], name=R4, color=blue
pseudoatom reactant, pos=[38.020046, 30.446459, 49.341152], name=R5, color=blue
pseudoatom reactant, pos=[38.1949, 30.447695, 49.237392], name=R6, color=blue
pseudoatom reactant, pos=[37.156162, 30.746124, 48.565517], name=R7, color=blue
pseudoatom reactant, pos=[37.43785, 30.646025, 48.699467], name=R8, color=blue
pseudoatom reactant, pos=[38.427288, 30.699585, 49.12562], name=R9, color=blue
pseudoatom reactant, pos=[38.29286, 30.76408, 49.679092], name=R10, color=blue
pseudoatom reactant, pos=[37.64492, 30.818428, 48.78962], name=R11, color=blue
pseudoatom reactant, pos=[38.294113, 30.07999, 49.195885], name=R12, color=blue
pseudoatom reactant, pos=[38.121944, 30.3505, 49.39873], name=R13, color=blue
pseudoatom reactant, pos=[38.103943, 30.347553, 49.232426], name=R14, color=blue
pseudoatom reactant, pos=[38.090004, 30.473003, 49.340508], name=R15, color=blue
pseudoatom reactant, pos=[37.622086, 30.386456, 49.16475], name=R16, color=blue
pseudoatom reactant, pos=[38.473347, 30.525652, 48.7401], name=R17, color=blue
pseudoatom reactant, pos=[38.54416, 30.25829, 49.378735], name=R18, color=blue
pseudoatom reactant, pos=[38.054882, 30.262663, 49.598194], name=R19, color=blue
pseudoatom reactant, pos=[37.303703, 30.998478, 48.94133], name=R20, color=blue
pseudoatom reactant, pos=[38.395676, 30.563738, 49.54009], name=R21, color=blue
pseudoatom reactant, pos=[38.62611, 30.2072, 49.378048], name=R22, color=blue
pseudoatom reactant, pos=[38.425884, 30.681011, 49.26977], name=R23, color=blue
pseudoatom reactant, pos=[38.305347, 30.701567, 49.495415], name=R24, color=blue
pseudoatom reactant, pos=[38.331192, 30.7837, 49.123665], name=R25, color=blue
pseudoatom reactant, pos=[38.272984, 30.335974, 49.045937], name=R26, color=blue
pseudoatom reactant, pos=[38.364674, 30.660141, 49.279583], name=R27, color=blue
pseudoatom reactant, pos=[38.099487, 30.45179, 49.22938], name=R28, color=blue
pseudoatom reactant, pos=[38.110497, 30.535126, 49.37653], name=R29, color=blue
pseudoatom reactant, pos=[37.89543, 30.563015, 49.142727], name=R30, color=blue
pseudoatom reactant, pos=[37.862133, 30.55597, 49.013462], name=R31, color=blue
pseudoatom reactant, pos=[38.732723, 30.65799, 49.322807], name=R32, color=blue
pseudoatom reactant, pos=[38.204315, 30.224487, 49.252], name=R33, color=blue
pseudoatom reactant, pos=[37.896244, 30.183527, 48.92082], name=R34, color=blue
pseudoatom reactant, pos=[38.2203, 30.778385, 49.662315], name=R35, color=blue
pseudoatom reactant, pos=[37.865215, 30.689867, 49.344204], name=R36, color=blue
pseudoatom reactant, pos=[38.561707, 30.457348, 49.388412], name=R37, color=blue
pseudoatom reactant, pos=[38.00765, 29.901648, 48.937286], name=R38, color=blue
pseudoatom reactant, pos=[36.85421, 30.618683, 48.146107], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[39.772003, 29.715876, 50.827854], name=P0, color=red
pseudoatom product, pos=[39.348576, 29.92045, 51.076824], name=P1, color=red
pseudoatom product, pos=[39.459972, 29.43592, 50.873573], name=P2, color=red
pseudoatom product, pos=[39.372787, 30.415123, 50.93203], name=P3, color=red
pseudoatom product, pos=[39.095707, 30.479933, 50.072094], name=P4, color=red
pseudoatom product, pos=[39.31698, 30.960655, 50.71205], name=P5, color=red
pseudoatom product, pos=[39.38313, 30.502892, 50.13287], name=P6, color=red
pseudoatom product, pos=[39.40139, 30.389723, 50.631058], name=P7, color=red
pseudoatom product, pos=[39.580383, 29.638882, 50.849644], name=P8, color=red
pseudoatom product, pos=[39.283825, 30.32871, 50.272747], name=P9, color=red
pseudoatom product, pos=[39.50892, 29.738686, 50.699722], name=P10, color=red
pseudoatom product, pos=[39.089634, 29.634209, 50.99624], name=P11, color=red
pseudoatom product, pos=[39.38072, 30.033255, 50.884007], name=P12, color=red
pseudoatom product, pos=[39.05978, 30.214123, 50.540474], name=P13, color=red
pseudoatom product, pos=[39.778778, 30.059452, 50.88163], name=P14, color=red
pseudoatom product, pos=[39.657528, 30.016388, 50.889908], name=P15, color=red
pseudoatom product, pos=[39.80397, 29.765543, 50.018467], name=P16, color=red
pseudoatom product, pos=[39.09988, 29.494478, 50.544304], name=P17, color=red
pseudoatom product, pos=[39.135883, 29.758951, 51.001873], name=P18, color=red
pseudoatom product, pos=[38.987732, 29.971739, 50.430016], name=P19, color=red
pseudoatom product, pos=[39.51513, 30.340574, 50.96224], name=P20, color=red
pseudoatom product, pos=[39.298912, 30.542198, 50.314816], name=P21, color=red
pseudoatom product, pos=[39.584763, 29.972805, 51.056374], name=P22, color=red
pseudoatom product, pos=[39.532814, 29.742676, 50.314445], name=P23, color=red
pseudoatom product, pos=[39.23028, 30.228088, 50.397907], name=P24, color=red
pseudoatom product, pos=[39.250214, 30.238516, 50.529564], name=P25, color=red
pseudoatom product, pos=[39.23977, 30.024866, 50.134426], name=P26, color=red
pseudoatom product, pos=[39.359207, 29.92854, 50.848286], name=P27, color=red
pseudoatom product, pos=[39.782444, 30.322767, 50.47947], name=P28, color=red
pseudoatom product, pos=[39.330826, 29.045395, 49.880497], name=P29, color=red
pseudoatom product, pos=[39.207264, 29.601452, 50.474308], name=P30, color=red
pseudoatom product, pos=[40.026985, 29.776339, 50.4532], name=P31, color=red
pseudoatom product, pos=[39.479107, 30.368473, 50.591564], name=P32, color=red
pseudoatom product, pos=[39.455845, 30.180836, 50.7232], name=P33, color=red
pseudoatom product, pos=[39.18843, 30.615812, 50.352673], name=P34, color=red
pseudoatom product, pos=[39.53624, 30.04914, 51.444843], name=P35, color=red
pseudoatom product, pos=[39.499645, 30.045528, 51.355934], name=P36, color=red
pseudoatom product, pos=[39.069363, 30.051157, 50.42592], name=P37, color=red
pseudoatom product, pos=[39.478897, 29.902056, 50.376774], name=P38, color=red
pseudoatom product, pos=[39.226112, 30.00449, 50.336483], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.045593675, 30.520991199999997, 49.14464982500001], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[39.395245649999985, 30.036417399999998, 50.618007725000005], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[38.72041966249999, 30.278704299999998, 49.88132877500001], label=Magnitude:2.06

# Set display options
hide everything
show spheres, reactant
show spheres, product
set sphere_scale, 0.5
show spheres, reactant_com
show spheres, product_com
set sphere_scale, 1.5, reactant_com
set sphere_scale, 1.5, product_com
hide label, dipole_label
show label, reactant_com
show label, product_com

# Set initial view
orient
zoom dipole, 5

# Create session states for different visualization options
scene dipole_only, store
disable reactant
disable product
scene coms_only, store
enable reactant
enable product
scene all_atoms, store
scene all_atoms, recall

# Print help message
cmd.wizard("message", "Dipole Magnitude: 2.0560\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
