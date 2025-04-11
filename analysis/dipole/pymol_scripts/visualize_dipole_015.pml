# PyMOL script for visualizing dipole in replica 015
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[39.365894, 30.565332, 49.429714], name=R0, color=blue
pseudoatom reactant, pos=[38.70072, 30.605688, 49.654625], name=R1, color=blue
pseudoatom reactant, pos=[38.76334, 31.03475, 49.126247], name=R2, color=blue
pseudoatom reactant, pos=[38.01323, 30.304613, 49.369316], name=R3, color=blue
pseudoatom reactant, pos=[38.293133, 30.323002, 49.30005], name=R4, color=blue
pseudoatom reactant, pos=[38.025074, 31.027445, 49.210094], name=R5, color=blue
pseudoatom reactant, pos=[38.572365, 30.972706, 49.48375], name=R6, color=blue
pseudoatom reactant, pos=[38.469418, 30.964844, 49.450916], name=R7, color=blue
pseudoatom reactant, pos=[38.27365, 30.55482, 49.343193], name=R8, color=blue
pseudoatom reactant, pos=[37.98721, 30.961851, 49.24136], name=R9, color=blue
pseudoatom reactant, pos=[38.280514, 31.23079, 49.520912], name=R10, color=blue
pseudoatom reactant, pos=[38.287, 31.054901, 48.926388], name=R11, color=blue
pseudoatom reactant, pos=[38.57603, 31.008091, 49.131645], name=R12, color=blue
pseudoatom reactant, pos=[38.01737, 30.850441, 49.277702], name=R13, color=blue
pseudoatom reactant, pos=[38.06811, 31.153082, 49.577744], name=R14, color=blue
pseudoatom reactant, pos=[38.194527, 30.758303, 49.30055], name=R15, color=blue
pseudoatom reactant, pos=[37.90435, 31.045322, 49.434433], name=R16, color=blue
pseudoatom reactant, pos=[38.017536, 30.752695, 49.14783], name=R17, color=blue
pseudoatom reactant, pos=[37.622948, 31.332497, 48.94617], name=R18, color=blue
pseudoatom reactant, pos=[37.668995, 30.908998, 49.30324], name=R19, color=blue
pseudoatom reactant, pos=[38.250885, 31.109402, 48.992744], name=R20, color=blue
pseudoatom reactant, pos=[38.349113, 31.557112, 49.144173], name=R21, color=blue
pseudoatom reactant, pos=[38.46308, 31.036617, 49.292007], name=R22, color=blue
pseudoatom reactant, pos=[38.236816, 30.966293, 49.30322], name=R23, color=blue
pseudoatom reactant, pos=[38.214493, 30.802885, 49.033424], name=R24, color=blue
pseudoatom reactant, pos=[37.38492, 31.291933, 49.33928], name=R25, color=blue
pseudoatom reactant, pos=[38.390686, 30.790024, 49.211708], name=R26, color=blue
pseudoatom reactant, pos=[38.04958, 30.922659, 49.570312], name=R27, color=blue
pseudoatom reactant, pos=[38.02183, 30.742474, 49.75643], name=R28, color=blue
pseudoatom reactant, pos=[38.208553, 30.74818, 49.221672], name=R29, color=blue
pseudoatom reactant, pos=[38.327515, 31.023718, 49.601665], name=R30, color=blue
pseudoatom reactant, pos=[38.0431, 31.010803, 49.32441], name=R31, color=blue
pseudoatom reactant, pos=[38.381542, 30.803976, 49.800346], name=R32, color=blue
pseudoatom reactant, pos=[38.606018, 31.231497, 49.48496], name=R33, color=blue
pseudoatom reactant, pos=[38.184303, 31.176428, 49.4198], name=R34, color=blue
pseudoatom reactant, pos=[37.78275, 30.667774, 49.16573], name=R35, color=blue
pseudoatom reactant, pos=[37.863255, 30.861162, 49.706894], name=R36, color=blue
pseudoatom reactant, pos=[38.23011, 31.270927, 49.406414], name=R37, color=blue
pseudoatom reactant, pos=[38.01535, 30.822077, 49.260494], name=R38, color=blue
pseudoatom reactant, pos=[38.335087, 31.194012, 49.68625], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[37.82891, 29.906605, 51.30451], name=P0, color=red
pseudoatom product, pos=[38.367992, 29.717386, 51.423485], name=P1, color=red
pseudoatom product, pos=[38.0306, 29.644138, 51.67339], name=P2, color=red
pseudoatom product, pos=[38.10336, 30.831467, 50.582294], name=P3, color=red
pseudoatom product, pos=[38.22116, 30.316828, 51.23703], name=P4, color=red
pseudoatom product, pos=[38.388466, 30.258917, 51.383736], name=P5, color=red
pseudoatom product, pos=[37.609985, 30.042227, 51.16877], name=P6, color=red
pseudoatom product, pos=[37.804695, 30.251066, 51.64065], name=P7, color=red
pseudoatom product, pos=[38.193886, 30.39168, 51.07368], name=P8, color=red
pseudoatom product, pos=[38.543457, 30.113691, 51.231476], name=P9, color=red
pseudoatom product, pos=[37.807434, 30.289621, 50.834927], name=P10, color=red
pseudoatom product, pos=[37.898647, 30.123802, 51.163208], name=P11, color=red
pseudoatom product, pos=[37.952374, 30.307037, 51.164352], name=P12, color=red
pseudoatom product, pos=[38.139343, 30.102165, 51.37348], name=P13, color=red
pseudoatom product, pos=[38.100822, 30.109024, 51.339935], name=P14, color=red
pseudoatom product, pos=[38.065704, 30.035112, 51.29191], name=P15, color=red
pseudoatom product, pos=[37.482292, 29.996664, 51.310856], name=P16, color=red
pseudoatom product, pos=[37.314438, 30.3848, 51.24618], name=P17, color=red
pseudoatom product, pos=[37.94418, 30.118006, 51.53105], name=P18, color=red
pseudoatom product, pos=[37.64958, 29.96431, 51.370876], name=P19, color=red
pseudoatom product, pos=[37.919823, 30.317938, 51.210964], name=P20, color=red
pseudoatom product, pos=[37.50181, 30.18197, 51.47823], name=P21, color=red
pseudoatom product, pos=[37.489876, 30.142357, 51.260967], name=P22, color=red
pseudoatom product, pos=[37.622475, 29.97836, 51.566223], name=P23, color=red
pseudoatom product, pos=[37.683025, 29.735239, 51.721428], name=P24, color=red
pseudoatom product, pos=[37.0292, 30.120897, 51.162422], name=P25, color=red
pseudoatom product, pos=[37.644894, 30.28207, 51.6825], name=P26, color=red
pseudoatom product, pos=[37.69724, 30.01967, 51.352474], name=P27, color=red
pseudoatom product, pos=[38.379494, 29.708164, 51.281467], name=P28, color=red
pseudoatom product, pos=[38.07802, 29.34271, 51.311844], name=P29, color=red
pseudoatom product, pos=[37.45775, 30.158966, 51.264004], name=P30, color=red
pseudoatom product, pos=[38.26935, 29.863379, 51.00606], name=P31, color=red
pseudoatom product, pos=[38.104465, 29.490076, 51.33658], name=P32, color=red
pseudoatom product, pos=[38.418915, 29.97269, 51.49203], name=P33, color=red
pseudoatom product, pos=[38.03479, 29.814852, 51.29938], name=P34, color=red
pseudoatom product, pos=[37.77984, 30.496113, 51.13312], name=P35, color=red
pseudoatom product, pos=[37.687206, 30.348146, 51.13992], name=P36, color=red
pseudoatom product, pos=[38.31729, 30.405378, 50.96539], name=P37, color=red
pseudoatom product, pos=[38.336685, 30.34055, 51.41122], name=P38, color=red
pseudoatom product, pos=[38.20784, 30.151154, 51.565166], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.211009999999995, 30.936003099999994, 49.34744529999999], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[37.92768282499999, 30.094380624999992, 51.299679600000005], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[38.069346412499996, 30.515191862499993, 50.32356245], label=Magnitude:2.14

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
cmd.wizard("message", "Dipole Magnitude: 2.1447\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
