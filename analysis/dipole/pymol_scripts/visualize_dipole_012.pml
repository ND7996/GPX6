# PyMOL script for visualizing dipole in replica 012
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[37.784695, 30.46321, 49.139843], name=R0, color=blue
pseudoatom reactant, pos=[39.390133, 31.286928, 49.5653], name=R1, color=blue
pseudoatom reactant, pos=[37.186306, 31.511019, 48.50825], name=R2, color=blue
pseudoatom reactant, pos=[38.35589, 30.906973, 49.418713], name=R3, color=blue
pseudoatom reactant, pos=[37.96305, 30.679842, 49.713516], name=R4, color=blue
pseudoatom reactant, pos=[38.488304, 30.800753, 49.578365], name=R5, color=blue
pseudoatom reactant, pos=[38.696438, 30.759161, 49.741783], name=R6, color=blue
pseudoatom reactant, pos=[38.751423, 30.697304, 49.79205], name=R7, color=blue
pseudoatom reactant, pos=[38.60231, 30.674335, 49.362003], name=R8, color=blue
pseudoatom reactant, pos=[38.08093, 30.731709, 48.361828], name=R9, color=blue
pseudoatom reactant, pos=[38.662064, 30.783072, 49.672585], name=R10, color=blue
pseudoatom reactant, pos=[38.113747, 30.66371, 49.595894], name=R11, color=blue
pseudoatom reactant, pos=[37.98565, 30.510994, 49.620914], name=R12, color=blue
pseudoatom reactant, pos=[38.109684, 30.910498, 49.59094], name=R13, color=blue
pseudoatom reactant, pos=[38.311947, 30.343521, 49.04988], name=R14, color=blue
pseudoatom reactant, pos=[39.05827, 30.573307, 49.228096], name=R15, color=blue
pseudoatom reactant, pos=[38.02476, 30.626596, 48.91999], name=R16, color=blue
pseudoatom reactant, pos=[38.557148, 30.866308, 49.285606], name=R17, color=blue
pseudoatom reactant, pos=[38.322105, 30.761318, 49.1317], name=R18, color=blue
pseudoatom reactant, pos=[38.21162, 30.694408, 49.194733], name=R19, color=blue
pseudoatom reactant, pos=[38.246544, 30.55546, 49.04586], name=R20, color=blue
pseudoatom reactant, pos=[38.287243, 30.938042, 48.864998], name=R21, color=blue
pseudoatom reactant, pos=[37.953022, 30.716688, 49.48452], name=R22, color=blue
pseudoatom reactant, pos=[38.563007, 30.803024, 48.905987], name=R23, color=blue
pseudoatom reactant, pos=[38.1547, 30.976988, 49.21318], name=R24, color=blue
pseudoatom reactant, pos=[38.227203, 30.922709, 49.44398], name=R25, color=blue
pseudoatom reactant, pos=[38.339233, 30.865698, 49.29678], name=R26, color=blue
pseudoatom reactant, pos=[38.420853, 30.55729, 49.579292], name=R27, color=blue
pseudoatom reactant, pos=[38.373226, 30.292633, 49.58946], name=R28, color=blue
pseudoatom reactant, pos=[38.837524, 30.95478, 49.29388], name=R29, color=blue
pseudoatom reactant, pos=[38.528145, 30.773901, 49.274765], name=R30, color=blue
pseudoatom reactant, pos=[37.8513, 30.564217, 49.894104], name=R31, color=blue
pseudoatom reactant, pos=[37.934715, 30.793076, 49.285904], name=R32, color=blue
pseudoatom reactant, pos=[37.529507, 30.520731, 49.707706], name=R33, color=blue
pseudoatom reactant, pos=[38.510803, 30.974827, 49.593563], name=R34, color=blue
pseudoatom reactant, pos=[38.48573, 31.054613, 49.193966], name=R35, color=blue
pseudoatom reactant, pos=[38.335735, 30.931961, 49.349583], name=R36, color=blue
pseudoatom reactant, pos=[37.961086, 30.394905, 49.063427], name=R37, color=blue
pseudoatom reactant, pos=[38.415066, 31.1728, 49.82978], name=R38, color=blue
pseudoatom reactant, pos=[38.170235, 31.377655, 49.935616], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[38.777813, 30.409624, 51.0101], name=P0, color=red
pseudoatom product, pos=[39.097576, 30.230247, 51.257008], name=P1, color=red
pseudoatom product, pos=[38.195652, 30.830942, 50.93366], name=P2, color=red
pseudoatom product, pos=[38.606495, 30.883783, 50.2191], name=P3, color=red
pseudoatom product, pos=[39.23719, 30.553064, 50.748478], name=P4, color=red
pseudoatom product, pos=[39.19081, 30.113253, 50.703552], name=P5, color=red
pseudoatom product, pos=[38.973816, 30.994667, 50.325092], name=P6, color=red
pseudoatom product, pos=[38.86651, 29.930904, 50.990578], name=P7, color=red
pseudoatom product, pos=[39.046818, 30.28061, 50.558483], name=P8, color=red
pseudoatom product, pos=[38.826622, 30.471392, 50.30193], name=P9, color=red
pseudoatom product, pos=[38.811924, 31.324856, 49.873966], name=P10, color=red
pseudoatom product, pos=[39.54231, 31.579185, 49.771324], name=P11, color=red
pseudoatom product, pos=[39.388268, 30.632603, 50.359295], name=P12, color=red
pseudoatom product, pos=[39.520397, 30.35618, 50.138054], name=P13, color=red
pseudoatom product, pos=[38.900166, 30.417562, 49.97565], name=P14, color=red
pseudoatom product, pos=[38.86085, 30.791664, 50.254745], name=P15, color=red
pseudoatom product, pos=[38.718147, 30.217869, 50.25319], name=P16, color=red
pseudoatom product, pos=[38.632977, 30.027668, 50.777496], name=P17, color=red
pseudoatom product, pos=[39.490963, 30.350798, 50.48986], name=P18, color=red
pseudoatom product, pos=[39.357018, 30.264753, 50.719498], name=P19, color=red
pseudoatom product, pos=[39.132015, 30.371998, 50.70418], name=P20, color=red
pseudoatom product, pos=[39.29964, 30.000122, 50.099895], name=P21, color=red
pseudoatom product, pos=[39.009068, 30.146221, 50.74108], name=P22, color=red
pseudoatom product, pos=[38.58264, 30.045431, 51.251183], name=P23, color=red
pseudoatom product, pos=[39.32046, 29.864723, 50.417496], name=P24, color=red
pseudoatom product, pos=[38.893215, 30.304482, 50.37974], name=P25, color=red
pseudoatom product, pos=[38.935253, 30.18921, 50.589737], name=P26, color=red
pseudoatom product, pos=[38.64556, 30.144037, 50.211987], name=P27, color=red
pseudoatom product, pos=[38.973957, 30.136583, 50.94919], name=P28, color=red
pseudoatom product, pos=[38.82332, 30.027966, 50.93931], name=P29, color=red
pseudoatom product, pos=[38.96972, 30.54033, 50.858223], name=P30, color=red
pseudoatom product, pos=[38.83809, 30.710657, 50.853256], name=P31, color=red
pseudoatom product, pos=[39.34739, 30.8306, 50.901913], name=P32, color=red
pseudoatom product, pos=[39.342773, 30.224424, 50.508877], name=P33, color=red
pseudoatom product, pos=[39.14658, 30.417421, 50.93658], name=P34, color=red
pseudoatom product, pos=[39.277985, 30.360374, 51.129772], name=P35, color=red
pseudoatom product, pos=[38.944286, 30.295969, 50.45303], name=P36, color=red
pseudoatom product, pos=[39.09268, 29.515282, 51.627384], name=P37, color=red
pseudoatom product, pos=[39.208805, 30.391935, 50.729263], name=P38, color=red
pseudoatom product, pos=[38.73551, 30.163282, 50.88341], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.294533775000005, 30.784674100000007, 49.3579585], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[39.014031724999995, 30.383566775000013, 50.620664125], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[38.65428275, 30.58412043750001, 49.9893113125], label=Magnitude:1.51

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
cmd.wizard("message", "Dipole Magnitude: 1.5076\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
