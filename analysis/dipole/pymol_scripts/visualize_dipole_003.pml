# PyMOL script for visualizing dipole in replica 003
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[37.37926, 31.203756, 49.79027], name=R0, color=blue
pseudoatom reactant, pos=[38.600063, 31.199045, 49.725655], name=R1, color=blue
pseudoatom reactant, pos=[37.866207, 30.861483, 49.548187], name=R2, color=blue
pseudoatom reactant, pos=[38.191113, 31.085321, 49.168087], name=R3, color=blue
pseudoatom reactant, pos=[37.669754, 30.957285, 49.16494], name=R4, color=blue
pseudoatom reactant, pos=[38.190247, 31.341534, 49.764523], name=R5, color=blue
pseudoatom reactant, pos=[38.526386, 30.65303, 49.64173], name=R6, color=blue
pseudoatom reactant, pos=[38.455795, 30.481783, 49.966713], name=R7, color=blue
pseudoatom reactant, pos=[38.060677, 30.506239, 49.42967], name=R8, color=blue
pseudoatom reactant, pos=[38.675385, 30.866503, 49.330746], name=R9, color=blue
pseudoatom reactant, pos=[38.219696, 30.907833, 49.104347], name=R10, color=blue
pseudoatom reactant, pos=[38.24946, 31.070847, 49.26503], name=R11, color=blue
pseudoatom reactant, pos=[38.17329, 31.31348, 48.960136], name=R12, color=blue
pseudoatom reactant, pos=[37.865273, 31.260584, 49.07603], name=R13, color=blue
pseudoatom reactant, pos=[37.671642, 31.111359, 48.861984], name=R14, color=blue
pseudoatom reactant, pos=[38.419754, 30.546843, 49.477444], name=R15, color=blue
pseudoatom reactant, pos=[38.031124, 31.138168, 49.233467], name=R16, color=blue
pseudoatom reactant, pos=[38.30537, 31.248913, 49.30122], name=R17, color=blue
pseudoatom reactant, pos=[38.25045, 30.54323, 49.38722], name=R18, color=blue
pseudoatom reactant, pos=[38.414192, 31.115028, 49.41511], name=R19, color=blue
pseudoatom reactant, pos=[37.69944, 30.981928, 48.726547], name=R20, color=blue
pseudoatom reactant, pos=[38.73405, 30.929518, 49.43008], name=R21, color=blue
pseudoatom reactant, pos=[38.08489, 30.827148, 48.96449], name=R22, color=blue
pseudoatom reactant, pos=[37.88043, 31.172863, 49.214962], name=R23, color=blue
pseudoatom reactant, pos=[38.246326, 30.98539, 49.148663], name=R24, color=blue
pseudoatom reactant, pos=[37.973957, 30.949667, 49.44274], name=R25, color=blue
pseudoatom reactant, pos=[37.89374, 31.157137, 49.466545], name=R26, color=blue
pseudoatom reactant, pos=[38.328228, 31.249168, 48.95719], name=R27, color=blue
pseudoatom reactant, pos=[39.124092, 31.17901, 49.156513], name=R28, color=blue
pseudoatom reactant, pos=[38.654716, 30.786919, 49.623135], name=R29, color=blue
pseudoatom reactant, pos=[38.395245, 31.270687, 49.337593], name=R30, color=blue
pseudoatom reactant, pos=[38.289253, 30.651165, 49.38915], name=R31, color=blue
pseudoatom reactant, pos=[38.064182, 30.640265, 49.341343], name=R32, color=blue
pseudoatom reactant, pos=[38.455746, 31.07, 49.725075], name=R33, color=blue
pseudoatom reactant, pos=[38.17784, 30.686747, 49.394627], name=R34, color=blue
pseudoatom reactant, pos=[38.947266, 30.866814, 50.10541], name=R35, color=blue
pseudoatom reactant, pos=[37.914352, 31.130487, 49.242455], name=R36, color=blue
pseudoatom reactant, pos=[38.59901, 30.952993, 49.234398], name=R37, color=blue
pseudoatom reactant, pos=[38.337578, 30.444578, 49.25104], name=R38, color=blue
pseudoatom reactant, pos=[38.166145, 30.795645, 49.16624], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[35.898468, 30.721643, 49.642387], name=P0, color=red
pseudoatom product, pos=[36.144264, 31.48323, 49.706623], name=P1, color=red
pseudoatom product, pos=[36.05621, 31.20231, 49.638912], name=P2, color=red
pseudoatom product, pos=[36.19962, 31.220865, 49.514595], name=P3, color=red
pseudoatom product, pos=[35.927666, 30.98769, 50.078583], name=P4, color=red
pseudoatom product, pos=[36.48415, 30.984613, 50.130157], name=P5, color=red
pseudoatom product, pos=[35.917835, 31.124138, 49.259857], name=P6, color=red
pseudoatom product, pos=[35.87269, 31.406895, 50.055454], name=P7, color=red
pseudoatom product, pos=[36.24759, 30.861166, 49.913322], name=P8, color=red
pseudoatom product, pos=[36.139046, 31.025108, 49.758156], name=P9, color=red
pseudoatom product, pos=[36.39294, 30.859325, 49.487835], name=P10, color=red
pseudoatom product, pos=[35.861523, 31.272104, 50.02645], name=P11, color=red
pseudoatom product, pos=[36.154316, 30.986116, 49.748695], name=P12, color=red
pseudoatom product, pos=[36.119114, 31.097464, 50.075104], name=P13, color=red
pseudoatom product, pos=[35.86307, 30.958633, 50.21426], name=P14, color=red
pseudoatom product, pos=[36.032658, 31.228561, 49.81311], name=P15, color=red
pseudoatom product, pos=[35.67242, 31.850954, 50.400986], name=P16, color=red
pseudoatom product, pos=[35.67955, 31.60719, 49.89012], name=P17, color=red
pseudoatom product, pos=[35.813526, 31.048777, 50.015812], name=P18, color=red
pseudoatom product, pos=[36.12247, 31.558487, 49.80683], name=P19, color=red
pseudoatom product, pos=[36.637615, 31.17108, 49.610786], name=P20, color=red
pseudoatom product, pos=[36.936104, 30.371653, 49.384506], name=P21, color=red
pseudoatom product, pos=[36.310143, 31.285517, 49.74152], name=P22, color=red
pseudoatom product, pos=[36.453403, 31.345446, 49.689907], name=P23, color=red
pseudoatom product, pos=[36.235672, 31.111265, 49.838345], name=P24, color=red
pseudoatom product, pos=[35.727375, 31.659754, 49.808018], name=P25, color=red
pseudoatom product, pos=[36.468704, 31.29173, 49.80609], name=P26, color=red
pseudoatom product, pos=[36.637257, 31.578295, 49.52232], name=P27, color=red
pseudoatom product, pos=[36.45105, 31.1863, 50.27248], name=P28, color=red
pseudoatom product, pos=[36.523125, 31.551308, 49.394382], name=P29, color=red
pseudoatom product, pos=[37.057735, 31.960217, 49.812386], name=P30, color=red
pseudoatom product, pos=[36.38001, 31.1586, 50.51324], name=P31, color=red
pseudoatom product, pos=[37.120018, 31.072306, 49.470497], name=P32, color=red
pseudoatom product, pos=[36.396713, 31.066748, 49.640205], name=P33, color=red
pseudoatom product, pos=[36.52685, 31.394352, 49.80795], name=P34, color=red
pseudoatom product, pos=[36.915054, 31.799337, 49.28323], name=P35, color=red
pseudoatom product, pos=[36.3488, 31.417437, 49.757183], name=P36, color=red
pseudoatom product, pos=[36.22548, 31.385092, 49.96879], name=P37, color=red
pseudoatom product, pos=[36.288696, 31.401558, 49.735546], name=P38, color=red
pseudoatom product, pos=[36.29093, 31.33312, 50.25857], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.22954059999999, 30.953509825000005, 49.348267625000005], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[36.263246500000015, 31.2506596, 49.81232997499999], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[37.24639355000001, 31.1020847125, 49.580298799999994], label=Magnitude:2.04

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
cmd.wizard("message", "Dipole Magnitude: 2.0420\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
