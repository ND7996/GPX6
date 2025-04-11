# PyMOL script for visualizing dipole in replica 014
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.33345, 30.765846, 49.10782], name=R0, color=blue
pseudoatom reactant, pos=[38.521637, 30.31032, 49.312004], name=R1, color=blue
pseudoatom reactant, pos=[38.00655, 30.566874, 49.645844], name=R2, color=blue
pseudoatom reactant, pos=[38.2808, 30.882772, 49.015274], name=R3, color=blue
pseudoatom reactant, pos=[37.72955, 30.639902, 49.447914], name=R4, color=blue
pseudoatom reactant, pos=[37.7044, 30.48652, 48.85786], name=R5, color=blue
pseudoatom reactant, pos=[37.687584, 30.35008, 49.553795], name=R6, color=blue
pseudoatom reactant, pos=[38.41013, 30.70482, 49.541264], name=R7, color=blue
pseudoatom reactant, pos=[37.675076, 30.362171, 49.221497], name=R8, color=blue
pseudoatom reactant, pos=[37.83414, 31.905458, 47.3216], name=R9, color=blue
pseudoatom reactant, pos=[37.670868, 30.383259, 49.013042], name=R10, color=blue
pseudoatom reactant, pos=[37.835365, 31.241385, 48.81586], name=R11, color=blue
pseudoatom reactant, pos=[38.05603, 30.521128, 48.19756], name=R12, color=blue
pseudoatom reactant, pos=[38.540627, 31.124418, 49.109108], name=R13, color=blue
pseudoatom reactant, pos=[37.735115, 31.001291, 48.894318], name=R14, color=blue
pseudoatom reactant, pos=[38.1881, 30.869883, 49.380676], name=R15, color=blue
pseudoatom reactant, pos=[38.51487, 30.772722, 49.41417], name=R16, color=blue
pseudoatom reactant, pos=[37.797623, 30.797298, 49.748165], name=R17, color=blue
pseudoatom reactant, pos=[37.446312, 30.865484, 48.53219], name=R18, color=blue
pseudoatom reactant, pos=[38.276463, 31.449757, 48.693546], name=R19, color=blue
pseudoatom reactant, pos=[38.03638, 30.758133, 48.934612], name=R20, color=blue
pseudoatom reactant, pos=[37.695744, 30.647446, 48.953587], name=R21, color=blue
pseudoatom reactant, pos=[38.02668, 30.800928, 49.261665], name=R22, color=blue
pseudoatom reactant, pos=[37.85678, 30.612114, 48.73121], name=R23, color=blue
pseudoatom reactant, pos=[38.42903, 30.890497, 49.40612], name=R24, color=blue
pseudoatom reactant, pos=[38.015705, 30.586573, 48.784374], name=R25, color=blue
pseudoatom reactant, pos=[38.300823, 30.739635, 48.84036], name=R26, color=blue
pseudoatom reactant, pos=[38.469822, 30.95976, 48.96094], name=R27, color=blue
pseudoatom reactant, pos=[37.798904, 30.469484, 49.425064], name=R28, color=blue
pseudoatom reactant, pos=[38.512005, 30.645947, 49.130875], name=R29, color=blue
pseudoatom reactant, pos=[38.06936, 30.40063, 49.246506], name=R30, color=blue
pseudoatom reactant, pos=[38.25557, 30.959492, 49.373486], name=R31, color=blue
pseudoatom reactant, pos=[37.63952, 30.695444, 48.85769], name=R32, color=blue
pseudoatom reactant, pos=[37.930992, 30.615099, 49.840733], name=R33, color=blue
pseudoatom reactant, pos=[38.15762, 30.828896, 49.007214], name=R34, color=blue
pseudoatom reactant, pos=[38.465355, 30.523912, 49.48161], name=R35, color=blue
pseudoatom reactant, pos=[37.63612, 30.33455, 48.57382], name=R36, color=blue
pseudoatom reactant, pos=[38.567432, 31.015663, 49.488693], name=R37, color=blue
pseudoatom reactant, pos=[38.104115, 30.966003, 48.92817], name=R38, color=blue
pseudoatom reactant, pos=[38.30638, 30.744333, 49.23293], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[37.894806, 30.016571, 51.362183], name=P0, color=red
pseudoatom product, pos=[37.528294, 29.83754, 51.379696], name=P1, color=red
pseudoatom product, pos=[38.247738, 29.431406, 51.41859], name=P2, color=red
pseudoatom product, pos=[37.543545, 29.959776, 51.54291], name=P3, color=red
pseudoatom product, pos=[37.701378, 29.338394, 51.036808], name=P4, color=red
pseudoatom product, pos=[38.386986, 29.392199, 51.35794], name=P5, color=red
pseudoatom product, pos=[37.86298, 29.709936, 51.493584], name=P6, color=red
pseudoatom product, pos=[37.86436, 30.03058, 51.834187], name=P7, color=red
pseudoatom product, pos=[38.012596, 29.82544, 51.439693], name=P8, color=red
pseudoatom product, pos=[38.051044, 29.926937, 51.376648], name=P9, color=red
pseudoatom product, pos=[38.089355, 30.121859, 51.69454], name=P10, color=red
pseudoatom product, pos=[38.064194, 29.614733, 51.34563], name=P11, color=red
pseudoatom product, pos=[37.763912, 29.820923, 51.250874], name=P12, color=red
pseudoatom product, pos=[37.319332, 30.18537, 51.328636], name=P13, color=red
pseudoatom product, pos=[37.68278, 29.749609, 51.468872], name=P14, color=red
pseudoatom product, pos=[37.932076, 30.047731, 51.764435], name=P15, color=red
pseudoatom product, pos=[37.80721, 30.040432, 51.344364], name=P16, color=red
pseudoatom product, pos=[37.417156, 30.10989, 51.56943], name=P17, color=red
pseudoatom product, pos=[37.45711, 30.061123, 51.52255], name=P18, color=red
pseudoatom product, pos=[38.047676, 30.079313, 51.48034], name=P19, color=red
pseudoatom product, pos=[37.84272, 29.844374, 50.9118], name=P20, color=red
pseudoatom product, pos=[37.829113, 29.653631, 51.31262], name=P21, color=red
pseudoatom product, pos=[37.89995, 29.68871, 51.538456], name=P22, color=red
pseudoatom product, pos=[37.35875, 29.607542, 51.741333], name=P23, color=red
pseudoatom product, pos=[37.498283, 29.812134, 51.509865], name=P24, color=red
pseudoatom product, pos=[37.742737, 29.591259, 51.237484], name=P25, color=red
pseudoatom product, pos=[37.681293, 29.507814, 51.092056], name=P26, color=red
pseudoatom product, pos=[37.606636, 29.566622, 51.530224], name=P27, color=red
pseudoatom product, pos=[38.297417, 29.689434, 51.149994], name=P28, color=red
pseudoatom product, pos=[38.04168, 30.067982, 51.2609], name=P29, color=red
pseudoatom product, pos=[37.41015, 29.818615, 51.08882], name=P30, color=red
pseudoatom product, pos=[37.8054, 29.040512, 51.229366], name=P31, color=red
pseudoatom product, pos=[38.11445, 29.861034, 50.844387], name=P32, color=red
pseudoatom product, pos=[38.46723, 29.620739, 51.11978], name=P33, color=red
pseudoatom product, pos=[38.639423, 29.879381, 51.5042], name=P34, color=red
pseudoatom product, pos=[38.717896, 29.879913, 51.481865], name=P35, color=red
pseudoatom product, pos=[38.45252, 29.716118, 51.426594], name=P36, color=red
pseudoatom product, pos=[37.518562, 29.690977, 50.84805], name=P37, color=red
pseudoatom product, pos=[38.09253, 29.785997, 51.694748], name=P38, color=red
pseudoatom product, pos=[37.879047, 29.56992, 51.69887], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.062975675000004, 30.754898175, 49.08207914999999], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[37.889257875000006, 29.77981175, 51.38083305], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[37.97611677500001, 30.2673549625, 50.231456099999996], label=Magnitude:2.50

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
cmd.wizard("message", "Dipole Magnitude: 2.5030\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
