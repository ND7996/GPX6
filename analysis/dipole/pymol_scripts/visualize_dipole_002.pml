# PyMOL script for visualizing dipole in replica 002
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[37.913006, 30.664553, 49.138042], name=R0, color=blue
pseudoatom reactant, pos=[38.43093, 30.67015, 49.65391], name=R1, color=blue
pseudoatom reactant, pos=[38.348686, 30.438662, 49.318035], name=R2, color=blue
pseudoatom reactant, pos=[38.05116, 30.771368, 49.18137], name=R3, color=blue
pseudoatom reactant, pos=[37.738064, 30.559319, 48.898148], name=R4, color=blue
pseudoatom reactant, pos=[38.166008, 30.192247, 49.50433], name=R5, color=blue
pseudoatom reactant, pos=[37.884663, 30.640797, 49.457413], name=R6, color=blue
pseudoatom reactant, pos=[37.808315, 30.42348, 49.10775], name=R7, color=blue
pseudoatom reactant, pos=[37.899567, 30.434874, 49.49118], name=R8, color=blue
pseudoatom reactant, pos=[38.360416, 30.80329, 49.57016], name=R9, color=blue
pseudoatom reactant, pos=[38.394382, 30.814035, 48.96115], name=R10, color=blue
pseudoatom reactant, pos=[38.91242, 30.514023, 49.46881], name=R11, color=blue
pseudoatom reactant, pos=[38.510277, 30.80201, 49.400852], name=R12, color=blue
pseudoatom reactant, pos=[38.34893, 30.306711, 49.610947], name=R13, color=blue
pseudoatom reactant, pos=[38.217106, 30.656807, 49.22202], name=R14, color=blue
pseudoatom reactant, pos=[37.92665, 30.833134, 49.00449], name=R15, color=blue
pseudoatom reactant, pos=[38.937977, 30.986008, 49.218002], name=R16, color=blue
pseudoatom reactant, pos=[38.135895, 31.009586, 49.19914], name=R17, color=blue
pseudoatom reactant, pos=[37.851204, 30.81899, 48.979645], name=R18, color=blue
pseudoatom reactant, pos=[37.687008, 30.669641, 48.832344], name=R19, color=blue
pseudoatom reactant, pos=[38.283752, 30.80073, 49.24549], name=R20, color=blue
pseudoatom reactant, pos=[37.943577, 30.689829, 49.26317], name=R21, color=blue
pseudoatom reactant, pos=[37.870766, 30.621426, 48.964466], name=R22, color=blue
pseudoatom reactant, pos=[37.72693, 30.650665, 49.405613], name=R23, color=blue
pseudoatom reactant, pos=[38.468704, 30.795618, 49.641182], name=R24, color=blue
pseudoatom reactant, pos=[38.102154, 31.006361, 49.496136], name=R25, color=blue
pseudoatom reactant, pos=[38.525707, 30.71907, 49.35999], name=R26, color=blue
pseudoatom reactant, pos=[38.251232, 30.83414, 49.37472], name=R27, color=blue
pseudoatom reactant, pos=[38.1752, 30.304298, 49.54523], name=R28, color=blue
pseudoatom reactant, pos=[37.892994, 30.670345, 49.672462], name=R29, color=blue
pseudoatom reactant, pos=[38.00095, 30.656082, 49.43092], name=R30, color=blue
pseudoatom reactant, pos=[38.148148, 30.67135, 49.203793], name=R31, color=blue
pseudoatom reactant, pos=[38.27986, 30.480608, 49.476498], name=R32, color=blue
pseudoatom reactant, pos=[37.984272, 30.373592, 48.921703], name=R33, color=blue
pseudoatom reactant, pos=[38.50561, 30.912367, 49.63383], name=R34, color=blue
pseudoatom reactant, pos=[38.70152, 30.824781, 48.833508], name=R35, color=blue
pseudoatom reactant, pos=[37.806763, 30.606945, 48.88198], name=R36, color=blue
pseudoatom reactant, pos=[38.535538, 30.857384, 49.1958], name=R37, color=blue
pseudoatom reactant, pos=[37.918163, 30.521852, 49.156746], name=R38, color=blue
pseudoatom reactant, pos=[38.215508, 30.598183, 49.06374], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[38.42226, 30.942625, 50.839565], name=P0, color=red
pseudoatom product, pos=[37.29394, 29.82954, 51.49665], name=P1, color=red
pseudoatom product, pos=[37.982998, 29.670181, 51.67075], name=P2, color=red
pseudoatom product, pos=[37.740086, 30.309036, 51.33651], name=P3, color=red
pseudoatom product, pos=[37.663925, 30.186878, 51.265434], name=P4, color=red
pseudoatom product, pos=[38.03981, 29.974138, 51.16969], name=P5, color=red
pseudoatom product, pos=[37.80834, 29.99724, 51.31582], name=P6, color=red
pseudoatom product, pos=[37.756607, 30.377647, 51.11071], name=P7, color=red
pseudoatom product, pos=[37.880627, 30.43444, 51.490368], name=P8, color=red
pseudoatom product, pos=[37.295074, 30.357607, 51.64443], name=P9, color=red
pseudoatom product, pos=[37.92353, 30.165785, 51.09391], name=P10, color=red
pseudoatom product, pos=[37.082226, 29.870743, 51.032154], name=P11, color=red
pseudoatom product, pos=[37.213024, 30.104572, 51.209396], name=P12, color=red
pseudoatom product, pos=[36.948586, 30.589777, 50.86288], name=P13, color=red
pseudoatom product, pos=[37.729713, 30.644903, 50.661728], name=P14, color=red
pseudoatom product, pos=[37.626385, 30.386223, 50.385685], name=P15, color=red
pseudoatom product, pos=[38.305145, 29.892399, 51.212784], name=P16, color=red
pseudoatom product, pos=[37.833046, 30.081713, 51.317944], name=P17, color=red
pseudoatom product, pos=[37.989655, 30.04725, 51.536453], name=P18, color=red
pseudoatom product, pos=[38.119907, 30.097818, 51.396202], name=P19, color=red
pseudoatom product, pos=[37.888943, 29.728914, 51.420372], name=P20, color=red
pseudoatom product, pos=[38.170033, 29.703526, 51.24639], name=P21, color=red
pseudoatom product, pos=[38.20969, 29.98237, 51.532024], name=P22, color=red
pseudoatom product, pos=[38.511845, 30.28164, 50.871986], name=P23, color=red
pseudoatom product, pos=[38.418503, 29.880016, 51.754593], name=P24, color=red
pseudoatom product, pos=[38.753635, 30.04768, 51.679092], name=P25, color=red
pseudoatom product, pos=[38.59809, 30.188976, 51.400608], name=P26, color=red
pseudoatom product, pos=[38.828575, 30.05627, 50.74045], name=P27, color=red
pseudoatom product, pos=[37.996803, 30.026815, 51.572624], name=P28, color=red
pseudoatom product, pos=[38.656376, 29.856474, 51.158455], name=P29, color=red
pseudoatom product, pos=[38.215527, 30.192879, 51.256084], name=P30, color=red
pseudoatom product, pos=[38.083897, 30.306475, 51.045975], name=P31, color=red
pseudoatom product, pos=[38.154472, 30.571606, 50.876415], name=P32, color=red
pseudoatom product, pos=[38.36886, 30.400604, 51.18026], name=P33, color=red
pseudoatom product, pos=[38.261997, 30.360325, 51.472836], name=P34, color=red
pseudoatom product, pos=[38.104435, 30.013918, 51.89286], name=P35, color=red
pseudoatom product, pos=[38.520924, 30.164722, 51.101685], name=P36, color=red
pseudoatom product, pos=[37.945267, 30.625286, 50.60908], name=P37, color=red
pseudoatom product, pos=[38.22494, 30.236584, 51.215782], name=P38, color=red
pseudoatom product, pos=[38.34097, 30.010715, 51.64218], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.17150030000001, 30.665132775, 49.27461787499999], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[38.02271664999999, 30.164907749999998, 51.24297034999999], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[38.097108475, 30.4150202625, 50.25879411249999], label=Magnitude:2.04

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
cmd.wizard("message", "Dipole Magnitude: 2.0364\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
