# PyMOL script for visualizing dipole in replica 005
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.298218, 30.773178, 49.1234], name=R0, color=blue
pseudoatom reactant, pos=[37.825916, 30.70045, 49.32433], name=R1, color=blue
pseudoatom reactant, pos=[37.901943, 30.932592, 49.100735], name=R2, color=blue
pseudoatom reactant, pos=[37.850273, 31.021008, 49.088017], name=R3, color=blue
pseudoatom reactant, pos=[37.86727, 31.107283, 48.83209], name=R4, color=blue
pseudoatom reactant, pos=[37.988037, 31.351336, 49.40514], name=R5, color=blue
pseudoatom reactant, pos=[38.22545, 30.947071, 49.547115], name=R6, color=blue
pseudoatom reactant, pos=[37.995712, 30.80623, 49.590897], name=R7, color=blue
pseudoatom reactant, pos=[37.851727, 31.040808, 49.527966], name=R8, color=blue
pseudoatom reactant, pos=[38.399826, 31.10557, 49.486538], name=R9, color=blue
pseudoatom reactant, pos=[38.120274, 30.948942, 49.31127], name=R10, color=blue
pseudoatom reactant, pos=[38.289867, 30.930416, 49.447582], name=R11, color=blue
pseudoatom reactant, pos=[37.486183, 30.614185, 49.254353], name=R12, color=blue
pseudoatom reactant, pos=[38.11418, 31.087847, 49.285187], name=R13, color=blue
pseudoatom reactant, pos=[37.867897, 30.220095, 49.43576], name=R14, color=blue
pseudoatom reactant, pos=[37.71865, 30.660646, 49.20072], name=R15, color=blue
pseudoatom reactant, pos=[37.935654, 30.768677, 49.547565], name=R16, color=blue
pseudoatom reactant, pos=[38.31523, 30.7323, 49.358944], name=R17, color=blue
pseudoatom reactant, pos=[37.394836, 31.124659, 48.10555], name=R18, color=blue
pseudoatom reactant, pos=[38.215977, 30.986502, 48.835236], name=R19, color=blue
pseudoatom reactant, pos=[38.304775, 31.08134, 49.28345], name=R20, color=blue
pseudoatom reactant, pos=[37.901897, 30.791128, 49.33919], name=R21, color=blue
pseudoatom reactant, pos=[38.495815, 31.383205, 48.929234], name=R22, color=blue
pseudoatom reactant, pos=[38.56576, 30.8801, 49.057644], name=R23, color=blue
pseudoatom reactant, pos=[36.97109, 31.0142, 48.95528], name=R24, color=blue
pseudoatom reactant, pos=[37.636513, 30.675571, 49.30744], name=R25, color=blue
pseudoatom reactant, pos=[38.069233, 31.180218, 49.400715], name=R26, color=blue
pseudoatom reactant, pos=[37.465927, 30.799082, 49.43785], name=R27, color=blue
pseudoatom reactant, pos=[37.59907, 30.81842, 49.14873], name=R28, color=blue
pseudoatom reactant, pos=[37.37441, 30.910746, 49.272503], name=R29, color=blue
pseudoatom reactant, pos=[38.128105, 30.965063, 48.71263], name=R30, color=blue
pseudoatom reactant, pos=[37.663795, 30.95798, 49.048035], name=R31, color=blue
pseudoatom reactant, pos=[37.6191, 30.502146, 49.072662], name=R32, color=blue
pseudoatom reactant, pos=[37.923046, 30.375538, 49.49448], name=R33, color=blue
pseudoatom reactant, pos=[37.84759, 31.172012, 49.14049], name=R34, color=blue
pseudoatom reactant, pos=[38.10483, 30.888956, 49.496372], name=R35, color=blue
pseudoatom reactant, pos=[37.453796, 30.504482, 49.164726], name=R36, color=blue
pseudoatom reactant, pos=[37.886738, 31.092142, 48.62434], name=R37, color=blue
pseudoatom reactant, pos=[37.87527, 30.243824, 48.940132], name=R38, color=blue
pseudoatom reactant, pos=[37.91092, 30.863253, 49.408096], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[37.38163, 30.036842, 51.92012], name=P0, color=red
pseudoatom product, pos=[37.233883, 29.877031, 51.424118], name=P1, color=red
pseudoatom product, pos=[37.235218, 29.832668, 51.59191], name=P2, color=red
pseudoatom product, pos=[37.12944, 30.298157, 51.136467], name=P3, color=red
pseudoatom product, pos=[36.939503, 30.469576, 51.575203], name=P4, color=red
pseudoatom product, pos=[36.848286, 30.067171, 51.390816], name=P5, color=red
pseudoatom product, pos=[37.289593, 29.905828, 51.113518], name=P6, color=red
pseudoatom product, pos=[37.112446, 30.033167, 51.484703], name=P7, color=red
pseudoatom product, pos=[36.90292, 29.816357, 51.548], name=P8, color=red
pseudoatom product, pos=[37.299976, 29.781305, 51.406235], name=P9, color=red
pseudoatom product, pos=[37.02323, 29.89196, 51.756332], name=P10, color=red
pseudoatom product, pos=[36.705257, 29.97667, 51.5333], name=P11, color=red
pseudoatom product, pos=[37.230957, 29.670774, 51.424896], name=P12, color=red
pseudoatom product, pos=[37.465664, 30.018427, 51.427883], name=P13, color=red
pseudoatom product, pos=[37.56547, 30.130894, 51.698284], name=P14, color=red
pseudoatom product, pos=[37.673756, 30.012432, 51.708176], name=P15, color=red
pseudoatom product, pos=[38.71472, 29.95706, 51.700207], name=P16, color=red
pseudoatom product, pos=[37.532806, 30.715075, 51.00771], name=P17, color=red
pseudoatom product, pos=[37.690285, 30.096241, 51.480484], name=P18, color=red
pseudoatom product, pos=[36.453644, 30.200844, 51.1318], name=P19, color=red
pseudoatom product, pos=[36.933258, 29.577522, 51.47669], name=P20, color=red
pseudoatom product, pos=[36.625412, 29.859375, 51.49099], name=P21, color=red
pseudoatom product, pos=[36.81665, 29.829237, 51.79047], name=P22, color=red
pseudoatom product, pos=[37.250637, 29.743635, 51.414375], name=P23, color=red
pseudoatom product, pos=[37.546906, 29.955612, 51.50686], name=P24, color=red
pseudoatom product, pos=[37.284824, 29.64942, 51.156475], name=P25, color=red
pseudoatom product, pos=[37.17925, 29.658092, 51.445267], name=P26, color=red
pseudoatom product, pos=[36.873707, 30.122145, 51.10638], name=P27, color=red
pseudoatom product, pos=[36.97801, 29.774485, 51.189106], name=P28, color=red
pseudoatom product, pos=[37.247456, 30.0403, 51.49717], name=P29, color=red
pseudoatom product, pos=[37.43743, 30.260208, 51.543034], name=P30, color=red
pseudoatom product, pos=[37.555515, 30.266386, 51.391876], name=P31, color=red
pseudoatom product, pos=[37.39228, 29.433458, 51.521515], name=P32, color=red
pseudoatom product, pos=[37.233433, 30.23073, 51.171673], name=P33, color=red
pseudoatom product, pos=[36.833157, 29.79441, 51.60646], name=P34, color=red
pseudoatom product, pos=[36.515438, 29.36386, 51.41447], name=P35, color=red
pseudoatom product, pos=[36.52726, 29.873735, 51.35389], name=P36, color=red
pseudoatom product, pos=[36.759422, 29.91923, 51.072628], name=P37, color=red
pseudoatom product, pos=[37.36289, 30.127716, 51.356514], name=P38, color=red
pseudoatom product, pos=[37.540867, 30.055616, 51.374752], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[37.91151999999999, 30.87398002499999, 49.20105985], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[37.183062150000005, 29.958091274999997, 51.433518925], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[37.547291075, 30.416035649999994, 50.3172893875], label=Magnitude:2.52

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
cmd.wizard("message", "Dipole Magnitude: 2.5206\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
