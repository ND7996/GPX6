# PyMOL script for visualizing dipole in replica 011
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.12643, 30.781118, 49.38132], name=R0, color=blue
pseudoatom reactant, pos=[38.87394, 30.37039, 49.538876], name=R1, color=blue
pseudoatom reactant, pos=[38.19616, 30.414375, 49.519794], name=R2, color=blue
pseudoatom reactant, pos=[37.572964, 30.477009, 49.179783], name=R3, color=blue
pseudoatom reactant, pos=[37.767693, 30.029118, 49.303837], name=R4, color=blue
pseudoatom reactant, pos=[37.81207, 30.419636, 48.807545], name=R5, color=blue
pseudoatom reactant, pos=[38.303913, 30.32552, 49.092514], name=R6, color=blue
pseudoatom reactant, pos=[38.188034, 30.576458, 49.1711], name=R7, color=blue
pseudoatom reactant, pos=[38.76776, 30.801311, 49.184303], name=R8, color=blue
pseudoatom reactant, pos=[37.99143, 30.648129, 49.155704], name=R9, color=blue
pseudoatom reactant, pos=[38.925507, 30.76846, 49.136772], name=R10, color=blue
pseudoatom reactant, pos=[38.544548, 30.70274, 49.287098], name=R11, color=blue
pseudoatom reactant, pos=[38.35, 30.766531, 49.298306], name=R12, color=blue
pseudoatom reactant, pos=[38.52176, 30.43273, 49.033794], name=R13, color=blue
pseudoatom reactant, pos=[38.67948, 30.989573, 49.334232], name=R14, color=blue
pseudoatom reactant, pos=[38.45997, 30.414276, 49.24271], name=R15, color=blue
pseudoatom reactant, pos=[38.891422, 30.528732, 49.421814], name=R16, color=blue
pseudoatom reactant, pos=[38.199192, 30.838345, 48.859577], name=R17, color=blue
pseudoatom reactant, pos=[38.106182, 30.526133, 49.26248], name=R18, color=blue
pseudoatom reactant, pos=[38.831173, 30.833128, 48.81063], name=R19, color=blue
pseudoatom reactant, pos=[37.649124, 30.632116, 49.20602], name=R20, color=blue
pseudoatom reactant, pos=[37.963604, 30.607662, 49.09512], name=R21, color=blue
pseudoatom reactant, pos=[38.37905, 30.662695, 49.20963], name=R22, color=blue
pseudoatom reactant, pos=[38.268787, 30.498215, 49.668713], name=R23, color=blue
pseudoatom reactant, pos=[37.936783, 30.961203, 48.870987], name=R24, color=blue
pseudoatom reactant, pos=[38.455887, 30.40202, 49.533478], name=R25, color=blue
pseudoatom reactant, pos=[38.34503, 30.427273, 49.152477], name=R26, color=blue
pseudoatom reactant, pos=[37.7072, 30.458357, 49.350906], name=R27, color=blue
pseudoatom reactant, pos=[38.30124, 30.553822, 49.407394], name=R28, color=blue
pseudoatom reactant, pos=[38.493904, 30.472315, 49.331234], name=R29, color=blue
pseudoatom reactant, pos=[38.66031, 30.450613, 49.039425], name=R30, color=blue
pseudoatom reactant, pos=[38.242455, 30.509024, 49.320625], name=R31, color=blue
pseudoatom reactant, pos=[38.282364, 30.304462, 48.930584], name=R32, color=blue
pseudoatom reactant, pos=[38.14144, 30.481287, 49.024902], name=R33, color=blue
pseudoatom reactant, pos=[38.552246, 30.412817, 49.545933], name=R34, color=blue
pseudoatom reactant, pos=[37.910793, 30.541336, 48.951733], name=R35, color=blue
pseudoatom reactant, pos=[38.58239, 30.47302, 49.22509], name=R36, color=blue
pseudoatom reactant, pos=[38.427635, 30.680458, 49.15346], name=R37, color=blue
pseudoatom reactant, pos=[37.90656, 30.406013, 48.66226], name=R38, color=blue
pseudoatom reactant, pos=[38.14781, 30.396873, 49.26593], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[36.71872, 29.690063, 48.363712], name=P0, color=red
pseudoatom product, pos=[36.718357, 29.513304, 48.397205], name=P1, color=red
pseudoatom product, pos=[36.53168, 29.371843, 48.595905], name=P2, color=red
pseudoatom product, pos=[35.880806, 29.376635, 48.338852], name=P3, color=red
pseudoatom product, pos=[36.01556, 29.369318, 48.37847], name=P4, color=red
pseudoatom product, pos=[36.52124, 29.41151, 48.659904], name=P5, color=red
pseudoatom product, pos=[36.426785, 29.392183, 48.656445], name=P6, color=red
pseudoatom product, pos=[35.965, 29.397533, 48.544014], name=P7, color=red
pseudoatom product, pos=[36.360188, 29.754831, 48.669823], name=P8, color=red
pseudoatom product, pos=[35.97406, 29.012066, 48.675426], name=P9, color=red
pseudoatom product, pos=[36.017986, 29.334425, 48.764137], name=P10, color=red
pseudoatom product, pos=[35.647995, 29.315199, 48.537117], name=P11, color=red
pseudoatom product, pos=[36.066185, 28.910276, 48.578297], name=P12, color=red
pseudoatom product, pos=[36.431725, 29.851065, 48.701313], name=P13, color=red
pseudoatom product, pos=[35.904106, 29.107998, 48.300495], name=P14, color=red
pseudoatom product, pos=[36.550552, 29.254902, 48.708134], name=P15, color=red
pseudoatom product, pos=[36.07851, 29.448622, 48.534966], name=P16, color=red
pseudoatom product, pos=[36.116207, 29.480856, 48.651608], name=P17, color=red
pseudoatom product, pos=[35.951145, 29.267624, 48.8887], name=P18, color=red
pseudoatom product, pos=[35.935226, 29.461653, 48.743305], name=P19, color=red
pseudoatom product, pos=[36.72111, 28.748463, 48.96548], name=P20, color=red
pseudoatom product, pos=[36.098656, 28.979216, 48.906467], name=P21, color=red
pseudoatom product, pos=[36.45932, 29.351171, 48.606373], name=P22, color=red
pseudoatom product, pos=[35.504997, 29.284016, 48.700348], name=P23, color=red
pseudoatom product, pos=[36.574276, 29.121576, 48.78837], name=P24, color=red
pseudoatom product, pos=[36.115795, 29.443766, 48.58137], name=P25, color=red
pseudoatom product, pos=[35.856617, 29.424932, 48.367584], name=P26, color=red
pseudoatom product, pos=[36.148865, 29.405458, 48.835384], name=P27, color=red
pseudoatom product, pos=[36.55433, 29.804127, 48.650623], name=P28, color=red
pseudoatom product, pos=[36.64525, 29.445618, 48.49666], name=P29, color=red
pseudoatom product, pos=[36.40638, 29.567062, 48.44532], name=P30, color=red
pseudoatom product, pos=[35.89955, 29.182894, 48.766468], name=P31, color=red
pseudoatom product, pos=[35.974915, 29.487846, 48.415855], name=P32, color=red
pseudoatom product, pos=[36.214066, 29.059067, 48.930923], name=P33, color=red
pseudoatom product, pos=[36.365368, 29.171854, 48.28765], name=P34, color=red
pseudoatom product, pos=[36.13765, 28.939827, 48.430325], name=P35, color=red
pseudoatom product, pos=[36.66744, 29.548695, 48.738094], name=P36, color=red
pseudoatom product, pos=[36.4741, 29.31725, 48.58763], name=P37, color=red
pseudoatom product, pos=[37.121193, 29.589828, 48.964878], name=P38, color=red
pseudoatom product, pos=[36.75014, 29.30466, 48.56651], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.286606, 30.54938232499999, 49.19920225000002], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[36.26255127500001, 29.347480799999992, 48.61800350000001], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[37.2745786375, 29.94843156249999, 48.908602875000014], label=Magnitude:2.42

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
cmd.wizard("message", "Dipole Magnitude: 2.4247\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
