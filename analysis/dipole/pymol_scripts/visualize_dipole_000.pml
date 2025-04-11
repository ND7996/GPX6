# PyMOL script for visualizing dipole in replica 000
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.78232, 31.175106, 49.68994], name=R0, color=blue
pseudoatom reactant, pos=[38.78481, 31.459213, 49.446434], name=R1, color=blue
pseudoatom reactant, pos=[37.79894, 30.249958, 48.954224], name=R2, color=blue
pseudoatom reactant, pos=[38.083855, 30.75732, 49.81625], name=R3, color=blue
pseudoatom reactant, pos=[37.954235, 31.002577, 49.2827], name=R4, color=blue
pseudoatom reactant, pos=[38.71813, 30.997564, 49.401894], name=R5, color=blue
pseudoatom reactant, pos=[37.603775, 31.019289, 47.892082], name=R6, color=blue
pseudoatom reactant, pos=[37.673676, 30.77483, 48.04654], name=R7, color=blue
pseudoatom reactant, pos=[38.4756, 30.352613, 50.08604], name=R8, color=blue
pseudoatom reactant, pos=[38.677876, 31.054153, 49.792118], name=R9, color=blue
pseudoatom reactant, pos=[38.261753, 30.76364, 49.91998], name=R10, color=blue
pseudoatom reactant, pos=[38.18845, 30.667368, 49.63562], name=R11, color=blue
pseudoatom reactant, pos=[38.265316, 30.30331, 49.947655], name=R12, color=blue
pseudoatom reactant, pos=[38.32737, 30.63191, 48.80061], name=R13, color=blue
pseudoatom reactant, pos=[38.661697, 30.479229, 49.61397], name=R14, color=blue
pseudoatom reactant, pos=[38.634365, 30.751133, 49.150066], name=R15, color=blue
pseudoatom reactant, pos=[38.955814, 30.67723, 49.604145], name=R16, color=blue
pseudoatom reactant, pos=[38.181923, 30.525095, 49.211834], name=R17, color=blue
pseudoatom reactant, pos=[38.198383, 31.19769, 49.143787], name=R18, color=blue
pseudoatom reactant, pos=[38.637833, 31.337574, 49.314686], name=R19, color=blue
pseudoatom reactant, pos=[37.87968, 30.337992, 49.429142], name=R20, color=blue
pseudoatom reactant, pos=[39.33326, 31.6471, 49.0915], name=R21, color=blue
pseudoatom reactant, pos=[38.507694, 30.977373, 49.617935], name=R22, color=blue
pseudoatom reactant, pos=[37.957447, 31.19189, 49.571896], name=R23, color=blue
pseudoatom reactant, pos=[37.976933, 30.900331, 49.826443], name=R24, color=blue
pseudoatom reactant, pos=[38.200706, 31.061583, 49.79951], name=R25, color=blue
pseudoatom reactant, pos=[38.119934, 30.801632, 49.863117], name=R26, color=blue
pseudoatom reactant, pos=[38.42302, 31.109352, 49.70342], name=R27, color=blue
pseudoatom reactant, pos=[38.401985, 30.780527, 49.74099], name=R28, color=blue
pseudoatom reactant, pos=[38.74483, 31.476896, 49.54057], name=R29, color=blue
pseudoatom reactant, pos=[37.52998, 31.381178, 49.60243], name=R30, color=blue
pseudoatom reactant, pos=[38.120766, 31.02139, 49.846184], name=R31, color=blue
pseudoatom reactant, pos=[37.881657, 31.022377, 49.501053], name=R32, color=blue
pseudoatom reactant, pos=[38.231773, 30.42694, 50.077156], name=R33, color=blue
pseudoatom reactant, pos=[37.946987, 30.618425, 49.520588], name=R34, color=blue
pseudoatom reactant, pos=[37.759583, 30.915474, 49.14055], name=R35, color=blue
pseudoatom reactant, pos=[38.268475, 30.99651, 49.304176], name=R36, color=blue
pseudoatom reactant, pos=[38.719383, 30.794079, 49.753277], name=R37, color=blue
pseudoatom reactant, pos=[38.490242, 30.58299, 49.310642], name=R38, color=blue
pseudoatom reactant, pos=[38.288273, 30.634953, 49.061653], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[39.215252, 30.071198, 50.789642], name=P0, color=red
pseudoatom product, pos=[39.37917, 30.076767, 51.622417], name=P1, color=red
pseudoatom product, pos=[39.117374, 29.66883, 50.882576], name=P2, color=red
pseudoatom product, pos=[39.080025, 29.901352, 50.600967], name=P3, color=red
pseudoatom product, pos=[39.010826, 29.820255, 51.05012], name=P4, color=red
pseudoatom product, pos=[38.73045, 29.59415, 50.971783], name=P5, color=red
pseudoatom product, pos=[39.148888, 29.925545, 50.545967], name=P6, color=red
pseudoatom product, pos=[39.143173, 29.862835, 51.23013], name=P7, color=red
pseudoatom product, pos=[39.047253, 30.045675, 50.94167], name=P8, color=red
pseudoatom product, pos=[38.91447, 29.834867, 51.300327], name=P9, color=red
pseudoatom product, pos=[39.01992, 29.579908, 50.692883], name=P10, color=red
pseudoatom product, pos=[39.017616, 29.74107, 51.128017], name=P11, color=red
pseudoatom product, pos=[38.67815, 29.228012, 51.380226], name=P12, color=red
pseudoatom product, pos=[39.03225, 29.729488, 51.610752], name=P13, color=red
pseudoatom product, pos=[39.120934, 29.992235, 51.50172], name=P14, color=red
pseudoatom product, pos=[39.45696, 29.800013, 51.3212], name=P15, color=red
pseudoatom product, pos=[39.072033, 29.775108, 51.06761], name=P16, color=red
pseudoatom product, pos=[39.372692, 30.122272, 50.965763], name=P17, color=red
pseudoatom product, pos=[39.232956, 29.748789, 51.58871], name=P18, color=red
pseudoatom product, pos=[39.046516, 29.666258, 51.452923], name=P19, color=red
pseudoatom product, pos=[38.732857, 29.830904, 51.39323], name=P20, color=red
pseudoatom product, pos=[38.998577, 29.735025, 51.832165], name=P21, color=red
pseudoatom product, pos=[39.201447, 29.516804, 51.553905], name=P22, color=red
pseudoatom product, pos=[38.489758, 28.912697, 51.538017], name=P23, color=red
pseudoatom product, pos=[39.105934, 29.362524, 51.774616], name=P24, color=red
pseudoatom product, pos=[39.17962, 29.75064, 51.240135], name=P25, color=red
pseudoatom product, pos=[38.93145, 30.214676, 51.443615], name=P26, color=red
pseudoatom product, pos=[38.984524, 29.8966, 51.46807], name=P27, color=red
pseudoatom product, pos=[39.019714, 29.934605, 52.072422], name=P28, color=red
pseudoatom product, pos=[38.874203, 30.048235, 51.596264], name=P29, color=red
pseudoatom product, pos=[39.69769, 29.555115, 51.448193], name=P30, color=red
pseudoatom product, pos=[39.522346, 30.100323, 51.309296], name=P31, color=red
pseudoatom product, pos=[39.066334, 29.133135, 51.220703], name=P32, color=red
pseudoatom product, pos=[38.93568, 29.725697, 51.334877], name=P33, color=red
pseudoatom product, pos=[39.017715, 29.055643, 50.992256], name=P34, color=red
pseudoatom product, pos=[39.273235, 29.480581, 51.67565], name=P35, color=red
pseudoatom product, pos=[38.941994, 29.747158, 51.372803], name=P36, color=red
pseudoatom product, pos=[39.524822, 29.637465, 51.540478], name=P37, color=red
pseudoatom product, pos=[39.182285, 29.102236, 51.22256], name=P38, color=red
pseudoatom product, pos=[39.20517, 28.94341, 51.45153], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.29121822500001, 30.87139485, 49.45132017500001], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[39.09305657500001, 29.696702499999994, 51.30315469999998], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[38.69213740000001, 30.284048674999998, 50.3772374375], label=Magnitude:2.33

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
cmd.wizard("message", "Dipole Magnitude: 2.3350\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
