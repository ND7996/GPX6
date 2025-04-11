# PyMOL script for visualizing dipole in replica 004
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.000935, 31.071796, 49.482323], name=R0, color=blue
pseudoatom reactant, pos=[37.168262, 32.029934, 48.071518], name=R1, color=blue
pseudoatom reactant, pos=[37.73692, 31.151392, 49.447163], name=R2, color=blue
pseudoatom reactant, pos=[38.732746, 31.172424, 49.953636], name=R3, color=blue
pseudoatom reactant, pos=[38.533363, 30.682304, 49.341232], name=R4, color=blue
pseudoatom reactant, pos=[37.639477, 30.505205, 49.505676], name=R5, color=blue
pseudoatom reactant, pos=[37.919765, 30.360046, 49.91607], name=R6, color=blue
pseudoatom reactant, pos=[37.520123, 30.73802, 49.30033], name=R7, color=blue
pseudoatom reactant, pos=[36.7018, 31.241695, 48.503258], name=R8, color=blue
pseudoatom reactant, pos=[38.068764, 30.668325, 49.62619], name=R9, color=blue
pseudoatom reactant, pos=[38.07701, 31.375938, 49.54886], name=R10, color=blue
pseudoatom reactant, pos=[37.02999, 31.863657, 47.85163], name=R11, color=blue
pseudoatom reactant, pos=[38.254128, 31.927734, 47.421017], name=R12, color=blue
pseudoatom reactant, pos=[37.77199, 30.426449, 48.858612], name=R13, color=blue
pseudoatom reactant, pos=[38.416077, 31.063038, 49.409454], name=R14, color=blue
pseudoatom reactant, pos=[37.279068, 30.471024, 48.299686], name=R15, color=blue
pseudoatom reactant, pos=[38.13566, 31.15653, 49.63805], name=R16, color=blue
pseudoatom reactant, pos=[38.189, 30.924398, 49.806385], name=R17, color=blue
pseudoatom reactant, pos=[37.880337, 30.93031, 49.478764], name=R18, color=blue
pseudoatom reactant, pos=[38.158012, 31.20111, 49.545113], name=R19, color=blue
pseudoatom reactant, pos=[37.535286, 30.707493, 50.285], name=R20, color=blue
pseudoatom reactant, pos=[37.661854, 30.79771, 49.849697], name=R21, color=blue
pseudoatom reactant, pos=[38.25792, 30.77627, 49.831505], name=R22, color=blue
pseudoatom reactant, pos=[37.79531, 30.541632, 49.994457], name=R23, color=blue
pseudoatom reactant, pos=[37.970554, 30.689539, 49.8581], name=R24, color=blue
pseudoatom reactant, pos=[37.681847, 31.115505, 49.08028], name=R25, color=blue
pseudoatom reactant, pos=[37.732418, 30.65882, 49.472588], name=R26, color=blue
pseudoatom reactant, pos=[38.2335, 30.850603, 49.317024], name=R27, color=blue
pseudoatom reactant, pos=[38.25718, 31.45758, 50.16092], name=R28, color=blue
pseudoatom reactant, pos=[38.32304, 31.12996, 49.8006], name=R29, color=blue
pseudoatom reactant, pos=[37.88379, 30.9989, 49.670795], name=R30, color=blue
pseudoatom reactant, pos=[37.70143, 31.044846, 49.081818], name=R31, color=blue
pseudoatom reactant, pos=[37.636055, 31.230396, 49.5355], name=R32, color=blue
pseudoatom reactant, pos=[37.846222, 30.615583, 49.127777], name=R33, color=blue
pseudoatom reactant, pos=[37.69662, 31.27712, 50.246746], name=R34, color=blue
pseudoatom reactant, pos=[37.850544, 31.393843, 49.703197], name=R35, color=blue
pseudoatom reactant, pos=[37.940853, 31.641409, 49.97079], name=R36, color=blue
pseudoatom reactant, pos=[38.044506, 31.437843, 49.766937], name=R37, color=blue
pseudoatom reactant, pos=[37.99942, 31.100027, 50.72903], name=R38, color=blue
pseudoatom reactant, pos=[37.478546, 31.626963, 49.557568], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[36.511932, 31.047117, 49.671936], name=P0, color=red
pseudoatom product, pos=[36.124268, 30.160414, 50.164917], name=P1, color=red
pseudoatom product, pos=[36.563545, 30.428148, 49.892445], name=P2, color=red
pseudoatom product, pos=[36.44515, 30.455704, 50.200314], name=P3, color=red
pseudoatom product, pos=[36.072376, 29.925835, 50.41812], name=P4, color=red
pseudoatom product, pos=[36.192295, 30.575369, 50.056232], name=P5, color=red
pseudoatom product, pos=[35.592995, 30.451752, 49.988533], name=P6, color=red
pseudoatom product, pos=[35.764263, 30.795094, 49.507423], name=P7, color=red
pseudoatom product, pos=[36.091606, 31.145348, 49.673344], name=P8, color=red
pseudoatom product, pos=[35.68041, 30.861801, 50.052628], name=P9, color=red
pseudoatom product, pos=[36.449974, 30.598354, 49.384964], name=P10, color=red
pseudoatom product, pos=[36.23249, 30.212624, 49.39221], name=P11, color=red
pseudoatom product, pos=[35.908936, 30.613317, 49.930115], name=P12, color=red
pseudoatom product, pos=[35.715168, 30.808681, 49.770805], name=P13, color=red
pseudoatom product, pos=[36.06833, 31.147427, 49.52685], name=P14, color=red
pseudoatom product, pos=[36.18956, 30.523848, 49.498837], name=P15, color=red
pseudoatom product, pos=[36.54102, 30.204706, 49.623547], name=P16, color=red
pseudoatom product, pos=[36.678165, 30.554222, 49.23821], name=P17, color=red
pseudoatom product, pos=[36.389896, 30.734842, 49.420647], name=P18, color=red
pseudoatom product, pos=[35.97323, 30.633604, 49.342644], name=P19, color=red
pseudoatom product, pos=[36.02606, 31.085354, 50.137306], name=P20, color=red
pseudoatom product, pos=[36.095215, 31.234287, 49.96036], name=P21, color=red
pseudoatom product, pos=[36.48199, 30.955206, 50.23444], name=P22, color=red
pseudoatom product, pos=[36.486217, 30.728807, 49.309574], name=P23, color=red
pseudoatom product, pos=[36.66314, 30.781797, 49.784725], name=P24, color=red
pseudoatom product, pos=[36.07214, 30.776846, 49.57438], name=P25, color=red
pseudoatom product, pos=[36.244083, 31.151937, 49.824284], name=P26, color=red
pseudoatom product, pos=[36.411327, 31.104061, 49.502003], name=P27, color=red
pseudoatom product, pos=[36.517643, 30.95104, 49.509663], name=P28, color=red
pseudoatom product, pos=[36.257767, 30.524216, 49.85552], name=P29, color=red
pseudoatom product, pos=[36.34552, 31.16443, 50.1895], name=P30, color=red
pseudoatom product, pos=[36.483036, 31.04034, 49.460476], name=P31, color=red
pseudoatom product, pos=[36.234413, 30.69867, 49.59636], name=P32, color=red
pseudoatom product, pos=[35.63912, 30.57516, 49.584953], name=P33, color=red
pseudoatom product, pos=[36.343876, 30.673794, 49.49919], name=P34, color=red
pseudoatom product, pos=[35.756077, 30.996044, 50.09876], name=P35, color=red
pseudoatom product, pos=[36.215878, 31.082531, 49.80598], name=P36, color=red
pseudoatom product, pos=[35.880527, 30.776329, 49.546253], name=P37, color=red
pseudoatom product, pos=[36.088615, 30.621962, 49.64331], name=P38, color=red
pseudoatom product, pos=[36.55181, 31.070581, 49.556053], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[37.86850805, 31.051334274999995, 49.451132400000006], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[36.199501575000006, 30.746789974999995, 49.735695275000005], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[37.0340048125, 30.899062124999993, 49.59341383750001], label=Magnitude:1.72

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
cmd.wizard("message", "Dipole Magnitude: 1.7203\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
