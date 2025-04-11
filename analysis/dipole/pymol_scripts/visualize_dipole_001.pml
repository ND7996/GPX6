# PyMOL script for visualizing dipole in replica 001
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.16695, 30.908566, 49.12702], name=R0, color=blue
pseudoatom reactant, pos=[38.05617, 30.970783, 49.16213], name=R1, color=blue
pseudoatom reactant, pos=[38.132324, 30.93639, 49.407093], name=R2, color=blue
pseudoatom reactant, pos=[37.80614, 30.840202, 48.99602], name=R3, color=blue
pseudoatom reactant, pos=[37.416496, 30.802126, 48.138496], name=R4, color=blue
pseudoatom reactant, pos=[37.973362, 31.165873, 49.06769], name=R5, color=blue
pseudoatom reactant, pos=[36.95766, 31.429129, 48.69911], name=R6, color=blue
pseudoatom reactant, pos=[37.960323, 31.215584, 49.782627], name=R7, color=blue
pseudoatom reactant, pos=[37.662888, 31.034395, 49.12807], name=R8, color=blue
pseudoatom reactant, pos=[37.503536, 30.588947, 49.624466], name=R9, color=blue
pseudoatom reactant, pos=[37.4845, 30.888102, 49.377296], name=R10, color=blue
pseudoatom reactant, pos=[37.543804, 31.06289, 49.219402], name=R11, color=blue
pseudoatom reactant, pos=[36.770733, 31.675617, 48.302742], name=R12, color=blue
pseudoatom reactant, pos=[36.86116, 31.707432, 48.05078], name=R13, color=blue
pseudoatom reactant, pos=[36.836708, 31.892063, 47.84098], name=R14, color=blue
pseudoatom reactant, pos=[37.46015, 30.84465, 49.25534], name=R15, color=blue
pseudoatom reactant, pos=[38.082226, 30.70581, 49.353077], name=R16, color=blue
pseudoatom reactant, pos=[37.790092, 31.113907, 49.256603], name=R17, color=blue
pseudoatom reactant, pos=[38.293854, 31.177055, 49.30079], name=R18, color=blue
pseudoatom reactant, pos=[37.999584, 31.226868, 48.90853], name=R19, color=blue
pseudoatom reactant, pos=[37.686436, 31.141071, 49.12146], name=R20, color=blue
pseudoatom reactant, pos=[37.67387, 30.757801, 49.467693], name=R21, color=blue
pseudoatom reactant, pos=[37.785587, 30.320732, 49.22816], name=R22, color=blue
pseudoatom reactant, pos=[37.39771, 30.605268, 49.400234], name=R23, color=blue
pseudoatom reactant, pos=[37.072643, 30.841043, 48.478657], name=R24, color=blue
pseudoatom reactant, pos=[37.757225, 30.813507, 49.145424], name=R25, color=blue
pseudoatom reactant, pos=[38.018883, 31.070488, 49.888954], name=R26, color=blue
pseudoatom reactant, pos=[37.37218, 30.789204, 48.583706], name=R27, color=blue
pseudoatom reactant, pos=[38.27428, 31.007818, 48.929615], name=R28, color=blue
pseudoatom reactant, pos=[36.938168, 32.85447, 47.47676], name=R29, color=blue
pseudoatom reactant, pos=[35.97416, 32.75547, 48.21456], name=R30, color=blue
pseudoatom reactant, pos=[37.654705, 31.238365, 48.963566], name=R31, color=blue
pseudoatom reactant, pos=[37.238987, 30.783043, 49.244064], name=R32, color=blue
pseudoatom reactant, pos=[37.92869, 30.940786, 49.324654], name=R33, color=blue
pseudoatom reactant, pos=[37.76586, 31.064255, 48.735508], name=R34, color=blue
pseudoatom reactant, pos=[38.01696, 31.045322, 48.89439], name=R35, color=blue
pseudoatom reactant, pos=[37.90504, 30.771255, 48.66058], name=R36, color=blue
pseudoatom reactant, pos=[37.24399, 30.840998, 48.88347], name=R37, color=blue
pseudoatom reactant, pos=[37.912003, 30.706266, 49.10741], name=R38, color=blue
pseudoatom reactant, pos=[38.093174, 31.128963, 49.33745], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[39.222733, 30.137453, 51.562073], name=P0, color=red
pseudoatom product, pos=[39.271847, 30.5822, 50.721645], name=P1, color=red
pseudoatom product, pos=[39.121536, 30.484772, 50.84261], name=P2, color=red
pseudoatom product, pos=[38.901684, 30.020786, 50.817383], name=P3, color=red
pseudoatom product, pos=[39.44332, 30.210781, 50.137505], name=P4, color=red
pseudoatom product, pos=[39.06652, 30.270325, 50.050446], name=P5, color=red
pseudoatom product, pos=[39.267174, 30.428629, 51.051453], name=P6, color=red
pseudoatom product, pos=[38.91593, 30.25263, 50.879642], name=P7, color=red
pseudoatom product, pos=[39.200977, 29.886564, 50.469296], name=P8, color=red
pseudoatom product, pos=[38.910366, 30.728828, 50.68206], name=P9, color=red
pseudoatom product, pos=[38.520916, 30.543161, 50.371346], name=P10, color=red
pseudoatom product, pos=[38.868694, 30.3043, 50.260002], name=P11, color=red
pseudoatom product, pos=[39.17855, 30.52887, 50.46941], name=P12, color=red
pseudoatom product, pos=[39.136448, 30.13974, 50.443214], name=P13, color=red
pseudoatom product, pos=[38.99015, 30.588137, 50.999874], name=P14, color=red
pseudoatom product, pos=[38.495037, 30.363375, 50.464832], name=P15, color=red
pseudoatom product, pos=[39.170002, 30.51642, 50.755943], name=P16, color=red
pseudoatom product, pos=[39.01092, 30.341955, 50.760754], name=P17, color=red
pseudoatom product, pos=[38.447018, 30.662264, 50.713707], name=P18, color=red
pseudoatom product, pos=[38.5975, 30.25712, 50.67283], name=P19, color=red
pseudoatom product, pos=[38.815353, 30.638649, 50.470062], name=P20, color=red
pseudoatom product, pos=[38.49759, 30.450125, 50.36163], name=P21, color=red
pseudoatom product, pos=[38.5266, 30.621426, 50.3025], name=P22, color=red
pseudoatom product, pos=[38.351643, 30.665909, 50.032616], name=P23, color=red
pseudoatom product, pos=[38.40134, 30.51308, 50.339176], name=P24, color=red
pseudoatom product, pos=[38.373077, 30.383673, 50.864735], name=P25, color=red
pseudoatom product, pos=[38.116264, 30.26055, 50.76915], name=P26, color=red
pseudoatom product, pos=[38.246723, 30.198626, 50.72693], name=P27, color=red
pseudoatom product, pos=[38.723373, 29.564972, 50.86086], name=P28, color=red
pseudoatom product, pos=[38.49843, 29.665714, 50.641743], name=P29, color=red
pseudoatom product, pos=[38.75732, 30.710903, 50.487858], name=P30, color=red
pseudoatom product, pos=[38.583862, 29.96756, 50.502205], name=P31, color=red
pseudoatom product, pos=[38.662483, 30.339008, 50.791607], name=P32, color=red
pseudoatom product, pos=[39.09199, 30.406162, 50.897793], name=P33, color=red
pseudoatom product, pos=[38.719494, 30.201038, 50.519424], name=P34, color=red
pseudoatom product, pos=[38.946808, 30.611408, 50.317677], name=P35, color=red
pseudoatom product, pos=[39.040623, 30.794037, 50.71036], name=P36, color=red
pseudoatom product, pos=[38.916847, 30.469604, 50.92475], name=P37, color=red
pseudoatom product, pos=[38.779633, 30.52988, 50.421894], name=P38, color=red
pseudoatom product, pos=[39.061348, 30.54949, 50.87317], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[37.61173027499999, 31.091562850000003, 48.97711442499999], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[38.82120307499999, 30.369753100000004, 50.62355412500001], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[38.21646667499999, 30.730657975000003, 49.800334275], label=Magnitude:2.17

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
cmd.wizard("message", "Dipole Magnitude: 2.1667\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
