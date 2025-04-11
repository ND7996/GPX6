# PyMOL script for visualizing dipole in replica 009
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[37.95715, 30.596273, 49.26489], name=R0, color=blue
pseudoatom reactant, pos=[37.532814, 30.739843, 49.42573], name=R1, color=blue
pseudoatom reactant, pos=[37.967075, 31.030672, 49.554], name=R2, color=blue
pseudoatom reactant, pos=[36.808723, 31.75833, 48.32093], name=R3, color=blue
pseudoatom reactant, pos=[38.35744, 30.624205, 49.775425], name=R4, color=blue
pseudoatom reactant, pos=[38.09245, 30.308516, 49.68684], name=R5, color=blue
pseudoatom reactant, pos=[37.57955, 30.795868, 49.656544], name=R6, color=blue
pseudoatom reactant, pos=[38.156116, 30.96132, 49.66751], name=R7, color=blue
pseudoatom reactant, pos=[38.356453, 31.09227, 49.558876], name=R8, color=blue
pseudoatom reactant, pos=[38.044548, 30.609976, 49.362488], name=R9, color=blue
pseudoatom reactant, pos=[38.021503, 30.96028, 49.65995], name=R10, color=blue
pseudoatom reactant, pos=[37.825184, 31.081617, 49.28896], name=R11, color=blue
pseudoatom reactant, pos=[37.66465, 30.972366, 49.383064], name=R12, color=blue
pseudoatom reactant, pos=[38.01014, 30.852766, 49.53111], name=R13, color=blue
pseudoatom reactant, pos=[37.98594, 30.521112, 49.83437], name=R14, color=blue
pseudoatom reactant, pos=[37.937717, 30.926052, 49.53985], name=R15, color=blue
pseudoatom reactant, pos=[38.05535, 30.575909, 49.637405], name=R16, color=blue
pseudoatom reactant, pos=[38.21303, 30.793161, 49.62142], name=R17, color=blue
pseudoatom reactant, pos=[37.445766, 30.633028, 48.909206], name=R18, color=blue
pseudoatom reactant, pos=[38.307327, 30.73889, 49.620544], name=R19, color=blue
pseudoatom reactant, pos=[37.72099, 30.623528, 50.049538], name=R20, color=blue
pseudoatom reactant, pos=[38.098766, 30.897139, 50.038876], name=R21, color=blue
pseudoatom reactant, pos=[37.99344, 30.751469, 50.301548], name=R22, color=blue
pseudoatom reactant, pos=[37.34663, 31.36886, 48.309086], name=R23, color=blue
pseudoatom reactant, pos=[37.278095, 31.415527, 48.1868], name=R24, color=blue
pseudoatom reactant, pos=[38.24427, 30.706211, 50.044914], name=R25, color=blue
pseudoatom reactant, pos=[38.131134, 30.315922, 49.765118], name=R26, color=blue
pseudoatom reactant, pos=[37.818314, 30.594769, 50.39194], name=R27, color=blue
pseudoatom reactant, pos=[38.228218, 30.886467, 50.077408], name=R28, color=blue
pseudoatom reactant, pos=[38.213326, 30.92928, 49.920795], name=R29, color=blue
pseudoatom reactant, pos=[38.160225, 30.50442, 49.71879], name=R30, color=blue
pseudoatom reactant, pos=[38.330406, 30.823582, 49.389187], name=R31, color=blue
pseudoatom reactant, pos=[37.895496, 30.651758, 49.545513], name=R32, color=blue
pseudoatom reactant, pos=[39.144657, 30.526566, 50.147404], name=R33, color=blue
pseudoatom reactant, pos=[37.671085, 32.013405, 47.822647], name=R34, color=blue
pseudoatom reactant, pos=[36.980865, 31.563972, 48.000866], name=R35, color=blue
pseudoatom reactant, pos=[37.756577, 30.285433, 50.324944], name=R36, color=blue
pseudoatom reactant, pos=[38.22724, 30.894703, 49.74195], name=R37, color=blue
pseudoatom reactant, pos=[38.026936, 30.599274, 49.23302], name=R38, color=blue
pseudoatom reactant, pos=[38.387928, 31.065615, 49.754414], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[37.6376, 29.904272, 51.055782], name=P0, color=red
pseudoatom product, pos=[38.10081, 29.740126, 51.218777], name=P1, color=red
pseudoatom product, pos=[37.870728, 30.203835, 51.681667], name=P2, color=red
pseudoatom product, pos=[37.609535, 29.757599, 50.89556], name=P3, color=red
pseudoatom product, pos=[37.70133, 29.96502, 51.05641], name=P4, color=red
pseudoatom product, pos=[38.098747, 29.501375, 51.14122], name=P5, color=red
pseudoatom product, pos=[38.382378, 29.725492, 51.6797], name=P6, color=red
pseudoatom product, pos=[37.846962, 29.558052, 51.62823], name=P7, color=red
pseudoatom product, pos=[37.368603, 29.548925, 51.5064], name=P8, color=red
pseudoatom product, pos=[37.071404, 29.815617, 51.29952], name=P9, color=red
pseudoatom product, pos=[37.720894, 29.807251, 51.125114], name=P10, color=red
pseudoatom product, pos=[37.355938, 29.73494, 51.0355], name=P11, color=red
pseudoatom product, pos=[37.47155, 30.039011, 51.455433], name=P12, color=red
pseudoatom product, pos=[37.81222, 29.550127, 51.449097], name=P13, color=red
pseudoatom product, pos=[37.832203, 29.934685, 51.331398], name=P14, color=red
pseudoatom product, pos=[37.485645, 29.347898, 51.189224], name=P15, color=red
pseudoatom product, pos=[37.79029, 29.684202, 51.065567], name=P16, color=red
pseudoatom product, pos=[37.369526, 29.929893, 51.300198], name=P17, color=red
pseudoatom product, pos=[37.69261, 29.86882, 51.568962], name=P18, color=red
pseudoatom product, pos=[37.54938, 29.904188, 51.01091], name=P19, color=red
pseudoatom product, pos=[37.701702, 30.154037, 50.787273], name=P20, color=red
pseudoatom product, pos=[37.604557, 30.09719, 51.198967], name=P21, color=red
pseudoatom product, pos=[37.615566, 29.775028, 51.1799], name=P22, color=red
pseudoatom product, pos=[37.660084, 29.860712, 51.200912], name=P23, color=red
pseudoatom product, pos=[37.259262, 30.029297, 51.291622], name=P24, color=red
pseudoatom product, pos=[37.70925, 29.748787, 51.34359], name=P25, color=red
pseudoatom product, pos=[37.641205, 30.068544, 50.9193], name=P26, color=red
pseudoatom product, pos=[38.09487, 29.731546, 51.14493], name=P27, color=red
pseudoatom product, pos=[37.965633, 29.891443, 51.440727], name=P28, color=red
pseudoatom product, pos=[37.630226, 30.089544, 51.4156], name=P29, color=red
pseudoatom product, pos=[37.754143, 29.758825, 51.496716], name=P30, color=red
pseudoatom product, pos=[37.816246, 29.806728, 51.58608], name=P31, color=red
pseudoatom product, pos=[37.07856, 29.87843, 51.021896], name=P32, color=red
pseudoatom product, pos=[37.121002, 29.868544, 51.08359], name=P33, color=red
pseudoatom product, pos=[37.108734, 29.534174, 51.622074], name=P34, color=red
pseudoatom product, pos=[37.11891, 30.042192, 51.610474], name=P35, color=red
pseudoatom product, pos=[37.195686, 29.182291, 51.399708], name=P36, color=red
pseudoatom product, pos=[37.207386, 29.835493, 51.445187], name=P37, color=red
pseudoatom product, pos=[37.771137, 29.680538, 51.59391], name=P38, color=red
pseudoatom product, pos=[38.11336, 30.159292, 51.34015], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[37.94933810000001, 30.849758849999994, 49.501596750000004], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[37.62339680000001, 29.81784907500001, 51.295431875000006], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[37.786367450000014, 30.3338039625, 50.398514312500005], label=Magnitude:2.09

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
cmd.wizard("message", "Dipole Magnitude: 2.0950\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
