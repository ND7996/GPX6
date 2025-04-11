# PyMOL script for visualizing dipole in replica 010
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.55397, 30.58946, 49.360878], name=R0, color=blue
pseudoatom reactant, pos=[37.52232, 31.127272, 49.19856], name=R1, color=blue
pseudoatom reactant, pos=[37.661514, 30.480192, 48.674942], name=R2, color=blue
pseudoatom reactant, pos=[38.00585, 31.11481, 49.66373], name=R3, color=blue
pseudoatom reactant, pos=[38.44134, 30.656658, 49.14857], name=R4, color=blue
pseudoatom reactant, pos=[38.63764, 30.52961, 49.68055], name=R5, color=blue
pseudoatom reactant, pos=[37.711357, 30.61147, 48.75602], name=R6, color=blue
pseudoatom reactant, pos=[37.362885, 31.101591, 48.673405], name=R7, color=blue
pseudoatom reactant, pos=[37.277367, 31.351162, 48.54171], name=R8, color=blue
pseudoatom reactant, pos=[37.153797, 31.278814, 49.435566], name=R9, color=blue
pseudoatom reactant, pos=[37.66067, 31.204134, 49.475544], name=R10, color=blue
pseudoatom reactant, pos=[38.03445, 30.954456, 49.142178], name=R11, color=blue
pseudoatom reactant, pos=[37.90709, 30.916107, 48.70639], name=R12, color=blue
pseudoatom reactant, pos=[37.68187, 30.887177, 49.22406], name=R13, color=blue
pseudoatom reactant, pos=[38.676914, 30.890556, 49.07731], name=R14, color=blue
pseudoatom reactant, pos=[37.83258, 30.85198, 49.347168], name=R15, color=blue
pseudoatom reactant, pos=[37.560307, 30.824509, 48.788074], name=R16, color=blue
pseudoatom reactant, pos=[38.372517, 31.061502, 49.568962], name=R17, color=blue
pseudoatom reactant, pos=[37.68429, 31.406588, 48.650787], name=R18, color=blue
pseudoatom reactant, pos=[37.945457, 30.666084, 48.51189], name=R19, color=blue
pseudoatom reactant, pos=[37.96809, 30.649714, 49.1047], name=R20, color=blue
pseudoatom reactant, pos=[38.114227, 30.710196, 49.567963], name=R21, color=blue
pseudoatom reactant, pos=[38.365578, 30.854208, 49.106018], name=R22, color=blue
pseudoatom reactant, pos=[38.25529, 30.82629, 49.63688], name=R23, color=blue
pseudoatom reactant, pos=[37.881084, 30.621067, 48.987343], name=R24, color=blue
pseudoatom reactant, pos=[38.230698, 30.520596, 48.93934], name=R25, color=blue
pseudoatom reactant, pos=[37.696102, 31.42943, 48.501663], name=R26, color=blue
pseudoatom reactant, pos=[38.406677, 30.449322, 49.08948], name=R27, color=blue
pseudoatom reactant, pos=[38.10904, 30.949024, 49.056488], name=R28, color=blue
pseudoatom reactant, pos=[38.9874, 31.19928, 49.15166], name=R29, color=blue
pseudoatom reactant, pos=[38.10606, 30.61989, 49.301052], name=R30, color=blue
pseudoatom reactant, pos=[38.258224, 30.78955, 48.718468], name=R31, color=blue
pseudoatom reactant, pos=[38.510223, 30.2557, 49.068928], name=R32, color=blue
pseudoatom reactant, pos=[38.283745, 30.821266, 49.22324], name=R33, color=blue
pseudoatom reactant, pos=[39.12501, 31.478401, 48.975765], name=R34, color=blue
pseudoatom reactant, pos=[38.158184, 30.018896, 48.65285], name=R35, color=blue
pseudoatom reactant, pos=[38.501762, 30.334202, 48.939735], name=R36, color=blue
pseudoatom reactant, pos=[37.692223, 30.345112, 48.83314], name=R37, color=blue
pseudoatom reactant, pos=[38.192547, 30.46393, 48.396477], name=R38, color=blue
pseudoatom reactant, pos=[38.162067, 30.45423, 49.04619], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[38.47367, 30.084454, 51.384415], name=P0, color=red
pseudoatom product, pos=[37.779934, 30.240047, 51.989826], name=P1, color=red
pseudoatom product, pos=[38.17433, 29.822758, 51.382465], name=P2, color=red
pseudoatom product, pos=[37.52744, 30.365961, 51.192623], name=P3, color=red
pseudoatom product, pos=[37.952663, 30.200426, 51.26926], name=P4, color=red
pseudoatom product, pos=[38.170708, 29.886896, 52.037895], name=P5, color=red
pseudoatom product, pos=[38.144096, 29.996183, 51.013115], name=P6, color=red
pseudoatom product, pos=[37.691044, 30.392422, 51.378536], name=P7, color=red
pseudoatom product, pos=[38.172092, 30.04097, 51.568737], name=P8, color=red
pseudoatom product, pos=[37.731915, 29.990837, 51.4992], name=P9, color=red
pseudoatom product, pos=[37.844463, 30.252419, 51.55441], name=P10, color=red
pseudoatom product, pos=[37.437042, 30.655558, 51.49414], name=P11, color=red
pseudoatom product, pos=[37.684135, 29.529247, 51.725513], name=P12, color=red
pseudoatom product, pos=[37.88011, 29.811743, 51.10797], name=P13, color=red
pseudoatom product, pos=[37.849277, 29.97769, 51.030647], name=P14, color=red
pseudoatom product, pos=[37.088783, 30.228546, 51.699978], name=P15, color=red
pseudoatom product, pos=[37.68191, 29.973392, 51.44074], name=P16, color=red
pseudoatom product, pos=[37.57153, 30.221102, 51.487778], name=P17, color=red
pseudoatom product, pos=[37.911716, 29.814425, 51.15419], name=P18, color=red
pseudoatom product, pos=[38.31633, 30.148508, 51.18354], name=P19, color=red
pseudoatom product, pos=[38.55692, 29.885761, 51.49653], name=P20, color=red
pseudoatom product, pos=[38.684666, 30.194729, 51.54273], name=P21, color=red
pseudoatom product, pos=[38.545322, 29.923222, 51.07809], name=P22, color=red
pseudoatom product, pos=[38.77359, 29.968538, 51.083714], name=P23, color=red
pseudoatom product, pos=[37.85723, 30.180271, 51.120117], name=P24, color=red
pseudoatom product, pos=[38.602425, 30.321674, 51.43475], name=P25, color=red
pseudoatom product, pos=[38.265175, 30.384714, 51.451073], name=P26, color=red
pseudoatom product, pos=[38.27436, 30.218018, 51.07903], name=P27, color=red
pseudoatom product, pos=[38.129906, 29.75484, 51.024326], name=P28, color=red
pseudoatom product, pos=[38.26876, 30.39455, 50.882843], name=P29, color=red
pseudoatom product, pos=[39.06894, 30.408886, 50.373615], name=P30, color=red
pseudoatom product, pos=[38.738922, 30.092924, 50.9521], name=P31, color=red
pseudoatom product, pos=[38.1539, 30.080149, 50.889923], name=P32, color=red
pseudoatom product, pos=[38.601536, 29.75968, 50.800335], name=P33, color=red
pseudoatom product, pos=[38.47598, 30.009378, 50.960453], name=P34, color=red
pseudoatom product, pos=[37.962173, 30.27192, 50.597645], name=P35, color=red
pseudoatom product, pos=[38.64496, 29.836815, 51.10586], name=P36, color=red
pseudoatom product, pos=[38.39844, 30.454596, 51.338398], name=P37, color=red
pseudoatom product, pos=[37.757423, 30.460693, 51.1231], name=P38, color=red
pseudoatom product, pos=[37.584972, 29.667543, 51.57661], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.06721040000001, 30.807360900000003, 49.04809185], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[38.11071969999999, 30.097562125000003, 51.262655499999994], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[38.08896505, 30.4524615125, 50.15537367499999], label=Magnitude:2.33

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
cmd.wizard("message", "Dipole Magnitude: 2.3259\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
