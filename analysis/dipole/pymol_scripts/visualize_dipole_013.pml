# PyMOL script for visualizing dipole in replica 013
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.489525, 30.572184, 48.875412], name=R0, color=blue
pseudoatom reactant, pos=[38.770195, 30.34363, 49.155754], name=R1, color=blue
pseudoatom reactant, pos=[37.698963, 30.651264, 48.49222], name=R2, color=blue
pseudoatom reactant, pos=[38.481834, 30.450808, 48.927002], name=R3, color=blue
pseudoatom reactant, pos=[38.086086, 30.339516, 49.65275], name=R4, color=blue
pseudoatom reactant, pos=[38.32587, 30.242384, 49.336544], name=R5, color=blue
pseudoatom reactant, pos=[38.263966, 30.060091, 49.270283], name=R6, color=blue
pseudoatom reactant, pos=[37.850258, 30.588036, 49.697895], name=R7, color=blue
pseudoatom reactant, pos=[38.300774, 30.238665, 49.77116], name=R8, color=blue
pseudoatom reactant, pos=[38.453087, 30.367111, 49.394028], name=R9, color=blue
pseudoatom reactant, pos=[38.16737, 30.45452, 48.60214], name=R10, color=blue
pseudoatom reactant, pos=[38.36863, 30.134043, 49.39139], name=R11, color=blue
pseudoatom reactant, pos=[38.22715, 30.540535, 49.342197], name=R12, color=blue
pseudoatom reactant, pos=[38.055973, 30.05183, 49.214554], name=R13, color=blue
pseudoatom reactant, pos=[38.178413, 30.449438, 49.498363], name=R14, color=blue
pseudoatom reactant, pos=[37.89457, 30.051826, 49.772846], name=R15, color=blue
pseudoatom reactant, pos=[38.015083, 30.192587, 49.304886], name=R16, color=blue
pseudoatom reactant, pos=[37.589176, 30.578402, 48.17154], name=R17, color=blue
pseudoatom reactant, pos=[37.35981, 30.509283, 48.795567], name=R18, color=blue
pseudoatom reactant, pos=[38.66623, 30.456598, 49.41316], name=R19, color=blue
pseudoatom reactant, pos=[38.652874, 30.403727, 49.285717], name=R20, color=blue
pseudoatom reactant, pos=[38.72054, 31.063457, 49.461403], name=R21, color=blue
pseudoatom reactant, pos=[38.685337, 30.801357, 49.18436], name=R22, color=blue
pseudoatom reactant, pos=[37.967655, 30.777288, 48.70208], name=R23, color=blue
pseudoatom reactant, pos=[38.183163, 30.94185, 48.816036], name=R24, color=blue
pseudoatom reactant, pos=[38.771202, 30.922604, 48.98492], name=R25, color=blue
pseudoatom reactant, pos=[38.401806, 30.713123, 49.10316], name=R26, color=blue
pseudoatom reactant, pos=[37.771645, 30.870441, 48.634857], name=R27, color=blue
pseudoatom reactant, pos=[38.49717, 30.983736, 48.718124], name=R28, color=blue
pseudoatom reactant, pos=[37.840347, 31.000994, 48.53662], name=R29, color=blue
pseudoatom reactant, pos=[37.626213, 31.110313, 48.441322], name=R30, color=blue
pseudoatom reactant, pos=[37.89433, 30.379454, 49.30984], name=R31, color=blue
pseudoatom reactant, pos=[38.483147, 30.747923, 48.891357], name=R32, color=blue
pseudoatom reactant, pos=[38.56212, 31.63558, 48.912193], name=R33, color=blue
pseudoatom reactant, pos=[38.167747, 30.611242, 49.08755], name=R34, color=blue
pseudoatom reactant, pos=[38.408207, 30.442295, 49.239964], name=R35, color=blue
pseudoatom reactant, pos=[38.232418, 30.260479, 49.214417], name=R36, color=blue
pseudoatom reactant, pos=[38.47597, 30.63898, 49.44402], name=R37, color=blue
pseudoatom reactant, pos=[38.987385, 30.43764, 49.746376], name=R38, color=blue
pseudoatom reactant, pos=[38.07279, 30.660751, 49.379616], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[36.115105, 30.430122, 49.473614], name=P0, color=red
pseudoatom product, pos=[36.179157, 30.46011, 49.527725], name=P1, color=red
pseudoatom product, pos=[36.006348, 30.173428, 49.938538], name=P2, color=red
pseudoatom product, pos=[35.8369, 30.446009, 49.5541], name=P3, color=red
pseudoatom product, pos=[35.54105, 30.935535, 49.158318], name=P4, color=red
pseudoatom product, pos=[35.98579, 30.531792, 49.45809], name=P5, color=red
pseudoatom product, pos=[35.87089, 30.518425, 49.73501], name=P6, color=red
pseudoatom product, pos=[36.317245, 30.457321, 49.062645], name=P7, color=red
pseudoatom product, pos=[36.156628, 30.341024, 49.456467], name=P8, color=red
pseudoatom product, pos=[35.598988, 30.488758, 49.312717], name=P9, color=red
pseudoatom product, pos=[35.581, 29.986366, 49.466484], name=P10, color=red
pseudoatom product, pos=[35.675934, 30.850367, 49.24903], name=P11, color=red
pseudoatom product, pos=[36.01731, 30.576927, 49.377094], name=P12, color=red
pseudoatom product, pos=[35.770054, 30.59103, 49.643764], name=P13, color=red
pseudoatom product, pos=[35.96324, 30.562311, 49.272343], name=P14, color=red
pseudoatom product, pos=[35.94452, 30.551886, 49.620758], name=P15, color=red
pseudoatom product, pos=[35.75704, 30.549154, 49.057552], name=P16, color=red
pseudoatom product, pos=[35.464523, 30.836956, 49.269466], name=P17, color=red
pseudoatom product, pos=[35.413242, 30.36154, 50.043015], name=P18, color=red
pseudoatom product, pos=[35.791046, 30.739225, 48.907867], name=P19, color=red
pseudoatom product, pos=[36.082844, 30.66125, 49.6485], name=P20, color=red
pseudoatom product, pos=[35.94341, 30.335648, 49.54274], name=P21, color=red
pseudoatom product, pos=[36.22105, 30.21426, 49.54006], name=P22, color=red
pseudoatom product, pos=[35.563118, 30.550482, 49.433086], name=P23, color=red
pseudoatom product, pos=[35.729954, 30.602577, 49.213352], name=P24, color=red
pseudoatom product, pos=[35.825928, 30.32249, 49.27223], name=P25, color=red
pseudoatom product, pos=[35.793606, 30.69254, 49.07179], name=P26, color=red
pseudoatom product, pos=[35.655052, 30.831427, 49.05364], name=P27, color=red
pseudoatom product, pos=[35.660225, 30.438742, 49.197945], name=P28, color=red
pseudoatom product, pos=[35.91336, 30.443323, 49.597454], name=P29, color=red
pseudoatom product, pos=[35.858917, 30.811426, 48.919487], name=P30, color=red
pseudoatom product, pos=[36.06854, 31.065504, 49.240166], name=P31, color=red
pseudoatom product, pos=[35.85179, 30.55267, 49.26949], name=P32, color=red
pseudoatom product, pos=[36.14857, 30.823305, 49.3339], name=P33, color=red
pseudoatom product, pos=[35.90038, 30.53927, 49.283035], name=P34, color=red
pseudoatom product, pos=[36.12645, 30.503819, 49.762463], name=P35, color=red
pseudoatom product, pos=[35.79439, 30.393366, 49.075733], name=P36, color=red
pseudoatom product, pos=[35.86774, 30.629118, 49.37518], name=P37, color=red
pseudoatom product, pos=[35.786198, 31.178364, 49.545296], name=P38, color=red
pseudoatom product, pos=[35.646362, 30.4977, 49.472137], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.241125725, 30.566899625, 49.129340575], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[35.86059735, 30.561889175000005, 49.38580702499999], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[37.0508615375, 30.564394400000005, 49.257573799999996], label=Magnitude:2.39

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
cmd.wizard("message", "Dipole Magnitude: 2.3943\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
