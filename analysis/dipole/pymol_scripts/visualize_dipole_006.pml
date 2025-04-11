# PyMOL script for visualizing dipole in replica 006
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.462803, 31.460754, 49.51355], name=R0, color=blue
pseudoatom reactant, pos=[37.448776, 30.795912, 48.918182], name=R1, color=blue
pseudoatom reactant, pos=[38.01334, 30.853092, 49.42219], name=R2, color=blue
pseudoatom reactant, pos=[37.9911, 31.01036, 49.28425], name=R3, color=blue
pseudoatom reactant, pos=[37.89256, 31.45475, 48.991512], name=R4, color=blue
pseudoatom reactant, pos=[37.937023, 30.384705, 49.76645], name=R5, color=blue
pseudoatom reactant, pos=[38.023838, 30.977156, 49.472576], name=R6, color=blue
pseudoatom reactant, pos=[37.717876, 30.889536, 49.325356], name=R7, color=blue
pseudoatom reactant, pos=[37.339794, 30.671526, 48.900604], name=R8, color=blue
pseudoatom reactant, pos=[37.447815, 31.039852, 49.164574], name=R9, color=blue
pseudoatom reactant, pos=[38.009033, 31.565876, 49.57087], name=R10, color=blue
pseudoatom reactant, pos=[38.20981, 31.661068, 49.44266], name=R11, color=blue
pseudoatom reactant, pos=[37.916355, 31.05902, 49.44399], name=R12, color=blue
pseudoatom reactant, pos=[38.533882, 31.04341, 49.559948], name=R13, color=blue
pseudoatom reactant, pos=[37.425694, 30.808012, 49.127293], name=R14, color=blue
pseudoatom reactant, pos=[37.84782, 31.111845, 49.385017], name=R15, color=blue
pseudoatom reactant, pos=[37.843704, 31.635372, 49.821262], name=R16, color=blue
pseudoatom reactant, pos=[38.06398, 31.365892, 49.891464], name=R17, color=blue
pseudoatom reactant, pos=[36.06709, 31.791101, 48.757133], name=R18, color=blue
pseudoatom reactant, pos=[36.65872, 31.56038, 48.79913], name=R19, color=blue
pseudoatom reactant, pos=[37.371456, 31.346636, 49.766026], name=R20, color=blue
pseudoatom reactant, pos=[37.41585, 30.664585, 49.561573], name=R21, color=blue
pseudoatom reactant, pos=[37.4474, 30.464357, 49.37139], name=R22, color=blue
pseudoatom reactant, pos=[38.39994, 30.627111, 49.559288], name=R23, color=blue
pseudoatom reactant, pos=[38.010956, 30.908049, 49.490425], name=R24, color=blue
pseudoatom reactant, pos=[37.328663, 30.715221, 48.932293], name=R25, color=blue
pseudoatom reactant, pos=[38.30758, 30.342487, 49.8806], name=R26, color=blue
pseudoatom reactant, pos=[38.503426, 31.160976, 49.560955], name=R27, color=blue
pseudoatom reactant, pos=[38.546944, 31.243242, 49.506947], name=R28, color=blue
pseudoatom reactant, pos=[37.798733, 30.553743, 49.38765], name=R29, color=blue
pseudoatom reactant, pos=[38.12447, 30.715565, 49.36793], name=R30, color=blue
pseudoatom reactant, pos=[38.544987, 30.9635, 49.753105], name=R31, color=blue
pseudoatom reactant, pos=[38.018513, 30.774284, 49.499203], name=R32, color=blue
pseudoatom reactant, pos=[38.11967, 30.970062, 49.699577], name=R33, color=blue
pseudoatom reactant, pos=[37.884716, 30.83789, 49.30565], name=R34, color=blue
pseudoatom reactant, pos=[37.71873, 30.832518, 48.778675], name=R35, color=blue
pseudoatom reactant, pos=[37.779194, 31.421679, 49.005623], name=R36, color=blue
pseudoatom reactant, pos=[37.906124, 31.457342, 49.021027], name=R37, color=blue
pseudoatom reactant, pos=[37.66226, 31.435661, 48.723057], name=R38, color=blue
pseudoatom reactant, pos=[37.97522, 31.009773, 48.75599], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[37.37148, 33.15791, 50.31975], name=P0, color=red
pseudoatom product, pos=[37.04934, 33.27363, 50.949863], name=P1, color=red
pseudoatom product, pos=[36.909187, 33.010002, 50.300823], name=P2, color=red
pseudoatom product, pos=[37.36718, 32.838665, 50.62434], name=P3, color=red
pseudoatom product, pos=[36.9977, 32.595, 49.81738], name=P4, color=red
pseudoatom product, pos=[36.98727, 33.20036, 50.04547], name=P5, color=red
pseudoatom product, pos=[37.008232, 32.76895, 50.54514], name=P6, color=red
pseudoatom product, pos=[37.408386, 33.195293, 50.669598], name=P7, color=red
pseudoatom product, pos=[37.044918, 33.46789, 50.36432], name=P8, color=red
pseudoatom product, pos=[37.198235, 33.559494, 50.95501], name=P9, color=red
pseudoatom product, pos=[36.801304, 33.327095, 50.669422], name=P10, color=red
pseudoatom product, pos=[36.99452, 33.420067, 50.952927], name=P11, color=red
pseudoatom product, pos=[36.87534, 33.21098, 50.241806], name=P12, color=red
pseudoatom product, pos=[37.111843, 33.37781, 51.31032], name=P13, color=red
pseudoatom product, pos=[36.96526, 33.457092, 50.67555], name=P14, color=red
pseudoatom product, pos=[37.035522, 33.3375, 50.929535], name=P15, color=red
pseudoatom product, pos=[37.402046, 33.186306, 50.7846], name=P16, color=red
pseudoatom product, pos=[37.217335, 33.095615, 50.66327], name=P17, color=red
pseudoatom product, pos=[37.272903, 33.48462, 50.41388], name=P18, color=red
pseudoatom product, pos=[37.20854, 33.485363, 50.368656], name=P19, color=red
pseudoatom product, pos=[37.151382, 33.26464, 50.574463], name=P20, color=red
pseudoatom product, pos=[37.236702, 33.246593, 50.61195], name=P21, color=red
pseudoatom product, pos=[37.45352, 32.977562, 50.932663], name=P22, color=red
pseudoatom product, pos=[37.451595, 33.27767, 51.07821], name=P23, color=red
pseudoatom product, pos=[37.202366, 33.41841, 51.079456], name=P24, color=red
pseudoatom product, pos=[37.435207, 33.29243, 51.530933], name=P25, color=red
pseudoatom product, pos=[38.073742, 32.76668, 51.04835], name=P26, color=red
pseudoatom product, pos=[37.030136, 33.232018, 50.69393], name=P27, color=red
pseudoatom product, pos=[36.921696, 33.166405, 50.534016], name=P28, color=red
pseudoatom product, pos=[37.150116, 33.176567, 51.040386], name=P29, color=red
pseudoatom product, pos=[37.19164, 33.13375, 50.70294], name=P30, color=red
pseudoatom product, pos=[37.64849, 33.323227, 50.778915], name=P31, color=red
pseudoatom product, pos=[37.97198, 33.373722, 51.265957], name=P32, color=red
pseudoatom product, pos=[37.894512, 33.2715, 50.86142], name=P33, color=red
pseudoatom product, pos=[37.23281, 33.173748, 51.228382], name=P34, color=red
pseudoatom product, pos=[36.814663, 33.488598, 50.76814], name=P35, color=red
pseudoatom product, pos=[37.19968, 33.362938, 50.66787], name=P36, color=red
pseudoatom product, pos=[37.528164, 34.009506, 50.85819], name=P37, color=red
pseudoatom product, pos=[37.84329, 33.310158, 50.7586], name=P38, color=red
pseudoatom product, pos=[37.189896, 33.24433, 50.46613], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[37.842896125, 31.0396075, 49.33712487499999], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[37.24620320000002, 33.24900235, 50.727064025], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[37.54454966250001, 32.144304925, 50.032094449999995], label=Magnitude:2.68

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
cmd.wizard("message", "Dipole Magnitude: 2.6776\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
