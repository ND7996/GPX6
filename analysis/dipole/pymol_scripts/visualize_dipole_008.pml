# PyMOL script for visualizing dipole in replica 008
reinitialize
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Create reactant structure
create reactant, none
pseudoatom reactant, pos=[38.81455, 30.819065, 49.729774], name=R0, color=blue
pseudoatom reactant, pos=[38.800777, 30.767605, 49.43851], name=R1, color=blue
pseudoatom reactant, pos=[38.65864, 30.642754, 49.265427], name=R2, color=blue
pseudoatom reactant, pos=[37.6219, 32.538128, 47.63701], name=R3, color=blue
pseudoatom reactant, pos=[38.428432, 30.674963, 49.204792], name=R4, color=blue
pseudoatom reactant, pos=[37.488613, 30.596693, 48.628525], name=R5, color=blue
pseudoatom reactant, pos=[38.384663, 30.515444, 49.892906], name=R6, color=blue
pseudoatom reactant, pos=[38.11486, 30.684269, 49.287952], name=R7, color=blue
pseudoatom reactant, pos=[38.21935, 30.419737, 49.554768], name=R8, color=blue
pseudoatom reactant, pos=[38.215145, 30.660631, 49.433975], name=R9, color=blue
pseudoatom reactant, pos=[38.074314, 30.317686, 49.004242], name=R10, color=blue
pseudoatom reactant, pos=[38.456825, 30.62644, 49.197327], name=R11, color=blue
pseudoatom reactant, pos=[38.042522, 30.822464, 49.60382], name=R12, color=blue
pseudoatom reactant, pos=[38.10754, 30.406057, 49.193302], name=R13, color=blue
pseudoatom reactant, pos=[38.44547, 30.800482, 49.283627], name=R14, color=blue
pseudoatom reactant, pos=[38.231346, 30.867727, 49.401024], name=R15, color=blue
pseudoatom reactant, pos=[38.00746, 30.851738, 49.286312], name=R16, color=blue
pseudoatom reactant, pos=[37.96716, 30.457506, 49.34063], name=R17, color=blue
pseudoatom reactant, pos=[38.14187, 30.857553, 49.15967], name=R18, color=blue
pseudoatom reactant, pos=[38.16832, 30.56785, 49.365704], name=R19, color=blue
pseudoatom reactant, pos=[38.28119, 30.547005, 48.824642], name=R20, color=blue
pseudoatom reactant, pos=[38.61727, 30.57191, 48.97508], name=R21, color=blue
pseudoatom reactant, pos=[38.322575, 30.592064, 49.33449], name=R22, color=blue
pseudoatom reactant, pos=[37.71164, 30.806583, 49.706383], name=R23, color=blue
pseudoatom reactant, pos=[38.52396, 30.79941, 49.183064], name=R24, color=blue
pseudoatom reactant, pos=[38.18102, 30.503607, 49.22789], name=R25, color=blue
pseudoatom reactant, pos=[37.649418, 30.206192, 49.265835], name=R26, color=blue
pseudoatom reactant, pos=[38.166527, 30.64301, 49.389656], name=R27, color=blue
pseudoatom reactant, pos=[38.246517, 30.7234, 49.86076], name=R28, color=blue
pseudoatom reactant, pos=[37.521545, 31.070705, 49.181465], name=R29, color=blue
pseudoatom reactant, pos=[37.89596, 30.316399, 49.292027], name=R30, color=blue
pseudoatom reactant, pos=[37.954247, 30.40287, 49.470524], name=R31, color=blue
pseudoatom reactant, pos=[37.99317, 30.7511, 49.308094], name=R32, color=blue
pseudoatom reactant, pos=[37.916317, 30.477678, 48.987003], name=R33, color=blue
pseudoatom reactant, pos=[37.693813, 30.842363, 48.705635], name=R34, color=blue
pseudoatom reactant, pos=[37.57725, 31.068287, 49.33621], name=R35, color=blue
pseudoatom reactant, pos=[38.016052, 30.638832, 49.470802], name=R36, color=blue
pseudoatom reactant, pos=[37.549713, 30.909666, 49.454475], name=R37, color=blue
pseudoatom reactant, pos=[38.694088, 30.66523, 49.107883], name=R38, color=blue
pseudoatom reactant, pos=[38.208683, 30.701317, 48.90867], name=R39, color=blue

# Create product structure
create product, none
pseudoatom product, pos=[39.504475, 30.282135, 51.509716], name=P0, color=red
pseudoatom product, pos=[39.827164, 29.642492, 52.092426], name=P1, color=red
pseudoatom product, pos=[39.43913, 29.708513, 51.88353], name=P2, color=red
pseudoatom product, pos=[39.683067, 30.203825, 51.228615], name=P3, color=red
pseudoatom product, pos=[39.35104, 30.144402, 51.63331], name=P4, color=red
pseudoatom product, pos=[39.216476, 29.942543, 51.73613], name=P5, color=red
pseudoatom product, pos=[39.159164, 30.068348, 50.77395], name=P6, color=red
pseudoatom product, pos=[39.644943, 29.899647, 50.850178], name=P7, color=red
pseudoatom product, pos=[39.415134, 29.059912, 51.432602], name=P8, color=red
pseudoatom product, pos=[39.289692, 30.029016, 51.46838], name=P9, color=red
pseudoatom product, pos=[39.23103, 29.528112, 51.55154], name=P10, color=red
pseudoatom product, pos=[39.331367, 29.91004, 50.909443], name=P11, color=red
pseudoatom product, pos=[39.27048, 29.761244, 51.100624], name=P12, color=red
pseudoatom product, pos=[38.91991, 29.721485, 50.733234], name=P13, color=red
pseudoatom product, pos=[38.935844, 29.727325, 51.02478], name=P14, color=red
pseudoatom product, pos=[39.202732, 29.651861, 50.954227], name=P15, color=red
pseudoatom product, pos=[39.111694, 29.886648, 50.94669], name=P16, color=red
pseudoatom product, pos=[38.991314, 29.517847, 51.169827], name=P17, color=red
pseudoatom product, pos=[39.40876, 29.827787, 51.155354], name=P18, color=red
pseudoatom product, pos=[39.11412, 29.72651, 51.125668], name=P19, color=red
pseudoatom product, pos=[39.764206, 30.10769, 51.043095], name=P20, color=red
pseudoatom product, pos=[39.264782, 29.749569, 50.71012], name=P21, color=red
pseudoatom product, pos=[38.86272, 29.617887, 50.820515], name=P22, color=red
pseudoatom product, pos=[39.01426, 29.898298, 50.304146], name=P23, color=red
pseudoatom product, pos=[39.146717, 30.003021, 49.89936], name=P24, color=red
pseudoatom product, pos=[39.206326, 29.88705, 50.300274], name=P25, color=red
pseudoatom product, pos=[39.802635, 29.799376, 49.987247], name=P26, color=red
pseudoatom product, pos=[39.535034, 30.008133, 51.064888], name=P27, color=red
pseudoatom product, pos=[39.118942, 29.674091, 51.275925], name=P28, color=red
pseudoatom product, pos=[38.73904, 29.503117, 51.04344], name=P29, color=red
pseudoatom product, pos=[39.544914, 29.931297, 51.073776], name=P30, color=red
pseudoatom product, pos=[39.38113, 30.14845, 51.280598], name=P31, color=red
pseudoatom product, pos=[38.857132, 29.974274, 50.81546], name=P32, color=red
pseudoatom product, pos=[38.66209, 30.121815, 51.159885], name=P33, color=red
pseudoatom product, pos=[38.73644, 29.631205, 51.1823], name=P34, color=red
pseudoatom product, pos=[38.809723, 29.767029, 51.485237], name=P35, color=red
pseudoatom product, pos=[39.18874, 30.158663, 51.624847], name=P36, color=red
pseudoatom product, pos=[39.169704, 29.66533, 50.767284], name=P37, color=red
pseudoatom product, pos=[39.469612, 29.738344, 50.487736], name=P38, color=red
pseudoatom product, pos=[38.907585, 29.954529, 51.030045], name=P39, color=red

# Create center of mass markers
pseudoatom reactant_com, pos=[38.1277678, 30.703310500000008, 49.24749712500001], color=marine, label=Reactant_COM
pseudoatom product_com, pos=[39.23073169999999, 29.839471500000002, 51.06591004999999], color=firebrick, label=Product_COM

# Create dipole vector
distance dipole, /reactant_com, /product_com
set dash_gap, 0, dipole
set dash_width, 6, dipole
set dash_color, purple, dipole

# Add dipole magnitude label
pseudoatom dipole_label, pos=[38.67924975, 30.271391000000005, 50.156703587500004], label=Magnitude:2.30

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
cmd.wizard("message", "Dipole Magnitude: 2.2955\nPress F1-F3 to switch between views\nF1: All atoms\nF2: Centers of mass only\nF3: Dipole vector only")

# Set function keys for scenes
cmd.set_key('F1', lambda: cmd.scene('all_atoms', 'recall'))
cmd.set_key('F2', lambda: cmd.scene('coms_only', 'recall'))
cmd.set_key('F3', lambda: cmd.scene('dipole_only', 'recall'))
