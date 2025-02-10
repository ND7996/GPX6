from pymol import cmd

# Define a function to export distances to a file
def export_distances(filename):
    with open(filename, 'w') as f:
        # Define atom selections
        atom_pairs = [
            ("SEC`49/SE", "PRX`196/O1"),
            ("GLN`83/OE1", "PRX`196/O1"),
            ("SEC`49/SE", "SEC`49/HG1"),
            ("PRX`196/O1", "PRX`196/H1"),
            ("GLN`83/OE1", "PRX`196/H1"),
            ("PRX`196/O1", "SEC`49/HG1")
        ]

        # Calculate distance for each frame
        for frame in range(1, cmd.count_states() + 1):
            f.write(f'Frame {frame}:\n')
            for atom1, atom2 in atom_pairs:
                # Set the frame
                cmd.frame(frame)
                # Calculate distance
                distance_value = cmd.distance(atom1, atom2)
                # Write distance information to file
                f.write(f'Distance between {atom1} and {atom2}: {distance_value}\n')

# Call the function to export distances to a file
export_distances('distances.txt')