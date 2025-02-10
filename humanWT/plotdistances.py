with open('distances.txt', 'r') as file:
    lines = file.readlines()

# Initialize a dictionary to store distances for each frame
distances = {}

# Parse the data
for line in lines:
    if line.startswith('Frame'):
        frame_number = line.split()[1].strip(':')
        distances[frame_number] = {}
    elif line.startswith('Distance between'):
        try:
            # Split the line by ':', then extract atom names and distance
            parts = line.split(':')
            atoms = parts[0].split()[3:5]
            atom1_name = atoms[0].split('`')[-1]  # Take the last part after splitting by '`'
            atom2_name = atoms[1].split('`')[-1]  # Take the last part after splitting by '`'
            distance = float(parts[1])
            distances[frame_number][(atom1_name, atom2_name)] = distance
        except (IndexError, ValueError) as e:
            print("Error parsing line:", line)

# Print the distances for each frame
for frame, frame_distances in distances.items():
    print(f"Frame {frame}:")
    for atom_pair, distance in frame_distances.items():
        print(f"    {atom_pair[0]} - {atom_pair[1]}: {distance}")
    print()

    # Print a text-based plot
    print("Distance Plot:")
    for atom_pair, distance in frame_distances.items():
        print(f"    {atom_pair[0]} - {atom_pair[1]}: {'=' * int(distance)}")
    print()