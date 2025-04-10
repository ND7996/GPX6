import numpy as np
import csv
import os

# Load the CSV files and handle missing or invalid data
def load_csv(csv_file):
    if not os.path.exists(csv_file):
        print(f"ERROR: CSV file not found: {csv_file}")
        return None, None
    
    coordinates = []
    charges = []  # We'll add atomic charges if available
    
    try:
        with open(csv_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)  # Get header row
            print(f"Header: {header}")  # Debugging: Print header
            
            # Check if charges are included (typically column 4 or 5)
            charge_col = None
            for i, col_name in enumerate(header):
                if col_name.lower() in ['charge', 'q']:
                    charge_col = i
                    break
            
            for i, row in enumerate(reader):
                print(f"Row {i}: {row}")  # Debugging: Print rows
                
                # Check if row has at least 4 fields with non-empty values
                if len(row) < 4:
                    print(f"Skipping row {i} due to insufficient number of columns: {row}")
                    continue
                
                # Check if the coordinates are missing or empty
                if not row[1] or not row[2] or not row[3]:
                    print(f"Skipping row {i} due to missing coordinate values: {row}")
                    continue
                
                try:
                    x, y, z = map(float, row[1:4])  # Assumes x, y, z are in columns 2, 3, 4
                except ValueError as ve:
                    print(f"Skipping row {i} due to invalid coordinates (cannot convert to float): {row} ({ve})")
                    continue  # Skip rows where coordinates can't be converted to float
                
                coordinates.append([x, y, z])
                
                # Add charge if available, otherwise use 1.0 as placeholder
                if charge_col is not None and len(row) > charge_col:
                    try:
                        charge = float(row[charge_col])
                    except ValueError:
                        charge = 1.0  # Default if charge can't be converted
                else:
                    charge = 1.0  # Default if charge column not found
                
                charges.append(charge)
                
    except Exception as e:
        print(f"ERROR reading CSV file {csv_file}: {str(e)}")
        return None, None
        
    return np.array(coordinates), np.array(charges)

# Center of mass calculation
def center_of_mass(coords):
    if coords.size == 0:
        print("ERROR: Coordinates array is empty.")
        return None
    return np.mean(coords, axis=0)

# Main function to process CSV files and calculate the dipole
def process_dipole(start_replica, end_replica):
    for replica_num in range(start_replica, end_replica + 1):  # Loop through user-defined range
        replica_str = f"{replica_num:03d}"  # Format as 000, 001, ..., 015
        
        # Construct file paths for reactant and product CSV files
        reactant_csv = f"/home/hp/results/MOUSE/level3/F48Y/replica{replica_str}/reactant.csv"
        product_csv = f"/home/hp/results/MOUSE/level3/F48Y/replica{replica_str}/product.csv"
        
        # Debugging: Print the file paths being processed
        print(f"Processing replica {replica_str}...")
        print(f"Reactant file path: {reactant_csv}")
        print(f"Product file path: {product_csv}")
        
        # Load reactant and product data
        reactant_coords, reactant_charges = load_csv(reactant_csv)
        product_coords, product_charges = load_csv(product_csv)
        
        # Check if data was loaded successfully
        if reactant_coords is None or product_coords is None:
            print(f"ERROR: Failed to load coordinate data for replica {replica_str}.")
            continue
        
        if reactant_coords.size == 0 or product_coords.size == 0:
            print(f"ERROR: No valid coordinate data found in CSV files for replica {replica_str}.")
            continue
        
        # Calculate the center of mass for both reactant and product
        reactant_com = center_of_mass(reactant_coords)
        product_com = center_of_mass(product_coords)
        
        if reactant_com is None or product_com is None:
            print(f"ERROR: Could not compute center of mass for replica {replica_str}.")
            continue
        
        print(f"Replica {replica_str}: Reactant Center of Mass: {reactant_com}")
        print(f"Replica {replica_str}: Product Center of Mass: {product_com}")

        # Additional calculations for dipoles, charges, etc. would go here
        # For now, just a placeholder for dipole-related calculations
        dipole_moment = np.sum((product_coords - reactant_coords) * (product_charges - reactant_charges)[:, None], axis=0)
        print(f"Replica {replica_str}: Dipole Moment: {dipole_moment}")

# User-defined range for replicas (can be changed manually)
start_replica = 0   # Change this to the starting replica number
end_replica = 15    # Change this to the ending replica number (for testing with replica 000, 001, 002)

# Run the main function with user-defined range
if __name__ == "__main__":
    process_dipole(start_replica, end_replica)
