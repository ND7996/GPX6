import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Paths to the reactant and product CSV files (change this to your actual directory paths)
reactant_file = '/home/hp/results/MOUSE/level3/F48Y/replica000/reactant.csv'
product_file = 'home/hp/results/MOUSE/level3/F48Y/replica000/product.csv'

# Constants (Adjust as necessary for your system)
charges_used_reactant = {1: 0.1, 2: -0.1, 3: 0.2}  # Example: Atom index -> charge mapping for reactant
charges_used_product = {1: 0.1, 2: -0.1, 3: 0.2}  # Example: Atom index -> charge mapping for product

# PDB file (for loading and visualization)
pdb_file = '/home/hp/results/MOUSE/level3/F48Y/minim/minim.pdb'

def load_csv_data(csv_file):
    try:
        data = pd.read_csv(csv_file)
        print(f"Loaded {csv_file}, Columns: {data.columns}")
        print(data.head())  # Print first few rows for debugging
        return data
    except Exception as e:
        print(f"Error loading {csv_file}: {e}")
        return pd.DataFrame()

def calculate_com(coordinates):
    if len(coordinates) == 0:
        print("Warning: Empty coordinates for COM calculation.")
        return np.array([np.nan, np.nan, np.nan])
    return np.mean(coordinates, axis=0)

def calculate_dipole(coordinates, charges, com):
    if len(coordinates) == 0:
        print("Warning: No coordinates for dipole calculation.")
        return np.array([np.nan, np.nan, np.nan])
    dipole = np.zeros(3)
    for i in range(len(coordinates)):
        r_i = coordinates[i] - com  # Relative position vector
        dipole += charges.get(i + 1, 0) * r_i  # Ensure charge index matches
    return dipole

def process_replica(reactant_file, product_file):
    reactant_data = load_csv_data(reactant_file)
    product_data = load_csv_data(product_file)

    if {'x', 'y', 'z'}.issubset(reactant_data.columns):
        reactant_coord = reactant_data[['x', 'y', 'z']].values
    else:
        print(f"Warning: Missing coordinates in {reactant_file}")
        reactant_coord = np.array([])

    if {'x', 'y', 'z'}.issubset(product_data.columns):
        product_coord = product_data[['x', 'y', 'z']].values
    else:
        print(f"Warning: Missing coordinates in {product_file}")
        product_coord = np.array([])

    reactant_com = calculate_com(reactant_coord)
    product_com = calculate_com(product_coord)

    reactant_dipole = calculate_dipole(reactant_coord, charges_used_reactant, reactant_com)
    product_dipole = calculate_dipole(product_coord, charges_used_product, product_com)

    print(f"Reactant COM: {reactant_com}")
    print(f"Product COM: {product_com}")
    print(f"Reactant Dipole: {reactant_dipole}")
    print(f"Product Dipole: {product_dipole}")

# Example usage
process_replica('/home/hp/results/MOUSE/level3/F48Y/replica000/reactant.csv',
                '/home/hp/results/MOUSE/level3/F48Y/replica000/product.csv')
