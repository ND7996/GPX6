# Master PyMOL script for dipole visualization
# This script allows you to load any replica's visualization

def load_replica(replica_num):
    cmd.reinitialize()
    replica = '{:03d}'.format(int(replica_num))
    script_path = cmd.exp_path('visualize_dipole_' + replica + '.pml')
    cmd.do('run ' + script_path)
    print('Loaded replica ' + replica)

# Create menu for loading replicas
import os
pymol_dir = cmd.exp_path('')
replica_scripts = [f for f in os.listdir(pymol_dir) if f.startswith('visualize_dipole_') and f.endswith('.pml')]
replica_nums = sorted([int(f.split('_')[2].split('.')[0]) for f in replica_scripts])

# Add menu items
for num in replica_nums:
    cmd.extend('load_replica_' + str(num), lambda n=num: load_replica(n))
    cmd.set_key('F' + str(min(num+4, 12)), lambda n=num: load_replica(n))

print('Master dipole visualization script loaded.')
print('Available replicas: ' + ', '.join([str(n) for n in replica_nums]))
print('Use load_replica(N) to load replica N')
print('Or use function keys F4-F12 to load replicas 0-8')
