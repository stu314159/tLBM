#!/usr/bin/env python
"""
generate a set of restart files for pyNFC based on
contents of a user-specified dump number.

Usage:
>>aprun ./pyNFC_create_restart.py dN

dN = # of velocity and pressure data files.
e.g. dN = 4 --> ux4.b_dat, uy4.b_dat, uz4.b_dat, density4.b_dat
will be loaded to create
restart_x, restart_y, restart_z, restart_rho


"""

import sys
sys.path.insert(1,'.')

import numpy as np
import argparse
import h5py

parser = argparse.ArgumentParser(prog="pyNFC_create_restart.py", 
                                description="restart file creation script.")
parser.add_argument('dumpNum',type=int)

# parse input arguments
args = parser.parse_args()

# assign to the variables
dN = args.dumpNum

rho_fn = 'density'+str(dN)+'.b_dat'
ux_fn = 'ux'+str(dN)+'.b_dat'
uy_fn = 'uy'+str(dN)+'.b_dat'
uz_fn = 'uz'+str(dN)+'.b_dat'

# load the node ordering

# Create numpy array from the binary data files
ux_i = np.fromfile(ux_fn,dtype=np.float32)
uy_i = np.fromfile(uy_fn,dtype=np.float32)
uz_i = np.fromfile(uz_fn,dtype=np.float32)
pressure_i = np.fromfile(rho_fn,dtype=np.float32)

order_map = np.fromfile('ordering.b_dat',dtype=np.int32).astype(np.int32)
# re-order per order_map
ux = np.zeros_like(ux_i); uy = np.zeros_like(uy_i); uz = np.zeros_like(uz_i);
pressure = np.zeros_like(pressure_i)
ux[order_map] = ux_i
uy[order_map] = uy_i
uz[order_map] = uz_i
pressure[order_map] = pressure_i

# write back out to files in correct order
f = h5py.File('restart.h5')
# Store velocity data into the velo_group of h5 file
velo_group = f.create_group("velocity")
x_velo = velo_group.create_dataset("x",data=ux)
y_velo = velo_group.create_dataset("y",data=uy)
z_velo = velo_group.create_dataset("z",data=uz)

# Store velocity data into the velo_group of h5 file
pres_group = f.create_group("density")
presmag = pres_group.create_dataset("rho",data=pressure)

f.close()
