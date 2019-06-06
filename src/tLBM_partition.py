#!/usr/bin/env python3
##!/home/users/sblair/anaconda2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 14:23:52 2017

@author: stu
"""

import sys
sys.path.insert(1,'.')

import pyPartition as pp
#from pymetis import part_graph #<-- requires that the PrgEnv-intel module be selected
import numpy as np
import scipy.io
import math
import argparse

parser = argparse.ArgumentParser(prog='pyNFC_partition.py',
                                 description='lattice partitioning script for pyNFC')

parser.add_argument('geom_filename',type=str)
parser.add_argument('lattice_type',type=str)
parser.add_argument('partition_style',type=str)
parser.add_argument('numProcs',type=int)

# parse input arguments
args = parser.parse_args()

# assign to required variables
geom_filename = args.geom_filename
lattice_type = args.lattice_type
partition_style = args.partition_style
numProcs = args.numProcs

geom_input = scipy.io.loadmat(geom_filename)
# overall domain dimensions
Lx_p = float(geom_input['Lx_p'])
Ly_p = float(geom_input['Ly_p'])
Lz_p = float(geom_input['Lz_p'])
Lo = float(geom_input['Lo'])
Ny_divs = int(geom_input['Ny_divs'])
rho_p = float(geom_input['rho_p'])
nu_p = float(geom_input['nu_p'])


Ny = math.ceil((Ny_divs-1)*(Ly_p/Lo))+1
Nx = math.ceil((Ny_divs-1)*(Lx_p/Lo))+1
Nz = math.ceil((Ny_divs-1)*(Lz_p/Lo))+1
nnodes = Nx*Ny*Nz

# compute geometric data only once
x = np.linspace(0.,Lx_p,Nx).astype(np.float32);
y = np.linspace(0.,Ly_p,Ny).astype(np.float32);
z = np.linspace(0.,Lz_p,Nz).astype(np.float32);
numEl = Nx*Ny*Nz
Y,Z,X = np.meshgrid(y,z,x);

XX = np.reshape(X,int(numEl))
YY = np.reshape(Y,int(numEl))
ZZ = np.reshape(Z,int(numEl))




if lattice_type == 'D3Q15':
   lat = pp.D3Q15Lattice(int(Nx),int(Ny),int(Nz))
elif lattice_type == 'D3Q19':
   lat = pp.D3Q19Lattice(int(Nx),int(Ny),int(Nz))
else:
   lat = pp.D3Q27Lattice(int(Nx),int(Ny),int(Nz))


print("initializing the adjacency list")
lat.initialize_adjDict();
print("creating %s partition for %d processes" % (partition_style, numProcs))
lat.set_Partition(numParts= numProcs, style = partition_style)
lat.compute_cutSize()
print("cut size for %s partition = %g" % (partition_style, lat.get_cutSize()))
#print "writing vtk file for %s partition" % partition_style
#partition_vtk_filename = "partition_%s.vtk" % partition_style
#lat.partition.write_vtk(partition_vtk_filename)
print("writing %s partition to disk" % partition_style)
lat.partition.write_partition()
