#!/usr/bin/env python
##!/home/users/sblair/anaconda2/bin/python
"""
Data processing script for binary data files.

Call in the directory of LBM input file and *.b_dat files.

Usage:
>>python processNFC.py

Will produce *.h5 storage files as well as *.xmf files to 
be read by Paraview.
"""

import sys
sys.path.insert(1,'.')

import math
import os
import numpy as np
from hdf5Helper import *

# Read data from params.lbm
input_file_name = 'params.lbm'
input_data = open(input_file_name,'r')

latticeType = str(input_data.readline())
dyamics = int(input_data.readline())
Num_ts = int(input_data.readline())
ts_rep_freq = int(input_data.readline())
Warmup_ts = int(input_data.readline())
plot_freq = int(input_data.readline())
Cs = float(input_data.readline())
rho_lbm = float(input_data.readline())
u_lbm = float(input_data.readline())
omega = float(input_data.readline())
Nx = int(input_data.readline())
Ny = int(input_data.readline())
Nz = int(input_data.readline())
Restart_flag = int(input_data.readline())
TimeAvg_flag = int(input_data.readline())
Lx_p = float(input_data.readline())
Ly_p = float(input_data.readline())
Lz_p = float(input_data.readline())
t_conv_fact = float(input_data.readline())
l_conv_fact = float(input_data.readline())
p_conv_fact = float(input_data.readline())
pRef_idx = int(input_data.readline())

input_data.close()

u_conv_fact = t_conv_fact/l_conv_fact;
nnodes = Nx*Ny*Nz
x = np.linspace(0.,Lx_p,Nx).astype(np.float64);
y = np.linspace(0.,Ly_p,Ny).astype(np.float64);
z = np.linspace(0.,Lz_p,Nz).astype(np.float64);
numEl = Nx*Ny*Nz
Y,Z,X = np.meshgrid(y,z,x);
XX = np.reshape(X,numEl)
YY = np.reshape(Y,numEl)
ZZ = np.reshape(Z,numEl)
dx = x[1] - x[0]

# compute the number of data dumps I expect to process
nDumps = (Num_ts-Warmup_ts)/plot_freq 

order_map = np.fromfile('ordering.b_dat',dtype=np.int32).astype(np.int32)

# Runs each data dump in serial
for i in range(int(nDumps)):
  rho_fn = 'density'+str(i)+'.b_dat'
  ux_fn = 'ux'+str(i)+'.b_dat'
  uy_fn = 'uy'+str(i)+'.b_dat'
  uz_fn = 'uz'+str(i)+'.b_dat'

  # Create numpy array from the binary data files
  ux_i = np.fromfile(ux_fn,dtype=np.float64)
  uy_i = np.fromfile(uy_fn,dtype=np.float64)
  uz_i = np.fromfile(uz_fn,dtype=np.float64)
  pressure_i = np.fromfile(rho_fn,dtype=np.float64)

  # Convert to physical units
  ux_i /= u_conv_fact
  uy_i /= u_conv_fact
  uz_i /= u_conv_fact
  pressure_i *= p_conv_fact 
    
# re-order per order_map
  ux = np.zeros_like(ux_i); uy = np.zeros_like(uy_i); uz = np.zeros_like(uz_i);
  pressure = np.zeros_like(pressure_i)
  ux[order_map] = ux_i
  uy[order_map] = uy_i
  uz[order_map] = uz_i
  pressure[order_map] = pressure_i
  pRef = pressure[pRef_idx];
  pressure -= pRef; # adjust for reference pressure
  
  velmag = np.sqrt(ux**2+uy**2+uz**2)

  # Create dimensions tuple for pressure reshape and XMF writer
  dims = (Nz,Ny,Nx)
  
  # Reshape pressure array
  # MUST BE DONE TO RUN PARAVIEW IN PARALLEL
  pressure = pressure.reshape(dims)
  velmag = velmag.reshape(dims)

  if i==0:
    print("Dimensions are",dims)
  print("Processing data dump #",i)

  # Write output files
  h5_file = 'out'+str(i)+'.h5'
  xmf_file = 'data'+str(i)+'.xmf'
  writeH5(pressure,ux,uy,uz,velmag,h5_file)
  writeXdmf(dims,dx,xmf_file,h5_file)

if TimeAvg_flag == 1:
  """
  process time average velocity and pressure data files so they can be visualized
  """
  rho_fn = 'rhoAvg.b_dat'
  ux_fn = 'uAvg.b_dat'
  uy_fn = 'vAvg.b_dat'
  uz_fn = 'wAvg.b_dat'

  # Create numpy array from the binary data files
  ux_i = np.fromfile(ux_fn,dtype=np.float64)
  uy_i = np.fromfile(uy_fn,dtype=np.float64)
  uz_i = np.fromfile(uz_fn,dtype=np.float64)
  pressure_i = np.fromfile(rho_fn,dtype=np.float64)
  
    
  ux_i /= Num_ts;
  uy_i /= Num_ts;
  uz_i /= Num_ts;
  pressure_i /= Num_ts;
  
  # Convert to physical units
  ux_i /= u_conv_fact
  uy_i /= u_conv_fact
  uz_i /= u_conv_fact
  pressure_i *= p_conv_fact 
    
# re-order per order_map
  ux = np.zeros_like(ux_i); uy = np.zeros_like(uy_i); uz = np.zeros_like(uz_i);
  pressure = np.zeros_like(pressure_i)
  ux[order_map] = ux_i
  uy[order_map] = uy_i
  uz[order_map] = uz_i
  pressure[order_map] = pressure_i
  pRef = pressure[pRef_idx];
  pressure -= pRef; # adjust for reference pressure
  
  velmag = np.sqrt(ux**2+uy**2+uz**2)

  # Create dimensions tuple for pressure reshape and XMF writer
  dims = (Nz,Ny,Nx)
  
  # Reshape pressure array
  # MUST BE DONE TO RUN PARAVIEW IN PARALLEL
  pressure = pressure.reshape(dims)
  velmag = velmag.reshape(dims)

  
  # Write output files
  h5_file = 'timeAvg.h5'
  xmf_file = 'timeAvg.xmf'
  writeH5(pressure,ux,uy,uz,velmag,h5_file)
  writeXdmf(dims,dx,xmf_file,h5_file)
