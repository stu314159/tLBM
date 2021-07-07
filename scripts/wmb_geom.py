#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 08:54:21 2017

@author: sblair
"""
import sys
sys.path.insert(1,'.')


import FluidChannel as fc
import argparse

parser = argparse.ArgumentParser(prog='wmb_geom.py',
                                 description='create geometry files for cavity channel problem')
                                 
parser.add_argument('nDivs',type=int)

# parse input arguments
args = parser.parse_args()

#overall channel dimensions
aLx_p = 6.4
aLy_p = 3.0
aLz_p = 7.0
aNdivs = args.nDivs

# wall mounted brick parameters
x_c = 3.5;
z_c = 3.2;
W = 1.;
H = W;
L = W;

myObst = fc.WallMountedBrick(x_c,z_c,L,W,H);

myChan = fc.FluidChannel(Lx_p=aLx_p,Ly_p=aLy_p,Lz_p=aLz_p,obst=myObst,
                         N_divs=aNdivs)                         

# write the mat file
myChan.write_mat_file('wall_mounted_brick');

# write vtk of boundary conditions so you can visualize them
#myChan.write_bc_vtk();

myChan.set_pRef_indx(aLx_p/2.,aLy_p/2.,0.95*aLz_p);
#print "selected pressure reference is node number %d \n"%myChan.pRef_indx
#print "X = %g \n"%myChan.x[myChan.pRef_indx];
#print "Y = %g \n"%myChan.y[myChan.pRef_indx];
#print "Z = %g \n"%myChan.z[myChan.pRef_indx];
