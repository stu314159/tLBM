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
import numpy as np

parser = argparse.ArgumentParser(prog='psph_geom.py',
                                 description='create geometry files for prolate spheroid flow problem')
                                 
parser.add_argument('nDivs',type=int)

# parse input arguments
args = parser.parse_args()

#overall channel dimensions
aLx_p = 0.41
aLy_p = 0.41
aLz_p = 1.5
ab = 0.032
c = 0.216
aoa = np.pi/8.
aNdivs = args.nDivs

myObst = fc.ProlateSpheroid(aLx_p/2.,aLy_p/2.,aLz_p/2.,ab,c,aoa);

myChan = fc.FluidChannel(Lx_p=aLx_p,Ly_p=aLy_p,Lz_p=aLz_p,obst=myObst,
                         N_divs=aNdivs)                         

# write the mat file
myChan.write_mat_file('prolate_spheroid');

# write vtk of boundary conditions so you can visualize them
#myChan.write_bc_vtk();

myChan.set_pRef_indx(aLx_p/2.,aLy_p/2.,0.95*aLz_p);
#print "selected pressure reference is node number %d \n"%myChan.pRef_indx
#print "X = %g \n"%myChan.x[myChan.pRef_indx];
#print "Y = %g \n"%myChan.y[myChan.pRef_indx];
#print "Z = %g \n"%myChan.z[myChan.pRef_indx];
