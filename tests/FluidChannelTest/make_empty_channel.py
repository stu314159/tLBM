#!/usr/bin/env python3
"""
create a simple empty channel with FluidChannel, 
write the *.mat file and the *.vtk files

"""

# ajust python path so required modules are visible
import sys
sys.path.insert(1,'../../src')


import FluidChannel as fc
import numpy as np

# overall channel dimensions
aLx_p = 1.0;
aLy_p = 1.0;
aLz_p = 5.0;
aNdivs = 11;

emptyChan = fc.FluidChannel(Lx_p = aLx_p, Ly_p = aLy_p,
                            Lz_p = aLz_p, N_divs = aNdivs);

emptyChan.write_mat_file('empty_chan');
emptyChan.write_bc_vtk();
