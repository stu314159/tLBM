#!/usr/bin/env python3

import sys
sys.path.insert(1,'../../python')

import FluidChannel as fc
import numpy as np

aLx_p = 1.0;
aLy_p = 1.0;
aLz_p = 5.0;
aNdivs = 4;

emptyChan = fc.FluidChannel(Lx_p = aLx_p,Ly_p = aLy_p,
                           Lz_p = aLz_p, N_divs = aNdivs);

emptyChan.write_mat_file('test_geom');
