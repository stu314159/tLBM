#!/usr/bin/env python3

import argparse
import numpy as np
import numpy.linalg as la

parser = argparse.ArgumentParser(prog='validate.py',description='compare current data output to established gold standard')
parser.add_argument('dump_number',type=int)

args = parser.parse_args()

dump_num = args.dump_number

ux_fn = 'ux%d.b_dat'%dump_num
uy_fn = 'uy%d.b_dat'%dump_num
uz_fn = 'uz%d.b_dat'%dump_num

ux_i = np.fromfile(ux_fn, dtype=np.float32)
uy_i = np.fromfile(uy_fn, dtype=np.float32)
uz_i = np.fromfile(uz_fn, dtype=np.float32)

order_map = np.fromfile('ordering.b_dat',dtype=np.int32).astype(np.int32)
ux = np.zeros_like(ux_i)
uy = np.zeros_like(uy_i)
uz = np.zeros_like(uz_i)

ux[order_map] = ux_i
uy[order_map] = uy_i
uz[order_map] = uz_i



umag_tst = np.sqrt(ux**2+uy**2+uz**2)

umag_gold = np.load('gold_standard.npy')

dif_abs = la.norm((umag_tst - umag_gold),ord=1)
la_norm = la.norm(umag_gold,ord=1)
diff_rel= dif_abs / la_norm

print(la.norm(umag_tst,ord=1))
print('La norm %g.'%la_norm)
print('The relative 1-norm difference = %g.'%diff_rel)
print('the absolute 1-norm difference = %g.'%dif_abs)
