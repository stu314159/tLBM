#!/usr/bin/env python3

import argparse
import numpy as np

parser = argparse.ArgumentParser(prog='gen_gold_standard.py',description='generate gold standard for validation')
parser.add_argument('dump_number',type=int)

args = parser.parse_args()

dump_num = args.dump_number

ux_fn = 'ux%d.b_dat'%dump_num
uy_fn = 'uy%d.b_dat'%dump_num
uz_fn = 'uz%d.b_dat'%dump_num

ux_i = np.fromfile(ux_fn, dtype=np.float32)
uy_i = np.fromfile(uy_fn, dtype=np.float32)
uz_i = np.fromfile(uz_fn, dtype=np.float32)

ux = np.zeros_like(ux_i)
uy = np.zeros_like(uy_i)
uz = np.zeros_like(uz_i)

order_map = np.fromfile('ordering.b_dat',dtype=np.int32).astype(np.int32)

ux[order_map] = ux_i
uy[order_map] = uy_i
uz[order_map] = uz_i

umag = np.sqrt(ux**2+uy**2+uz**2)

filename = 'gold_standard'
np.save(filename,umag)

