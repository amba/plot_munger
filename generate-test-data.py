#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('output', help="output filename")
args = parser.parse_args()

def save_3d_file(output_file, data, header):
    fh = open(output_file, 'w')
    fh.write(header + "\n")
    shape = data.shape
    for i in range(shape[0]):
        block = data[i]
        np.savetxt(fh, block, fmt="%.17g", delimiter="\t")
        fh.write("\n")
    fh.close()
    
N_x = 500
N_y = 500
x_min, x_max, y_min, y_max=[0, 30, 30, 0]
x = np.linspace(x_min, x_max, N_x)
y = np.linspace(y_min, y_max, N_y)
xv, yv = np.meshgrid(x,y)

zv = np.sin(xv - yv)

data = np.stack([yv, xv, zv],axis=-1)
save_3d_file(args.output, data , "# x\ty\tz")
