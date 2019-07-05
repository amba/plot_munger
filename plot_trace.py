#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import io

def open_3d_file(file):
    print(np.__version__)
    contents = open(file, 'r').read().rstrip()
    
    list_of_blocks = contents.split("\n\n")
    print("number of blocks: ", len(list_of_blocks))
    arrays = []
    for block in (list_of_blocks):
        arrays.append(np.genfromtxt(io.StringIO(block)))
    first_shape = arrays[0].shape
    for i in range(len(arrays)-1, -1, -1):
        shape = arrays[i].shape
        if shape != first_shape:
            print("block ", i, " with first line", arrays[i][0], " does not match :", shape, " != ", first_shape)
            del arrays[i]
    return np.stack(arrays)


parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', help="basename for output data file")
parser.add_argument('-f', '--force', action="store_true", help="overwrite existing files")
parser.add_argument('-x', '--xrange', help="xrange for plots")
parser.add_argument('-y', '--yrange', help="yrange for plots")
parser.add_argument('filename', help="input file")
parser.add_argument('XYZ', help="xcol:ycol:zcol description")
parser.add_argument('trace', help="x=value or y=value")
args = parser.parse_args()

cols = [int(x)-1 for x in args.XYZ.split(':')]
if len(cols) != len(set(cols)):
    sys.exit("x, y, and z need to be unique")
    
trace_index, trace_value = args.trace.split('=')
if trace_index != 'x' and trace_index != 'y':
    sys.exit("trace argument need to be x=value or y=value")
    
trace_value = float(trace_value)

print("XYZ: ", cols)
print(trace_index,"=",trace_value)

data_3d = open_3d_file(args.filename)
print("input data shape: ", data_3d.shape)
data_3d = data_3d[..., cols]

print("3d_data.shape: ", data_3d.shape)

# assume regular data: x11 = x12 = x13, ...; y11 = y21 = y31, ...




