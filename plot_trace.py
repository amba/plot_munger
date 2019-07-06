#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import io
import code

if np.__version__ < '1.16.4':
    sys.exit("numpy version ", np.__version__, "is too old")
    
def open_3d_file(file):
    print(np.__version__)
    fh = open(file, 'r')
    header = fh.readline().rstrip()
    print("header: ", header)
    contents = fh.read().rstrip()
    
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
    return np.stack(arrays), header

def closest_index(array, value):
    if len(array.shape) != 1:
        sys.exit("argument of closest_index needs dimension 1")
        
    abs_array = np.abs(array - value)
    index = np.where(abs_array == np.amin(abs_array))
    if len(index) > 1:
        print("warning: multiple optimal values in closest_index")
    index = index[0][0]
    value = array[index]
    return index, value

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', help="basename for output data file")
parser.add_argument('-f', '--force', action="store_true", help="overwrite existing files")
parser.add_argument('-x', '--xrange', help="xrange for plots")
parser.add_argument('-y', '--yrange', help="yrange for plots")
parser.add_argument('filename', help="input file")
parser.add_argument('OIZ', help="outer:inner:zcol description")
parser.add_argument('trace', help="o=value or i=value")
args = parser.parse_args()
print(args)
cols = [int(x)-1 for x in args.OIZ.split(':')]
if len(cols) != len(set(cols)):
    sys.exit("outer, inner, and z columns need to be unique")

o_col, i_col, z_col = cols
trace_index, trace_value = args.trace.split('=')
if trace_index != 'o' and trace_index != 'i':
    sys.exit("trace argument need to be o=value or i=value")
    
trace_value = float(trace_value)

data_3d, header = open_3d_file(args.filename)
col_legends = header.split()[1:]
col_dict = {}
print("legends:", col_legends)
for i, val in enumerate(col_legends):
    col_dict[i] = val

print("input data shape: ", data_3d.shape)

# assume regular data: o11 = o12 = o13, ...; i11 = i21 = i31, ...

if trace_index == 'o':
    trace_col = o_col
    o_vals = data_3d[:,0,o_col]
    o_index, value = closest_index(o_vals, trace_value)
    output_block = data_3d[o_index,:,:]
    x_col = i_col
    
elif trace_index == 'i':
    trace_col = i_col
    i_vals = data_3d[0, :, i_col]
    i_index, value = closest_index(i_vals, trace_value)
    output_block = data_3d[:, i_index, :]
    x_col = o_col

print("output_block shape: ", output_block.shape)
x_vals = output_block[...,x_col] 
z_vals = output_block[...,z_col]

if args.output:
    np.savetxt(args.output, output_block, fmt="%.17g", header=header, comments='')
    
plt.plot(x_vals, z_vals, label="%s=%g" %( col_dict[trace_col], value))
plt.grid()
plt.xlabel(col_dict[x_col])
plt.ylabel(col_dict[z_col])
plt.legend()
plt.show(block=False)
code.interact()
