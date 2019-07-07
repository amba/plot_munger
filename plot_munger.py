#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import io
import code
import scipy.ndimage.filters
import os.path

if np.__version__ < '1.14.1':
    sys.exit("numpy version " + np.__version__ + " is too old")
    
def open_3d_file(file):
    fh = open(file, 'r')
    header = fh.readline().rstrip()
    contents = fh.read().rstrip()
    
    list_of_blocks = contents.split("\n\n")
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


def save_3d_file(output_file, data, header):
    fh = open(output_file, 'w')
    fh.write(header + "\n")
    shape = data.shape
    for i in range(shape[0]):
        block = data[i]
        np.savetxt(fh, block, fmt="%.17g", delimiter="\t")
        fh.write("\n")
    fh.close()
    




    
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--command', action='append', default=[], help=" transformation command, see below; can be provided multiple times")
parser.add_argument('-o', '--output', help="basename for output data file")
parser.add_argument('-f', '--force', action="store_true", help="overwrite existing files")
parser.add_argument('-x', '--xrange', help="xrange for plots")
parser.add_argument('-y', '--yrange', help="yrange for plots")
parser.add_argument('filename', help="input file")
parser.add_argument('OIZ', help="outer:inner:zcol description")
parser.add_argument('-g', '--grid', action="store_true", help="add grid in plot")
parser.add_argument('--cmap', help='name of color map (default: seismic)', default='seismic')
args = parser.parse_args()
cols = [int(x)-1 for x in args.OIZ.split(':')]
if len(cols) != len(set(cols)):
    sys.exit("outer, inner, and z columns need to be unique")
o_col, i_col, z_col = cols

data_3d, header = open_3d_file(args.filename)
print("input data shape: ", data_3d.shape)

col_legends = header.split()[1:]
col_dict = {}
for i, val in enumerate(col_legends):
    col_dict[i] = val

commands = []
for cmd in (args.command):
    if cmd == "diff_log":
        commands.extend('dzdi abs G0 log10'.split())
    else:
        commands.append(cmd)

i_block = data_3d[...,i_col]
o_block = data_3d[...,o_col]
z_block = data_3d[...,z_col]
z_label = col_dict[z_col]

def apply_commands(commands, z_label, blocks):
    for cmd in (commands):
        blocks, z_label = apply_command(cmd, z_label, blocks)
    return blocks[0], blocks[1], blocks[2], z_label

def apply_command(cmd, z_label, blocks):
    o_block = blocks[0]
    i_block = blocks[1]
    z_block = blocks[2]
    print("apply command ", cmd)
    if cmd in 'abs log log10'.split():
        z_block = getattr(np, cmd)(z_block)
        z_label = cmd + '(' + z_label + ')'
    elif cmd == 'dzdi':
        z_block = np.gradient(z_block,axis=1) / np.gradient(i_block,axis=1)
        z_label = 'dz/di(' + z_label + ')'
    elif cmd == 'dzdo':
        z_block = np.gradient(z_block,axis=0) / np.gradient(o_block,axis=o)
        z_label = 'dz/do(' + z_label + ')'
    elif cmd.startswith('_kernel=', 1):
        axis, kernel = cmd.split('_')
        tmp,kernel = cmd.split('=')
        kernel = np.array([float(x) for x in kernel.split(',')])
        kernel = kernel / kernel.sum()
        if axis == 'o':
            axis = 0
        elif axis == 'i':
            axis = 1
        else:
            sys.exit("unknown kernel axis " + axis)
        z_block = scipy.ndimage.filters.convolve1d(z_block, kernel, axis=axis)
    elif cmd == 'G0':
        z_block = z_block / 3.874045865410302e-05 # e**2 / h
        z_label = z_label + ' (e²/h)'
    elif cmd.startswith('min='):
        tmp,value = cmd.split('=')
        value = float(value)
        z_block = np.clip(z_block, value, None)
    elif cmd.startswith('max='):
        tmp,value = cmd.split('=')
        value = float(value)
        z_block = np.clip(z_block, None, value)
    elif cmd.startswith('add='):
        tmp,value = cmd.split('=')
        value = float(value)
        z_block = z_block + value
        z_label = z_label + ("%g" % value)
    elif cmd.startswith('factor='):
        tmp,value = cmd.split('=')
        value = float(value)
        z_block = value * z_block
        z_label = ("%g • " % value) + z_label
    else:
        sys.exit("unknown command " + cmd)
    return [o_block, i_block, z_block], z_label
        
o_block, i_block, z_block, z_label = apply_commands(commands, z_label, [o_block, i_block, z_block])

if args.output:
    output_header = "# %s\t%s\t%s" % (col_dict[o_col], col_dict[i_col], "transformed")
    output_filename = args.output + '_' + '_'.join(commands) + '.dat'
    if not args.force:
        if os.path.isfile(output_filename):
            sys.exit("file %s already exists. Use -f option to overwrite" % output_filename)
    print("writing output to ", output_filename)
    data = np.stack([o_block, i_block, z_block], axis=-1)
    save_3d_file(output_filename, data, output_header)



# need to flip the dimensions for imshow
o_min = o_block[0,0]
o_max = o_block[-1,0]
if o_min > o_max:
    o_max, o_min = o_min, o_max
    z_block = np.flip(z_block, axis=0)
    
i_min = i_block[0,0]
i_max = i_block[0,-1]
if i_min > i_max:
    i_max, i_min = i_min, i_max
    z_block = np.flip(z_block, axis=1)



# show outer axis left to right, inner axis bottom to top
z_block = np.flip(z_block, axis=1) # imshow plots the first axis top to bottom
z_block = np.swapaxes(z_block, 0, 1)

    
plt.imshow(z_block, aspect='auto', extent=[o_min, o_max, i_min, i_max],
           interpolation='none', cmap=args.cmap)
plt.xlabel(col_dict[o_col])
plt.ylabel(col_dict[i_col])
if args.grid:
    plt.grid()
plt.colorbar(format="%.1e", label=z_label)
plt.show(block=False)
code.interact()
