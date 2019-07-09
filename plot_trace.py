#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import io
import code
import os.path
import scipy.optimize

if np.__version__ < '1.14.1':
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

# mandatory arguments
parser.add_argument('filename', help="input file")
parser.add_argument('OIZ', help="outer:inner:zcol description")
parser.add_argument('trace', help="o=value or i=value")

# output arguments
parser.add_argument('-o', '--output', help="basename for output data file")
parser.add_argument('-f', '--force', action="store_true", help="overwrite existing files")

# commands
parser.add_argument('-c', '--command', action='append', default=[], help=" transformation command, see below; can be provided multiple times")

# plot options
parser.add_argument('-s', '--save-plot', help='save plot to filename. Suffix determines the format')
parser.add_argument('--line', action='store_true')

# fit functions
parser.add_argument('--fit', help="fit type: 'linear', 'gaussian', 'lorentzian'") 

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
    x_col = i_col
    x_vals = data_3d[o_index,:,x_col]
    z_vals = data_3d[o_index,:,z_col]
    
elif trace_index == 'i':
    trace_col = i_col
    i_vals = data_3d[0, :, i_col]
    i_index, value = closest_index(i_vals, trace_value)
    x_col = o_col
    x_vals = data_3d[:,i_index, x_col]
    z_vals = data_3d[:,i_index, z_col]

x_label = col_dict[x_col]
z_label = col_dict[z_col]

def apply_commands(commands, x_vals, z_vals, x_label, z_label):
    for cmd in (commands):
        x_vals, z_vals, x_label, z_label = apply_command(cmd, x_vals, z_vals, x_label,
                                                z_label)
    return x_vals, z_vals, x_label, z_label

def apply_command(cmd, x_vals, z_vals, x_label, z_label):
    print("apply command", cmd)
    if cmd in 'abs log log10'.split():
        z_vals = getattr(np, cmd)(z_vals)
        z_label = cmd + '(' + z_label + ')'
    elif cmd.startswith('xmin='):
        tmp,value = cmd.split('=')
        value = float(value)
        mask = np.where(x_vals < value)
        x_vals = np.delete(x_vals, mask)
        z_vals = np.delete(z_vals, mask)
    elif cmd.startswith('xmax='):
        tmp,value = cmd.split('=')
        value = float(value)
        mask = np.where(x_vals > value)
        x_vals = np.delete(x_vals, mask)
        z_vals = np.delete(z_vals, mask)
    elif cmd.startswith('zmin='):
        tmp,value = cmd.split('=')
        value = float(value)
        mask = np.where(z_vals < value)
        x_vals = np.delete(x_vals, mask)
        z_vals = np.delete(z_vals, mask)
    elif cmd.startswith('zmax='):
        tmp,value = cmd.split('=')
        value = float(value)
        mask = np.where(z_vals > value)
        x_vals = np.delete(x_vals, mask)
        z_vals = np.delete(z_vals, mask)
    elif cmd.startswith('add='):
        tmp,value = cmd.split('=')
        value = float(value)
        z_vals = z_vals + value
        z_label = z_label + ("%g" % value)
    elif cmd.startswith('factor='):
        tmp,value = cmd.split('=')
        value = float(value)
        z_vals = value * z_vals
        z_label = ("%g • " % value) + z_label
    elif cmd == 'fft':
        z_vals = np.fft.rfft(z_vals)
        x_vals = np.fft.rfftfreq(x_vals.shape[0], np.abs(x_vals[1]-x_vals[0]))
        z_label = "|fft(%s)|" % z_label
        x_label = "freq(%s)" % x_label
    
    else:
        sys.exit("unknown command " + cmd)
    return x_vals, z_vals, x_label, z_label
    
x_vals, z_vals, x_label, z_label = apply_commands(args.command, x_vals, z_vals,
                                                  x_label, z_label)

if args.output:
    if not args.force and os.path.isfile(args.output):
        sys.exit("file %s already exists. Use -f option to overwrite" % args.output)
    output_block = np.stack([x_vals, z_vals], axis=-1)
    header = "# %s\t%s" % (x_label, z_label)
    np.savetxt(args.output, output_block, fmt="%.17g", header=header, comments='')


linestyle = "-" if args.line else ""
plt.plot(x_vals, z_vals, marker="x", linestyle=linestyle, label="%s=%g" %( col_dict[trace_col], value))

def gaussian(x, *p):
    x0, w, A, a, b = p
    return A * np.exp(-1/2 * ((x-x0)/w)**2) + a*x + b

def lorentzian(x, *p):
    x0, w, A, a, b = p
    return A /(w**2 + (x-x0)**2) + a*x + b

if args.fit:
    if args.fit == 'linear':
        coeff, V = np.polyfit(x_vals, z_vals, 1, cov=True)
        print("coeffs of linear fit: ", coeff)
        p = np.poly1d(coeff)
        cov = np.sqrt(np.diag(V))
        print("standard deviations: ", cov)
        label = "%.3g(±%.2g) • %s %+.3g(±%.2g)" % (coeff[0], cov[0], col_dict[x_col], coeff[1], cov[1])
        plt.plot(x_vals, p(x_vals), label=label)
    elif args.fit == 'gaussian':
        p0 = [(x_vals[-1] + x_vals[0])/2, 1, 1, 0, 0]
        popt, pcov = scipy.optimize.curve_fit(gaussian, x_vals, z_vals, p0=p0)
        print("fit parameters: ", popt)
        z_plot = gaussian(x_vals, *popt)
        label = 'gaussian(x_0 = %.4g, σ = %.3g)' % (popt[0], popt[1])
        plt.plot(x_vals, z_plot, label=label)
    elif args.fit == 'lorentzian':
        p0 = [(x_vals[-1] + x_vals[0])/2, 1, 1, 0, 0]
        popt, pcov = scipy.optimize.curve_fit(lorentzian, x_vals, z_vals, p0=p0)
        print("fit parameters: ", popt)
        z_plot = lorentzian(x_vals, *popt)
        label = 'lorentzian(x_0 = %.4g, w = %.3g)' % (popt[0], popt[1])
        plt.plot(x_vals, z_plot, label=label)
    else:
        sys.exit("unknown fit command %s" % args.fit)
plt.grid()
plt.xlabel(x_label)
plt.ylabel(z_label)
plt.legend()
plt.ticklabel_format(style='sci', axis='both')

if args.save_plot:
    if not args.force and os.path.isfile(args.save_plot):
        sys.exit("file %s already exists. Use -f option to overwrite" % args.save_plot)
    plt.savefig(args.save_plot, bbox_inches='tight')
plt.show(block=False)
code.interact()
