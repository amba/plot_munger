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


parser = argparse.ArgumentParser()

# mandatory arguments
parser.add_argument('filename', help="input file")
parser.add_argument('xy', help="x:y description")

# output arguments
parser.add_argument('-o', '--output', help="basename for output data file")
parser.add_argument('-f', '--force', action="store_true", help="overwrite existing files")

# commands
parser.add_argument('-c', '--command', action='append', default=[], help=" transformation commands: abs, log, log10, xmin=..., xmax=..., ymin=..., ymax=..., add=..., factor=..., fft, diff")

# plot options
parser.add_argument('-s', '--save-plot', help='save plot to filename. Suffix determines the format')
parser.add_argument('--line', action='store_true')

# fit functions
parser.add_argument('--fit', help="fit type: 'linear', 'gaussian', 'RLC', 'TIA'") 

args = parser.parse_args()
print(args)
cols = [int(x)-1 for x in args.xy.split(':')]
if len(cols) != len(set(cols)):
    sys.exit("x and y columns need to be unique")

x_col, y_col = cols

data_3d, header = open_3d_file(args.filename)
data_3d = data_3d[0,...]
col_legends = header.split()[1:]
col_dict = {}
print("legends:", col_legends)
for i, val in enumerate(col_legends):
    col_dict[i] = val

print("input data shape: ", data_3d.shape)

# assume regular data: o11 = o12 = o13, ...; i11 = i21 = i31, ...

x_vals = data_3d[...,x_col]
y_vals = data_3d[...,y_col]

x_label = col_dict[x_col]
y_label = col_dict[y_col]

def apply_commands(commands, x_vals, y_vals, x_label, y_label):
    for cmd in (commands):
        x_vals, y_vals, x_label, y_label = apply_command(cmd, x_vals, y_vals, x_label, y_label)
    return x_vals, y_vals, x_label, y_label

def apply_command(cmd, x_vals, y_vals, x_label, y_label):
    print("apply command", cmd)
    if cmd in 'abs log log10'.split():
        y_vals = getattr(np, cmd)(y_vals)
        y_label = cmd + '(' + y_label + ')'
    elif cmd == 'db_to_s':
        y_vals = 10**(y_vals / 20)
        y_label = 'db_to_s(' + y_label + ')'
    elif cmd.startswith('xmin='):
        tmp,value = cmd.split('=')
        value = float(value)
        mask = np.where(x_vals < value)
        x_vals = np.delete(x_vals, mask)
        y_vals = np.delete(y_vals, mask)
    elif cmd.startswith('xmax='):
        tmp,value = cmd.split('=')
        value = float(value)
        mask = np.where(x_vals > value)
        x_vals = np.delete(x_vals, mask)
        y_vals = np.delete(y_vals, mask)
    elif cmd.startswith('ymin='):
        tmp,value = cmd.split('=')
        value = float(value)
        mask = np.where(y_vals < value)
        x_vals = np.delete(x_vals, mask)
        y_vals = np.delete(y_vals, mask)
    elif cmd.startswith('ymax='):
        tmp,value = cmd.split('=')
        value = float(value)
        mask = np.where(y_vals > value)
        x_vals = np.delete(x_vals, mask)
        y_vals = np.delete(y_vals, mask)
    elif cmd.startswith('add='):
        tmp,value = cmd.split('=')
        value = float(value)
        y_vals = y_vals + value
        y_label = y_label + ("%g" % value)
    elif cmd.startswith('factor='):
        tmp,value = cmd.split('=')
        value = float(value)
        y_vals = value * y_vals
        y_label = ("%g • " % value) + y_label
    elif cmd == 'fft':
        y_vals = np.abs(np.fft.rfft(y_vals))
        x_vals = np.fft.rfftfreq(x_vals.shape[0], np.abs(x_vals[1]-x_vals[0]))
        y_label = "|fft(%s)|" % y_label
        x_label = "freq(%s)" % x_label
    elif cmd == 'diff':
        y_vals = np.gradient(y_vals) / np.gradient(x_vals)
        y_label = "d %s / d %s" % (y_label, x_label)
    else:
        sys.exit("unknown command " + cmd)
    return x_vals, y_vals, x_label, y_label
    
x_vals, y_vals, x_label, y_label = apply_commands(args.command, x_vals, y_vals,
                                                  x_label, y_label)

if args.output:
    if not args.force and os.path.isfile(args.output):
        sys.exit("file %s already exists. Use -f option to overwrite" % args.output)
    output_block = np.stack([x_vals, y_vals], axis=-1)
    header = "# %s\t%s" % (x_label, y_label)


linestyle = "-" if args.line else ""
plt.plot(x_vals, y_vals, marker="x", linestyle=linestyle)

def gaussian(x, *p):
    x0, w, A, a, b = p
    return A * np.exp(-1/2 * ((x-x0)/w)**2) + a*x + b


# 1/Z of R-L-C circuit Y = 1/(R + iωL + 1/(iωC))
#
# omega_0 = 1/sqrt(L * C)
# Q = 1/R * sqrt(L/C)
# max_val = 1/R

def RLC_Y(omega, omega_0, Q, max_val):
    return (max_val / Q) / (1/Q + 1j * (omega**2 - omega_0**2) / (omega * omega_0))

def damped_resonator(x, *p):
    x0, Q, max_val = p
    omega = 2 * np.pi * x
    omega_0 = 2 * np.pi * x0
    return np.abs(RLC_Y(omega, omega_0, Q, max_val))

# Butterworth van-Dyke model
def RLC_Cp(x, *p):
    x0, Q, max_val, p_val = p
    omega = 2 * np.pi * x
    omega_0 = 2 * np.pi * x0
    Y1 = RLC_Y(omega, omega_0, Q, max_val)
    Y2 = 1j * omega / omega_0 * p_val
    return np.abs(Y1 + Y2)


# return sensitvity V_0 / I_D
def TIA(x, *p):
    Rf, Cf, Rg, Cg, A_OL, omega_A = p
    omega = 2 * np.pi * x
    s = 1j * omega
    A = A_OL * omega_A / (s + omega_A)
    Zf = 1/(1/Rf + s * Cf)
    Zg = 1/(1/Rg + s* Cg)
    
    # V_0 / I_D
    return np.abs(-Zf / (1 + ((1 + Zf / Zg) / A)))
    
    


if args.fit:
    if args.fit == 'linear':
        coeff, V = np.polyfit(x_vals, y_vals, 1, cov=True)
        print("coeffs of linear fit: ", coeff)
        p = np.poly1d(coeff)
        cov = np.sqrt(np.diag(V))
        print("standard deviations: ", cov)
        label = "%.3g(±%.2g) • %s %+.3g(±%.2g)" % (coeff[0], cov[0], col_dict[x_col], coeff[1], cov[1])
        fit_vals = p(x_vals)
        plt.plot(x_vals, fit_vals, label=label)
    elif args.fit == 'gaussian':
        p0 = [(x_vals[-1] + x_vals[0])/2, 1, 1, 0, 0]
        popt, pcov = scipy.optimize.curve_fit(gaussian, x_vals, y_vals, p0=p0)
        print("fit parameters: ", popt)
        fit_vals = gaussian(x_vals, *popt)
        label = 'gaussian(x_0 = %.4g, σ = %.3g)' % (popt[0], popt[1])
        plt.plot(x_vals, fit_vals, label=label)
    elif args.fit == 'RLC':
        i_max = np.argmax(y_vals)
        f0 = x_vals[i_max]
        max_val = y_vals[i_max]
        Q = np.abs(f0 / (x_vals[-1] - x_vals[0]))
        print("Q = ", Q)
        p_val = (y_vals[-1] + y_vals[0]) / 2
        p0 = [f0, Q, max_val, p_val]
        popt, pcov = scipy.optimize.curve_fit(RLC_Cp, x_vals, y_vals, p0=p0)
        print("fit parameters: ", popt)
        fit_vals = RLC_Cp(x_vals, *popt)
        label = (
    'RLC(f_0 = %.1f Hz, Q = %.1f, res-gain-max = %.1f, parallel = %.1f)'
        % (popt[0], popt[1], popt[2], popt[3]))
        plt.plot(x_vals, fit_vals, label=label)
        ylabel = 'gain'
    elif args.fit == 'TIA':
        f0 = 1e6
        sigma = 0.01
        Rf = 1.6
        p0 = [Rf, f0, sigma]
        popt, pcov = scipy.optimize.curve_fit(TIA, x_vals, y_vals, p0=p0, maxfev=2000)
        fit_vals = TIA(x_vals, *popt)
        label = "TIA(Rf = %g, f0 = %g, sigma = %g, f0*sigma = %g, f0/sigma = %g)" % (popt[0], popt[1], popt[2], popt[1]*popt[2], popt[1] / popt[2])
        
        plt.plot(x_vals, fit_vals, label=label)
    else:
        sys.exit("unknown fit command %s" % args.fit)
    if args.output:
        fit_vals = np.expand_dims(fit_vals, axis=1)
        output_block = np.concatenate([output_block, fit_vals], axis=1)
        header += "\tfit"
        
plt.grid()
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.legend()
plt.ticklabel_format(style='sci', axis='both')

if args.output:
   np.savetxt(args.output, output_block, fmt="%.17g", header=header, comments='')
   
if args.save_plot:
    if not args.force and os.path.isfile(args.save_plot):
        sys.exit("file %s already exists. Use -f option to overwrite" % args.save_plot)
    plt.savefig(args.save_plot, bbox_inches='tight')
plt.show(block=False)
code.interact()
