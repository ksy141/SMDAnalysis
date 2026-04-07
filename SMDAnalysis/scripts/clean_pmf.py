import numpy as np
import pandas as pd
#import sys
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
kj_to_kcal = 0.239006 

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument("-f",    nargs='+', help="filenames e.g. fes_*.dat")
parser.add_argument("-zero", nargs='*', help="where to make it zero: default=-0.5")
parser.add_argument("-kcal", nargs='*', help="convert it to kcal e.g. True True False default: True True True")
parser.add_argument("-x",    nargs='*', help='xrange e.g. -0.5 0.5')
parser.add_argument("-y",    nargs='*', help='yrange e.g. -0.5 0.5')

parser.set_defaults(zero=True)
parser.set_defaults(kcal=True)
parser.set_defaults(x=True)
parser.set_defaults(y=True)

args = parser.parse_args()


f = args.f
zero = args.zero
if zero == True:
    zeros = [-0.5] * len(f)
else:
    zeros = [float(val) for val in zero]

kcal = args.kcal
if kcal == True:
    kcals = ['True'] * len(f)
else:
    kcals = kcal


pmf = {}
#files = sys.argv[1:]
files = f
print(files)
fig, ax = plt.subplots()

for ifile, zero, kcal in zip(files, zeros, kcals):
    a = pd.read_csv(ifile, 
                    delim_whitespace=True, 
                    header=None, 
                    usecols=[0,1], 
                    comment='#',
                    engine='python')
    
    pmf[ifile] = {}
    pmf[ifile]['x'] = a[0].values
    pmf[ifile]['y'] = a[1].values
    idx = (np.abs(pmf[ifile]['x'] - zero)).argmin()
    pmf[ifile]['y'] -= pmf[ifile]['y'][idx]
    if kcal == 'True':
        pmf[ifile]['y'] *= kj_to_kcal

    ax.plot(pmf[ifile]['x'], pmf[ifile]['y'], label=ifile, lw=3)

if args.x != True:
    ax.set_xlim([float(args.x[0]), float(args.x[1])])

if args.y != True:
    ax.set_ylim([float(args.y[0]), float(args.y[1])])

ax.legend()
plt.show()

