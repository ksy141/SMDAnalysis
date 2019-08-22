import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

kj_to_kcal = 0.239006

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument("-f",    nargs='+', help="filenames e.g. fes_*.dat")
parser.add_argument("-o",    help="output filename default=pmf")
parser.add_argument("-zero", nargs='*', help="where to make it zero: default=-0.5")
parser.add_argument("-kcal", nargs='*', help="convert it to kcal e.g. True True False default: True True True")
parser.add_argument("-std",  help="whether include std or not default=True")
parser.add_argument("-sym",  help="symmetrize pmf default=False")
        
parser.set_defaults(zero=True)
parser.set_defaults(kcal=True)
parser.set_defaults(x=True)
parser.set_defaults(y=True)
parser.set_defaults(std=True)
parser.set_defaults(sym=False)
parser.set_defaults(o='pmf')

args = parser.parse_args()
# print(type(args.sym))
output = args.o
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


pmf  = {}
pmfs = []
files = f
print(files)

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

    if (pmf[ifile]['x'][0] > 0) and (pmf[ifile]['x'][-1]):
        pmf[ifile]['x'] *= -1

    idx = (np.abs(pmf[ifile]['x'] - zero)).argmin()
    pmf[ifile]['y'] -= pmf[ifile]['y'][idx]
    if kcal == 'True':
        pmf[ifile]['y'] *= kj_to_kcal

    x = pmf[ifile]['x']
    pmfs.append(pmf[ifile]['y'])


average = np.average(pmfs, axis=0)
std     = np.std(pmfs, axis=0)

if args.sym == 'True':
    average = (average + average[::-1])/2
    args.std = 'False'

np.savetxt(output, np.transpose([x, average, std, average-std, average+std]), fmt='%7.3f')
d = np.loadtxt(output)
fig, ax = plt.subplots()
ax.plot(d[:,0], d[:,1], color='C0')
ax.fill_between(d[:,0], d[:,3], d[:,4], color='C0', alpha=0.4)
fig.tight_layout()
fig.savefig(output + '.pdf')


if args.std != 'True':
    np.savetxt(output, np.transpose([x, average]), fmt='%7.3f')

