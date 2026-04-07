import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
mpl.rcParams['lines.linewidth'] = 5
mpl.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt
font = {'family'    : 'sans-serif',
        'sans-serif': ['Tahoma'],
        'weight'    : 'normal',
        'size'      : 16}
mpl.rc('font', **font)

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument("-f",    help="filenames")
parser.add_argument("-col",  help="column index (1-based)")
parser.add_argument("-frac", help="the first fraction will be ignored e.g.0.2")
parser.add_argument("-draw", help="draw average, FWHM")
parser.add_argument("-nbins", help="nbins")
parser.set_defaults(frac=0)
parser.set_defaults(draw=False)
parser.set_defaults(nbins=200)

args = parser.parse_args()

cv = np.loadtxt(args.f)[:, int(args.col) - 1]
cv_equ = cv[int(len(cv) * float(args.frac)):]
cv_avg = np.average(cv_equ)

y, x = np.histogram(cv_equ, density=True, bins=int(args.nbins))
x = x[1:]

idx = np.abs(x -  np.average(cv_equ)).argmin()
y_avg = y[idx]

iidx = np.abs(y - y_avg/2).argmin()

FWHM = np.abs(cv_avg - x[iidx])

print("FWHM: {:.3f}".format(FWHM * 2))
print("X1={:.5f}  X2={:.5f}".format(cv_avg - FWHM, cv_avg + FWHM))

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(x, y, color='black')

if args.draw: 
    ax.axvline(cv_avg, color='C0', lw=2, linestyle='--')
    ax.axvline(cv_avg - FWHM, color='C0', lw=2, linestyle='--')
    ax.axvline(cv_avg + FWHM, color='C0', lw=2, linestyle='--')
plt.show()
#fig.savefig('histo.pdf')

