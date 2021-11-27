import numpy as np
import matplotlib
matplotlib.use('PDF')
params = {'font.size': 14, 'font.family': 'Times New Roman', 'mathtext.fontset': 'stix'}
matplotlib.rcParams.update(params)
import matplotlib.pyplot as plt

s = 8
m = '^'

fig, ax = plt.subplots(figsize=(3.5, 2.5))
ax.set_yscale('log')

d = np.loadtxt('PLacyl.dat')
ax.scatter(d[:,0], d[:,1], color='red',    s=s, marker=m, label='PLacyl')

d = np.loadtxt('TGacyl.dat')
ax.scatter(d[:,0], d[:,1], color='C1',     s=s, marker=m, label='TGacyl')

d = np.loadtxt('TGglyc.dat')
ax.scatter(d[:,0], d[:,1], color='purple', s=s, marker=m, label='TGglyc')

d = np.loadtxt('Deep.dat')
ax.scatter(d[:,0], d[:,1], color='pink',   s=s, marker=m, label='Deep')

ax.legend(frameon=False, fontsize=10)

ax.set_xlabel('defect size ($\AA$)')
ax.set_ylabel('probability')

fig.tight_layout()
fig.savefig('plot.pdf')


