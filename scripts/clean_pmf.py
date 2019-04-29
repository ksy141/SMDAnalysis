import numpy as np
import pandas as pd
import sys
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
kj_to_kcal = 0.239006 

pmf = {}
files = sys.argv[1:]
print(files)
fig, ax = plt.subplots()

for ifile in files:
    a = pd.read_csv(ifile, 
                    delim_whitespace=True, 
                    header=None, 
                    usecols=[0,1], 
                    comment='#')
    
    pmf[ifile] = {}
    pmf[ifile]['x'] = a[0].values
    pmf[ifile]['y'] = a[1].values
    pmf[ifile]['y'] -= pmf[ifile]['y'][0]
    pmf[ifile]['y'] *= kj_to_kcal
    ax.plot(pmf[ifile]['x'], pmf[ifile]['y'], label=ifile)

ax.legend()
plt.show()

