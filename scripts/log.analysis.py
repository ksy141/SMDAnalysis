import sys
import numpy as np

log = open(sys.argv[1], 'r')
var = sys.argv[2]

read = False
data = False

for line in log:
    if line.startswith('Step'):
        indices = line.split()
        read = True

        if not data:
            data = np.empty((0, len(indices)))
        continue

    elif line.startswith('Loop'):
        read = False

    if read:
        data = np.append(data, np.array([line.split()]).astype(np.float), axis=0)


if var in indices:
    index = indices.index(var)
else:
    print("No thermo variable called " + var)


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.plot(data[:,1], data[:,index])
plt.show(block=True)


