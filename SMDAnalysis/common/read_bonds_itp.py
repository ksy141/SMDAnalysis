import numpy as np
import os

class BBItp:
    def __init__(self):
        pass

    def read(self, fpath=''):
        if not os.path.exists(fpath):
            print('please provide a valid file path')
            return

        bonds = []
        readbonds = False
        with open(fpath) as W:
            for line in W.readlines():
                if line.startswith(';') or \
                        line == '\n' or \
                        line.startswith('#'):
                    continue
        
                if line.startswith('['):
                    if line.startswith('[ bonds') or \
                            line.startswith('[ cons'):
                        readbonds = True
                    else:
                        readbonds = False
                    continue
        
                if readbonds:
                    sl = line.split()
                    i  = int(sl[0])
                    j  = int(sl[1])
                    bonds.append([i, j])
        
        bonds = np.array(bonds)
        return bonds - 1

