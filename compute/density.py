from __future__ import absolute_import, division, print_function
from MDAnalysis import Universe
from ..common.block import Block
from ..common.frame import Frame
import numpy as np
import density_help
import time

mass = {'H': 1.008, 'O': 15.999, 'C': 12.011, 'N': 14.007, 
        'P': 30.976, 'S':32.06, 'SOD':22.99, 'CLA':35.45}

class Density:
    def __init__(self):
        pass
    
    def density(self, u, selection=False, nbins=100, 
                nblocks=5, b=0, e=10000):

        if not selection:
            print("provide selection i.e. 'resname TIP3'")
        
        start_time = time.time()
        group = u.select_atoms(selection)

        masses = []
        for atom in group:
            if atom.name == 'SOD' or atom.name == 'CLA':
                name = atom.name
            else:
                name = atom.name[0]
            masses.append(mass[name])
        masses = np.array(masses)
    
        
        ### Get the closest begin and end frame corresponding to b, e
        bframe, eframe = Frame().frame(u, b, e)
        print("frame starts at: %d" %bframe)
        print("frame ends   at: %d" %eframe)
        
        zs = []
        densities = []
        for ts in u.trajectory[bframe:eframe]:
        
            density = np.zeros(nbins, dtype=np.double)
            
            x, y, z = u.dimensions[0:3]
            dz = z/nbins
            zs.append(z)
        
            vecs = group.positions[:,2]/z
            vecs -= np.floor(vecs)
            indices = np.digitize(vecs, np.linspace(0, 1, nbins))
            
            density = density_help.density_help(indices, masses, density)
            density[0] = density[1]
            
            density /= x * y * dz * 0.602214
            densities.append(density)
        
        X = np.linspace(0, np.average(zs), num=nbins)
        X /= 10
        
        average, std = Block().block(densities, nblocks)
        
        print("--- %s seconds ---" % (time.time() - start_time))
        print("Z (nm), density (kg/m^3), std (kg/m^3)")
        
        return X, average*1000, std*1000 # nm, kg/m3 kg/m3



