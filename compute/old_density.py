from __future__ import absolute_import, division, print_function
from MDAnalysis import Universe
from multiprocessing import Queue, Process, cpu_count
from ..common.block import Block
from ..common.frame import Frame
from ..common.dataclass import DataClass
import multiprocessing
import numpy as np
import time

mass = {'H': 1.008, 'O': 15.999, 'C': 12.011, 'N': 14.007, 
        'P': 30.976, 'S':32.06, 'SOD':22.99, 'CLA':35.45}

class OldDensity:
    def __init__(self):
        pass
    
    def density_mp(self, u, selection=False, nbins=100, 
                   nblocks=5, b=0, e=10000, cores=cpu_count()):

        if not selection:
            print("provide selection i.e. 'resname TIP3'")
        
        ### Analyze a fragment of trajectory in each core
        def worker(u, bb, ee, selection, nbins, queue, pos):
            #print(multiprocessing.current_process().name, 'Starting')
            
            zs = []
            densities = []
            group = u.select_atoms(selection)
        
            for ts in u.trajectory[bb:ee]:
                density = np.zeros(nbins)
                x, y, z = u.dimensions[0:3]
                zs.append(z)
                dz = z/nbins
            
                for atom in group:
                    index = int(atom.position[2]/dz) % nbins
                    if atom.name == 'SOD' or atom.name == 'CLA':
                        name = atom.name
                    else:
                        name = atom.name[0]
                        
                    density[index] += mass[name]
                density /= x*y*dz*0.602
                densities.append(density)
            
            dataclass = DataClass(pos, zs, densities)
            queue.put(dataclass)
            #queue.put(pos, zs, densities)
            #print(multiprocessing.current_process().name, 'Finishing')
        
        
        ### Get the closest begin and end frame corresponding to b, e
        bframe, eframe = Frame().frame(u, b, e)
        #print(bframe, eframe)
        eeframe = bframe + ((eframe - bframe)//cores) * cores
        fs = np.linspace(bframe, eeframe, cores+1, dtype=int)
        print(fs)
        us = []
        for process_number in range(cores):
            us.append(u.copy())
        
        queue = Queue()
        ps = [Process(target=worker, 
                      args=(us[i], fs[i], fs[i+1], selection, nbins, queue, i)) 
                      for i in range(cores)]
        for p in ps: p.start()
        ### first-finished, first-appeared
        ### need to sort so that result is in order
        

        results = []
        for p in ps:
            dc = queue.get()
            results.append((dc.data[0], dc.data[1], dc.data[2]))
        #results = [queue.get() for p in ps]
        results.sort()
        
        for p in ps: p.join()
        
        X = np.array([r[1] for r in results]).reshape(-1)
        # print(X)
        X = np.average(X)
        Y = np.array([r[2] for r in results]).reshape(-1, nbins)

        average, std = Block().block(Y, nblocks)
        print("Z (nm), density (kg/m^3), std (kg/m^3)")
        
        return np.linspace(0, X/10, num=nbins), average*1000, std*1000 # nm, kg/m3 kg/m3



    def density_frame(self, u, selection=False, nbins=100):
        if not selection:
            print("provide selection i.e. 'resname TIP3'")




    def density(self, u, selection=False, nbins=100, nblocks=5, b=0, e=10000):
        if not selection:
            print("provide selection i.e. 'resname TIP3'")
        
        group = u.select_atoms(selection)
        
        zs = []
        blocks = []
        
        bframe, eframe = Frame().frame(u, b, e)
        numframes = eframe - bframe
        size = numframes/nblocks

        for block in range(nblocks):
            densities = []
            start = bframe + int(block*size)
            end = bframe + int((block+1)*size)
            for ts in u.trajectory[start:end]:
                ptime = int(ts.time/1000)
                if ptime % 100 == 0: print("Now %4d ns " %ptime)
                density = np.zeros(nbins)
                x, y, z = u.dimensions[0:3]
                zs.append(z)
                dz = z/nbins

                for atom in group:
                    index = int(atom.position[2]/dz) % nbins
                    if atom.name == 'SOD' or atom.name == 'CLA':
                        name = atom.name
                    else:
                        name = atom.name[0]
                    density[index] += mass[name]
                density /= x*y*dz*0.602
                densities.append(density)
            
            blocks.append(np.average(densities, axis=0))

        blockaverage = np.average(blocks, axis=0)
        blockstd = np.std(blocks, axis=0)
        
        return np.linspace(0, np.average(zs), num=nbins), blockaverage, blockstd

