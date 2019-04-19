# calculate coordination number
import numpy as np
from ..common.frame import Frame
from MDAnalysis.analysis.distances import distance_array

class Coord:
    def __init__(self, u, sel1, sel2):
        self.u = u
        self.sel1 = sel1
        self.sel2 = sel2
        
    def run(self, nn=6, mm=12, d0=0, r0=0.25, b=0, e=10000, skip=1, backend='serial'):
        bframe, eframe = Frame().frame(self.u, b, e)
        ag1 = self.u.select_atoms(self.sel1)
        ag2 = self.u.select_atoms(self.sel2)
        times = []
        coords = []
        for ts in self.u.trajectory[bframe:eframe+1:skip]:
            times.append(ts.time/1000)
            d = distance_array(ag1.positions, ag2.positions, box=self.u.dimensions, backend=backend)/10
            D = (d - d0)/r0
            sij = (1 - np.power(D, nn))/(1 - np.power(D, mm))
            coords.append(np.sum(sij))
        
        output = np.array([times, coords])
        return np.transpose(output)


