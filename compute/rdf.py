from __future__ import absolute_import, division, print_function
import numpy as np
from MDAnalysis.core.groups import AtomGroup
# from MDAnalysis.core.universe import Universe
from MDAnalysis.analysis.distances import distance_array
from ..common.block import Block
from ..common.frame import Frame

class RDF:
    def __init__(self, g1, g2,
                 nbins=75, limits=(0.0, 1.5),
                 b=0, e=100000, skip=1,
                 serial=True,
                 nblocks = 5):
        
        self.self_rdf = False
        if g1 and g2:
            assert isinstance(g1, AtomGroup)
            assert isinstance(g2, AtomGroup)
        
            if g1 == g2:
                self.self_rdf = True

            self.g1 = g1
            self.g2 = g2
            self.u  = g1.universe
        
            if serial:
                bframe, eframe = Frame().frame(self.u, b, e)
            else:
                bframe = 0
                eframe = -1

            self.bframe  = bframe
            self.eframe  = eframe
            self.skip    = skip
            self.nblocks = nblocks
    

        self.rdf_settings = {'bins': nbins, 'range': limits}
        self.rmax = limits[1]
        
        _, edges = np.histogram([-1], **self.rdf_settings)
        self.bins  = 0.5 * (edges[1:] + edges[:-1])
        self.shell_vol  = 4.0/3.0 * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
        self.shell_area = np.pi * (np.power(edges[1:], 2) - np.power(edges[:-1], 2))

    def run(self, D = 3):
        print("frame starts at: %d" %self.bframe)
        print("frame ends   at: %d" %self.eframe)

        nframes = 0
        rdfs = []

        for ts in self.u.trajectory[self.bframe:self.eframe:self.skip]:
            if D == 3:
                rdf = self.run_frame(ts)
            elif D == 2:
                rdf = self.run2d_frame(ts)
            
            rdfs.append(rdf)
            nframes += 1
        
        print("total %d frames" %nframes)
        avg, std = Block().block(rdfs, self.nblocks)
        return np.transpose([self.bins, avg, std])


    def run_frame(self, ts, *args):
        ts
        if args:
            g1_pos = args[0]
            g2_pos = args[1]
            if g1_pos == g2_pos:
                self.self_rdf = True
        
        else:
            g1_pos = self.g1.positions
            g2_pos = self.g2.positions
        
        nA = len(g1_pos)
        nB = len(g2_pos)
        N = nA * nB
        
        vol = ts.volume / np.power(10, 3)
        density = N / vol

        d = distance_array(g1_pos, g2_pos, box=ts.dimensions)/10
        if self.self_rdf:
            np.fill_diagonal(d, self.rmax + 1)
        
        count = np.histogram(d, **self.rdf_settings)[0]  
        count = count.astype(np.float64)
        rdf = count / density / self.shell_vol
        return rdf
    

    def run2d_frame(self, ts, *args):
        ts
        if args:
            g1_pos = args[0]
            g2_pos = args[1]
            if g1_pos == g2_pos:
                self.self_rdf = True
 
        else:
            g1_pos = self.g1.positions
            g2_pos = self.g2.positions

        nA = len(g1_pos)
        nB = len(g2_pos)
        N = nA * nB

        area = (ts.dimensions[0] * ts.dimensions[1])/100
        density = N / area
        
        g1_pos[:,2] = 0.0
        g2_pos[:,2] = 0.0
            
        d = distance_array(g1_pos, g2_pos, box=ts.dimensions)/10
        if self.self_rdf:
            np.fill_diagonal(d, self.rmax + 1)
        
        count = np.histogram(d, **self.rdf_settings)[0]  
        count = count.astype(np.float64)
        rdf = count / density / self.shell_area
        return rdf

