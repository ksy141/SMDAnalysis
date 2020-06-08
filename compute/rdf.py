from __future__ import absolute_import, division, print_function
import numpy as np
from MDAnalysis.analysis.distances import distance_array, self_distance_array
from MDAnalysis.core.groups import AtomGroup
from ..common.block import Block
from ..common.frame import Frame

class RDF:
    def __init__(self, nbins=100, limits=(0.0, 15.0)):
        self.rdf_settings = {'bins': nbins, 'range': limits}
        self.rmax = limits[1]
        
        _, edges = np.histogram([-1], **self.rdf_settings)
        self.bins  = 0.5 * (edges[1:] + edges[:-1])
        self.shell_vol  = 4.0/3.0 * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
        self.shell_area = np.pi * (np.power(edges[1:], 2) - np.power(edges[:-1], 2))

    def run(self, ag1, ag2, D = None, b=0, e=1000000, nblocks=1):
        assert isinstance(ag1, AtomGroup)
        assert isinstance(ag2, AtomGroup)
        u = ag1.universe
        
        bframe, eframe = Frame().frame(u, b, e)
        print("frame starts at %d" %bframe)
        print("frame ends   at %d" %eframe)

        rdfs = []
        for nframes, ts in enumerate(u.trajectory[bframe:eframe], 1):
            if D == 3:
                rdf = self.run3d_frame(ag1.positions, ag2.positions, u.dimensions)
            elif D == 2:
                rdf = self.run2d_frame(ag1.positions, ag2.positions, u.dimensions)
            else:
                assert 1==0, 'Specify D=2 or D=3'
            rdfs.append(rdf)
        
        print("total %d frames" %nframes)
        avg, std = Block().block(rdfs, nblocks)
        return np.transpose([self.bins, avg, std])

        
    def run3d_frame(self, g1_pos, g2_pos, dimensions):
        if np.all(g1_pos == g2_pos):
            self_rdf = True
        else:
            self_rdf = False

        N = len(g1_pos) * len(g2_pos)
        if N == 0:
            return np.zeros(len(self.bins))

        vol = dimensions[0] * dimensions[1] * dimensions[2]
        density = N / vol
        
        d = distance_array(g1_pos, g2_pos, box=dimensions)
        
        count = np.histogram(d[d!=0], **self.rdf_settings)[0]  
        count = count.astype(np.float64)
        rdf = count / density / self.shell_vol
        return rdf
    

    def run2d_frame(self, g1_pos, g2_pos, dimensions):
        if np.all(g1_pos == g2_pos):
            self_rdf = True
        else:
            self_rdf = False

        N = len(g1_pos) * len(g2_pos)
        if N == 0:
            return np.zeros(len(self.bins))

        g1_pos[:,2] = 0.0
        g2_pos[:,2] = 0.0

        area = dimensions[0] * dimensions[1]
        density = N / area
        
        d = distance_array(g1_pos, g2_pos, box=dimensions)

        count = np.histogram(d[d!=0], **self.rdf_settings)[0]  
        count = count.astype(np.float64)
        rdf = count / density / self.shell_area
        return rdf



