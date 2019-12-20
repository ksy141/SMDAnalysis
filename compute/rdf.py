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

        assert isinstance(g1, AtomGroup)
        assert isinstance(g2, AtomGroup)
        
        if g1 == g2:
            self.self_rdf = True
        else:
            self.self_rdf = False

        self.g1 = g1
        self.g2 = g2
        self.u  = g1.universe

        self.rdf_settings = {'bins': nbins, 'range': limits}
        self.rmax = limits[1]
        
        if serial:
            bframe, eframe = Frame().frame(self.u, b, e)
            print("frame starts at: %d" %bframe)
            print("frame ends   at: %d" %eframe)

        else:
            bframe = 0
            eframe = -1

        self.bframe  = bframe
        self.eframe  = eframe
        self.skip    = skip
        self.nblocks = nblocks


    def run(self):
        _, edges = np.histogram([-1], **self.rdf_settings)
        bins  = 0.5 * (edges[1:] + edges[:-1])

        shell_vol = 4.0/3.0 * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
        nframes = 0
        rdfs = []

        for ts in self.u.trajectory[self.bframe:self.eframe:self.skip]:
            nA = len(self.g1)
            nB = len(self.g2)
            N  = nA*nB

            vol = ts.volume / np.power(10, 3)
            density = N / vol
            
            if self.self_rdf:
                d = distance_array(self.g1.positions, self.g2.positions, box=self.u.dimensions)/10
                np.fill_diagonal(d, self.rmax + 1)
            else:
                d = distance_array(self.g1.positions, self.g2.positions, box=self.u.dimensions)/10
            
            count = np.histogram(d, **self.rdf_settings)[0]  
            count = count.astype(np.float64)
            rdf = count / density / shell_vol

            rdfs.append(rdf)
            nframes += 1
        
        print("total %d frames" %nframes)
        avg, std = Block().block(rdfs, self.nblocks)
        
        return np.transpose([bins, avg, std])


    def run2d(self):
        _, edges = np.histogram([-1], **self.rdf_settings)
        bins  = 0.5 * (edges[1:] + edges[:-1])

        shell_vol = np.pi * (np.power(edges[1:], 2) - np.power(edges[:-1], 2))
        nframes = 0
        rdfs = []
        
        for ts in self.u.trajectory[self.bframe:self.eframe:self.skip]:
            nA = len(self.g1)
            nB = len(self.g2)
            N  = nA*nB
            
            area = (self.u.dimensions[0] * self.u.dimensions[1])/100
            density = N / area 

            g1_pos = self.g1.positions
            g1_pos[:,2] = 0.0

            g2_pos = self.g2.positions
            g2_pos[:,2] = 0.0
            
            if self.self_rdf:
                d = distance_array(g1_pos, g2_pos, box=self.u.dimensions)/10
                np.fill_diagonal(d, self.rmax + 1)
            else:
                d = distance_array(g1_pos, g2_pos, box=self.u.dimensions)/10
             
            count = np.histogram(d, **self.rdf_settings)[0]
            count = count.astype(np.float64)
            rdf = count / density / shell_vol

            rdfs.append(rdf)
            nframes += 1
        
        print("total %d frames" %nframes)
        avg, std = Block().block(rdfs, self.nblocks)

        return np.transpose([bins, avg, std])



        
