from __future__ import absolute_import, division, print_function
import numpy as np
from MDAnalysis.analysis.distances import distance_array, self_distance_array
from MDAnalysis.core.groups import AtomGroup
from ..common.block import Block
from ..common.frame import Frame

class RDF:
    """
    Compute the RDF
    r.bins, r.shell_vol, r.shell_area available
    >>> r = smda.RDF(nbins=100, limits=(0.0, 15.0)[A])

    1) For static atomic groups
    >>> rdf = r.run(ag1, ag2, D=2/3, 
    ...             b=0[ns], e=1e10[ns], nblocks=1)
    >>> rdf[:,0] = bins [A]
    >>> rdf[:,1] = RDF (avg) [unitless]
    >>> rdf[:,2] = RDF (std) [unitless]

    2) For dynamic atomic groups
    >>> data = []
    >>> for ts in u.trajectory:
    ...    g1_pos = g1.positions[SELECT]
    ...    g2_pos = g2.positions[SELECT]
    ...    data.append(r.run2d_frame(g1_pos, g2_pos, u.dimensions))
    ...    data.append(r.run3d_frame(g1_pos, g2_pos, u.dimensions))
    >>> rdf = np.transpose([r.bins, np.average(data, axis=0)])
    >>> rdf[:,0] = bins [A]
    >>> rdf[:,1] = RDF  [unitless]
    """
    
    def __init__(self, nbins=100, limits=(0.0, 15.0)):
        """
        Set up a RDF calculation.
        Define the number of bins and limits.

        Parameter
        ---------
        nbins  = 100  [int]
        limits = (0.0, 15.0) [A]
        """
        
        self.rdf_settings = {'bins': nbins, 'range': limits}
        self.rmax = limits[1]
        
        _, edges = np.histogram([-1], **self.rdf_settings)
        self.bins  = 0.5 * (edges[1:] + edges[:-1])
        self.shell_vol  = 4.0/3.0 * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
        self.shell_area = np.pi * (np.power(edges[1:], 2) - np.power(edges[:-1], 2))


    def run(self, ag1, ag2, D = None, b=0, e=1e10, skip=1, nblocks=1):
        """
        Run a RDF calculation for static atomic groups
        
        Parameter
        ---------
        ag1, ag2:   atomic groups
        D  = None   [int]
        b  = 0      [ns]
        e  = 1e10   [ns]
        skip = 1    [int]
        nblocks = 1 [int]

        Output
        ------
        [:,0] = bins [A]
        [:,1] = RDF (avg) [unitless]
        [:,2] = RDF (std) [unitless]
        """

        assert isinstance(ag1, AtomGroup)
        assert isinstance(ag2, AtomGroup)
        u = ag1.universe
        
        bframe, eframe = Frame().frame(u, b, e)
        print("frame starts at %d" %bframe)
        print("frame ends   at %d" %eframe)

        rdfs = []
        for ts in u.trajectory[bframe:eframe+1:skip]:
            if D == 3:
                rdf = self.run3d_frame(ag1.positions, ag2.positions, u.dimensions)
            elif D == 2:
                rdf = self.run2d_frame(ag1.positions, ag2.positions, u.dimensions)
            else:
                assert 1==0, 'Specify D=2 or D=3'
            rdfs.append(rdf)
        
        avg, std = Block().block(rdfs, nblocks)
        return np.transpose([self.bins, avg, std])

        
    def run3d_frame(self, g1_pos, g2_pos, dimensions):
        """
        Run a RDF calculation for one frame.
        Great flexibility as it takes atomic positions

        Parameters
        ----------
        g1_pos = g1.positions [A]
        g2_pos = g2.positions [A]
        dimensions = u.dimensions
        
        Outputs
        -------
        RDF array [A]
        """

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
        """
        Run a RDF calculation for one frame.
        Great flexibility as it takes atomic positions

        Parameters
        ----------
        g1_pos = g1.positions [A]
        g2_pos = g2.positions [A]
        dimensions = u.dimensions
        
        Outputs
        -------
        RDF array [A]
        """

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



