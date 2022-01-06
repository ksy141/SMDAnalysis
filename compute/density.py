from __future__ import absolute_import, division, print_function
import numpy as np
from MDAnalysis.core.groups import AtomGroup
from ..common.block import Block

class Density:
    """
    Compute the density with respect to Z dimension
    >>> d = smda.Density().run(ag)
    """

    def __init__(self):
        pass

    def run(self, ag, nbins=100, b=0, e=None, skip=1, nblocks=1, center=False):
        """
        Run a density calculation

        Parameters
        ----------
        ag              Atomic Group
        b = 0           frame begins at b [int]
        e = None        frame ends   at e [int or None]
        skip = 1        every skip frame  [int]
        nblocks = 1     the number of blocks [int]
        center = False  centering based on this atomic group 
                        (e.g. Phosphorus atoms)

        Output
        ------
        [:,0] = bins [A]
        [:,1] = avg  [g/cm3]
        [:,2] = std  [g/cm3]
        """

        assert isinstance(ag, AtomGroup)
        if center: assert isinstance(center, AtomGroup)

        u = ag.universe

        print("frame begins at %d" %b)
        if e == None:
            nframes = u.trajectory.n_frames - 1
            print("frame ends at %d" %nframes)
        else:
            print("frame ends at %d" %e)
        
        z = np.zeros(nbins)
        ds = []
        nframes = 0

        for ts in u.trajectory[b:e:skip]:
            pos  = ag.positions[:,2]
            mass = ag.masses
            pbc  = u.dimensions[0:3]

            if center:
                ### CENTERING
                pos -= center.center_of_mass()[2]

                ### WRAP from -pbcz/2 to pbcz/2
                pos -= pbc[2] * np.around(pos/pbc[2])
                bins = np.linspace(-pbc[2]/2, pbc[2]/2, nbins + 1)
            
            else:
                ### WRAP from 0 to pbcz
                pos -= pbc[2] * np.floor(pos/pbc[2])
                bins = np.linspace(0, pbc[2], nbins + 1)

            z += 0.5 * (bins[1:] + bins[:-1])
            nframes += 1
            
            d = self.run_frame(pos, mass, pbc, bins)
            ds.append(d)

        z /= nframes
        avg, std = Block().block(ds, nblocks)
        return np.transpose([z, avg, std])


    def run_frame(self, pos, mass, pbc, bins):
        dz = bins[1] - bins[0]
        h, _ = np.histogram(pos, weights=mass, bins=bins)
        h /= pbc[0] * pbc[1] * dz * 0.602214
        return h



