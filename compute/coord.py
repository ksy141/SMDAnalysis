from __future__ import absolute_import, division, print_function
import numpy as np
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.core.groups import AtomGroup

class Coord:
    """
    Compute the coordination numbers of atoms
    >>> c = smda.Coord().run(ag1, ag2, nn=6, mm=12, 
    ...                      d0=0[A], r0=2.5[A], 
    ...                      density=False, 
    ...                      b=0, e=None, skip=1)
    >>> c[:,0] = frame
    >>> c[:,1] = coordination number [unitless]

    Compute the contact map
    >>> c = smda.Coord().contact(ag1, ag2, rcut=5[A],
    ...                          density=False,
    ...                          b=0, e=None, skip=1)
    >>> c[:,0] = frame
    >>> c[:,1] = contact number [unitless]
    """

    def __init__(self):
        pass
        
    def run(self, ag1, ag2, nn=6, mm=12, d0=0, r0=2.5, density=False, b=0, e=None, skip=1):
        """
        Compute a coordination number
        s = [1 - ((r-d0)/r0)**n] / [1 - ((r-d0)/r0)**m] 

        Parameter
        ---------
        ag1, ag2: atomic groups
        density:  False [bool]
        nn = 6    [int]
        mm = 12   [int]
        d0 = 0    [A]
        r0 = 2.5  [A]
        b  = 0
        e  = None
        skip = 1

        Output
        ------
        [:,0] = frame
        [:,1] = coordination number [unitless]
        """

        assert isinstance(ag1, AtomGroup)
        assert isinstance(ag2, AtomGroup)
        u = ag1.universe

        times = []; coords = []
        for i, ts in enumerate(u.trajectory[b:e:skip]):
            times.append(i)
            d = distance_array(ag1.positions, ag2.positions, box=u.dimensions)
            d[ d < d0 ] = 1
            D = (d - d0)/r0
            sij = (1 - np.power(D, nn))/(1 - np.power(D, mm))
            coords.append(np.sum(sij))
        
        if density:
            coords = np.array(coords)
            coords /= (ag1.n_atoms * ag2.n_atoms)
        
        return np.transpose([times, coords])


    def contact(self, ag1, ag2, rcut, density=False, b=0, e=None, skip=1):
        """
        Compute the contact map
        s = 1 if d <= rcut
        s = 0 if d >  rcut

        Parameter
        ---------
        ag1, ag2: atomic groups
        density:  False [bool]
        rcut      [A]
        b  = 0    
        e  = None
        skip = 1

        Output
        ------
        [:,0] = frame
        [:,1] = contact number [unitless]
        """

        assert isinstance(ag1, AtomGroup)
        assert isinstance(ag2, AtomGroup)
        u = ag1.universe

        times = []; coords = []
        for i, ts in enumerate(u.trajectory[b:e:skip]):
            times.append(i)
            d = distance_array(ag1.positions, ag2.positions, box=u.dimensions)
            d[ d <= rcut ] = 1
            d[ d >  rcut ] = 0
            coords.append(np.sum(d))
        
        if density:
            coords = np.array(coords)
            coords /= (ag1.n_atoms * ag2.n_atoms)
        
        return np.transpose([times, coords])


