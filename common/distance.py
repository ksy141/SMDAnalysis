from __future__ import absolute_import, division, print_function

import numpy as np
# from .wholemolecules import Wholemolecules
from MDAnalysis import AtomGroup

class Distance:
    def __init__(self, ag1, ag2):
        assert isinstance(ag1, AtomGroup)
        assert isinstance(ag2, AtomGroup)
        assert all(ag1.dimensions == ag2.dimensions)
        #Wholemolecules(ag1)
        #Wholemolecules(ag2)
        self.pbc  = ag1.dimensions[0:3]
        self.cag1 = ag1.center_of_geometry()
        self.cag2 = ag2.center_of_geometry()

    def actual_distance(self, pbc=None):
        dr = self.cag2 - self.cag1
        if pbc:
            dr -= self.pbc * np.round(dr/self.pbc)
            return dr
        else:
            return dr

    def scale_distance(self, pbc=None):
        self.cag1 -= self.pbc/2
        self.cag2 -= self.pbc/2
        dr = (self.cag2 - self.cag1)/self.pbc
        if pbc:
            dr -= np.round(dr)
            return dr
        else:
            return dr






