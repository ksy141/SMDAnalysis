from __future__ import absolute_import, division, print_function
from MDAnalysis.core.universe import Universe
import numpy as np


class PBC:
    '''
    Obtain PBC in A unit
    PBC.get_pbc will return dimensions
    '''
    def __init__(self, u):
        assert (u, Universe)
        pbcx, pbcy, pbcz, alpha, beta, gamma = u.dimensions
        self.pbc = np.array(pbcx, pbcy, pbcz)

    @property
    def get_pbc(self):
        return self.pbc
