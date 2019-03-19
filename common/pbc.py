from __future__ import absolute_import, division, print_function
from MDAnalysis import Universe
from MDAnalysis import AtomGroup
import numpy as np


class PBC:
    '''
    Obtain PBC in A unit
    PBC.get_pbc will return dimensions
    '''
    def __init__(self, *argv):
        arg = argv[0]
        assert (isinstance(arg, Universe) or isinstance(arg, AtomGroup))
        self.pbc = arg.dimensions[0:3]

    @property
    def get_pbc(self):
        return self.pbc
