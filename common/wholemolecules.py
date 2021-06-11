from __future__ import absolute_import, division, print_function

import numpy as np
from MDAnalysis import AtomGroup

class Wholemolecules:
    '''
    Make wholemolecules of selection
    Wholemolecules(AtomGroup)
    '''
    def __init__(self, atomgroup, point=None):
        '''
        Make wholemolecules based on the first atom or point
        '''
        assert isinstance(atomgroup, AtomGroup)
        pbc = atomgroup.dimensions[0:3]

        if point is None:
            point = atomgroup.positions[0]

        dr  = atomgroup.positions - point
        dr -= pbc * np.round(dr/pbc)
        atomgroup.positions = point + dr


